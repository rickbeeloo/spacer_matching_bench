"""
Results comparison and validation command (compares the tool outputs against the ground truth).
"""

import os
import logging
from typing import Dict, Optional
import polars as pl

# import polars_bio as pb
from rich.console import Console

from bench.utils.functions import (
    read_fasta,
    read_results,
    read_hyperfine_results,
    test_alignment,
    populate_pldf_withseqs_needletail,
    prettify_alignment,
    prettify_alignment_hamming,
    get_seq_from_fastx,
    vstack_easy,
    match_intervals_with_tolerance,
)
from bench.utils.tool_commands import load_tool_configs

# Logger is configured by cli.py with RichHandler
logger = logging.getLogger(__name__)
console = Console(soft_wrap=False)

# so printing to stdout doesn't break, line wrap, or truncate.
pl.Config.set_tbl_rows(123123)
pl.Config.set_tbl_cols(123123)  # should be large enough
pl.Config.set_fmt_str_lengths(2100)
pl.Config.set_tbl_width_chars(2100)


def merge_overlapping_intervals(
    df: pl.DataFrame, tolerance: int = 5
) -> pl.DataFrame:
    """
    Merge overlapping or nearby intervals within tolerance.
    
    Groups by spacer_id, contig_id, strand and merges intervals that overlap
    or are within tolerance bp of each other.
    
    Args:
        df: DataFrame with columns [spacer_id, contig_id, start, end, strand]
        tolerance: Maximum gap between intervals to merge (default: 5bp)
        
    Returns:
        DataFrame with merged intervals
    """
    if df.height == 0:
        return df
    
    merged_intervals = []
    
    # Group by spacer_id, contig_id, strand
    for group_key, group_df in df.group_by(["spacer_id", "contig_id", "strand"]):
        spacer_id, contig_id, strand = group_key
        
        # Sort by start position
        intervals = group_df.sort("start").select(["start", "end"]).to_numpy()
        
        if len(intervals) == 0:
            continue
            
        # Merge overlapping/nearby intervals
        current_start = intervals[0][0]
        current_end = intervals[0][1]
        
        for start, end in intervals[1:]:
            # Check if this interval overlaps or is within tolerance of current
            if start <= current_end + tolerance:
                # Merge: extend current interval
                current_end = max(current_end, end)
            else:
                # No overlap: save current and start new interval
                merged_intervals.append({
                    "spacer_id": spacer_id,
                    "contig_id": contig_id,
                    "start": current_start,
                    "end": current_end,
                    "strand": strand,
                })
                current_start = start
                current_end = end
        
        # Don't forget the last interval
        merged_intervals.append({
            "spacer_id": spacer_id,
            "contig_id": contig_id,
            "start": current_start,
            "end": current_end,
            "strand": strand,
        })
    
    return pl.DataFrame(merged_intervals)


def classify_unique_alignments_across_tools(
    ground_truth: pl.DataFrame,
    all_tool_results: pl.DataFrame,
    contigs_file: str,
    spacers_file: str,
    max_distance: int = 5,
    distance_metric: str = "hamming",
    gap_open_penalty: int = 5,
    gap_extend_penalty: int = 5,
    coordinate_tolerance: int = 5,
) -> pl.DataFrame:
    """
    Validate all unique alignments across all tools once.

    This consolidates validation to avoid redundant checking of the same alignment
    reported by multiple tools.

    Args:
        ground_truth: Ground truth DataFrame with planned insertions
        all_tool_results: Combined results from all tools
        contigs_file: Path to contigs FASTA file
        spacers_file: Path to spacers FASTA file
        max_distance: Maximum distance (mismatches/edit distance) to consider valid
        distance_metric: Distance metric for validation: 'hamming' (substitutions only) or 'edit' (substitutions + indels)
        gap_open_penalty: Gap open penalty for alignment validation, will be used if distance_metric is 'edit' (default: 5)
        gap_extend_penalty: Gap extension penalty for alignment validation, will be used if distance_metric is 'edit' (default: 5)
        coordinate_tolerance: Maximum bases difference in coordinates to consider as same alignment (default: 5)

    Returns:
        dataframe with alignment classifications.

    Note:
        this doens't add the ground truth results that have no matches in any tool results (these should be considered false negatives, and counted somewhere).
    """
    logger.debug("Validating unique alignments across all tools...")

    # Get unique alignments across all tools (by coordinates)
    unique_alignments = all_tool_results.select(
        ["spacer_id", "contig_id", "start", "end", "strand", "mismatches"]
    ).unique()

    logger.debug(f"Total unique alignments to validate: {unique_alignments.height}")

    # Add index for tracking
    unique_alignments = unique_alignments.with_row_index("alignment_idx")

    # Match tool results to ground truth allowing coordinate tolerance
    # This handles cases where tools report coordinates slightly differently
    # (e.g., 0-based vs 1-based, or boundary variations due to indels)
    matches_with_tolerance = match_intervals_with_tolerance(
        ground_truth=ground_truth,
        tool_results=unique_alignments,
        tolerance=coordinate_tolerance,
    )

    logger.debug(
        f"Matches within tolerance ({coordinate_tolerance}bp): {matches_with_tolerance.height}"
    )

    # Get alignments that didn't match ground truth within tolerance
    need_verify_df = unique_alignments.join(
        matches_with_tolerance.select(
            ["spacer_id", "contig_id", "start", "end", "strand"]
        ),
        on=["spacer_id", "contig_id", "start", "end", "strand"],
        how="anti",
    )
    # Categorize results not in ground truth
    logger.debug("Loading sequences for FP verification...")
    # Load spacer sequences
    need_verify_df = populate_pldf_withseqs_needletail(
        need_verify_df,
        seqfile=spacers_file,
        idcol="spacer_id",
        seqcol="spacer_seq",
        trim_to_region=False,
        reverse_by_strand_col=False,
    )

    # Load contig sequences
    need_verify_df = populate_pldf_withseqs_needletail(
        need_verify_df,
        seqfile=contigs_file,
        idcol="contig_id",
        seqcol="contig_seq",
        trim_to_region=True,
        reverse_by_strand_col=True,
    )

    logger.debug("Verifying false positive alignments...")
    logger.debug("Algorithm: Needleman-Wunsch (parasail)")
    logger.debug(f"Gap open penalty: {gap_open_penalty}")
    logger.debug(f"Gap extend penalty: {gap_extend_penalty}")
    logger.debug(f"Distance metric: {distance_metric})")
    logger.debug(f"Max allowed distance: {max_distance}")
    logger.debug(f"Alignments to verify: {need_verify_df.height}")

    # Test each non-planned alignment
    # contig_seq is already trimmed, so do not pass start/end
    need_verify_df = need_verify_df.with_columns(
        pl.struct(pl.col("spacer_seq"), pl.col("contig_seq"), pl.col("strand"))
        .map_elements(
            lambda x: test_alignment(
                contig_seq=x["contig_seq"],
                spacer_seq=x["spacer_seq"],
                strand=x["strand"],
                distance_metric=distance_metric,
                gap_cost=gap_open_penalty,
                extend_cost=gap_extend_penalty,
            ),
            return_dtype=pl.Int64,
        )
        .alias("recalculated_distnace")
    )
    # add classifications based on verification
    need_verify_df = need_verify_df.with_columns(
        pl.when(pl.col("recalculated_distnace") <= max_distance)
        .then(pl.lit("positive_not_in_plan"))
        .otherwise(pl.lit("invalid_alignment"))
        .alias("classification")
    )
    logger.debug("verification completed:")
    logger.debug(f"{need_verify_df['classification'].value_counts(sort=True)}")
    logger.debug("Combining matches with tolerance and verified alignments")

    # Prepare matched results - these already have classification from match_intervals_with_tolerance
    # They also have start_gt and end_gt columns
    matches_out = matches_with_tolerance.with_columns(
        pl.col("mismatches").alias("recalculated_distnace")
    ).select(
        [
            "spacer_id",
            "contig_id",
            "start",
            "end",
            "strand",
            "start_gt",
            "end_gt",
            "classification",
            "recalculated_distnace",
        ]
    )

    # For non-matches, there are no ground truth coordinates (they don't match any GT entry)
    need_verify_out = need_verify_df.with_columns(
        [
            pl.lit(None).cast(pl.Int64).alias("start_gt"),
            pl.lit(None).cast(pl.Int64).alias("end_gt"),
        ]
    ).select(
        [
            "spacer_id",
            "contig_id",
            "start",
            "end",
            "strand",
            "start_gt",
            "end_gt",
            "classification",
            "recalculated_distnace",
        ]
    )

    alignment_classifications = vstack_easy(matches_out, need_verify_out)

    alignment_classifications = (
        all_tool_results.join(
            alignment_classifications.select(
                [
                    "spacer_id",
                    "contig_id",
                    "start",
                    "end",
                    "strand",
                    "start_gt",
                    "end_gt",
                    "classification",
                    "recalculated_distnace",
                ]
            ),
            on=["spacer_id", "contig_id", "start", "end", "strand"],
            how="left",
        )
        .select(
            [
                "spacer_id",
                "contig_id",
                "start",
                "end",
                "strand",
                "start_gt",
                "end_gt",
                "classification",
                "recalculated_distnace",
            ]
        )
        .unique()
    )

    logger.debug("Validation of unique alignments completed.")
    return alignment_classifications


def display_example_alignments(
    alignment_classifications: pl.DataFrame,
    tools_results: pl.DataFrame,
    contigs_file: str,
    spacers_file: str,
    max_distance: int = 5,
    num_examples: int = 3,
    distance_metric: str = "hamming",
    gap_open_penalty: int = 5,
    gap_extend_penalty: int = 5,
):
    """
    Display example alignments for each classification type.

    Args:
        alignment_classifications: DataFrame with classified alignments
        tools_results: DataFrame with tool results to show which tools found each alignment
        contigs_file: Path to contigs FASTA file
        spacers_file: Path to spacers FASTA file
        max_mismatches: Maximum allowed mismatches
        num_examples: Number of examples to show per classification
        distance_metric: Distance metric for validation: 'hamming' (substitutions only) or 'edit' (substitutions + indels)
        gap_open_penalty: Gap open penalty for alignment
        gap_extend_penalty: Gap extend penalty for alignment
    Note:
        if the query region (spacer) is of a different length than the target region (area on the contig where it matched), the hamming distance is not meaningful.
    """
    logger.debug("[bold cyan]EXAMPLE ALIGNMENTS BY CLASSIFICATION[/bold cyan]")

    # Get examples for each classification type
    classifications_to_show = [
        "positive_not_in_plan",
        "invalid_alignment",
        "positive_in_plan",
    ]

    for classification in classifications_to_show:
        examples = alignment_classifications.filter(
            pl.col("classification") == classification
        ).head(num_examples)

        total_count = alignment_classifications.filter(
            pl.col("classification") == classification
        ).height

        if examples.height == 0:
            # Still log that this classification exists but has no examples
            classification_name = classification.upper().replace("_", " ")
            logger.debug(f"\n[dim]{classification_name} (0 found)[/dim]")
            continue

        logger.debug(
            f"\n[bold yellow]{classification.upper().replace('_', ' ')} ({total_count} total, showing {min(num_examples, examples.height)}):[/bold yellow]"
        )

        for idx, row in enumerate(examples.iter_rows(named=True), 1):
            spacer_id = row["spacer_id"]
            contig_id = row["contig_id"]
            start = row["start"]
            end = row["end"]
            strand = row["strand"]

            # Get sequences (only fetch what we need!)
            try:
                spacer_df = get_seq_from_fastx(
                    spacers_file,
                    [spacer_id],
                    return_df=True,
                    idcol="spacer_id",
                    seqcol="spacer_seq",
                )
                contig_df = get_seq_from_fastx(
                    contigs_file,
                    [contig_id],
                    return_df=True,
                    idcol="contig_id",
                    seqcol="contig_seq",
                )

                if spacer_df.height == 0 or contig_df.height == 0:
                    logger.warning(
                        f"Could not find sequences for {spacer_id} or {contig_id}"
                    )
                    continue

                spacer_seq = spacer_df["spacer_seq"][0]
                contig_seq = contig_df["contig_seq"][0]
            except Exception as e:
                logger.error(f"Error fetching sequences: {e}")
                continue

            # Find which tools reported this alignment
            tools_with_this = tools_results.filter(
                (pl.col("spacer_id") == spacer_id)
                & (pl.col("contig_id") == contig_id)
                & (pl.col("start") == start)
                & (pl.col("end") == end)
                & (pl.col("strand") == strand)
            )
            tool_names = (
                tools_with_this["tool"].unique().sort().to_list()
                if "tool" in tools_with_this.columns
                else []
            )

            # When fetching from FASTA, contig_seq is full-length, so pass start/end
            edit_distance = test_alignment(
                spacer_seq,
                contig_seq,
                strand=strand,
                start=start,
                end=end,
                distance_metric="edit",
                gap_cost=gap_open_penalty,
                extend_cost=gap_extend_penalty,
            )
            # Calculate hamming distance (substitutions only, no indels)
            hamming_distance = test_alignment(
                spacer_seq,
                contig_seq,
                strand=strand,
                start=start,
                end=end,
                distance_metric="hamming",
                gap_cost=gap_open_penalty,
                extend_cost=gap_extend_penalty,
            )

            # Determine if there are indels (gaps)
            has_indels = edit_distance != hamming_distance

            # Generate alignment visualization
            # Use hamming (ungapped) for hamming distance metric, NW (gapped) for edit distance
            if distance_metric == "hamming":
                alignment_str = prettify_alignment_hamming(
                    spacer_seq,
                    contig_seq,
                    strand=strand,
                    start=start,
                    end=end,
                )
            else:
                alignment_str = prettify_alignment(
                    spacer_seq,
                    contig_seq,
                    strand=strand,
                    start=start,
                    end=end,
                    gap_cost=gap_open_penalty,
                    extend_cost=gap_extend_penalty,
                )

            # Add labels to the alignment
            strand_str = "(-)" if strand else "(+)"
            location_str = f"{contig_id}:{start}-{end} {strand_str}"

            lines = alignment_str.split("\n")
            if len(lines) >= 3:
                # Add sequence IDs on the right side
                max_len = max(len(line) for line in lines)

                logger.debug(f"[bold]Example {idx}:[/bold]")
                logger.debug(f"  {lines[0]:<{max_len}}  [cyan]{spacer_id}[/cyan]")
                logger.debug(f"  {lines[1]:<{max_len}}  [dim](alignment)[/dim]")
                logger.debug(f"  {lines[2]:<{max_len}}  [cyan]{location_str}[/cyan]")

                # Determine which distance was used for validation
                used_distance = (
                    edit_distance if distance_metric == "edit" else hamming_distance
                )
                used_metric_name = "Edit" if distance_metric == "edit" else "Hamming"

                # Color code the distances based on validity
                edit_color = (
                    "bold green" if edit_distance <= max_distance else "bold red"
                )
                hamming_color = (
                    "bold green" if hamming_distance <= max_distance else "bold red"
                )
                used_color = (
                    "bold green" if used_distance <= max_distance else "bold red"
                )

                # Always show edit distance
                logger.debug(
                    f"  Edit distance (substitutions + indels): [{edit_color}]{edit_distance}[/{edit_color}]"
                )

                # Show hamming distance, but note if it's not meaningful due to indels
                if has_indels:
                    logger.debug(
                        f"  Hamming distance (substitutions only): [{hamming_color}]{hamming_distance}[/{hamming_color}] [dim](indels present)[/dim]"
                    )
                else:
                    logger.debug(
                        f"  Hamming distance (substitutions only): [{hamming_color}]{hamming_distance}[/{hamming_color}]"
                    )

                logger.debug(f"  Max allowed mismatches: [bold]{max_distance}[/bold]")
                logger.debug(
                    f"  Metric used for validation: [bold]{used_metric_name} distance[/bold] ([{used_color}]{used_distance}[/{used_color}])"
                )

                if has_indels:
                    num_indels = edit_distance - hamming_distance
                    logger.debug(
                        f"  Indels (gap positions): [yellow]{num_indels}[/yellow]"
                    )
                else:
                    logger.debug("  Indels (gap positions): [green]0[/green]")

                if tool_names:
                    logger.debug(
                        f"  Found by tools: [magenta]{', '.join(tool_names)}[/magenta]"
                    )


def display_false_negatives(
    ground_truth: pl.DataFrame,
    alignment_classifications: pl.DataFrame,
    tools_results: pl.DataFrame,
    contigs_file: str,
    spacers_file: str,
    num_examples: int = 5,
    distance_metric: str = "hamming",
    gap_open_penalty: int = 5,
    gap_extend_penalty: int = 5,
    show_all_tools_missed: bool = False,
    specific_tool: Optional[str] = None,
    all_tools: Optional[Dict] = None,
):
    """
    Display example false negative alignments (in ground truth but not found by tools).

    Args:
        ground_truth: Ground truth DataFrame
        alignment_classifications: DataFrame with classified alignments
        tools_results: DataFrame with tool results
        contigs_file: Path to contigs FASTA file
        spacers_file: Path to spacers FASTA file
        num_examples: Number of examples to show per tool
        distance_metric: Distance metric for validation
        gap_open_penalty: Gap open penalty for alignment
        gap_extend_penalty: Gap extend penalty for alignment
        show_all_tools_missed: If True, show FNs missed by ALL tools
        specific_tool: If provided, show FNs for this tool only
        all_tools: Dictionary of all tools (used when showing per-tool FNs)
    """

    if show_all_tools_missed:
        # Show only FNs missed by ALL tools
        # Use ground truth coordinates (start_gt, end_gt) from matches
        found_gt = (
            alignment_classifications.filter(
                pl.col("classification") == "positive_in_plan"
            )
            .select(["spacer_id", "contig_id", "start_gt", "end_gt", "strand"])
            .unique()
        )

        # Rename to match ground truth column names for join
        found_gt = found_gt.rename({"start_gt": "start", "end_gt": "end"})

        false_negatives = ground_truth.join(
            found_gt,
            on=["spacer_id", "contig_id", "start", "end", "strand"],
            how="anti",
        )

        if false_negatives.height == 0:
            logger.debug("\n[green]No false negatives missed by ALL tools![/green]")
            return

        logger.debug(
            f"\n[bold red]FALSE NEGATIVES MISSED BY ALL TOOLS ({false_negatives.height} total, showing {min(num_examples, false_negatives.height)})[/bold red]"
        )
        _display_fn_examples(
            false_negatives.head(num_examples),
            tools_results,
            contigs_file,
            spacers_file,
            distance_metric,
            gap_open_penalty,
            gap_extend_penalty,
        )
        return

    # Show per-tool FNs
    if specific_tool:
        tools_to_check = [specific_tool]
    elif all_tools:
        tools_to_check = sorted(all_tools.keys())
    else:
        tools_to_check = tools_results["tool"].unique().sort().to_list()

    for tool_name in tools_to_check:
        # Get what this tool found (using ground truth coordinates)
        tool_found = (
            tools_results.filter(pl.col("tool") == tool_name)
            .join(
                alignment_classifications.filter(
                    pl.col("classification") == "positive_in_plan"
                ),
                on=["spacer_id", "contig_id", "start", "end", "strand"],
                how="inner",
            )
            .select(["spacer_id", "contig_id", "start_gt", "end_gt", "strand"])
            .unique()
        )

        # Rename to match ground truth column names for join
        tool_found = tool_found.rename({"start_gt": "start", "end_gt": "end"})

        # Find what this tool missed
        false_negatives = ground_truth.join(
            tool_found,
            on=["spacer_id", "contig_id", "start", "end", "strand"],
            how="anti",
        )

        if false_negatives.height == 0:
            logger.debug(f"\n[green]{tool_name.upper()}: No false negatives![/green]")
            continue

        pct_missed = 100 * false_negatives.height / ground_truth.height
        logger.debug(
            f"\n[bold yellow]{tool_name.upper()} FALSE NEGATIVES ({false_negatives.height:,} total = {pct_missed:.2f}%, showing {min(num_examples, false_negatives.height)})[/bold yellow]"
        )

        _display_fn_examples(
            false_negatives.head(num_examples),
            tools_results,
            contigs_file,
            spacers_file,
            distance_metric,
            gap_open_penalty,
            gap_extend_penalty,
            highlight_tool=tool_name,
        )


def _display_fn_examples(
    examples: pl.DataFrame,
    tools_results: pl.DataFrame,
    contigs_file: str,
    spacers_file: str,
    distance_metric: str,
    gap_open_penalty: int,
    gap_extend_penalty: int,
    highlight_tool: Optional[str] = None,
):
    """Helper to display FN examples."""
    for idx, row in enumerate(examples.iter_rows(named=True), 1):
        spacer_id = row["spacer_id"]
        contig_id = row["contig_id"]
        start = row["start"]
        end = row["end"]
        strand = row["strand"]
        expected_mismatches = row.get("mismatches", "unknown")

        # Get sequences
        try:
            spacer_df = get_seq_from_fastx(
                spacers_file,
                [spacer_id],
                return_df=True,
                idcol="spacer_id",
                seqcol="spacer_seq",
            )
            contig_df = get_seq_from_fastx(
                contigs_file,
                [contig_id],
                return_df=True,
                idcol="contig_id",
                seqcol="contig_seq",
            )

            if spacer_df.height == 0 or contig_df.height == 0:
                logger.warning(
                    f"Could not find sequences for {spacer_id} or {contig_id}"
                )
                continue

            spacer_seq = spacer_df["spacer_seq"][0]
            contig_seq = contig_df["contig_seq"][0]
        except Exception as e:
            logger.error(f"Error fetching sequences: {e}")
            continue

        # Calculate actual distances
        edit_distance = test_alignment(
            spacer_seq,
            contig_seq,
            strand=strand,
            start=start,
            end=end,
            distance_metric="edit",
            gap_cost=gap_open_penalty,
            extend_cost=gap_extend_penalty,
        )
        hamming_distance = test_alignment(
            spacer_seq,
            contig_seq,
            strand=strand,
            start=start,
            end=end,
            distance_metric="hamming",
            gap_cost=gap_open_penalty,
            extend_cost=gap_extend_penalty,
        )

        # Generate alignment visualization
        alignment_str = prettify_alignment(
            spacer_seq,
            contig_seq,
            strand=strand,
            start=start,
            end=end,
            gap_cost=gap_open_penalty,
            extend_cost=gap_extend_penalty,
        )

        # Check which tools reported this location (even if wrongly classified)
        tools_that_reported = tools_results.filter(
            (pl.col("spacer_id") == spacer_id)
            & (pl.col("contig_id") == contig_id)
            & (pl.col("start") == start)
            & (pl.col("end") == end)
            & (pl.col("strand") == strand)
        )
        tools_reported_names = (
            tools_that_reported["tool"].unique().to_list()
            if tools_that_reported.height > 0
            else []
        )

        logger.debug(f"\n[bold]False Negative #{idx}:[/bold]")
        logger.debug(f"\n  Spacer:            {spacer_id}")
        logger.debug(f"\n  Contig:            {contig_id}")
        logger.debug(
            f"\n  Position:          {start}-{end} ({'reverse' if strand else 'forward'})"
        )
        logger.debug(f"\n  Expected distance: {expected_mismatches}")
        logger.debug(f"\n  Actual hamming:    {hamming_distance}")
        logger.debug(f"\n  Actual edit:       {edit_distance}")

        if highlight_tool:
            # Highlight which tools DID find this that the highlight_tool missed
            if tools_reported_names:
                logger.debug(
                    f"\n  [green]Found by: {', '.join(tools_reported_names)}[/green]"
                )
            else:
                logger.debug("\n  [red]Not reported by any tool[/red]")
        else:
            if tools_reported_names:
                logger.debug(
                    f"\n  [yellow]Reported by (but not classified as TP): {', '.join(tools_reported_names)}[/yellow]"
                )
            else:
                logger.debug("\n  [red]Not reported by any tool[/red]")

        # Display alignment
        # lines = alignment_str.split("\n )")
        # if len(lines) >= 3:
        #     logger.debug("\n  Alignment:")
        #     for line in lines:
        #         logger.debug(f"    {line}")

        logger.debug(f"\n  Alignment::\n{alignment_str}")


def calculate_all_tool_performance(
    tools_results: pl.DataFrame,
    alignment_classifications: pl.DataFrame,
    ground_truth: pl.DataFrame,
) -> pl.DataFrame:
    """
    Calculate performance metrics for all tools at once using group_by aggregations.

    Always calculates both planned-only and augmented (including non-planned) metrics.

    Args:
        tools_results: DataFrame with all tool results (must have 'tool' column)
        alignment_classifications: DataFrame with classifications for unique alignments
                                  (columns: spacer_id, contig_id, start, end, strand, classification, recalculated_distance)
        ground_truth: Original ground truth DataFrame (needed to identify false negatives)

    Returns:
        DataFrame with performance metrics for all tools (both planned-only and augmented)
    """
    # Join tool results with classifications
    classified_results = tools_results.join(
        alignment_classifications,
        on=["spacer_id", "contig_id", "start", "end", "strand"],
        how="left",
    )

    # Handle cases where some results don't have classifications (shouldn't happen, but defensive)
    # This can occur if alignment_classifications is incomplete
    if classified_results.filter(pl.col("classification").is_null()).height > 0:
        logger.warning(
            f"Found {classified_results.filter(pl.col('classification').is_null()).height} tool results without classifications - skipping them"
        )
        classified_results = classified_results.filter(
            pl.col("classification").is_not_null()
        )

    # If no classified results remain, return empty performance metrics
    if classified_results.height == 0:
        logger.warning(
            "No classified results found - returning empty performance metrics"
        )
        return pl.DataFrame(
            {
                "tool": [],
                "recall_planned": [],
                "precision_planned": [],
                "recall_augmented": [],
                "precision_augmented": [],
                "all_true_positives": [],
                "planned_true_positives": [],
                "positives_not_in_plan": [],
                "invalid_alignments": [],
                "false_negatives_planned": [],
            }
        )

    # Calculate ground truth counts
    planned_gt_count = ground_truth.height

    # For augmented GT, count unique genomic locations (not tool-reported coordinates)
    # For positive_not_in_plan, merge overlapping intervals to avoid counting
    # the same spurious match multiple times when tools report slightly different boundaries
    not_in_plan_df = alignment_classifications.filter(
        pl.col("classification") == "positive_not_in_plan"
    ).select(["spacer_id", "contig_id", "start", "end", "strand"])
    
    # Merge overlapping intervals within tolerance
    merged_not_in_plan = merge_overlapping_intervals(not_in_plan_df, tolerance=3)
    unique_positives_not_in_plan = merged_not_in_plan.height

    augmented_gt_count = planned_gt_count + unique_positives_not_in_plan

    logger.debug(f"Ground truth (planned only): {planned_gt_count}")
    logger.debug(
        f"Verified non-planned alignments (before merging overlaps): {not_in_plan_df.height}"
    )
    logger.debug(
        f"Verified non-planned alignments (after merging overlaps within 3bp): {unique_positives_not_in_plan}"
    )
    logger.debug(f"Ground truth (augmented): {augmented_gt_count}")

    # Aggregate by tool to count each classification type
    # Important: We need to count unique genomic locations, not duplicate tool reports
    # For positive_in_plan, use GROUND TRUTH coordinates (start_gt, end_gt) to avoid counting
    # the same GT location multiple times when a tool reports it with slightly different coordinates
    
    # For positives_not_in_plan, we need to merge overlapping intervals per tool
    # to avoid counting the same spurious match multiple times when reported with slight boundary differences
    tool_not_in_plan_counts = []
    for tool_name in classified_results["tool"].unique().to_list():
        tool_not_in_plan = (
            classified_results.filter(
                (pl.col("tool") == tool_name)
                & (pl.col("classification") == "positive_not_in_plan")
            )
            .select(["spacer_id", "contig_id", "start", "end", "strand"])
        )
        merged = merge_overlapping_intervals(tool_not_in_plan, tolerance=3)
        tool_not_in_plan_counts.append({
            "tool": tool_name,
            "positives_not_in_plan": merged.height,
        })
    
    not_in_plan_df = pl.DataFrame(tool_not_in_plan_counts)
    
    performance = classified_results.group_by("tool").agg(
        [
            # Count unique GROUND TRUTH locations for positive_in_plan (use start_gt, end_gt)
            pl.struct(
                [
                    "spacer_id",
                    "contig_id",
                    "start_gt",
                    "end_gt",
                    "strand",
                    "classification",
                ]
            )
            .filter(pl.col("classification") == "positive_in_plan")
            .n_unique()
            .alias("planned_true_positives"),
            # For invalid alignments, use tool-reported coordinates (exact count is fine here)
            pl.struct(
                ["spacer_id", "contig_id", "start", "end", "strand", "classification"]
            )
            .filter(pl.col("classification") == "invalid_alignment")
            .n_unique()
            .alias("invalid_alignments"),
        ]
    )
    
    # Join the merged positives_not_in_plan counts
    performance = performance.join(not_in_plan_df, on="tool", how="left")

    # Calculate false negatives per tool
    # FN = ground truth entries not found by the tool (within the distance threshold)
    fn_counts = []
    for tool_name in performance["tool"].to_list():
        tool_alignments = classified_results.filter(pl.col("tool") == tool_name)

        # GT entries found by this tool (positive_in_plan)
        # Use ground truth coordinates (start_gt, end_gt) instead of tool coordinates
        gt_found = (
            tool_alignments.filter(pl.col("classification") == "positive_in_plan")
            .select(["spacer_id", "contig_id", "start_gt", "end_gt", "strand"])
            .unique()
        )

        # Rename to match ground truth column names for join
        gt_found = gt_found.rename({"start_gt": "start", "end_gt": "end"})

        # False negatives = GT entries not in the found set
        fn = ground_truth.join(
            gt_found,
            on=["spacer_id", "contig_id", "start", "end", "strand"],
            how="anti",
        ).height

        fn_counts.append({"tool": tool_name, "false_negatives_planned": fn})

    fn_df = pl.DataFrame(fn_counts)
    performance = performance.join(fn_df, on="tool", how="left")

    # Add ground truth counts
    performance = performance.with_columns(
        [
            pl.lit(planned_gt_count).alias("ground_truth_planned"),
            pl.lit(augmented_gt_count).alias("ground_truth_augmented"),
        ]
    )

    # Calculate all true positives (planned + non-planned)
    performance = performance.with_columns(
        [
            (pl.col("planned_true_positives") + pl.col("positives_not_in_plan")).alias(
                "all_true_positives"
            ),
        ]
    )

    # Calculate metrics for PLANNED ONLY mode (strict)
    performance = performance.with_columns(
        [
            # FP (planned only) = verified non-planned + invalid alignments
            (pl.col("positives_not_in_plan") + pl.col("invalid_alignments")).alias(
                "false_positives_planned"
            ),
            # Precision (planned only) = planned TP / (planned TP + all FP)
            (
                pl.col("planned_true_positives")
                / (
                    pl.col("planned_true_positives")
                    + pl.col("positives_not_in_plan")
                    + pl.col("invalid_alignments")
                )
            )
            .fill_null(0.0)
            .alias("precision_planned"),
            # Recall (planned only) = planned TP / (planned TP + FN)
            (
                pl.col("planned_true_positives")
                / (pl.col("planned_true_positives") + pl.col("false_negatives_planned"))
            )
            .fill_null(0.0)
            .alias("recall_planned"),
        ]
    )

    # Calculate metrics for AUGMENTED mode (including non-planned)
    performance = performance.with_columns(
        [
            # FP (augmented) = only invalid alignments (non-planned positives are now TPs)
            pl.col("invalid_alignments").alias("false_positives_augmented"),
            # FN (augmented) = augmented GT - all TPs found by this tool
            (pl.col("ground_truth_augmented") - pl.col("all_true_positives")).alias(
                "false_negatives_augmented"
            ),
            # Precision (augmented) = all TP / (all TP + invalid only)
            (
                pl.col("all_true_positives")
                / (pl.col("all_true_positives") + pl.col("invalid_alignments"))
            )
            .fill_null(0.0)
            .alias("precision_augmented"),
            # Recall (augmented) = all TP / augmented GT (shared across all tools)
            (pl.col("all_true_positives") / pl.col("ground_truth_augmented"))
            .fill_null(0.0)
            .alias("recall_augmented"),
        ]
    )

    # Add percentage columns for easier interpretation (planned mode)
    performance = performance.with_columns(
        [
            (pl.col("recall_planned") * 100).alias("identified_pct_planned"),
            ((1 - pl.col("recall_planned")) * 100).alias("missed_pct_planned"),
            (
                (
                    pl.col("false_positives_planned")
                    / (
                        pl.col("planned_true_positives")
                        + pl.col("false_positives_planned")
                    )
                )
                * 100
            )
            .fill_null(0.0)
            .alias("false_positive_pct_planned"),
        ]
    )

    return performance


def compare_all_tools(
    tools: Dict,
    ground_truth: pl.DataFrame,
    tools_results: pl.DataFrame,
    contigs_file: Optional[str] = None,
    spacers_file: Optional[str] = None,
    verify_false_positives: bool = True,
    estimate_spurious: bool = False,
    hyperfine_results: Optional[pl.DataFrame] = None,
    max_mismatches: int = 5,
    distance_metric: str = "hamming",
    gap_open_penalty: int = 5,
    gap_extend_penalty: int = 5,
) -> pl.DataFrame:
    """
    Compare all tool results against ground truth with consolidated validation.

    This function:
    1. Validates unique alignments once across all tools
    2. Calculates per-tool performance metrics using efficient group_by aggregations
    3. Always computes both planned-only and augmented metrics

    This approach is much more efficient than validating each tool separately.

    Args:
        tools: Dictionary of tool configurations
        ground_truth: Ground truth DataFrame
        tools_results: Combined results from all tools
        contigs_file: Path to contigs FASTA file (for FP verification and spurious estimation)
        spacers_file: Path to spacers FASTA file (for FP verification and spurious estimation)
        verify_false_positives: Whether to verify false positives by checking actual alignments
        estimate_spurious: Whether to estimate expected spurious alignments
        hyperfine_results: Optional DataFrame with hyperfine benchmark results
        max_mismatches: Maximum number of mismatches to consider valid
        distance_metric: Distance metric for validation: 'hamming' (substitutions only) or 'edit' (substitutions + indels)
        gap_open_penalty: Gap open penalty for alignment validation (default: 5)
        gap_extend_penalty: Gap extension penalty for alignment validation (default: 5)

    Returns:
        DataFrame with comparison results for all tools (includes both planned and augmented metrics)
    """
    logger.debug("Starting consolidated tool comparison...")

    # Step 1: Classify unique alignments across all tools
    alignment_classifications = classify_unique_alignments_across_tools(
        ground_truth=ground_truth,
        all_tool_results=tools_results,
        contigs_file=contigs_file,
        spacers_file=spacers_file,
        max_distance=max_mismatches,
        distance_metric=distance_metric,
        gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty,
        coordinate_tolerance=5,  # Allow 5bp tolerance for coordinate matching
    )

    # Step 1.5: Create synthetic "all_tools_combined" entry
    # This combines all unique alignments from all tools into a single tool-independent view
    logger.debug("Creating tool-independent combined results...")
    
    # Get unique alignments by coordinates (drop the tool column, deduplicate, then add new tool name)
    all_tools_combined = (
        tools_results
        .select([col for col in tools_results.columns if col != "tool"])
        .unique(subset=["spacer_id", "contig_id", "start", "end", "strand"])
        .with_columns(pl.lit("all_tools_combined").alias("tool"))
    )
    
    # Add the combined results to tools_results for unified processing
    tools_results_with_combined = vstack_easy(tools_results, all_tools_combined)
    logger.debug(f"Added 'all_tools_combined' with {all_tools_combined.height} unique alignments")

    # Step 2: Display example alignments for each classification
    if contigs_file and spacers_file:
        try:
            display_example_alignments(
                alignment_classifications=alignment_classifications,
                tools_results=tools_results_with_combined,
                contigs_file=contigs_file,
                spacers_file=spacers_file,
                max_distance=max_mismatches,
                num_examples=1,
                distance_metric=distance_metric,
                gap_open_penalty=gap_open_penalty,
                gap_extend_penalty=gap_extend_penalty,
            )
        except Exception as e:
            logger.warning(f"Failed to display example alignments: {e}")

        # Display false negatives missed by all tools
        try:
            # First, show how many FNs each tool has
            logger.debug("\n[bold cyan]FALSE NEGATIVE COUNTS BY TOOL:[/bold cyan]")
            for tool_name in sorted(tools.keys()):
                tool_alignments = tools_results_with_combined.filter(pl.col("tool") == tool_name)
                # Use ground truth coordinates (start_gt, end_gt) from matches
                tool_found_gt = (
                    tool_alignments.join(
                        alignment_classifications.filter(
                            pl.col("classification") == "positive_in_plan"
                        ),
                        on=["spacer_id", "contig_id", "start", "end", "strand"],
                        how="inner",
                    )
                    .select(["spacer_id", "contig_id", "start_gt", "end_gt", "strand"])
                    .unique()
                )

                # Rename to match ground truth column names for join
                tool_found_gt = tool_found_gt.rename(
                    {"start_gt": "start", "end_gt": "end"}
                )

                tool_fn_count = ground_truth.join(
                    tool_found_gt,
                    on=["spacer_id", "contig_id", "start", "end", "strand"],
                    how="anti",
                ).height

                if tool_fn_count == 0:
                    logger.debug(f"  [green]{tool_name}: {tool_fn_count} FNs[/green]")
                else:
                    logger.debug(
                        f"  [yellow]{tool_name}: {tool_fn_count:,} FNs ({100 * tool_fn_count / ground_truth.height:.2f}%)[/yellow]"
                    )

            # Show FNs missed by ALL tools
            display_false_negatives(
                ground_truth=ground_truth,
                alignment_classifications=alignment_classifications,
                tools_results=tools_results_with_combined,
                contigs_file=contigs_file,
                spacers_file=spacers_file,
                num_examples=1,
                distance_metric=distance_metric,
                gap_open_penalty=gap_open_penalty,
                gap_extend_penalty=gap_extend_penalty,
                show_all_tools_missed=True,
            )

            # Show per-tool FN examples (2 examples per tool)
            display_false_negatives(
                ground_truth=ground_truth,
                alignment_classifications=alignment_classifications,
                tools_results=tools_results_with_combined,
                contigs_file=contigs_file,
                spacers_file=spacers_file,
                num_examples=1,
                distance_metric=distance_metric,
                gap_open_penalty=gap_open_penalty,
                gap_extend_penalty=gap_extend_penalty,
                show_all_tools_missed=False,
                all_tools=tools,
            )

            # Show per-tool FN examples (3 examples per tool)
            display_false_negatives(
                ground_truth=ground_truth,
                alignment_classifications=alignment_classifications,
                tools_results=tools_results_with_combined,
                contigs_file=contigs_file,
                spacers_file=spacers_file,
                num_examples=1,
                distance_metric=distance_metric,
                gap_open_penalty=gap_open_penalty,
                gap_extend_penalty=gap_extend_penalty,
                show_all_tools_missed=False,
                all_tools=tools,
            )
        except Exception as e:
            logger.warning(f"Failed to display false negatives: {e}")

    # Step 3: Calculate performance for all tools at once using group_by
    logger.info(
        "Calculating performance for all tools (both planned-only and augmented metrics)..."
    )
    try:
        results_df = calculate_all_tool_performance(
            tools_results=tools_results_with_combined,
            alignment_classifications=alignment_classifications,
            ground_truth=ground_truth,
        )

        # # Add expected spurious estimates if available (both methods)
        # if expected_spurious_hamming is not None:
        #     results_df = results_df.with_columns([
        #         pl.lit(expected_spurious_hamming).alias("expected_spurious_hamming"),
        #     ])

        logger.debug(
            f"Successfully calculated performance for {results_df.height} tools"
        )

    except Exception as e:
        logger.error(f"Failed to calculate performance metrics: {e}")
        import traceback

        logger.error(traceback.format_exc())
        return pl.DataFrame()

    # Step 4: Merge with hyperfine results if available
    if hyperfine_results is not None and hyperfine_results.height > 0:
        hyperfine_df = hyperfine_results.select(
            [pl.col("tool"), pl.col("mean_time").alias("avg_runtime_seconds")]
        )
        results_df = results_df.join(hyperfine_df, on="tool", how="left")

    logger.debug("Comparison completed successfully")
    return results_df


def run_compare_results(
    input_dir,
    max_mismatches=5,
    output_file=None,
    threads=4,
    skip_tools="",
    only_tools=None,
    contigs=None,
    spacers=None,
    verify_false_positives=True,
    distance_metric="hamming",
    gap_open_penalty=5,
    gap_extend_penalty=5,
    logfile=None,
    skip_hyperfine=False,
):
    """
    This function reads alignment tool outputs, compares them against ground truth,
    and calculates performance metrics including precision, recall...

    Always calculates both planned-only and augmented (including non-planned) metrics.

    Args:
        input_dir: Directory containing tool outputs and ground truth
        max_mismatches: Maximum number of mismatches to consider valid
        output_file: Output file for comparison results (None = stdout)
        threads: Number of threads for processing
        skip_tools: Comma-separated list of tools to skip
        only_tools: Comma-separated list of tools to process (overrides skip_tools)
        contigs: Path to custom contigs file (optional)
        spacers: Path to custom spacers file (optional)
        verify_false_positives: Whether to verify false positives by alignment checking (default: True)
        distance_metric: Distance metric for validation: 'hamming' (substitutions only) or 'edit' (substitutions + indels). Default: 'hamming'
        gap_open_penalty: Gap open penalty for alignment validation (default: 5)
        gap_extend_penalty: Gap extension penalty for alignment validation (default: 5)
        logfile: Path to log file for DEBUG output (handled by CLI, included here for documentation)
        skip_hyperfine: If True, skip reading hyperfine results (default: False)

    Returns:
        Tuple of (tools_results, hyperfine_results, performance_results)
    """
    logger.debug(f"Processing tool results from {input_dir}")

    # Ensure directory exists
    if not os.path.exists(input_dir):
        logger.error(f"Input directory {input_dir} does not exist")
        raise FileNotFoundError(f"Input directory {input_dir} does not exist")

    # Determine file paths
    contigs_path = (
        contigs if contigs else f"{input_dir}/simulated_data/simulated_contigs.fa"
    )
    spacers_path = (
        spacers if spacers else f"{input_dir}/simulated_data/simulated_spacers.fa"
    )

    logger.debug(f"Contigs: {contigs_path}")
    logger.debug(f"Spacers: {spacers_path}")

    # Load tool configurations using the new clean function
    logger.debug("Loading tool configurations...")
    tools = load_tool_configs(
        results_dir=input_dir,
        threads=threads,
        contigs_file=contigs,
        spacers_file=spacers,
    )
    logger.debug(f"Loaded {len(tools)} tool configurations")

    # Filter tools based on skip_tools
    if skip_tools:
        logger.debug(f"Skipping tools: {skip_tools}")
        skip_list = skip_tools.split(",")
        tools = {k: v for k, v in tools.items() if k not in skip_list}
        logger.debug(f"Remaining tools after skip: {len(tools)}")

    # Filter tools based on only_tools
    if only_tools:
        logger.debug(f"Only processing tools: {only_tools}")
        only_list = only_tools.split(",")
        tools = {k: v for k, v in tools.items() if k in only_list}
        logger.debug(f"Tools to process: {len(tools)}")

    # Create spacer length dataframe
    logger.debug("Reading spacer sequences...")
    spacers_dict = read_fasta(spacers_path)
    spacer_lendf = pl.DataFrame(
        {
            "spacer_id": spacers_dict.keys(),
            "length": [len(seq) for seq in spacers_dict.values()],
        }
    )
    logger.debug(f"Loaded {len(spacers_dict)} spacers")

    # Read tool results
    logger.debug("Reading tool alignment results...")
    tools_results = read_results(
        tools,
        max_mismatches=max_mismatches,
        spacer_lendf=spacer_lendf,
        ref_file=contigs_path,
    )
    logger.info(f"Read {tools_results.height} total alignment results")

    # Write tools results
    tools_output = f"{input_dir}/tools_results.tsv"
    tools_results.write_csv(tools_output, separator="\t")
    logger.debug(f"Wrote tool results to {tools_output}")

    if not skip_hyperfine:
        # Read hyperfine benchmarking results
        logger.debug("Reading hyperfine benchmark results...")
        hyperfine_results = read_hyperfine_results(tools, input_dir)
        hyperfine_output = f"{input_dir}/hyperfine_results.tsv"
        hyperfine_results.write_csv(hyperfine_output, separator="\t")
        logger.debug(f"Wrote hyperfine results to {hyperfine_output}")

    # Check for ground truth and run performance comparison
    ground_truth_file = f"{input_dir}/simulated_data/planned_ground_truth.tsv"
    if not os.path.exists(ground_truth_file):
        # Try alternate name
        ground_truth_file = f"{input_dir}/simulated_data/ground_truth.tsv"

    # Load ground truth or create empty DataFrame if not found
    if os.path.exists(ground_truth_file):
        logger.debug(f"Reading ground truth from {ground_truth_file}...")
        ground_truth = pl.read_csv(ground_truth_file, separator="\t")
        logger.info(f"Loaded {ground_truth.height} ground truth annotations")

        # Only filter if ground truth is not empty (avoid type inference issues with empty files)
        if ground_truth.height > 0:
            ground_truth = ground_truth.filter(pl.col("mismatches") <= max_mismatches)
            logger.info(
                f"{ground_truth.height} ground truth annotations within max mismatches of {max_mismatches}"
            )

            # Check for duplicates in ground truth
            unique_gt = ground_truth.unique(
                subset=["spacer_id", "contig_id", "start", "end", "strand"]
            )
            if unique_gt.height != ground_truth.height:
                logger.warning(
                    f"Ground truth contains {ground_truth.height - unique_gt.height} duplicate entries!"
                )
                logger.debug(
                    f"Using deduplicated ground truth: {unique_gt.height} unique entries"
                )
                ground_truth = unique_gt
        else:
            # Empty ground truth file - ensure correct schema to avoid type issues downstream
            logger.info("Ground truth file is empty (0 annotations)")
            ground_truth = pl.DataFrame(
                {
                    "spacer_id": [],
                    "contig_id": [],
                    "start": [],
                    "end": [],
                    "strand": [],
                    "mismatches": [],
                },
                schema={
                    "spacer_id": pl.Utf8,
                    "contig_id": pl.Utf8,
                    "start": pl.Int64,
                    "end": pl.Int64,
                    "strand": pl.Boolean,
                    "mismatches": pl.Int64,
                },
            )
    else:
        logger.warning("No ground truth file found, using empty ground truth")
        logger.debug(f"Expected file at: {ground_truth_file}")
        # Create empty ground truth DataFrame with expected schema
        ground_truth = pl.DataFrame(
            {
                "spacer_id": [],
                "contig_id": [],
                "start": [],
                "end": [],
                "strand": [],
                "mismatches": [],
            },
            schema={
                "spacer_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "start": pl.Int64,
                "end": pl.Int64,
                "strand": pl.Boolean,
                "mismatches": pl.Int64,
            },
        )
        logger.info("Created empty ground truth (0 annotations)")

    logger.debug(
        "Comparing results against ground truth using consolidated validation..."
    )
    logger.debug("Computing both planned-only and augmented metrics")

    performance_results = compare_all_tools(
        tools=tools,
        ground_truth=ground_truth,
        tools_results=tools_results,
        contigs_file=contigs_path,
        spacers_file=spacers_path,
        verify_false_positives=verify_false_positives,
        estimate_spurious=False,  # estimate_spurious,
        hyperfine_results=None if skip_hyperfine else hyperfine_results,
        max_mismatches=max_mismatches,
        distance_metric=distance_metric,
        gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty,
    )

    # Output performance results
    if output_file:
        performance_results.write_csv(output_file, separator="\t")
        logger.debug(f"Wrote performance results to {output_file}")
    # else:
    #     logger.debug("\n=== Performance Results ===")
    #     print(performance_results)

    # Also save to standard location with distance metric and threshold in filename
    perf_output = f"{input_dir}/performance_results_{distance_metric}_mm{max_mismatches}.tsv"
    performance_results.write_csv(perf_output, separator="\t")
    logger.info(f"Saved performance results to {perf_output}")

    # Print a sorted, slimmer summary table to screen

    logger.info("PERFORMANCE SUMMARY (sorted by recall_planned)")

    # Handle empty results (e.g., when ground truth is empty and no classifications succeeded)
    if performance_results.height == 0:
        logger.warning(
            "No performance results to display - all tool results may have been unclassified or ground truth was empty"
        )
        logger.info("Results processing completed (no valid comparisons)")
        return (
            tools_results,
            None if skip_hyperfine else hyperfine_results,
            performance_results,
        )

    # Select key columns and sort by recall_planned
    summary_cols = [
        "tool",
        "recall_planned",
        "precision_planned",
        "recall_augmented",
        "precision_augmented",
        "all_true_positives",
        "planned_true_positives",
        "positives_not_in_plan",
        "invalid_alignments",
        "false_negatives_planned",
    ]
    if "expected_spurious_hamming" in performance_results.columns:
        summary_cols.extend(["expected_spurious_hamming"])
    if "avg_runtime_seconds" in performance_results.columns:
        summary_cols.append("avg_runtime_seconds")

    summary_table = performance_results.select(summary_cols).sort(
        "recall_planned", descending=True
    )

    # Format and print
    pl.Config.set_tbl_cols(-1)
    pl.Config.set_tbl_rows(15)
    logger.info("")  # Add newline
    print(summary_table)  # Use regular print to avoid Rich padding/whitespace

    logger.info("Results processing completed successfully")

    return (
        tools_results,
        None if skip_hyperfine else hyperfine_results,
        performance_results,
    )
