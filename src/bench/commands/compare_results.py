"""
Results comparison and validation command (compares the tool outpus against the ground truth).
"""
import os
import logging
from typing import Dict, Optional, Union

import polars as pl
import polars_bio as pb
from rich.console import Console

from bench.utils.functions import (
    read_fasta, read_results, read_hyperfine_results,
    validate_intervals_with_polars_bio, test_alignment,
    populate_pldf_withseqs_needletail, reverse_complement,
    estimate_expected_spurious_alignments_fast, prettify_alignment,
    get_seq_from_fastx
)
from bench.utils.tool_commands import load_tool_configs

# Logger is configured by cli.py with RichHandler
logger = logging.getLogger(__name__)
console = Console()


def validate_unique_alignments_across_tools(
    ground_truth: pl.DataFrame,
    all_tool_results: pl.DataFrame,
    contigs_file: str,
    spacers_file: str,
    max_mismatches: int = 5,
    verify_false_positives: bool = True,
    distance_metric: str = 'hamming',
    gap_open_penalty: int = 5,
    gap_extend_penalty: int = 5
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
        max_mismatches: Maximum mismatches to consider valid
        verify_false_positives: Whether to verify FPs by checking actual alignments
        distance_metric: Distance metric for validation: 'hamming' (substitutions only) or 'edit' (substitutions + indels)
        gap_open_penalty: Gap open penalty for alignment validation (default: 5)
        gap_extend_penalty: Gap extension penalty for alignment validation (default: 5)
    
    Returns:
        Tuple of (classified_alignments, exact_matches_with_gt):
        - classified_alignments: DataFrame with alignment_idx and classification
        - exact_matches_with_gt: DataFrame with GT info for recall calculation (from polars-bio join)
    """
    logger.info("Validating unique alignments across all tools...")
    
    # Get unique alignments across all tools (by coordinates)
    unique_alignments = all_tool_results.select([
        "spacer_id", "contig_id", "start", "end", "strand", "mismatches"
    ]).unique()
    
    logger.info(f"Total unique alignments to validate: {unique_alignments.height}")
    
    # Add index for tracking
    unique_alignments = unique_alignments.with_row_index("alignment_idx")
    
    # Use interval-based validation with polars-bio
    try:
        validation_results = validate_intervals_with_polars_bio(
            ground_truth, unique_alignments, max_mismatches
        )
    except Exception as e:
        logger.error(f"Error in interval validation: {e}")
        logger.warning("Falling back to simple classification")
        # Fallback: mark everything as unclassified
        unique_alignments = unique_alignments.with_columns(
            pl.lit("unclassified").alias("classification")
        )
        return unique_alignments
    
    # Process validation results to categorize matches
    exact_matches = validation_results.filter(pl.col("overlap_type") == "exact_match")
    partial_overlaps = validation_results.filter(pl.col("overlap_type") == "partial_overlap")
    
    logger.info(f"Exact matches from validation: {exact_matches.height}")
    logger.info(f"Partial overlaps: {partial_overlaps.height} (will be verified)")
    
    # Debug: Check if validation created duplicate matches
    # Each unique_alignment should match at most one ground truth entry
    index_col_temp = "alignment_idx_2" if "alignment_idx_2" in exact_matches.columns else "alignment_idx"
    if exact_matches.height > 0:
        unique_exact_alignments = exact_matches[index_col_temp].n_unique()
        if unique_exact_alignments != exact_matches.height:
            logger.warning(f"[yellow]Validation created duplicate exact matches:[/yellow]")
            logger.warning(f"  {exact_matches.height} exact match rows, but only {unique_exact_alignments} unique alignments")
            logger.warning(f"  This means some alignments matched multiple GT entries (shouldn't happen!)")
    
    # The polars-bio overlap adds "_2" suffix to the second dataframe (unique_alignments) columns
    index_col = "alignment_idx_2" if "alignment_idx_2" in validation_results.columns else "alignment_idx"
    
    # Only exact matches are true positives - partial overlaps need verification
    # Keep it simple: just track alignment_idx and classification
    tp_classified = pl.DataFrame({"alignment_idx": [], "classification": []})
    
    if exact_matches.height > 0:
        tp_classified = (
            exact_matches.select(index_col)
            .rename({index_col: "alignment_idx"})
            .unique()  # One classification per alignment
            .with_columns(pl.lit("exact_match").alias("classification"))
        )
    
    # Partial overlaps and non-overlaps both need verification as potential FPs
    # Combine them into the fps set to be verified
    fps_indices = []
    
    if partial_overlaps.height > 0:
        partial_indices = partial_overlaps.select(index_col).rename({index_col: "alignment_idx"})
        fps_indices.append(partial_indices)
    
    # Also get alignments with no overlap at all
    if tp_classified.height > 0:
        no_overlap = unique_alignments.filter(
            ~pl.col("alignment_idx").is_in(tp_classified["alignment_idx"])
        )
    else:
        no_overlap = unique_alignments
    
    if no_overlap.height > 0:
        fps_indices.append(no_overlap.select("alignment_idx"))
    
    # Combine all FPs to verify
    if fps_indices:
        fps = unique_alignments.filter(
            pl.col("alignment_idx").is_in(
                pl.concat(fps_indices, how="vertical_relaxed")["alignment_idx"]
            )
        )
    else:
        fps = pl.DataFrame(schema=unique_alignments.schema)
    
    logger.info(f"False positives to categorize: {fps.height}")
    
    # Categorize false positives
    fp_classifications = []
    
    if fps.height > 0 and verify_false_positives:
        logger.info("Loading sequences for FP verification...")
        try:
            # Load sequences for false positives
            fps_with_seqs = populate_pldf_withseqs_needletail(
                fps,
                seqfile=contigs_file,
                idcol="contig_id",
                seqcol="contig_seq",
                trim_to_region=True,
                reverse_by_strand_col=True,
            )
            fps_with_seqs = populate_pldf_withseqs_needletail(
                fps_with_seqs,
                seqfile=spacers_file,
                idcol="spacer_id",
                seqcol="spacer_seq",
                trim_to_region=False,
                reverse_by_strand_col=False,
            )
            
            logger.info("Verifying false positive alignments...")
            logger.info(f"  Algorithm: [cyan]Needleman-Wunsch (parasail)[/cyan]")
            logger.info(f"  Gap open penalty: [bold]{gap_open_penalty}[/bold]")
            logger.info(f"  Gap extend penalty: [bold]{gap_extend_penalty}[/bold]")
            logger.info(f"  Distance metric: [bold]{distance_metric}[/bold] ({'substitutions + indels' if distance_metric == 'edit' else 'substitutions only'})")
            logger.info(f"  Max allowed mismatches: [bold]{max_mismatches}[/bold]")
            logger.info(f"  Alignments to verify: [bold]{fps_with_seqs.height}[/bold]")
            
            # Test each FP alignment
            gaps_as_mismatches = (distance_metric == 'edit')
            for row in fps_with_seqs.iter_rows(named=True):
                spacer_seq = row.get("spacer_seq")
                contig_seq = row.get("contig_seq")
                alignment_idx = row["alignment_idx"]
                
                if spacer_seq and contig_seq:
                    # Test alignment quality using provided gap parameters
                    actual_mismatches = test_alignment(
                        spacer_seq, contig_seq,
                        strand=False,  # Already reverse complemented if needed
                        start=None, end=None,
                        gaps_as_mismatch=gaps_as_mismatches,
                        gap_cost=gap_open_penalty,
                        extend_cost=gap_extend_penalty,
                    )
                    
                    if actual_mismatches <= max_mismatches:
                        # Valid alignment, just not in the plan
                        classification = "positive_not_in_plan"
                    else:
                        # Poor alignment quality
                        classification = "true_false_positive"
                else:
                    # Missing sequences
                    classification = "true_false_positive"
                
                fp_classifications.append({
                    "alignment_idx": alignment_idx,
                    "classification": classification
                })
        
        except Exception as e:
            logger.warning(f"Could not verify false positives: {e}")
            # Without verification, mark all as positives_not_in_plan
            fp_classifications = [
                {"alignment_idx": row["alignment_idx"], "classification": "positive_not_in_plan"}
                for row in fps.iter_rows(named=True)
            ]
    elif fps.height > 0:
        # Without verification, mark all as positives_not_in_plan
        fp_classifications = [
            {"alignment_idx": row["alignment_idx"], "classification": "positive_not_in_plan"}
            for row in fps.iter_rows(named=True)
        ]
    
    # Create FP classification dataframe
    if fp_classifications:
        fp_classified = pl.DataFrame(fp_classifications)
    else:
        fp_classified = pl.DataFrame({"alignment_idx": [], "classification": []})
    
    # Combine all classifications (simple schema: alignment_idx + classification)
    all_classified = pl.concat([tp_classified, fp_classified], how="vertical_relaxed") if tp_classified.height > 0 or fp_classified.height > 0 else pl.DataFrame({
        "alignment_idx": [], 
        "classification": []
    })
    
    # Join back to unique alignments
    result = unique_alignments.join(
        all_classified,
        on="alignment_idx",
        how="left"
    )
    
    # Fill any missing classifications (shouldn't happen, but just in case)
    result = result.with_columns(
        pl.col("classification").fill_null("unclassified")
    )
    
    # Log classification summary
    classification_counts = result.group_by("classification").agg(pl.count()).sort("count", descending=True)
    logger.info("Classification summary:")
    for row in classification_counts.iter_rows(named=True):
        logger.info(f"  {row['classification']}: {row['count']}")
    
    # Add interpretation note
    true_fps = result.filter(pl.col("classification") == "true_false_positive").height
    pos_not_in_plan = result.filter(pl.col("classification") == "positive_not_in_plan").height
    
    if true_fps == 0 and pos_not_in_plan > 0:
        logger.info(f"  [green]✓ All {pos_not_in_plan} non-planned alignments passed validation[/green]")
    elif true_fps > 0:
        logger.info(f"  [yellow]⚠ {true_fps} alignments failed validation (spurious matches)[/yellow]")
    
    # Return both classifications and exact matches (with GT info for recall calculation)
    return result, exact_matches


def display_example_alignments(
    alignment_classifications: pl.DataFrame,
    tools_results: pl.DataFrame,
    contigs_file: str,
    spacers_file: str,
    max_mismatches: int = 5,
    num_examples: int = 3,
    distance_metric: str = 'hamming',
    gap_open_penalty: int = 5,
    gap_extend_penalty: int = 5
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
    logger.info("[bold]" + "="*80 + "[/bold]")
    logger.info("[bold cyan]EXAMPLE ALIGNMENTS BY CLASSIFICATION[/bold cyan]")
    
    # Get examples for each classification type (except exact_match)
    classifications_to_show = ["positive_not_in_plan", "true_false_positive"]
    
    for classification in classifications_to_show:
        examples = alignment_classifications.filter(
            pl.col("classification") == classification
        ).head(num_examples)
        
        total_count = alignment_classifications.filter(
            pl.col("classification") == classification
        ).height
        
        if examples.height == 0:
            # Still log that this classification exists but has no examples
            classification_name = classification.upper().replace('_', ' ')
            logger.info(f"\n[dim]{classification_name} (0 found)[/dim]")
            continue
        
        logger.info(f"\n[bold yellow]{classification.upper().replace('_', ' ')} ({total_count} total, showing {min(num_examples, examples.height)}):[/bold yellow]")
        logger.info("[dim]" + "-" * 80 + "[/dim]")
        
        for idx, row in enumerate(examples.iter_rows(named=True), 1):
            spacer_id = row["spacer_id"]
            contig_id = row["contig_id"]
            start = row["start"]
            end = row["end"]
            strand = row["strand"]
            
            # Get sequences (only fetch what we need!)
            try:
                spacer_df = get_seq_from_fastx(spacers_file, [spacer_id], return_df=True, 
                                              idcol="spacer_id", seqcol="spacer_seq")
                contig_df = get_seq_from_fastx(contigs_file, [contig_id], return_df=True,
                                              idcol="contig_id", seqcol="contig_seq")
                
                if spacer_df.height == 0 or contig_df.height == 0:
                    logger.warning(f"Could not find sequences for {spacer_id} or {contig_id}")
                    continue
                    
                spacer_seq = spacer_df["spacer_seq"][0]
                contig_seq = contig_df["contig_seq"][0]
            except Exception as e:
                logger.error(f"Error fetching sequences: {e}")
                continue
            
            # Find which tools reported this alignment
            tools_with_this = tools_results.filter(
                (pl.col("spacer_id") == spacer_id) &
                (pl.col("contig_id") == contig_id) &
                (pl.col("start") == start) &
                (pl.col("end") == end) &
                (pl.col("strand") == strand)
            )
            tool_names = tools_with_this["tool"].unique().sort().to_list() if "tool" in tools_with_this.columns else []
            
            # Calculate edit distance (with indels/gaps counted)
            edit_distance = test_alignment(
                spacer_seq, contig_seq,
                strand=strand, start=start, end=end,
                gaps_as_mismatch=True,
                gap_cost=gap_open_penalty,
                extend_cost=gap_extend_penalty
            )
            
            # Calculate hamming distance (substitutions only, no indels)
            hamming_distance = test_alignment(
                spacer_seq, contig_seq,
                strand=strand, start=start, end=end,
                gaps_as_mismatch=False,
                gap_cost=gap_open_penalty,
                extend_cost=gap_extend_penalty
            )
            
            # Determine if there are indels (gaps)
            has_indels = (edit_distance != hamming_distance)
            
            # Generate alignment visualization
            alignment_str = prettify_alignment(
                spacer_seq,
                contig_seq,
                strand=strand,
                start=start,
                end=end,
                gap_cost=gap_open_penalty,
                extend_cost=gap_extend_penalty
            )
            
            # Add labels to the alignment
            strand_str = "(-)" if strand else "(+)"
            location_str = f"{contig_id}:{start}-{end} {strand_str}"
            
            lines = alignment_str.split("\n")
            if len(lines) >= 3:
                # Add sequence IDs on the right side
                max_len = max(len(line) for line in lines)
                
                logger.info(f"[bold]Example {idx}:[/bold]")
                logger.info(f"  {lines[0]:<{max_len}}  [cyan]{spacer_id}[/cyan]")
                logger.info(f"  {lines[1]:<{max_len}}  [dim](alignment)[/dim]")
                logger.info(f"  {lines[2]:<{max_len}}  [cyan]{location_str}[/cyan]")
                
                # Determine which distance was used for validation
                used_distance = edit_distance if distance_metric == 'edit' else hamming_distance
                used_metric_name = "Edit" if distance_metric == 'edit' else "Hamming"
                
                # Color code the distances based on validity
                edit_color = "bold green" if edit_distance <= max_mismatches else "bold red"
                hamming_color = "bold green" if hamming_distance <= max_mismatches else "bold red"
                used_color = "bold green" if used_distance <= max_mismatches else "bold red"
                
                # Always show edit distance
                logger.info(f"  Edit distance (substitutions + indels): [{edit_color}]{edit_distance}[/{edit_color}]")
                
                # Show hamming distance, but note if it's not meaningful due to indels
                if has_indels:
                    logger.info(f"  Hamming distance (substitutions only): [{hamming_color}]{hamming_distance}[/{hamming_color}] [dim](indels present)[/dim]")
                else:
                    logger.info(f"  Hamming distance (substitutions only): [{hamming_color}]{hamming_distance}[/{hamming_color}]")
                
                logger.info(f"  Max allowed mismatches: [bold]{max_mismatches}[/bold]")
                logger.info(f"  Metric used for validation: [bold]{used_metric_name} distance[/bold] ([{used_color}]{used_distance}[/{used_color}])")
                
                if has_indels:
                    num_indels = edit_distance - hamming_distance
                    logger.info(f"  Indels (gap positions): [yellow]{num_indels}[/yellow]")
                else:
                    logger.info(f"  Indels (gap positions): [green]0[/green]")
                
                if tool_names:
                    logger.info(f"  Found by tools: [magenta]{', '.join(tool_names)}[/magenta]")
            
    logger.info("[bold]" + "="*80 + "[/bold]")


def calculate_all_tool_performance(
    tools_results: pl.DataFrame,
    alignment_classifications: pl.DataFrame,
    exact_matches_with_gt: pl.DataFrame,
    ground_truth: pl.DataFrame,
    augment_ground_truth: bool = False
) -> pl.DataFrame:
    """
    Calculate performance metrics for all tools at once using group_by aggregations.
    
    This is more efficient and elegant than looping through tools individually.
    
    Args:
        tools_results: Combined results from all tools (must have 'tool' column)
        alignment_classifications: Pre-computed classifications for all unique alignments
        exact_matches_with_gt: Validation results with GT info (from polars-bio join)
        ground_truth: Original ground truth DataFrame
        augment_ground_truth: If True, count verified "positives_not_in_plan" as TPs.
                             Ground truth becomes: planned + unique verified non-planned alignments
    
    Returns:
        DataFrame with performance metrics for all tools
    """
    # Deduplicate tools_results before joining to avoid double-counting
    # Some tools may report the same alignment with different mismatch counts
    original_count = tools_results.height
    unique_tool_results = tools_results.unique(
        subset=["tool", "spacer_id", "contig_id", "start", "end", "strand"]
    )
    deduplicated_count = unique_tool_results.height
    
    if original_count != deduplicated_count:
        logger.info(f"Deduplicated tool results: {original_count} → {deduplicated_count} (removed {original_count - deduplicated_count} duplicates)")
    
    # Join deduplicated tool results with classifications
    classified_results = unique_tool_results.join(
        alignment_classifications,
        on=["spacer_id", "contig_id", "start", "end", "strand"],
        how="left"
    )
    
    # Determine ground truth count
    if augment_ground_truth:
        # Augmented GT = planned + unique verified non-planned alignments (shared across all tools)
        unique_positives_not_in_plan = alignment_classifications.filter(
            pl.col("classification") == "positive_not_in_plan"
        ).height
        gt_count = ground_truth.height + unique_positives_not_in_plan
        logger.info(f"Augmented ground truth: {ground_truth.height} (planned) + {unique_positives_not_in_plan} (verified non-planned) = {gt_count}")
    else:
        gt_count = ground_truth.height
        logger.info(f"Using standard ground truth count: {gt_count}")
    
    # Aggregate by tool to count each classification type
    performance = classified_results.group_by("tool").agg([
        # Count each classification type
        pl.col("classification").filter(pl.col("classification") == "exact_match").count().alias("planned_true_positives"),
        pl.col("classification").filter(pl.col("classification") == "positive_not_in_plan").count().alias("positives_not_in_plan"),
        pl.col("classification").filter(pl.col("classification") == "true_false_positive").count().alias("true_false_positives"),
    ])
    
    # For recall: count unique GT entries found per tool (not just total alignments)
    # Join tool results with exact match classifications, then with GT info
    if exact_matches_with_gt.height > 0:
        # Get alignment_idx for exact matches
        exact_match_alignments = classified_results.filter(
            pl.col("classification") == "exact_match"
        )
        
        # Join with exact_matches_with_gt using alignment_idx
        index_col = "alignment_idx_2" if "alignment_idx_2" in exact_matches_with_gt.columns else "alignment_idx"
        
        # First join tool results with classifications to get alignment_idx
        tool_exact_matches = unique_tool_results.join(
            alignment_classifications.filter(pl.col("classification") == "exact_match"),
            on=["spacer_id", "contig_id", "start", "end", "strand"],
            how="inner"
        )
        
        # Then join with GT info using alignment_idx
        exact_tool_alignments_with_gt = tool_exact_matches.join(
            exact_matches_with_gt.select([index_col, "spacer_id_1", "contig_id_1", "start_1", "end_1", "strand_1"]),
            left_on="alignment_idx",
            right_on=index_col,
            how="inner"
        )
        
        # Count unique GT entries per tool
        gt_counts = exact_tool_alignments_with_gt.group_by("tool").agg([
            pl.struct(["spacer_id_1", "contig_id_1", "start_1", "end_1", "strand_1"]).n_unique().alias("unique_gt_entries_found")
        ])
        
        performance = performance.join(gt_counts, on="tool", how="left")
        # Fill nulls (tools with no exact matches)
        performance = performance.with_columns(
            pl.col("unique_gt_entries_found").fill_null(0)
        )
    else:
        # No exact matches found
        performance = performance.with_columns(
            pl.lit(0).alias("unique_gt_entries_found")
        )
    
    # Calculate derived metrics using polars expressions
    if augment_ground_truth:
        # Augmented mode: positives_not_in_plan also count as TPs
        performance = performance.with_columns([
            # True positives = planned TPs + verified non-planned positives
            (pl.col("planned_true_positives") + pl.col("positives_not_in_plan")).alias("true_positives"),
            # False positives = only true_false_positives (spurious alignments)
            pl.col("true_false_positives").alias("false_positives"),
            # Ground truth is the augmented count
            pl.lit(gt_count).alias("ground_truth_total"),
            pl.lit(True).alias("augmented_ground_truth"),
        ])
    else:
        # Standard mode: only planned TPs count as true positives
        performance = performance.with_columns([
            # True positives = only alignments matching the plan
            pl.col("planned_true_positives").alias("true_positives"),
            # False positives = verified non-planned + spurious alignments
            (pl.col("positives_not_in_plan") + pl.col("true_false_positives")).alias("false_positives"),
            # Ground truth is the original count
            pl.lit(gt_count).alias("ground_truth_total"),
            pl.lit(False).alias("augmented_ground_truth"),
        ])
    
    # Calculate precision, recall, and F1 score
    if augment_ground_truth:
        # Augmented mode: recall based on total unique valid alignments (planned + verified non-planned)
        # We already deduplicated tool results, so true_positives represents unique valid alignments found
        performance = performance.with_columns([
            # Precision = TP / (TP + FP)
            (pl.col("true_positives") / (pl.col("true_positives") + pl.col("false_positives")))
                .fill_null(0.0)
                .alias("precision"),
            # Recall = unique valid alignments found / augmented GT, capped at 1.0
            # Cap at 1.0 because a tool can report overlapping positions that map to the same biological targets
            pl.min_horizontal(
                pl.col("true_positives") / pl.col("ground_truth_total"),
                pl.lit(1.0)
            )
                .fill_null(0.0)
                .alias("recall"),
        ])
        
        # Log if any tool exceeded augmented GT (recall capped at 1.0)
        tools_over_gt = performance.filter(pl.col("true_positives") > pl.col("ground_truth_total"))
        if tools_over_gt.height > 0:
            for row in tools_over_gt.iter_rows(named=True):
                logger.info(f"  [yellow]{row['tool']}:[/yellow] {row['true_positives']} TPs > {row['ground_truth_total']} augmented GT (recall capped at 1.0)")
                logger.info(f"    This indicates the tool reported overlapping/redundant positions for some targets")
    else:
        # Standard mode: recall based only on planned GT entries found
        performance = performance.with_columns([
            # Precision = TP / (TP + FP)
            (pl.col("true_positives") / (pl.col("true_positives") + pl.col("false_positives")))
                .fill_null(0.0)
                .alias("precision"),
            # Recall = unique planned GT entries found / planned GT total
            (pl.col("unique_gt_entries_found") / pl.col("ground_truth_total"))
                .fill_null(0.0)
                .alias("recall"),
        ])
    
    # F1 score = 2 * (precision * recall) / (precision + recall)
    performance = performance.with_columns([
        (2 * pl.col("precision") * pl.col("recall") / (pl.col("precision") + pl.col("recall")))
            .fill_null(0.0)
            .alias("f1_score"),
    ])
    
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
    augment_ground_truth: bool = False,
    distance_metric: str = 'hamming',
    gap_open_penalty: int = 5,
    gap_extend_penalty: int = 5,
) -> pl.DataFrame:
    """
    Compare all tool results against ground truth with consolidated validation.
    
    This function:
    1. Validates unique alignments once across all tools
    2. Calculates per-tool performance metrics using efficient group_by aggregations
    
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
        augment_ground_truth: If True, count verified "positives_not_in_plan" as TPs
        distance_metric: Distance metric for validation: 'hamming' (substitutions only) or 'edit' (substitutions + indels)
        gap_open_penalty: Gap open penalty for alignment validation (default: 5)
        gap_extend_penalty: Gap extension penalty for alignment validation (default: 5)
    
    Returns:
        DataFrame with comparison results for all tools
    """
    logger.info("Starting consolidated tool comparison...")
    
    # Step 1: Validate unique alignments across all tools
    alignment_classifications, exact_matches_with_gt = validate_unique_alignments_across_tools(
        ground_truth=ground_truth,
        all_tool_results=tools_results,
        contigs_file=contigs_file,
        spacers_file=spacers_file,
        max_mismatches=max_mismatches,
        verify_false_positives=verify_false_positives,
        distance_metric=distance_metric,
        gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty
    )
    
    # Step 1.5: Display example alignments for each classification
    if contigs_file and spacers_file:
        try:
            display_example_alignments(
                alignment_classifications=alignment_classifications,
                tools_results=tools_results,
                contigs_file=contigs_file,
                spacers_file=spacers_file,
                max_mismatches=max_mismatches,
                num_examples=3,
                distance_metric=distance_metric,
                gap_open_penalty=gap_open_penalty,
                gap_extend_penalty=gap_extend_penalty
            )
        except Exception as e:
            logger.warning(f"Failed to display example alignments: {e}")
    
    # Step 2: Estimate expected spurious alignments if requested
    expected_spurious_per_tool = None
    if estimate_spurious and contigs_file and spacers_file:
        try:
            logger.info("Estimating expected spurious alignments...")
            contigs_dict = read_fasta(contigs_file)
            spacers_dict = read_fasta(spacers_file)
            
            spurious_est = estimate_expected_spurious_alignments_fast(
                contigs_dict, spacers_dict,
                max_mismatches=max_mismatches
            )
            expected_spurious_per_tool = spurious_est['mean_spurious']
            logger.info(f"Expected spurious: {expected_spurious_per_tool:.2f} ± {spurious_est['std_spurious']:.2f}")
        except Exception as e:
            logger.warning(f"Failed to estimate spurious alignments: {e}")
    
    # Step 3: Calculate performance for all tools at once using group_by
    logger.info("Calculating performance for all tools using aggregations...")
    try:
        results_df = calculate_all_tool_performance(
            tools_results=tools_results,
            alignment_classifications=alignment_classifications,
            exact_matches_with_gt=exact_matches_with_gt,
            ground_truth=ground_truth,
            augment_ground_truth=augment_ground_truth
        )
        
        # Add expected spurious if available
        if expected_spurious_per_tool is not None:
            results_df = results_df.with_columns(
                pl.lit(expected_spurious_per_tool).alias("expected_spurious")
            )
        
        logger.info(f"Successfully calculated performance for {results_df.height} tools")
        
    except Exception as e:
        logger.error(f"Failed to calculate performance metrics: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return pl.DataFrame()
    
    # Step 4: Merge with hyperfine results if available
    if hyperfine_results is not None and hyperfine_results.height > 0:
        hyperfine_df = hyperfine_results.select([
            pl.col("tool"),
            pl.col("mean_time").alias("avg_runtime_seconds")
        ])
        results_df = results_df.join(hyperfine_df, on="tool", how="left")
    
    logger.info("Comparison completed successfully")
    return results_df


def run_compare_results(input_dir, max_mismatches=5, output_file=None, threads=4,
                       skip_tools='', only_tools=None, contigs=None, spacers=None,
                       augment_ground_truth=False, verify_false_positives=True,
                       estimate_spurious=True, distance_metric='hamming',
                       gap_open_penalty=5, gap_extend_penalty=5):
    """
    This function reads alignment tool outputs, compares them against ground truth,
    and calculates performance metrics including precision, recall, and F1 scores.
    
    Args:
        input_dir: Directory containing tool outputs and ground truth
        max_mismatches: Maximum number of mismatches to consider valid
        output_file: Output file for comparison results (None = stdout)
        threads: Number of threads for processing
        skip_tools: Comma-separated list of tools to skip
        only_tools: Comma-separated list of tools to process (overrides skip_tools)
        contigs: Path to custom contigs file (optional)
        spacers: Path to custom spacers file (optional)
        augment_ground_truth: If True, verified "positives_not_in_plan" count as TPs (default: False)
                             The augmented ground truth = original plan + unique verified non-planned alignments
                             This shared ground truth is used for all tools to ensure recall ≤ 1.0
        verify_false_positives: Whether to verify false positives by alignment checking (default: True)
        estimate_spurious: Whether to estimate expected spurious alignments (default: True)
        distance_metric: Distance metric for validation: 'hamming' (substitutions only) or 'edit' (substitutions + indels). Default: 'hamming'
        gap_open_penalty: Gap open penalty for alignment validation (default: 5)
        gap_extend_penalty: Gap extension penalty for alignment validation (default: 5)
    
    Returns:
        Tuple of (tools_results, hyperfine_results, performance_results)
    """
    logger.info(f"Processing tool results from {input_dir}")
    
    # Ensure directory exists
    if not os.path.exists(input_dir):
        logger.error(f"Input directory {input_dir} does not exist")
        raise FileNotFoundError(f"Input directory {input_dir} does not exist")
    
    # Determine file paths
    contigs_path = contigs if contigs else f"{input_dir}/simulated_data/simulated_contigs.fa"
    spacers_path = spacers if spacers else f"{input_dir}/simulated_data/simulated_spacers.fa"
    
    logger.debug(f"Contigs: {contigs_path}")
    logger.debug(f"Spacers: {spacers_path}")
    
    
    # Load tool configurations using the new clean function
    logger.info("Loading tool configurations...")
    tools = load_tool_configs(
        results_dir=input_dir,
        threads=threads,
        contigs_file=contigs,
        spacers_file=spacers
    )
    logger.info(f"Loaded {len(tools)} tool configurations")
    
    # Filter tools based on skip_tools
    if skip_tools:
        logger.info(f"Skipping tools: {skip_tools}")
        skip_list = skip_tools.split(",")
        tools = {k: v for k, v in tools.items() if k not in skip_list}
        logger.info(f"Remaining tools after skip: {len(tools)}")
    
    # Filter tools based on only_tools
    if only_tools:
        logger.info(f"Only processing tools: {only_tools}")
        only_list = only_tools.split(",")
        tools = {k: v for k, v in tools.items() if k in only_list}
        logger.info(f"Tools to process: {len(tools)}")
    
    # Create spacer length dataframe
    logger.info("Reading spacer sequences...")
    spacers_dict = read_fasta(spacers_path)
    spacer_lendf = pl.DataFrame({
        "spacer_id": spacers_dict.keys(),
        "length": [len(seq) for seq in spacers_dict.values()]
    })
    logger.debug(f"Loaded {len(spacers_dict)} spacers")
    
    # Read tool results
    logger.info("Reading tool alignment results...")
    tools_results = read_results(
        tools,
        max_mismatches=max_mismatches,
        spacer_lendf=spacer_lendf,
        ref_file=contigs_path,
    )
    logger.info(f"Read {len(tools_results)} total alignment results")
    
    # Write tools results
    tools_output = f"{input_dir}/tools_results.tsv"
    tools_results.write_csv(tools_output, separator="\t")
    logger.info(f"Wrote tool results to {tools_output}")
    
    # Read hyperfine benchmarking results
    logger.info("Reading hyperfine benchmark results...")
    hyperfine_results = read_hyperfine_results(tools, input_dir)
    hyperfine_output = f"{input_dir}/hyperfine_results.tsv"
    hyperfine_results.write_csv(hyperfine_output, separator="\t")
    logger.info(f"Wrote hyperfine results to {hyperfine_output}")
    
    # Check for ground truth and run performance comparison
    ground_truth_file = f"{input_dir}/simulated_data/planned_ground_truth.tsv"
    if not os.path.exists(ground_truth_file):
        # Try alternate name
        ground_truth_file = f"{input_dir}/simulated_data/ground_truth.tsv"
    
    performance_results = None
    if os.path.exists(ground_truth_file):
        logger.info(f"Reading ground truth from {ground_truth_file}...")
        ground_truth = pl.read_csv(ground_truth_file, separator="\t")
        logger.info(f"Loaded {len(ground_truth)} ground truth annotations")
        
        # Check for duplicates in ground truth
        unique_gt = ground_truth.unique(subset=["spacer_id", "contig_id", "start", "end", "strand"])
        if unique_gt.height != ground_truth.height:
            logger.warning(f"Ground truth contains {ground_truth.height - unique_gt.height} duplicate entries!")
            logger.info(f"Using deduplicated ground truth: {unique_gt.height} unique entries")
            ground_truth = unique_gt
        
        logger.info("Comparing results against ground truth using consolidated validation...")
        if augment_ground_truth:
            logger.info("  augment_ground_truth=True: Verified non-planned alignments will count as TPs")
        
        performance_results = compare_all_tools(
            tools=tools,
            ground_truth=ground_truth,
            tools_results=tools_results,
            contigs_file=contigs_path,
            spacers_file=spacers_path,
            verify_false_positives=verify_false_positives,
            estimate_spurious=estimate_spurious,
            hyperfine_results=hyperfine_results,
            max_mismatches=max_mismatches,
            augment_ground_truth=augment_ground_truth,
            distance_metric=distance_metric,
            gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty
        )
        
        # Output performance results
        if output_file:
            performance_results.write_csv(output_file, separator="\t")
            logger.info(f"Wrote performance results to {output_file}")
        else:
            logger.info("\n=== Performance Results ===")
            print(performance_results)
        
        # Also save to standard location
        perf_output = f"{input_dir}/performance_results.tsv"
        performance_results.write_csv(perf_output, separator="\t")
        logger.info(f"Saved performance results to {perf_output}")
        
        # Print a sorted, slimmer summary table to screen
        print("\n" + "="*80)
        logger.info("PERFORMANCE SUMMARY (sorted by recall score)")
        if augment_ground_truth:
            logger.info("  Note: augment_ground_truth=True, verified non-planned alignments count as TPs")
        
        # Select key columns and sort by recall score
        summary_cols = ["tool",  "recall",  "true_false_positives","positives_not_in_plan","false_positives","true_positives"]
        if "expected_spurious" in performance_results.columns:
            summary_cols.append("expected_spurious")
        if "augmented_ground_truth" in performance_results.columns:
            summary_cols.append("augmented_ground_truth")
        if "avg_runtime_seconds" in performance_results.columns:
            summary_cols.append("avg_runtime_seconds")
        
        summary_table = performance_results.select(summary_cols).sort("recall", descending=True)
        
        # Format and print
        pl.Config.set_tbl_cols(-1)
        pl.Config.set_tbl_rows(15)
        logger.info("")  # Add newline
        console.print(summary_table)  # Use console for polars table formatting
        logger.info("="*80)
    else:
        logger.warning("No ground truth file found, skipping performance comparison")
        logger.warning(f"Expected file at: {ground_truth_file}")
    
    logger.info("Results processing completed successfully")
    
    return tools_results, hyperfine_results, performance_results


