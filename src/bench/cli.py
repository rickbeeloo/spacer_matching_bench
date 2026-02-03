"""
Click entry point for the spacer_bencher tool.
"""

import click
import logging
import sys

# from pathlib import Path
from rich.logging import RichHandler

# Configure global logger
logger = logging.getLogger("spacer_bencher")


def setup_logging(verbose: bool = False, logfile: str = None):
    """
    Configure global logging for the application using Rich for pretty output.

    Args:
        verbose: If True, set logging level to DEBUG, otherwise INFO
        logfile: Path to log file for DEBUG output (always DEBUG level, independent of verbose)
    """
    # Set root logger to DEBUG if either verbose or logfile is specified
    # This ensures DEBUG messages are generated and can be captured
    root_level = logging.DEBUG if (verbose or logfile) else logging.INFO
    
    # Console handler level - INFO unless verbose
    console_level = logging.DEBUG if verbose else logging.INFO

    # Get the root logger
    root_logger = logging.getLogger()

    # Remove any existing handlers to avoid duplicates
    root_logger.handlers.clear()

    # Configure Rich handler for pretty output
    # Show file paths and line numbers only in verbose mode
    # Disable soft wrapping to prevent line breaks
    from rich.console import Console
    console = Console(soft_wrap=False, width=380)
    
    rich_handler = RichHandler(
        console=console,
        rich_tracebacks=True,
        show_path=verbose,  # Show file:// paths in verbose mode
        show_time=True,
        markup=True,  # Enable rich markup in log messages
        omit_repeated_times=False,
    )
    rich_handler.setLevel(console_level)  # Set console handler level explicitly

    # Configure root logger (all loggers inherit from this)
    root_logger.setLevel(root_level)
    root_logger.addHandler(rich_handler)

    # Add file handler if logfile is specified (always DEBUG level)
    if logfile:
        file_handler = logging.FileHandler(logfile, mode="a")
        file_handler.setLevel(logging.DEBUG)  # Always DEBUG for logfile
        file_formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)
        logger.info(f"DEBUG logging enabled to file: {logfile}")


@click.group()
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose (DEBUG) logging")
@click.version_option(version="1.0.0", prog_name="spacer_bencher")
def cli(verbose):
    """
    Spacer Bencher -  benchmarking spacer-protospacer identification tools.

    This tool provides a complete pipeline for:

    \b
    1. Generating simulated CRISPR spacer and contig sequences (spacer_bencher simulate)
    2. Creating execution scripts for various alignment tools (spacer_bencher generate-scripts)
    3. Running alignment tools and collecting results (spacer_bencher execute-tools)
    4. Comparing and validating tool performance (spacer_bencher compare-results)

    Use 'spacer_bencher COMMAND --help' for detailed information on each command.
    """
    setup_logging(verbose)
    logger.debug("CLI initialized in verbose mode")


@cli.command(context_settings=dict(show_default=True))
@click.option(
    "--num-contigs", "-nc", type=int, default=1000, help="Number of contigs to generate"
)
@click.option(
    "--num-spacers", "-ns", type=int, default=200, help="Number of spacers to generate"
)
@click.option(
    "--contig-length",
    "-cl",
    type=(int, int),
    default=(2000, 15000),
    help="Range of contig lengths (min max)",
)
@click.option(
    "--spacer-length",
    "-ls",
    type=(int, int),
    default=(20, 60),
    help="Range of spacer lengths (min max)",
)
@click.option(
    "--mismatch-range",
    "-lm",
    type=(int, int),
    default=(0, 3),
    help="Range of substitution mismatches per spacer (min max)",
)
@click.option(
    "--spacer-insertions",
    "-ir",
    type=(int, int),
    default=(1, 4),
    help="Number of times to insert each spacer into contigs (min max)",
)
@click.option(
    "--indel-insertions",
    "-nir",
    type=(int, int),
    default=(0, 0),
    help="Number of insertion mutations (indels) to add within each spacer (min max)",
)
@click.option(
    "--indel-deletions",
    "-ndr",
    type=(int, int),
    default=(0, 0),
    help="Number of deletion mutations (indels) to add within each spacer (min max)",
)
@click.option(
    "--reverse-complement",
    "-prc",
    type=float,
    default=0.5,
    help="Proportion of spacers to reverse complement (0.0-1.0)",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=4,
    help="Number of threads for parallel processing",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(),
    required=True,
    help="Output directory for simulated data",
)
@click.option(
    "--id-prefix",
    "-id",
    type=str,
    default=None,
    help="Prefix for sequence IDs (default: auto-generated)",
)
@click.option(
    "--gc-content",
    type=float,
    default=None,
    help="GC content percentage (0-100) - DEPRECATED: use --contig-gc-content and --spacer-gc-content",
)
@click.option(
    "--contig-gc-content",
    type=float,
    default=None,
    help="GC content percentage for contigs (0-100)",
)
@click.option(
    "--spacer-gc-content",
    type=float,
    default=None,
    help="GC content percentage for spacers (0-100)",
)
@click.option(
    "--contig-distribution",
    type=click.Choice(["uniform", "normal", "bell"]),
    default="uniform",
    help="Distribution for contig lengths",
)
@click.option(
    "--spacer-distribution",
    type=click.Choice(["uniform", "normal", "bell"]),
    default="uniform",
    help="Distribution for spacer lengths",
)
@click.option("--verify", is_flag=True, help="Verify simulation after generation")
@click.option(
    "--contigs",
    type=click.Path(exists=True),
    default=None,
    help="Path to existing contigs FASTA file (optional, for read-and-insert mode)",
)
@click.option(
    "--spacers",
    type=click.Path(exists=True),
    default=None,
    help="Path to existing spacers FASTA file (optional, for read-and-insert mode)",
)
def simulate(
    num_contigs,
    num_spacers,
    contig_length,
    spacer_length,
    mismatch_range,
    spacer_insertions,
    indel_insertions,
    indel_deletions,
    reverse_complement,
    threads,
    output_dir,
    id_prefix,
    gc_content,
    contig_gc_content,
    spacer_gc_content,
    contig_distribution,
    spacer_distribution,
    verify,
    contigs,
    spacers,
):
    """
    Generate simulated CRISPR spacer and contig sequences.

    This command creates synthetic test data by:

    \b
    1. Generating random contig sequences (or reading from --contigs file)
    2. Generating random spacer sequences (or reading from --spacers file)
    3. Inserting spacers into contigs with optional mutations
    4. Recording ground truth positions and mutations

    \b
    Key Parameters Explained:

    --spacer-insertions (-ir): Controls how many TIMES each spacer is inserted
                                into the contigs. Example: -ir 1 2 means each
                                spacer will be inserted 1 or 2 times.

    --indel-insertions (-nir): Controls how many INSERTION MUTATIONS (adding
                                bases) are applied to the spacer sequence.
                                Example: -nir 0 2 means 0-2 insertion indels.

    --indel-deletions (-ndr): Controls how many DELETION MUTATIONS (removing
                               bases) are applied to the spacer sequence.
                               Example: -ndr 0 1 means 0-1 deletion indels.

    \b
    Read-and-Insert Mode:
      If --contigs and/or --spacers are provided, those sequences will be read from
      files instead of being generated. This is useful for:
      - Testing tools on real CRISPR sequences
      - Creating synthetic insertions in real metagenomic contigs
      - Validating tool performance on known sequences
      - Note: that the number of contigs/spacers parameters will be used to subsample that many sequences from the files, or fail if not enough are present.
      - Note2: the subsampled/non-synthetic files will still be called "simulated_contigs.fa" and "simulated_spacers.fa" in the output, for consistency.

    \b
    Perfect Match Example:
      spacer_bencher simulate -nc 100 -ns 50 -lm 0 0 -ir 1 1 -nir 0 0 -ndr 0 0
      This creates spacers with NO mutations inserted exactly ONCE.

    \b
    Read-and-Insert Example:
      spacer_bencher simulate --contigs real_contigs.fa --spacers real_spacers.fa -ir 2 5 -lm 0 2 -o output/
      This reads real sequences and inserts spacers with 0-2 mismatches, 2-5 times each.

    \b
    Output Files:
      - simulated_contigs.fa: FASTA file with contig sequences
      - simulated_spacers.fa: FASTA file with spacer sequences
      - planned_ground_truth.tsv: Tab-separated ground truth annotations
    """
    from bench.commands.simulate import run_simulate

    logger.info(
        f"Starting simulation with {num_spacers} spacers and {num_contigs} contigs"
    )
    logger.debug(f"Output directory: {output_dir}")

    try:
        run_simulate(
            num_contigs=num_contigs,
            num_spacers=num_spacers,
            contig_length_range=contig_length,
            spacer_length_range=spacer_length,
            mismatch_range=mismatch_range,
            spacer_insertions=spacer_insertions,
            indel_insertions=indel_insertions,
            indel_deletions=indel_deletions,
            reverse_complement=reverse_complement,
            threads=threads,
            output_dir=output_dir,
            id_prefix=id_prefix,
            gc_content=gc_content,
            contig_gc_content=contig_gc_content,
            spacer_gc_content=spacer_gc_content,
            contig_distribution=contig_distribution,
            spacer_distribution=spacer_distribution,
            verify=verify,
            contigs=contigs,
            spacers=spacers,
        )
        logger.info("Simulation completed successfully")
        click.echo(click.style("✓ Simulation completed successfully", fg="green"))
    except Exception as e:
        logger.exception(f"Simulation failed: {e}")
        click.echo(click.style(f"✗ Simulation failed: {e}", fg="red"), err=True)
        sys.exit(1)


@cli.command("generate-scripts", context_settings=dict(show_default=True))
@click.option(
    "--input-dir",
    "-i",
    type=click.Path(exists=True),
    required=False,
    help="Input directory containing simulated data",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(),
    required=False,
    help="Output directory for generated scripts (defaults to input-dir if specified)",
)
@click.option(
    "--contigs",
    type=click.Path(exists=True),
    required=False,
    help="Path to custom contigs file",
)
@click.option(
    "--spacers",
    type=click.Path(exists=True),
    required=False,
    help="Path to custom spacers file",
)
@click.option(
    "--threads", "-t", type=int, default=4, help="Number of threads for tool execution"
)
@click.option(
    "--max-runs",
    "-mr",
    type=int,
    default=1,
    help="Maximum number of benchmark runs (hyperfine)",
)
@click.option(
    "--warmups", "-w", type=int, default=0, help="Number of warmup runs (hyperfine)"
)
@click.option(
    "--skip-tools",
    "-st",
    type=str,
    default="",
    help="Comma-separated list of tools to skip",
)
@click.option(
    "--only-tools",
    "-rt",
    type=str,
    default=None,
    help="Comma-separated list of tools to run (overrides skip-tools)",
)
@click.option(
    "--max-mismatches",
    "-mm",
    type=int,
    default=5,
    help="Maximum number of mismatches/edit-distance parameter for tools that support it (i.e. subs in indelfree.sh, k in sassy, -v in bowtie1)",
)
def generate_scripts(
    input_dir,
    output_dir,
    contigs,
    spacers,
    threads,
    max_runs,
    warmups,
    skip_tools,
    only_tools,
    max_mismatches,
):
    """
    Generate execution scripts for alignment tools.

    This command creates bash scripts that run various alignment tools
    on the simulated (or real) data. Scripts are generated using tool
    configuration files and include proper parameters for benchmarking.

    \b
    Generated Files:
      - bash_scripts/TOOL_NAME.sh: Individual tool execution scripts
      - Each script uses hyperfine for benchmarking if max-runs > 1

    \b
    Example:
      spacer_bencher generate-scripts -i tests/validation -o output_dir -t 4 -mr 1
      spacer_bencher generate-scripts -o output_dir --contigs contigs.fa --spacers spacers.fa -t 4
    """
    from bench.commands.generate_scripts import run_generate_scripts

    logger.info(f"Generating scripts for tools in {input_dir}")

    try:
        run_generate_scripts(
            input_dir=input_dir,
            output_dir=output_dir,
            threads=threads,
            max_runs=max_runs,
            warmups=warmups,
            skip_tools=skip_tools,
            only_tools=only_tools,
            contigs=contigs,
            spacers=spacers,
            max_mismatches=max_mismatches,
        )
        logger.info("Script generation completed successfully")
        click.echo(click.style("✓ Scripts generated successfully", fg="green"))
    except Exception as e:
        logger.exception(f"Script generation failed: {e}")
        click.echo(click.style(f"✗ Script generation failed: {e}", fg="red"), err=True)
        sys.exit(1)


@cli.command("execute-tools", context_settings=dict(show_default=True))
@click.option(
    "--input-dir",
    "-i",
    type=click.Path(exists=True),
    required=True,
    help="Input directory containing scripts and data",
)
@click.option(
    "--skip-tools",
    "-st",
    type=str,
    default="",
    help="Comma-separated list of tools to skip",
)
@click.option(
    "--only-tools",
    "-rt",
    type=str,
    default=None,
    help="Comma-separated list of tools to run",
)
@click.option(
    "--debug",
    is_flag=True,
    help="Run commands directly without hyperfine for better error messages",
)
def execute_tools(input_dir, skip_tools, only_tools, debug):
    """
    Execute alignment tools using generated bash scripts.

    This command runs the bash scripts generated by 'generate-scripts',
    executing each alignment tool and collecting their output files.

    Note: Thread counts and benchmark settings are already baked into the
    generated scripts, so you don't need to specify them here.

    \b
    Output Files:
      - raw_outputs/TOOL_NAME_output.sam: Tool alignment results
      - raw_outputs/TOOL_NAME.json: Benchmarking results (if max-runs > 1)

    \b
    Example:
      spacer_bencher execute-tools -i tests/validation
    """
    from bench.commands.run_tools import run_execute_tools

    logger.info(f"Executing tools from {input_dir}")

    try:
        run_execute_tools(
            input_dir=input_dir,
            skip_tools=skip_tools,
            only_tools=only_tools,
            debug=debug,
        )
        logger.info("Tool execution completed successfully")
        click.echo(click.style("✓ Tools executed successfully", fg="green"))
    except Exception as e:
        logger.exception(f"Tool execution failed: {e}")
        click.echo(click.style(f"✗ Tool execution failed: {e}", fg="red"), err=True)
        sys.exit(1)


@cli.command("full-run", context_settings=dict(show_default=True))
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(),
    required=True,
    help="Output directory for the entire pipeline",
)
@click.option(
    "--num-contigs", "-nc", type=int, default=1000, help="Number of contigs to generate"
)
@click.option(
    "--num-spacers", "-ns", type=int, default=200, help="Number of spacers to generate"
)
@click.option(
    "--contig-length",
    "-cl",
    type=(int, int),
    default=(2000, 5000),
    help="Range of contig lengths (min max)",
)
@click.option(
    "--spacer-length",
    "-ls",
    type=(int, int),
    default=(20, 60),
    help="Range of spacer lengths (min max)",
)
@click.option(
    "--mismatch-range",
    "-lm",
    type=(int, int),
    default=(0, 5),
    help="Range of substitution mismatches per spacer (min max)",
)
@click.option(
    "--spacer-insertions",
    "-ir",
    type=(int, int),
    default=(1, 4),
    help="Number of times to insert each spacer into contigs (min max)",
)
@click.option(
    "--indel-insertions",
    "-nir",
    type=(int, int),
    default=(0, 0),
    help="Number of insertion mutations per spacer (min max)",
)
@click.option(
    "--indel-deletions",
    "-ndr",
    type=(int, int),
    default=(0, 0),
    help="Number of deletion mutations per spacer (min max)",
)
@click.option(
    "--reverse-complement",
    "-prc",
    type=float,
    default=0.5,
    help="Proportion of spacers to reverse complement (0.0-1.0)",
)
@click.option(
    "--threads", "-t", type=int, default=4, help="Number of threads for tool execution"
)
@click.option(
    "--max-runs", "-mr", type=int, default=1, help="Maximum number of benchmark runs"
)
@click.option("--warmups", "-w", type=int, default=0, help="Number of warmup runs")
@click.option(
    "--skip-tools",
    "-st",
    type=str,
    default="",
    help="Comma-separated list of tools to skip",
)
@click.option(
    "--only-tools",
    "-rt",
    type=str,
    default=None,
    help="Comma-separated list of tools to run",
)
@click.option(
    "--max-mismatches",
    "-mm",
    type=int,
    default=5,
    help="Maximum number of mismatches for comparison",
)
@click.option(
    "--augment-ground-truth",
    is_flag=True,
    default=False,
    help="Count verified non-planned alignments as true positives",
)
@click.option(
    "--distance",
    type=click.Choice(["hamming", "edit"]),
    default="hamming",
    help="Distance metric for validation: hamming (substitutions only) or edit (substitutions + indels). ",
)
@click.option(
    "--gap-open-penalty",
    type=int,
    default=5,
    help="Gap open penalty for alignment validation",
)
@click.option(
    "--gap-extend-penalty",
    type=int,
    default=5,
    help="Gap extension penalty for alignment validation",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    help="Enable verbose (DEBUG) logging with file paths and line numbers",
)
def full_run(
    output_dir,
    num_contigs,
    num_spacers,
    contig_length,
    spacer_length,
    mismatch_range,
    spacer_insertions,
    indel_insertions,
    indel_deletions,
    reverse_complement,
    threads,
    max_runs,
    warmups,
    skip_tools,
    only_tools,
    max_mismatches,
    augment_ground_truth,
    distance,
    gap_open_penalty,
    gap_extend_penalty,
    verbose,
):
    """
    Run the complete benchmarking pipeline.

    This command orchestrates all steps:
    1. Simulate data
    2. Generate scripts
    3. Execute tools
    4. Compare results

    \b
    Example:
      spacer_bencher full-run -o results/benchmark_run -nc 100 -ns 50 -t 4
    """
    from bench.commands.simulate import run_simulate
    from bench.commands.generate_scripts import run_generate_scripts
    from bench.commands.run_tools import run_execute_tools
    from bench.commands.compare_results import run_compare_results

    # Update logging if verbose flag is set
    if verbose:
        setup_logging(verbose=True)

    logger.info("STARTING FULL PIPELINE RUN")

    try:
        # Step 1: Simulate data
        logger.info("[1/4] SIMULATING DATA")
        run_simulate(
            num_contigs=num_contigs,
            num_spacers=num_spacers,
            contig_length_range=contig_length,
            spacer_length_range=spacer_length,
            mismatch_range=mismatch_range,
            spacer_insertions=spacer_insertions,
            indel_insertions=indel_insertions,
            indel_deletions=indel_deletions,
            reverse_complement=reverse_complement,
            threads=threads,
            output_dir=output_dir,
        )
        click.echo(click.style("✓ Step 1/4: Data simulation completed", fg="green"))

        # Step 2: Generate scripts
        logger.info("[2/4] GENERATING TOOL SCRIPTS")
        run_generate_scripts(
            input_dir=output_dir,
            threads=threads,
            max_runs=max_runs,
            warmups=warmups,
            skip_tools=skip_tools,
            only_tools=only_tools,
        )
        click.echo(click.style("✓ Step 2/4: Script generation completed", fg="green"))

        # Step 3: Execute tools
        logger.info("[3/4] EXECUTING TOOLS")
        run_execute_tools(
            input_dir=output_dir, skip_tools=skip_tools, only_tools=only_tools
        )
        click.echo(click.style("✓ Step 3/4: Tool execution completed", fg="green"))

        # Step 4: Compare results
        logger.info("[4/4] COMPARING RESULTS")
        run_compare_results(
            input_dir=output_dir,
            max_mismatches=max_mismatches,
            output_file=f"{output_dir}/comparison_results.tsv",
            threads=threads,
            skip_tools=skip_tools,
            only_tools=only_tools,
            augment_ground_truth=augment_ground_truth,
            distance_metric=distance,
            gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty,
        )
        click.echo(click.style("✓ Step 4/4: Results comparison completed", fg="green"))

        # Success summary
        logger.info("PIPELINE COMPLETED SUCCESSFULLY")
        logger.info(f"All results saved to: {output_dir}")
        # Normalize paths to avoid double slashes
        output_dir_clean = output_dir.rstrip("/")
        logger.info(f"  - Simulated data: {output_dir_clean}/simulated_data/")
        logger.info(f"  - Tool outputs: {output_dir_clean}/raw_outputs/")
        logger.info(
            f"  - Comparison results: {output_dir_clean}/comparison_results.tsv"
        )

        click.echo(
            click.style(
                f"\n✓ Full pipeline completed! Results in: {output_dir}",
                fg="green",
                bold=True,
            )
        )

    except Exception as e:
        logger.exception(f"Pipeline failed: {e}")
        click.echo(
            click.style(f"✗ Pipeline failed at current step: {e}", fg="red"), err=True
        )
        sys.exit(1)


@cli.command("compare-results")
@click.option(
    "--input-dir",
    "-i",
    type=click.Path(exists=True),
    required=True,
    help="Input directory containing tool outputs and ground truth",
)
@click.option(
    "--max-mismatches",
    "-mm",
    type=int,
    default=5,
    help="Maximum number of mismatches to consider",
)
@click.option(
    "--output-file",
    "-o",
    type=click.Path(),
    default=None,
    help="Output file for comparison results (default: stdout)",
)
@click.option(
    "--logfile",
    "-l",
    type=click.Path(),
    default=None,
    help="Path to log file for detailed DEBUG output (always DEBUG level, independent of --verbose)",
)
@click.option(
    "--threads", "-t", type=int, default=4, help="Number of threads for processing"
)
@click.option(
    "--verify-false-positives",
    is_flag=True,
    default=True,
    help="Count verified non-planned alignments as true positives",
)
@click.option(
    "--skip-tools",
    "-st",
    type=str,
    default="",
    help="Comma-separated list of tools to skip",
)
@click.option(
    "--only-tools",
    "-rt",
    type=str,
    default=None,
    help="Comma-separated list of tools to run",
)
@click.option(
    "--distance",
    type=click.Choice(["hamming", "edit"]),
    default="hamming",
    help="Distance metric for validation: hamming (substitutions only) or edit (substitutions + indels). Default: hamming",
)
@click.option(
    "--gap-open-penalty",
    type=int,
    default=5,
    help="Gap open penalty for alignment validation",
)
@click.option(
    "--gap-extend-penalty",
    type=int,
    default=5,
    help="Gap extension penalty for alignment validation",
)
@click.option(
    "--verbose", "-v", is_flag=True, help="Enable verbose (DEBUG) logging to console"
)
@click.option(
    "--skip-hyperfine", "-s", is_flag=True, help="Skip hyperfine benchmarking"
)
@click.option(
    "--contigs",
    type=click.Path(exists=True),
    required=False,
    default = None,
    help="Path to custom contigs file",
)
@click.option(
    "--spacers",
    type=click.Path(exists=True),
    required=False,
    default = None,
    help="Path to custom spacers file",
)
def compare_results(
    input_dir,
    max_mismatches,
    output_file,
    logfile,
    threads,
    skip_tools,
    only_tools,
    distance,
    gap_open_penalty,
    gap_extend_penalty,
    verbose,
    skip_hyperfine,
    contigs,
    spacers
):
    """
    Compare and validate alignment tool results.

    This command analyzes the output from alignment tools, compares them
    against ground truth data, and generates performance metrics including
    precision, recall, and F1 scores.

    \b
    Metrics Calculated:
      - True Positives: Correctly identified spacer locations (that were in the simulation plan)
      - False Positives: Incorrectly reported locations
      - False Positives: Incorrectly reported locations
      - False Negatives: Missed spacer locations
      - Precision: TP / (TP + FP)
      - Recall: TP / (TP + FN)
      - F1 Score: Harmonic mean of precision and recall

    \b
    Augmented Ground Truth Mode (--augment-ground-truth):
      When enabled, verified non-planned alignments (valid alignments found by tools
      but not in the original simulation plan) are counted as true positives rather
      than false positives. This gives a more realistic performance assessment when
      sequences contain naturally occurring similar regions.

    \b
    Gap Penalty Options:
      --gap-open-penalty: Penalty for opening a gap
      --gap-extend-penalty: Penalty for extending a gap

    \b
    Example:
      spacer_bencher compare-results -i tests/validation -mm 3 -o results.tsv
      spacer_bencher compare-results -i tests/validation --augment-ground-truth
    """
    from bench.commands.compare_results import run_compare_results

    # Update logging if verbose flag is set, and setup logfile if specified
    if verbose or logfile:
        setup_logging(verbose=verbose, logfile=logfile)

    logger.info(f"Comparing tool results in {input_dir}")

    try:
        run_compare_results(
            input_dir=input_dir,
            max_mismatches=max_mismatches,
            output_file=output_file,
            threads=threads,
            skip_tools=skip_tools,
            only_tools=only_tools,
            contigs=contigs,
            spacers=spacers,
            distance_metric=distance,
            gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty,
            logfile=logfile,
            skip_hyperfine=skip_hyperfine,
        )
    except Exception as e:
        logger.exception(f"Results comparison failed: {e}")
        click.echo(click.style(f"✗ Results comparison failed: {e}", fg="red"), err=True)
        sys.exit(1)


def main():
    """Entry point for the CLI"""
    cli()


if __name__ == "__main__":
    main()
