"""
Results comparison and validation command (compares the tool outpus against the ground truth).
"""
import os
import logging

import polars as pl

from bench.utils.functions import *
from bench.utils.tool_commands import load_tool_configs

logger = logging.getLogger(__name__)


def run_compare_results(input_dir, max_mismatches=5, output_file=None, threads=4,
                       skip_tools='', only_tools=None, contigs=None, spacers=None):
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
        
        logger.info("Comparing results against ground truth using interval-based validation...")
        performance_results = compare_results(
            tools, ground_truth, tools_results,
            contigs_file=contigs_path,
            spacers_file=spacers_path,
            verify_false_positives=True,
            estimate_spurious=True
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
    else:
        logger.warning("No ground truth file found, skipping performance comparison")
        logger.warning(f"Expected file at: {ground_truth_file}")
    
    logger.info("Results processing completed successfully")
    
    return tools_results, hyperfine_results, performance_results


# Old argparse-based main() removed - use Click CLI via spacer_bencher command
