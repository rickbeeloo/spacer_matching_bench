"""
Script generation command.
Handles the generation of bash execution scripts for tools based on their JSON configuration files.
"""
import os
import logging

from bench.utils.functions import *
from bench.utils.tool_commands import load_tool_configs

logger = logging.getLogger(__name__)


def run_generate_scripts(input_dir=None, output_dir=None, threads=4, max_runs=1, warmups=0, 
                        skip_tools='', only_tools=None, contigs=None, spacers=None, max_mismatches=5):
    """
    Core script generation function.
    
    Args:
        input_dir: Directory containing simulated data (optional if contigs/spacers provided)
        output_dir: Directory for generated scripts (defaults to input_dir if not specified)
        threads: Number of threads for tool execution
        max_runs: Maximum number of benchmark runs (hyperfine)
        warmups: Number of warmup runs (hyperfine)
        skip_tools: Comma-separated list of tools to skip
        only_tools: Comma-separated list of tools to run (overrides skip_tools)
        contigs: Path to custom contigs file (optional)
        spacers: Path to custom spacers file (optional)
        max_mismatches: Maximum number of mismatches/distance for tools that support it (default: 5)
    """
    # Determine output directory
    if output_dir is None and input_dir is None:
        raise ValueError("Must provide either --input-dir (-i) or --output-dir (-o)")
    
    # Use input_dir as output_dir if output_dir not specified
    output_dir = output_dir if output_dir is not None else input_dir
    
    logger.debug(f"Generating scripts in {output_dir}")
    
    # Create necessary subdirectories
    os.makedirs(f"{output_dir}/raw_outputs", exist_ok=True)
    os.makedirs(f"{output_dir}/bash_scripts", exist_ok=True)
    logger.debug("Created output directories")
    
    # Load tool configurations using the new clean function
    logger.debug("Loading tool configurations...")
    tools = load_tool_configs(
        results_dir=output_dir,
        threads=threads,
        contigs_file=contigs,
        spacers_file=spacers,
        max_mismatches=max_mismatches
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
        logger.debug(f"Only running tools: {only_tools}")
        only_list = only_tools.split(",")
        tools = {k: v for k, v in tools.items() if k in only_list}
        logger.debug(f"Tools to run: {len(tools)}")
    
    # Generate scripts for each tool
    for tool_name, tool in tools.items():
        logger.debug(f"Generating script for {tool['name']}...")
        try:
            create_bash_script(tool, input_dir, max_runs=max_runs, warmups=warmups)
            logger.debug(f"  ✓ Created script for {tool_name}")
        except Exception as e:
            logger.error(f"  ✗ Failed to create script for {tool_name}: {e}", exc_info=True)
            raise
    
    logger.debug(f"Successfully generated {len(tools)} tool scripts in {input_dir}/bash_scripts/")
