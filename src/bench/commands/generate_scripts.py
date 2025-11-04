"""
Script generation command.
Handles the generation of bash execution scripts for tools based on their JSON configuration files.
"""
import os
import logging

from bench.utils.functions import *
from bench.utils.tool_commands import load_tool_configs

logger = logging.getLogger(__name__)


def run_generate_scripts(input_dir, threads=4, max_runs=1, warmups=0, 
                        skip_tools='', only_tools=None, contigs=None, spacers=None):
    """
    Core script generation function.
    
    Args:
        input_dir: Directory containing simulated data
        threads: Number of threads for tool execution
        max_runs: Maximum number of benchmark runs (hyperfine)
        warmups: Number of warmup runs (hyperfine)
        skip_tools: Comma-separated list of tools to skip
        only_tools: Comma-separated list of tools to run (overrides skip_tools)
        contigs: Path to custom contigs file (optional)
        spacers: Path to custom spacers file (optional)
    """
    logger.info(f"Generating scripts in {input_dir}")
    
    # Ensure directories exist
    if not os.path.exists(input_dir):
        logger.warning(f"Input directory {input_dir} does not exist, creating it")
        os.makedirs(input_dir, exist_ok=True)
    
    # Create necessary subdirectories
    os.makedirs(f"{input_dir}/raw_outputs", exist_ok=True)
    os.makedirs(f"{input_dir}/bash_scripts", exist_ok=True)
    logger.debug("Created output directories")
    
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
        logger.info(f"Only running tools: {only_tools}")
        only_list = only_tools.split(",")
        tools = {k: v for k, v in tools.items() if k in only_list}
        logger.info(f"Tools to run: {len(tools)}")
    
    # Generate scripts for each tool
    for tool_name, tool in tools.items():
        logger.info(f"Generating script for {tool['name']}...")
        try:
            create_bash_script(tool, input_dir, max_runs=max_runs, warmups=warmups)
            logger.debug(f"  ✓ Created script for {tool_name}")
        except Exception as e:
            logger.error(f"  ✗ Failed to create script for {tool_name}: {e}", exc_info=True)
            raise
    
    logger.info(f"Successfully generated {len(tools)} tool scripts in {input_dir}/bash_scripts/")
