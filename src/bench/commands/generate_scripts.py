import argparse
import os

import polars as pl

from bench.utils.functions import *
from bench.utils.tool_commands import populate_tools
from bench.utils.arg_parsers import add_common_args, add_tool_args


def main():
    """Generate tool execution scripts"""
    parser = argparse.ArgumentParser(
        description="Generate tool execution scripts from configurations",
        add_help=True,
        allow_abbrev=True,
    )
    
    add_common_args(parser)
    add_tool_args(parser)
    parser.add_argument(
        "-i", "--input_dir",
        type=str,
        required=True,
        help="Input directory containing simulated data"
    )
    
    args = parser.parse_args()
    
    print("Generating tool execution scripts...")
    
    input_dir = args.input_dir
    if not os.path.exists(input_dir):
        print(f"Error: Input directory {input_dir} does not exist")
        return
    
    # Create necessary subdirectories
    os.makedirs(f"{input_dir}/raw_outputs", exist_ok=True)
    os.makedirs(f"{input_dir}/bash_scripts", exist_ok=True)
    
    # Set up args object for populate_tools
    args.results_dir = input_dir
    
    # Create spacer length dataframe from existing data
    spacers_file = f"{input_dir}/simulated_data/simulated_spacers.fa"
    if not os.path.exists(spacers_file):
        print(f"Error: Spacers file {spacers_file} not found")
        return
    
    spacers = read_fasta(spacers_file)
    spacer_lendf = pl.DataFrame(
        {"spacer_id": spacers.keys(), "length": [len(seq) for seq in spacers.values()]}
    )
    
    print("Initializing tools...")
    tools = populate_tools(args, spacer_lendf=spacer_lendf)
    print(f"Initialized {len(tools)} tools")

    if args.skip_tools:
        print(f"Disabling user-specified tools: {args.skip_tools}")
        skip_tools = args.skip_tools.split(",")
        tools = {k: v for k, v in tools.items() if k not in skip_tools}

    # Generate scripts for each tool
    for tool in tools.values():
        print(f"Generating script for {tool['name']}...")
        create_bash_script(tool, input_dir, max_runs=args.max_runs, warmups=args.warmups)
    
    print(f"Tool scripts generated successfully in {input_dir}/bash_scripts/")


if __name__ == "__main__":
    main()
