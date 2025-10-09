import argparse
import os

import polars as pl

from bench.utils.functions import *
from bench.utils.tool_commands import populate_tools
from bench.utils.arg_parsers import add_common_args, add_tool_args

def main():
    """Execute tools on simulated data"""
    parser = argparse.ArgumentParser(
        description="Execute tools on simulated data",
        add_help=True,
        allow_abbrev=True,
    )
    
    add_common_args(parser)
    add_tool_args(parser)
    parser.add_argument(
        "-i", "--input_dir",
        type=str,
        required=True,
        help="Input directory containing simulated data and tool scripts"
    )
    
    args = parser.parse_args()
    
    print("Executing tools...")
    
    input_dir = args.input_dir
    if not os.path.exists(input_dir):
        print(f"Error: Input directory {input_dir} does not exist")
        return
    
    # Create necessary subdirectories, just in case
    os.makedirs(f"{input_dir}/raw_outputs", exist_ok=True)
    os.makedirs(f"{input_dir}/bash_scripts", exist_ok=True)
    
    # Set up args object for populate_tools
    args.results_dir = input_dir
    
    # Determine file paths
    if args.spacers:
        spacers_file = args.spacers
    else:
        spacers_file = f"{input_dir}/simulated_data/simulated_spacers.fa"
    
    if args.contigs:
        contigs_file = args.contigs
    else:
        contigs_file = f"{input_dir}/simulated_data/simulated_contigs.fa"
    
    # Store file paths in args for populate_tools
    args.contigs_file = contigs_file
    args.spacers_file = spacers_file
    

    print("Initializing tools...")
    tools = populate_tools(args)
    print(f"Initialized {len(tools)} tools")

    if args.skip_tools:
        print(f"Disabling user-specified tools: {args.skip_tools}")
        skip_tools = args.skip_tools.split(",")
        tools = {k: v for k, v in tools.items() if k not in skip_tools}

    if args.tools_to_run:
        print(f"Running user-specified tools: {args.tools_to_run}")
        tools_to_run = args.tools_to_run.split(",")
        tools = {k: v for k, v in tools.items() if k in tools_to_run}

    # Run tools and get results
    print("Running tools...")
    run_tools(tools, input_dir, max_runs=args.max_runs, warmups=args.warmups)
    print("Finished running tools")
    # Create spacer length dataframe from existing data
    spacers = read_fasta(spacers_file)
    spacer_lendf = pl.DataFrame(
        {"spacer_id": spacers.keys(), "length": [len(seq) for seq in spacers.values()]}
    )
    
    print("Reading tool results...")
    tools_results = read_results(
        tools,
        max_mismatches=args.max_mismatches,
        spacer_lendf=spacer_lendf,
        ref_file=contigs_file,
    )
    tools_results.write_csv(f"{input_dir}/tools_results.tsv", separator="\t")
    print(f"Wrote tool results to {input_dir}/tools_results.tsv")

    print("Reading hyperfine results...")
    hyperfine_results = read_hyperfine_results(tools, input_dir)
    hyperfine_results.write_csv(
        file=f"{input_dir}/hyperfine_results.tsv", separator="\t"
    )
    print(f"Wrote hyperfine results to {input_dir}/hyperfine_results.tsv")

    # Check for ground truth and run performance comparison
    ground_truth_file = f"{input_dir}/simulated_data/ground_truth.tsv"
    if os.path.exists(ground_truth_file):
        print("Reading ground truth...")
        ground_truth = pl.read_csv(ground_truth_file, separator="\t")
        print("Comparing results against ground truth...")
        performance_results = compare_results(
            tools, ground_truth, tools_results, soft_false_positive_threshold=1
        )
        print(performance_results)
        performance_results.write_csv(
            f"{input_dir}/performance_results.tsv", separator="\t"
        )
        print(f"Wrote performance results to {input_dir}/performance_results.tsv")
    else:
        print("No ground truth found, skipping performance comparison")

    print("Tool execution completed successfully")


if __name__ == "__main__":
    main()
