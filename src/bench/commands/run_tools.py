import argparse
import os

import polars as pl

from bench.utils.functions import *
from bench.utils.tool_commands import populate_tools

pl.Config(tbl_cols=13, tbl_rows=8)


def add_common_args(parser):
    """Add common arguments shared across all commands"""
    parser.add_argument(
        "-cl",
        "--contig_length_range",
        nargs=2,
        type=int,
        default=[1000, 32000],
        help="Range of contig lengths",
    )
    parser.add_argument(
        "-nc",
        "--sample_size_contigs",
        type=int,
        default=1500,
        help="Number of contigs to generate",
    )
    parser.add_argument(
        "-ns",
        "--sample_size_spacers",
        type=int,
        default=50,
        help="Number of spacers to generate",
    )
    parser.add_argument(
        "-ls",
        "--spacer_length_range",
        nargs=2,
        type=int,
        default=[20, 60],
        help="Range of spacer lengths",
    )
    parser.add_argument(
        "-lm",
        "--n_mismatch_range",
        nargs=2,
        type=int,
        default=[0, 5],
        help="Range of number of mismatches",
    )
    parser.add_argument(
        "-ir",
        "--insertion_range",
        nargs=2,
        type=int,
        default=[1, 5],
        help="Range of number of insertions per contig",
    )
    parser.add_argument(
        "-nir",
        "--n_insertion_range",
        nargs=2,
        type=int,
        default=[0, 2],
        help="Range of number of insertions within spacer sequences",
    )
    parser.add_argument(
        "-ndr",
        "--n_deletion_range",
        nargs=2,
        type=int,
        default=[0, 1],
        help="Range of number of deletions within spacer sequences",
    )
    parser.add_argument(
        "-prc",
        "--prop_rc",
        type=float,
        default=0.5,
        help="Proportion of spacers to reverse complement",
    )
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number threads")
    parser.add_argument(
        "-c",
        "--contigs",
        type=str,
        default=None,
        help="Path to contigs file (if not provided, will generate simulated data)",
    )
    parser.add_argument(
        "-s",
        "--spacers",
        type=str,
        default=None,
        help="Path to spacers file (if not provided, will generate simulated data)",
    )
    parser.add_argument(
        "-id",
        "--id_prefix",
        type=str,
        default=None,
        help="Prefix for sequence IDs (default is hash of parameters)",
    )


def add_tool_args(parser):
    """Add tool-specific arguments"""
    parser.add_argument(
        "-mm",
        "--max_mismatches",
        type=int,
        default=5,
        help="Maximum number of mismatches to allow",
    )
    parser.add_argument(
        "-mr",
        "--max_runs",
        type=int,
        default=5,
        help="Maximum number of runs (hyperfine)",
    )
    parser.add_argument(
        "-w", "--warmups", type=int, default=0, help="Number of warmups (hyperfine)"
    )
    parser.add_argument(
        "-st",
        "--skip_tools",
        type=str,
        default="vsearch",
        help="Comma-separated list of tools to skip",
    )


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

    # Run tools and get results
    print("Running tools...")
    run_tools(tools, input_dir, max_runs=args.max_runs, warmups=args.warmups)
    print("Finished running tools")

    print("Reading tool results...")
    tools_results = read_results(
        tools,
        max_mismatches=args.max_mismatches,
        spacer_lendf=spacer_lendf,
        ref_file=f"{input_dir}/simulated_data/simulated_contigs.fa",
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
