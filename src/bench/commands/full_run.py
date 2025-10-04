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
    """Run the complete benchmarking pipeline (original functionality)"""
    parser = argparse.ArgumentParser(
        description="Run the complete benchmarking pipeline (simulate + generate scripts + run tools)",
        add_help=True,
        allow_abbrev=True,
    )
    
    add_common_args(parser)
    add_tool_args(parser)
    
    args = parser.parse_args()
    
    print("Running full benchmark pipeline...")
    
    # Determine the results directory path
    results_dir = f"results/run_t_{args.threads}_nc_{args.sample_size_contigs}_ns_{args.sample_size_spacers}_ir_{args.insertion_range[0]}_{args.insertion_range[1]}_lm_{args.n_mismatch_range[0]}_{args.n_mismatch_range[1]}_nir_{args.n_insertion_range[0]}_{args.n_insertion_range[1]}_ndr_{args.n_deletion_range[0]}_{args.n_deletion_range[1]}_prc_{args.prop_rc}"
    args.results_dir = results_dir  # add to the args object for simplicity
    print(f"Results will be saved to: {results_dir}")

    # Create directory structure BEFORE simulation
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(f"{results_dir}/simulated_data", exist_ok=True)
    os.makedirs(f"{results_dir}/raw_outputs", exist_ok=True)
    os.makedirs(f"{results_dir}/bash_scripts", exist_ok=True)
    print("Created directory structure")

    print("Simulating data...")
    contigs, spacers, ground_truth, myers_ground_truth = simulate_data_rust(
        contigs=args.contigs,
        spacers=args.spacers,
        contig_length_range=args.contig_length_range,
        spacer_length_range=args.spacer_length_range,
        n_mismatch_range=args.n_mismatch_range,
        sample_size_contigs=args.sample_size_contigs,
        sample_size_spacers=args.sample_size_spacers,
        insertion_range=args.insertion_range,
        n_insertion_range=args.n_insertion_range,
        n_deletion_range=args.n_deletion_range,
        prop_rc=args.prop_rc,
        debug=True,
        threads=args.threads,
        results_dir=results_dir,
        id_prefix=args.id_prefix,
    )
    print(f"Generated {len(contigs)} contigs and {len(spacers)} spacers")

    args.sample_size_contigs = len(contigs)
    args.sample_size_spacers = len(spacers)

    spacer_lendf = pl.DataFrame(
        {"spacer_id": spacers.keys(), "length": [len(seq) for seq in spacers.values()]}
    )
    print("Created spacer length dataframe")

    print("Initializing tools...")
    tools = populate_tools(args, spacer_lendf=spacer_lendf)
    print(f"Initialized {len(tools)} tools")

    if args.skip_tools:
        print(f"Disabling user-specified tools: {args.skip_tools}")
        skip_tools = args.skip_tools.split(",")
        tools = {k: v for k, v in tools.items() if k not in skip_tools}

    # Run tools and get results
    print("Running tools...")
    run_tools(tools, results_dir, max_runs=args.max_runs, warmups=args.warmups)
    print("Finished running tools")

    print("Reading tool results...")
    tools_results = read_results(
        tools,
        max_mismatches=args.max_mismatches,
        spacer_lendf=spacer_lendf,
        ref_file=f"{results_dir}/simulated_data/simulated_contigs.fa",
    )
    tools_results.write_csv(f"{results_dir}/tools_results.tsv", separator="\t")

    print(f"Wrote tool results to {results_dir}/tools_results.tsv")

    print("Reading hyperfine results...")
    hyperfine_results = read_hyperfine_results(tools, results_dir)
    hyperfine_results.write_csv(
        file=f"{results_dir}/hyperfine_results.tsv", separator="\t"
    )
    print(f"Wrote hyperfine results to {results_dir}/hyperfine_results.tsv")

    if myers_ground_truth is not None or len(myers_ground_truth) > 0:
        print("Comparing results against ground truth...")
        performance_results = compare_results(
            tools, myers_ground_truth, tools_results, soft_false_positive_threshold=1
        )
        print(performance_results)
        performance_results.write_csv(
            f"{results_dir}/performance_results.tsv", separator="\t"
        )
        print(f"Wrote performance results to {results_dir}/performance_results.tsv")
    else:
        print("No ground truth provided, skipping performance comparison")

    print("Benchmark completed successfully")


if __name__ == "__main__":
    main()
