import argparse
import os

import polars as pl

from bench.utils.functions import *

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


def main():
    """Run only the sequence simulation step"""
    parser = argparse.ArgumentParser(
        description="Generate simulated sequences (contigs and spacers) with ground truth",
        add_help=True,
        allow_abbrev=True,
    )
    
    add_common_args(parser)
    parser.add_argument(
        "-o", "--output_dir", 
        type=str, 
        default=None,
        help="Output directory for simulated data (default: auto-generated based on parameters)"
    )
    
    args = parser.parse_args()
    
    print("Running sequence simulation...")
    
    # Determine output directory
    if args.output_dir is None:
        output_dir = f"results/sim_t_{args.threads}_nc_{args.sample_size_contigs}_ns_{args.sample_size_spacers}_ir_{args.insertion_range[0]}_{args.insertion_range[1]}_lm_{args.n_mismatch_range[0]}_{args.n_mismatch_range[1]}_nir_{args.n_insertion_range[0]}_{args.n_insertion_range[1]}_ndr_{args.n_deletion_range[0]}_{args.n_deletion_range[1]}_prc_{args.prop_rc}"
    else:
        output_dir = args.output_dir
    
    print(f"Simulation output will be saved to: {output_dir}")
    
    # Create directory structure
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/simulated_data", exist_ok=True)
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
        results_dir=output_dir,
        id_prefix=args.id_prefix,
    )
    print(f"Generated {len(contigs)} contigs and {len(spacers)} spacers")
    print(f"Simulation completed successfully. Data saved to {output_dir}")


if __name__ == "__main__":
    main()
