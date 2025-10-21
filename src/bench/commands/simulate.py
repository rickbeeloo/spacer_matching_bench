import argparse
import os


from bench.utils.functions import *
from bench.utils.arg_parsers import add_common_args


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
        # New parameters
        contig_distribution=args.contig_distribution,
        spacer_distribution=args.spacer_distribution,
        gc_content=args.gc_content,
        a_frac=args.a_frac,
        t_frac=args.t_frac,
        c_frac=args.c_frac,
        g_frac=args.g_frac,
    )
    print(f"Generated {len(contigs)} contigs and {len(spacers)} spacers")
    print(f"Simulation completed successfully. Data saved to {output_dir}")


if __name__ == "__main__":
    main()
