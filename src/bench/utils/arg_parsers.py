"""Shared argument parser functions for bench commands."""



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
    
    # Base composition arguments
    parser.add_argument(
        "--gc_content",
        type=float,
        default=None,
        help="GC content percentage (0-100). Mutually exclusive with individual base fractions.",
    )
    parser.add_argument(
        "--a_frac",
        type=float,
        default=None,
        help="A fraction (0-1). Requires t_frac, c_frac, g_frac to be specified.",
    )
    parser.add_argument(
        "--t_frac",
        type=float,
        default=None,
        help="T fraction (0-1). Requires a_frac, c_frac, g_frac to be specified.",
    )
    parser.add_argument(
        "--c_frac",
        type=float,
        default=None,
        help="C fraction (0-1). Requires a_frac, t_frac, g_frac to be specified.",
    )
    parser.add_argument(
        "--g_frac",
        type=float,
        default=None,
        help="G fraction (0-1). Requires a_frac, t_frac, c_frac to be specified.",
    )
    
    # Distribution type arguments
    parser.add_argument(
        "--contig_distribution",
        type=str,
        choices=["uniform", "normal", "bell"],
        default=None,
        help="Distribution type for contig lengths (default: uniform)",
    )
    parser.add_argument(
        "--spacer_distribution",
        type=str,
        choices=["uniform", "normal", "bell"],
        default=None,
        help="Distribution type for spacer lengths (default: uniform)",
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
    parser.add_argument(
        "-rt",
        "--tools_to_run",
        type=str,
        default=None,
        help="Comma-separated list of tools to run (default is all tools - skip_tools, so use tools_to_run to specify a subset of tools to run)",
    )
