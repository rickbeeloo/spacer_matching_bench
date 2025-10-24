"""
Sequence simulation command. mostly a wrapper for the rust simulator.
"""
import os
import logging

from bench.utils.functions import *

logger = logging.getLogger(__name__)


def run_simulate(num_contigs, num_spacers, contig_length_range, spacer_length_range,
                 mismatch_range, spacer_insertions, indel_insertions, indel_deletions,
                 reverse_complement, threads, output_dir, id_prefix=None, gc_content=None,
                 contig_distribution='uniform', spacer_distribution='uniform', verify=False,
                 contigs=None, spacers=None, a_frac=None, t_frac=None, c_frac=None, g_frac=None):
    """
    Core simulation function that can be called from both CLI interfaces.
    
    Args:
        num_contigs: Number of contigs to generate
        num_spacers: Number of spacers to generate
        contig_length_range: Tuple of (min, max) contig lengths
        spacer_length_range: Tuple of (min, max) spacer lengths
        mismatch_range: Tuple of (min, max) substitution mismatches per spacer
        spacer_insertions: Tuple of (min, max) number of times to insert each spacer
        indel_insertions: Tuple of (min, max) insertion mutations within spacers
        indel_deletions: Tuple of (min, max) deletion mutations within spacers
        reverse_complement: Proportion of spacers to reverse complement (0.0-1.0)
        threads: Number of threads for parallel processing
        output_dir: Output directory for simulated data
        id_prefix: Prefix for sequence IDs (optional)
        gc_content: GC content percentage (optional)
        contig_distribution: Distribution type for contig lengths
        spacer_distribution: Distribution type for spacer lengths
        verify: Whether to verify simulation after generation
        contigs: Path to existing contigs file (optional)
        spacers: Path to existing spacers file (optional)
        a_frac: A base fraction (optional)
        t_frac: T base fraction (optional)
        c_frac: C base fraction (optional)
        g_frac: G base fraction (optional)
    
    Returns:
        Tuple of (contigs, spacers, ground_truth) dictionaries/lists
    """
    logger.info(f"Starting simulation: {num_contigs} contigs, {num_spacers} spacers")
    logger.debug(f"Output directory: {output_dir}")
    
    # Create directory structure
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/simulated_data", exist_ok=True)
    logger.debug("Created directory structure")
    
    logger.info("Running Rust simulator...")
    contigs_out, spacers_out, ground_truth = simulate_data_rust(
        contigs=contigs,
        spacers=spacers,
        contig_length_range=contig_length_range,
        spacer_length_range=spacer_length_range,
        n_mismatch_range=mismatch_range,
        sample_size_contigs=num_contigs,
        sample_size_spacers=num_spacers,
        insertion_range=spacer_insertions,
        n_insertion_range=indel_insertions,
        n_deletion_range=indel_deletions,
        prop_rc=reverse_complement,
        debug=False,
        threads=threads,
        results_dir=output_dir,
        id_prefix=id_prefix,
        contig_distribution=contig_distribution,
        spacer_distribution=spacer_distribution,
        gc_content=gc_content,
        a_frac=a_frac,
        t_frac=t_frac,
        c_frac=c_frac,
        g_frac=g_frac,
        verify=verify
    )
    
    logger.info(f"Generated {len(contigs_out)} contigs and {len(spacers_out)} spacers")
    logger.info(f"Ground truth contains {len(ground_truth)} spacer insertions")
    logger.info(f"Simulation completed. Data saved to {output_dir}/simulated_data/")
    
    return contigs_out, spacers_out, ground_truth


# Old argparse-based main() removed - use Click CLI via spacer_bencher command
