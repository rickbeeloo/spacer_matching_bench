import argparse
import os
import random
import logging
from typing import Union, List, Dict, Tuple, Optional

import polars as pl

from bench.utils.functions import *

pl.Config(tbl_cols=13, tbl_rows=8)


def parse_reduction_factor(reduce_val: str) -> Tuple[float, bool]:
    """Parse reduction factor and determine if it's a fraction or absolute number.
    
    Args:
        reduce_val: String representation of reduction factor
        
    Returns:
        Tuple of (value, is_fraction)
    """
    try:
        reduce_val = float(reduce_val)
        if 0 < reduce_val <= 1:
            return reduce_val, True
        elif reduce_val > 1:
            return int(reduce_val), False
        else:
            raise ValueError("Reduction factor must be > 0")
    except ValueError as e:
        raise ValueError(f"Invalid reduction factor: {reduce_val}. Must be fraction (0-1) or absolute number (>1)")



def perform_stratified_sampling(
    metadata: pl.DataFrame, 
    target_contigs: int, 
    taxonomic_rank: str,
    gc_bins: int = 10,
    length_bins: int = 10
) -> List[str]:
    """Perform stratified sampling to preserve taxonomic, GC content, and length diversity.
    
    Args:
        metadata: Metadata DataFrame with taxonomic rank columns
        target_contigs: Target number of contigs to select
        taxonomic_rank: Column name for taxonomic rank to use for stratification
        gc_bins: Number of GC content bins for stratification
        length_bins: Number of length bins for stratification
        
    Returns:
        List of selected contig IDs
    """
    selected_contigs = []
    
    # Create bins for GC content and length
    gc_min, gc_max = metadata['gc'].min(), metadata['gc'].max()
    length_min, length_max = metadata['length'].min(), metadata['length'].max()
    
    # Add bin columns
    metadata = metadata.with_columns([
        pl.col('gc').qcut(gc_bins, labels=[f"gc_bin_{i}" for i in range(gc_bins)]).alias('gc_bin'),
        pl.col('length').qcut(length_bins, labels=[f"length_bin_{i}" for i in range(length_bins)]).alias('length_bin')
    ])
    
    # Check if taxonomic rank column exists
    if taxonomic_rank not in metadata.columns:
        print(f"Warning: Taxonomic rank column '{taxonomic_rank}' not found. Using GC and length stratification only.")
        # Stratify by GC and length bins only
        unique_combinations = metadata.select(['gc_bin', 'length_bin']).unique()
        
        # Check if we have valid combinations
        if len(unique_combinations) == 0:
            print("No valid GC/length combinations found, using random sampling")
            available_contigs = metadata['seqid'].to_list()
            return random.sample(available_contigs, min(target_contigs, len(available_contigs)))
        
        contigs_per_combo = max(1, target_contigs // len(unique_combinations))
        
        for row in unique_combinations.iter_rows(named=True):
            combo_contigs = metadata.filter(
                (pl.col('gc_bin') == row['gc_bin']) & 
                (pl.col('length_bin') == row['length_bin'])
            )['seqid'].to_list()
            sample_size = min(contigs_per_combo, len(combo_contigs))
            selected = random.sample(combo_contigs, sample_size)
            selected_contigs.extend(selected)
            if len(selected) > 0:
                print(f"Selected {len(selected)} contigs from GC bin {row['gc_bin']}, length bin {row['length_bin']}")
    else:
        # Get unique taxa, filtering out null values
        unique_taxa = metadata.filter(pl.col(taxonomic_rank).is_not_null())[taxonomic_rank].unique().to_list()
        print(f"Found {len(unique_taxa)} unique taxa at rank {taxonomic_rank}")
        
        if len(unique_taxa) > 0:
            # For each taxon, stratify by GC and length bins
            contigs_per_taxon = max(1, target_contigs // len(unique_taxa))
            
            for taxon in unique_taxa:
                # Filter out null values to avoid the warning
                taxon_data = metadata.filter(
                    pl.col(taxonomic_rank).is_not_null() & 
                    (pl.col(taxonomic_rank) == taxon)
                )
                taxon_combinations = taxon_data.select(['gc_bin', 'length_bin']).unique()
                
                # Skip if no valid combinations for this taxon
                if len(taxon_combinations) == 0:
                    print(f"No valid GC/length combinations for {taxon}, skipping")
                    continue
                
                contigs_per_combo = max(1, contigs_per_taxon // len(taxon_combinations))
                
                taxon_selected = 0
                for row in taxon_combinations.iter_rows(named=True):
                    combo_contigs = taxon_data.filter(
                        (pl.col('gc_bin') == row['gc_bin']) & 
                        (pl.col('length_bin') == row['length_bin'])
                    )['seqid'].to_list()
                    sample_size = min(contigs_per_combo, len(combo_contigs))
                    selected = random.sample(combo_contigs, sample_size)
                    selected_contigs.extend(selected)
                    taxon_selected += len(selected)
                
                if taxon_selected > 0:
                    print(f"Selected {taxon_selected} contigs from {taxon}")
        else:
            # Random sampling if no taxonomic info
            available_contigs = metadata['seqid'].to_list()
            selected_contigs = random.sample(available_contigs, min(target_contigs, len(available_contigs)))
    
    # If we need more contigs, add more randomly from remaining
    if len(selected_contigs) < target_contigs:
        selected_set = set(selected_contigs)
        remaining_contigs = [seqid for seqid in metadata['seqid'].to_list() if seqid not in selected_set]
        additional_needed = target_contigs - len(selected_contigs)
        if len(remaining_contigs) > 0:
            additional = random.sample(remaining_contigs, min(additional_needed, len(remaining_contigs)))
            selected_contigs.extend(additional)
            print(f"Added {len(additional)} additional contigs to reach target")
        else:
            print(f"Warning: Could only select {len(selected_contigs)} contigs out of {target_contigs} requested")
    
    return selected_contigs


def write_subsampled_sequences(seq_file: str, selected_seqs: List[str], output_file: str, extract_method: str = "iter") -> int:
    """Write subsampled contigs to output file.
    
    Args:
        contigs_file: Path to input contigs FASTA file
        selected_contigs: List of selected contig IDs
        output_file: Path to output FASTA file
        
    Returns:
        Number of sequences written
    """
    if extract_method == "iter":
        selected_contig_set = set(selected_seqs)
        target_iter = iter(selected_seqs)
        try:
            next_target = next(target_iter)
        except StopIteration:
            next_target = None  # nothing to extract
        written_count = 0
        f_out = open(output_file, 'wt', encoding='utf-8')
        for record in parse_fastx_file(seq_file):
            # stop iterating once we passed the last needed record
            if next_target is None and record.id not in selected_contig_set:
                break
            # fast O(1) membership test
            if record.id in selected_contig_set:
                # print(record.id)
                f_out.write(f">{record.id}\n{record.seq}\n")
                written_count += 1
                if written_count % 100 == 0:
                    print(f"Written {written_count} contigs...")
                # advance to the next needed index
                try:
                    next_target = next(target_iter)
                except StopIteration:
                    next_target = None   # no more records needed
            # close the input file (handled automatically if using a context manager)
        f_out.close()
        return written_count
    elif extract_method == "pyfastx":
        return get_seq_from_fastx(seq_file, selected_seqs, return_df=False, output_file=output_file)
    else:
        raise ValueError(f"Invalid extract method: {extract_method}")


# def write_subsampled_spacers(spacers_file: str, output_file: str) -> int:
#     """Copy all spacers to output file (no subsampling).
    
#     Args:
#         spacers_file: Path to input spacers FASTA file
#         output_file: Path to output FASTA file
        
#     Returns:
#         Number of spacers written
#     """
#     written_spacers = 0
#     with open(output_file, 'w') as f_out:
#         for record in parse_fastx_file(spacers_file):
#             f_out.write(f">{record.id}\n{record.seq}\n")
#             written_spacers += 1
#             if written_spacers % 10000 == 0:
#                 print(f"Written {written_spacers} spacers...")
    
#     return written_spacers


def create_subsampling_report(
    output_dir: str,
    original_contigs: int,
    subsampled_contigs: int,
    length_min: float,
    length_max: float,
    gc_min: float,
    gc_max: float,
    unique_taxa_count: int,
    hq_fraction: float,
    taxonomic_rank: str
) -> str:
    """Create a summary report of the subsampling process.
    
    Args:
        output_dir: Output directory path
        original_contigs: Number of original contigs
        subsampled_contigs: Number of subsampled contigs
        length_min: Minimum length
        length_max: Maximum length
        gc_min: Minimum GC content
        gc_max: Maximum GC content
        unique_taxa_count: Number of unique taxa
        hq_fraction: High-quality fraction used
        taxonomic_rank: Taxonomic rank used
        
    Returns:
        Path to the report file
    """
    report_output = f"{output_dir}/subsampling_report.txt"
    with open(report_output, 'w') as f:
        f.write("Subsampling Report\n")
        f.write("==================\n\n")
        f.write(f"Original contigs: {original_contigs}\n")
        f.write(f"Subsampled contigs: {subsampled_contigs}\n")
        f.write(f"Reduction factor: {subsampled_contigs/original_contigs:.3f}\n\n")
        f.write(f"Length range: {length_min} - {length_max} bp\n")
        f.write(f"GC content range: {gc_min:.3f} - {gc_max:.3f}\n\n")
        f.write(f"Unique taxa at rank {taxonomic_rank}: {unique_taxa_count}\n")
        f.write(f"High-quality contigs selected: {hq_fraction}\n")
    
    return report_output


def subsample_dataset(
    contigs_file: str,
    metadata_file: str,
    output_dir: str,
    reduce_factor: Union[str, float, int],
    hq_fraction: float = 1.0,
    taxonomic_rank: str = 'p',
    gc_bins: int = 10,
    length_bins: int = 10,
    logger: Optional[logging.Logger] = None,
    extract_method: str = "iter" # "iter" or "pyfastx"
) -> pl.DataFrame:
    """Main function to subsample a dataset while preserving diversity.
    
    Args:
        contigs_file: Path to contigs FASTA file
        metadata_file: Path to metadata TSV file
        output_dir: Output directory for subsampled data
        reduce_factor: Reduction factor (fraction 0-1 or absolute number >1)
        hq_fraction: Fraction of high-quality contigs to select
        taxonomic_rank: Taxonomic rank for diversity preservation
        gc_bins: Number of GC content bins for stratification
        length_bins: Number of length bins for stratification
        logger: Optional logger instance    
        extract_method: Method to extract sequences from the contigs file, "iter" or "pyfastx"
        
        
    Returns:
        Polars DataFrame with subsampled metadata
    """
    if logger is None:
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    
    logger.info("Starting intelligent subsampling...")
    
    # Parse reduction factor
    reduce_val, is_fraction = parse_reduction_factor(str(reduce_factor))
    if is_fraction:
        logger.info(f"Using fraction reduction: {reduce_val}")
    else:
        logger.info(f"Using absolute reduction: {reduce_val} contigs")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read metadata
    logger.info(f"Reading metadata from {metadata_file}")
    try:
        if metadata_file.endswith(".parquet"):
            metadata = pl.read_parquet(metadata_file)
        elif metadata_file.endswith(".csv"):
            metadata = pl.read_csv(metadata_file, separator=",")
        elif metadata_file.endswith(".tsv"):
            metadata = pl.read_csv(metadata_file, separator="\t")
        else:
            raise ValueError("Metadata file must be in .parquet, .csv, or .tsv format")
            
        logger.info(f"Loaded metadata with {metadata.height} entries")
    except Exception as e:
        logger.error(f"Failed to read metadata: {e}")
        raise
    
    # Validate required columns
    if 'length' not in metadata.columns or 'gc' not in metadata.columns:
        raise ValueError("Required columns 'length' and 'gc' not found in metadata table")
    
    # Get length and GC statistics from metadata
    all_lengths = metadata['length'].to_list()
    all_gc = metadata['gc'].to_list()
    contig_count = len(metadata)
    logger.info(f"Using metadata for {contig_count} contigs")
    
    # Filter by quality if specified
    hq_contigs = metadata.filter(pl.col('hq') == True)
    hq_contigs = hq_contigs.sample(fraction=hq_fraction)

    logger.info(f"Selected {hq_contigs.height} high-quality contigs")
    
    # Determine target number of contigs
    if is_fraction:
        target_contigs = int(len(hq_contigs) * reduce_val)
    else:
        target_contigs = min(reduce_val, len(hq_contigs))
    
    logger.info(f"Target number of contigs: {target_contigs}")
    
    # Get statistics
    length_min, length_max = min(all_lengths), max(all_lengths)
    gc_min, gc_max = min(all_gc), max(all_gc)
    
    logger.info(f"Length range: {length_min} - {length_max} bp")
    logger.info(f"GC content range: {gc_min:.3f} - {gc_max:.3f}")
    
    # Perform stratified sampling
    logger.info("Performing stratified sampling...")
    selected_contigs = perform_stratified_sampling(hq_contigs, target_contigs, taxonomic_rank, gc_bins, length_bins)
    logger.info(f"Selected {len(selected_contigs)} contigs total")
    
    # Write subsampled contigs
    contigs_output = f"{output_dir}/subsampled_contigs.fa"
    logger.info(f"Writing selected contigs to {contigs_output}")
    written_contigs = write_subsampled_sequences(contigs_file, selected_contigs, contigs_output, extract_method)
    logger.info(f"Written {written_contigs} subsampled contigs")
    
    # Write subsampled metadata
    selected_contig_set = set(selected_contigs)
    subsampled_metadata = hq_contigs.filter(pl.col('seqid').is_in(selected_contig_set))
    metadata_output = f"{output_dir}/subsampled_metadata.tsv"
    logger.info(f"Writing subsampled metadata to {metadata_output}")
    subsampled_metadata.write_csv(metadata_output, separator="\t")
    
    logger.info(f"Subsampling completed successfully!")
    logger.info(f"Results saved to: {output_dir}")
    logger.info(f"Original contigs: {contig_count}")
    logger.info(f"Subsampled contigs: {subsampled_metadata.height}")
    logger.info(f"Reduction factor: {(metadata.height) / subsampled_metadata.height:.3f}")
    
    logger.info(f"Results saved to: {output_dir}")
    logger.info(f"Original contigs: {metadata.height}")
    logger.info(f"Subsampled contigs: {subsampled_metadata.height}")
    logger.info(f"Reduction factor: {(metadata.height) / subsampled_metadata.height:.3f}")
    logger.info(f"Sampled contigs stats:")
    for col in ["length", "gc"]:
        logger.info(f"{col} range: {subsampled_metadata[col].min()} - {subsampled_metadata[col].max()}")
    logger.info(f"{taxonomic_rank} value  counts: {subsampled_metadata[taxonomic_rank].value_counts(sort=True)}")
    
        
    return subsampled_metadata


def main():
    """Intelligently subsample real dataset while preserving diversity"""
    parser = argparse.ArgumentParser(
        description="Intelligently subsample real dataset while preserving diversity",
        add_help=True,
        allow_abbrev=True,
    )
    
    parser.add_argument(
        "--contigs",
        type=str,
        required=True,
        help="Path to contigs FASTA file"
    )
    parser.add_argument(
        "--metadata",
        type=str,
        required=True,
        help="Path to metadata TSV file"
    )
    parser.add_argument(
        "--reduce",
        type=str,
        required=True,
        help="Reduction factor: fraction (0-1) or absolute number (>1)"
    )
    parser.add_argument(
        "--hq",
        type=float,
        default=1.0,
        help="Fraction of high-quality contigs to select [default: 1.0]"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory for subsampled data"
    )
    parser.add_argument(
        "--taxonomic_rank",
        type=str,
        default="p",
        help="Taxonomic rank for diversity preservation [default: p (phylum)]"
    )
    parser.add_argument(
        "--gc_bins",
        type=int,
        default=10,
        help="Number of GC content bins [default: 10]"
    )
    parser.add_argument(
        "--length_bins",
        type=int,
        default=10,
        help="Number of length bins [default: 10]"
    )
    
    args = parser.parse_args()
    
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    
    # Use the functionalized subsample_dataset function
    try:
        results = subsample_dataset(
            contigs_file=args.contigs,
            metadata_file=args.metadata,
            output_dir=args.output_dir,
            reduce_factor=args.reduce,
            hq_fraction=args.hq,
            taxonomic_rank=args.taxonomic_rank,
            gc_bins=args.gc_bins,
            length_bins=args.length_bins,
            logger=logger
        )
        results['seqid'].to_list()
        # Print summary

    except Exception as e:
        logger.error(f"Subsampling failed: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    main()
