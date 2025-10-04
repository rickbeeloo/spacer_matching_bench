import argparse
import os
import random
import logging

import polars as pl

from bench.utils.functions import *

pl.Config(tbl_cols=13, tbl_rows=8)


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
        "--spacers", 
        type=str,
        required=True,
        help="Path to spacers FASTA file"
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
    
    logger.info("Starting intelligent subsampling...")
    
    # Parse reduction factor
    try:
        reduce_val = float(args.reduce)
        if 0 < reduce_val <= 1:
            is_fraction = True
            logger.info(f"Using fraction reduction: {reduce_val}")
        elif reduce_val > 1:
            is_fraction = False
            target_contigs = int(reduce_val)
            logger.info(f"Using absolute reduction: {target_contigs} contigs")
        else:
            raise ValueError("Reduction factor must be > 0")
    except ValueError as e:
        logger.error(f"Invalid reduction factor: {args.reduce}. Must be fraction (0-1) or absolute number (>1)")
        return
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(f"{args.output_dir}/subsampled_data", exist_ok=True)
    
    # Read metadata
    logger.info(f"Reading metadata from {args.metadata}")
    try:
        metadata = pl.read_csv(args.metadata, separator="\t")
        logger.info(f"Loaded metadata with {metadata.height} entries")
    except Exception as e:
        logger.error(f"Failed to read metadata: {e}")
        return
    
    # Use metadata table for length and GC statistics (no need to read sequences)
    logger.info("Using metadata table for length and GC statistics")
    
    # Get length and GC statistics from metadata
    if 'Length' in metadata.columns and 'GC' in metadata.columns:
        all_lengths = metadata['Length'].to_list()
        all_gc = metadata['GC'].to_list()
        contig_count = len(metadata)
        logger.info(f"Using metadata for {contig_count} contigs")
    else:
        logger.error("Required columns 'Length' and 'GC' not found in metadata table")
        return
    
    # Filter by quality if specified
    if args.hq < 1.0:
        logger.info(f"Filtering for high-quality contigs (fraction: {args.hq})")
        # Assuming 'MIUViG quality' column exists and 'High-quality' is the value we want
        if 'MIUViG quality' in metadata.columns:
            hq_contigs = metadata.filter(pl.col('MIUViG quality') == 'High-quality')
            hq_count = int(len(hq_contigs) * args.hq)
            hq_contigs = hq_contigs.head(hq_count)
            logger.info(f"Selected {hq_contigs.height} high-quality contigs")
        else:
            logger.warning("MIUViG quality column not found, skipping quality filtering")
            hq_contigs = metadata
    else:
        hq_contigs = metadata
    
    # Determine target number of contigs
    if is_fraction:
        target_contigs = int(len(hq_contigs) * reduce_val)
    else:
        target_contigs = min(target_contigs, len(hq_contigs))
    
    logger.info(f"Target number of contigs: {target_contigs}")
    
    # Analyze taxonomic diversity
    taxon_col = f"Taxonomic classification"
    if taxon_col in hq_contigs.columns:
        # Extract taxonomic rank
        rank_map = {'d': 0, 'k': 0, 'p': 1, 'c': 2, 'o': 3, 'f': 4, 'g': 5, 's': 6}
        rank_idx = rank_map.get(args.taxonomic_rank, 1)  # default to phylum
        
        # Extract taxa at specified rank
        taxa = []
        for row in hq_contigs.iter_rows(named=True):
            taxon_str = row.get(taxon_col, "")
            if taxon_str and taxon_str != "":
                # Split by semicolon and get the rank
                taxon_parts = [part.strip() for part in taxon_str.split(';')]
                if len(taxon_parts) > rank_idx:
                    taxon = taxon_parts[rank_idx]
                    if taxon and not taxon.startswith('__'):
                        taxa.append(taxon)
        
        unique_taxa = list(set(taxa))
        logger.info(f"Found {len(unique_taxa)} unique taxa at rank {args.taxonomic_rank}")
        
        if len(unique_taxa) > target_contigs:
            logger.warning(f"More unique taxa ({len(unique_taxa)}) than target contigs ({target_contigs}). "
                          f"Diversity may not be fully represented.")
    else:
        logger.warning(f"Taxonomic classification column not found: {taxon_col}")
        unique_taxa = []
    
    # Create bins for GC content and length (using metadata values)
    all_lengths = all_lengths  # Already extracted from metadata
    all_gc = all_gc  # Already extracted from metadata
    
    length_min, length_max = min(all_lengths), max(all_lengths)
    gc_min, gc_max = min(all_gc), max(all_gc)
    
    logger.info(f"Length range: {length_min} - {length_max} bp")
    logger.info(f"GC content range: {gc_min:.3f} - {gc_max:.3f}")
    
    # Create bins
    length_bins = [length_min + i * (length_max - length_min) / args.length_bins 
                  for i in range(args.length_bins + 1)]
    gc_bins = [gc_min + i * (gc_max - gc_min) / args.gc_bins 
               for i in range(args.gc_bins + 1)]
    
    # Stratified sampling
    logger.info("Performing stratified sampling...")
    selected_contigs = []
    
    # Group by taxonomic diversity first
    if unique_taxa:
        taxa_groups = {}
        for row in hq_contigs.iter_rows(named=True):
            contig_id = row.get('UVIG', '')
            taxon_str = row.get(taxon_col, "")
            if taxon_str:
                taxon_parts = [part.strip() for part in taxon_str.split(';')]
                if len(taxon_parts) > rank_idx:
                    taxon = taxon_parts[rank_idx]
                    if taxon and not taxon.startswith('__'):
                        if taxon not in taxa_groups:
                            taxa_groups[taxon] = []
                        taxa_groups[taxon].append(contig_id)
        
        # Sample from each taxonomic group
        contigs_per_taxon = max(1, target_contigs // len(unique_taxa))
        for taxon, contig_list in taxa_groups.items():
            sample_size = min(contigs_per_taxon, len(contig_list))
            selected = random.sample(contig_list, sample_size)
            selected_contigs.extend(selected)
            logger.info(f"Selected {len(selected)} contigs from {taxon}")
    else:
        # Random sampling if no taxonomic info
        available_contigs = [row['UVIG'] for row in hq_contigs.iter_rows(named=True)]
        selected_contigs = random.sample(available_contigs, min(target_contigs, len(available_contigs)))
    
    # If we need more contigs, add more randomly
    if len(selected_contigs) < target_contigs:
        remaining_contigs = [row['UVIG'] for row in hq_contigs.iter_rows(named=True) 
                           if row['UVIG'] not in selected_contigs]
        additional_needed = target_contigs - len(selected_contigs)
        additional = random.sample(remaining_contigs, min(additional_needed, len(remaining_contigs)))
        selected_contigs.extend(additional)
        logger.info(f"Added {len(additional)} additional contigs to reach target")
    
    logger.info(f"Selected {len(selected_contigs)} contigs total")
    
    # Second pass: stream and write only selected contigs using fast O(1) membership test
    contigs_output = f"{args.output_dir}/subsampled_data/subsampled_contigs.fa"
    logger.info(f"Second pass: writing selected contigs to {contigs_output}")
    
    # Create target set and iterator for fast lookup
    selected_contig_set = set(selected_contigs)
    target_iter = iter(selected_contigs)
    try:
        next_target = next(target_iter)
    except StopIteration:
        next_target = None
    
    written_count = 0
    
    # Fast streaming approach with O(1) membership test
    with open(contigs_output, 'w') as f_out:
        for record in parse_fastx_file(args.contigs):
            contig_id = record.id
            # Fast O(1) membership test
            if contig_id in selected_contig_set:
                f_out.write(f">{contig_id}\n{record.seq}\n")
                written_count += 1
                # Advance to next needed contig
                try:
                    next_target = next(target_iter)
                except StopIteration:
                    next_target = None  # no more contigs needed
                if written_count % 1000 == 0:
                    logger.info(f"Written {written_count} contigs...")
    
    logger.info(f"Second pass complete: written {written_count} subsampled contigs")
    
    # Copy all spacers (no subsampling)
    spacers_output = f"{args.output_dir}/subsampled_data/subsampled_spacers.fa"
    logger.info(f"Copying all spacers to {spacers_output}")
    
    written_spacers = 0
    with open(spacers_output, 'w') as f_out:
        for record in parse_fastx_file(args.spacers):
            f_out.write(f">{record.id}\n{record.seq}\n")
            written_spacers += 1
            if written_spacers % 10000 == 0:
                logger.info(f"Written {written_spacers} spacers...")
    
    logger.info(f"Copied all {written_spacers} spacers")
    
    # Write subsampled metadata
    selected_contig_set = set(selected_contigs)
    subsampled_metadata = hq_contigs.filter(pl.col('UVIG').is_in(selected_contig_set))
    metadata_output = f"{args.output_dir}/subsampled_data/subsampled_metadata.tsv"
    logger.info(f"Writing subsampled metadata to {metadata_output}")
    subsampled_metadata.write_csv(metadata_output, separator="\t")
    
    # Generate summary report
    report_output = f"{args.output_dir}/subsampling_report.txt"
    with open(report_output, 'w') as f:
        f.write("Subsampling Report\n")
        f.write("==================\n\n")
        f.write(f"Original contigs: {contig_count}\n")
        f.write(f"Subsampled contigs: {len(selected_contigs)}\n")
        f.write(f"Reduction factor: {len(selected_contigs)/contig_count:.3f}\n\n")
        f.write(f"Original spacers: {written_spacers}\n")
        f.write(f"Copied spacers: {written_spacers} (all spacers copied, no subsampling)\n\n")
        f.write(f"Length range: {length_min} - {length_max} bp\n")
        f.write(f"GC content range: {gc_min:.3f} - {gc_max:.3f}\n\n")
        f.write(f"Unique taxa at rank {args.taxonomic_rank}: {len(unique_taxa)}\n")
        f.write(f"High-quality contigs selected: {args.hq}\n")
    
    logger.info(f"Subsampling completed successfully!")
    logger.info(f"Results saved to: {args.output_dir}")
    logger.info(f"Summary report: {report_output}")


if __name__ == "__main__":
    main()
