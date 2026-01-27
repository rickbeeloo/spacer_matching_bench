# some utility functions for the benchmark.
import glob
import hashlib
import json
import os
import re
from typing import Any, Optional, Union#, Dict, List, Tuple

import subprocess
import tempfile

import parasail as ps
import polars as pl
import psutil
import pyfastx as pfx
import pysam
from needletail import parse_fastx_file
from tqdm import tqdm
from collections import Counter
import math
import time
import polars_bio as pb
import _duckdb as duckdb
from rust_simulator import Simulator, BaseComposition, DistributionType


# so printing to stdout doesn't break, line wrap, or truncate.
pl.Config.set_tbl_rows(123123)
pl.Config.set_tbl_cols(123123)  # should be large enough
pl.Config.set_fmt_str_lengths(2100)
pl.Config.set_tbl_width_chars(2100)


global tab
global tab_b
global tab_extended
global tab_extended_b
tab = bytes.maketrans(b"ACTG", b"TGAC")
tab_b = bytearray(tab)
tab_extended = bytes.maketrans(b"ACGTURYSWKMBDHVN", b"TGCAAYRSWMKVHDBN")
tab_extended_b = bytearray(tab_extended)

def read_fasta_needletail(fasta_file: str) -> tuple[list[str], list[str]]:
    """Read sequences from a FASTA file using needletail.

    Args:
        fasta_file(str): Path to the FASTA file

    Returns:
        tuple: (seq_ids, seqs) where seq_ids is list of sequence IDs and seqs is list of sequences
    """
    seqs = []
    seq_ids = []
    for record in parse_fastx_file(fasta_file):
        seqs.append(record.seq)
        seq_ids.append(record.id)
    return seq_ids, seqs

def calculate_shannon_entropy(s: str) -> float:
    """Calculate the Shannon entropy of a string.

    Args:
        s(str): The input string

    Returns:
        float: The Shannon entropy value
    """
    if not s:
        return 0.0  # here in case of NULLs so that applting over ovject of size n still returns a number.  TODO: think if 0.0 is the best place holder value...

    # Count character frequencies
    char_counts = Counter(s)
    total_chars = len(s)

    entropy = 0.0
    for count in char_counts.values():
        probability = count / total_chars
        # entropy -= probability * math.log10(probability) # TODO: check which base other people use.
        entropy -= probability * math.log2(probability)
    return entropy


def count_kmers_df_explicit(
    df: pl.DataFrame,
    seq_col: str = "seq",
    id_col: str = "seqid",
    k: int = 3,
    relative: bool = False,
) -> pl.DataFrame:
    """Calculate ALL k-mers counts for all sequences in a DataFrame. all possible k-mers are counted, not just the complete ones."""
    # Split sequences into characters
    import itertools

    all_kmers = ["".join(p) for p in itertools.product("ATCGN", repeat=k)]
    count_df = (
        df.with_columns(
            pl.col(seq_col)
            .str.extract_many(all_kmers, overlapping=True, ascii_case_insensitive=True)
            .alias("kmers")
        )
        .group_by(id_col)
        .agg(
            pl.col("kmers")
            .explode()
            .value_counts(normalize=relative)
            .alias(f"kmer_{k}_relative" if relative else f"kmer_{k}_counts"),
        )
    )
    return count_df


def count_kmers_df(
    df: pl.DataFrame,
    seq_col: str = "seq",
    id_col: str = "seqid",
    k: int = 3,
    relative: bool = False,
) -> pl.DataFrame:
    """Calculate k-mer counts for all sequences in a DataFrame"""
    # Split sequences into characters
    split_chars_expr = pl.col(seq_col).str.split("").alias("chars")

    # Create k-mers by shifting and concatenating # TODO: look if this can be down in one step or if there is some sliding window function.
    create_kmers_expr = pl.concat_str(
        [pl.col("chars").shift(-i).over(id_col) for i in range(k)]
    ).alias("substrings")

    # Filter for complete k-mers only
    filter_complete_kmers_expr = pl.col("substrings").str.len_chars() == k

    # Aggregate expressions
    agg_exprs = [
        pl.first(seq_col),  # Keep the original sequence
        pl.col("substrings")
        .value_counts(normalize=relative)
        .alias(f"kmer_{k}_relative" if relative else f"kmer_{k}_counts"),
        pl.exclude(
            seq_col, "chars", "substrings"
        ).first(),  # Keep all other original columns
    ]

    return (
        df.with_columns(split_chars_expr)
        .explode("chars")
        .with_columns(create_kmers_expr)
        .filter(filter_complete_kmers_expr)
        .group_by(id_col, maintain_order=True)
        .agg(*agg_exprs)
    )


def filter_repetitive_kmers(
    df: pl.DataFrame,
    seq_col: str = "seq",
    id_col: str = "seqid",
    k: int = 6,
    max_count: int = 4,
) -> pl.DataFrame:
    """Filter sequences that have any k-mer appearing more than max_count times.
        Note, this is not used currently in the package but an almost identical process is done in the spacer_inspection.ipynb notebook.
    """
    # First get k-mer counts
    df_with_kmers = count_kmers_df(df, seq_col, id_col, k, relative=False)

    # Filter for sequences without highly repetitive k-mers
    filter_repetitive_expr = (
        ~pl.col("kmer_counts")
        .list.eval(pl.element().struct.field("count") > max_count)
        .list.any()
    )

    return df_with_kmers.filter(filter_repetitive_expr)


# lcc_mult and lcc_simp are sourced from the biopython library - we don't need to depend on it here just for these two functions.
def lcc_mult(seq: str, wsize: int) -> list[float]:
    """Calculate Local Composition Complexity (LCC) values over sliding window.
    sourced from: https://github.com/biopython/biopython/blob/e451db211bdd855a5d0f1f6bba18985ffee12696/Bio/SeqUtils/lcc.py#L13

    Args:
        seq(str): An unambiguous DNA sequence (a string or Seq object)
        wsize(int): Window size, integer

    Returns:
        list[float]: The LCC values for a sliding window over the sequence.

    Note:
        The result is the same as applying lcc_simp multiple times, but this
        version is optimized for speed. The optimization works by using the
        value of previous window as a base to compute the next one.
    """
    l4 = math.log(4)
    seq = seq.upper()
    tamseq = len(seq)
    compone = [0]
    lccsal = []
    for i in range(wsize):
        compone.append(((i + 1) / wsize) * math.log((i + 1) / wsize) / l4)
    window = seq[0:wsize]
    cant_a = window.count("A")
    cant_c = window.count("C")
    cant_t = window.count("T")
    cant_g = window.count("G")
    term_a = compone[cant_a]
    term_c = compone[cant_c]
    term_t = compone[cant_t]
    term_g = compone[cant_g]
    lccsal.append(-(term_a + term_c + term_t + term_g))
    tail = seq[0]
    for x in range(tamseq - wsize):
        window = seq[x + 1 : wsize + x + 1]
        if tail == window[-1]:
            lccsal.append(lccsal[-1])
        elif tail == "A":
            cant_a -= 1
            if window.endswith("C"):
                cant_c += 1
                term_a = compone[cant_a]
                term_c = compone[cant_c]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("T"):
                cant_t += 1
                term_a = compone[cant_a]
                term_t = compone[cant_t]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("G"):
                cant_g += 1
                term_a = compone[cant_a]
                term_g = compone[cant_g]
                lccsal.append(-(term_a + term_c + term_t + term_g))
        elif tail == "C":
            cant_c -= 1
            if window.endswith("A"):
                cant_a += 1
                term_a = compone[cant_a]
                term_c = compone[cant_c]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("T"):
                cant_t += 1
                term_c = compone[cant_c]
                term_t = compone[cant_t]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("G"):
                cant_g += 1
                term_c = compone[cant_c]
                term_g = compone[cant_g]
                lccsal.append(-(term_a + term_c + term_t + term_g))
        elif tail == "T":
            cant_t -= 1
            if window.endswith("A"):
                cant_a += 1
                term_a = compone[cant_a]
                term_t = compone[cant_t]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("C"):
                cant_c += 1
                term_c = compone[cant_c]
                term_t = compone[cant_t]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("G"):
                cant_g += 1
                term_t = compone[cant_t]
                term_g = compone[cant_g]
                lccsal.append(-(term_a + term_c + term_t + term_g))
        elif tail == "G":
            cant_g -= 1
            if window.endswith("A"):
                cant_a += 1
                term_a = compone[cant_a]
                term_g = compone[cant_g]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("C"):
                cant_c += 1
                term_c = compone[cant_c]
                term_g = compone[cant_g]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("T"):
                cant_t += 1
                term_t = compone[cant_t]
                term_g = compone[cant_g]
                lccsal.append(-(term_a + term_c + term_t + term_g))
        tail = window[0]
    return lccsal


# lcc_mult and lcc_simp are sourced from the biopython library - we don't need to depend on it here just for these two functions.
def lcc_simp(seq: str) -> float:
    """Calculate Local Composition Complexity (LCC) for a sequence.
    sourced from: https://github.com/biopython/biopython/blob/e451db211bdd855a5d0f1f6bba18985ffee12696/Bio/SeqUtils/lcc.py#L120

    Args:
        seq(str): An unambiguous DNA sequence (a string or Seq object)

    Returns:
        float: The Local Composition Complexity (LCC) value for the entire sequence.

    Reference:
        Andrzej K Konopka (2005) Sequence Complexity and Composition
        https://doi.org/10.1038/npg.els.0005260
    """
    wsize = len(seq)
    seq = seq.upper()
    l4 = math.log(4)
    # Check to avoid calculating the log of 0.
    if "A" not in seq:
        term_a = 0
    else:
        term_a = (seq.count("A") / wsize) * math.log(seq.count("A") / wsize) / l4
    if "C" not in seq:
        term_c = 0
    else:
        term_c = (seq.count("C") / wsize) * math.log(seq.count("C") / wsize) / l4
    if "T" not in seq:
        term_t = 0
    else:
        term_t = (seq.count("T") / wsize) * math.log(seq.count("T") / wsize) / l4
    if "G" not in seq:
        term_g = 0
    else:
        term_g = (seq.count("G") / wsize) * math.log(seq.count("G") / wsize) / l4
    return -(float(term_a + term_c + term_t + term_g))


def reverse_complement(seq: Union[bytes, str], return_str: bool = True) -> Union[bytes, str]:
    """Compute the reverse complement of a DNA sequence in bytes, handling ambiguous bases.
    Args:
        seq(bytes|str): The input DNA sequence in bytes or string format (will be converted to bytes if string)
        return_str(bool): Whether to return the result as a string (True) or bytes (False)
    Returns:
        bytes|str: The reverse complement sequence in bytes or string format (based on return_str)
    Note:
        This function is adapted from Jack Aidley's Stack overflow answer implemented in https://github.com/conchoecia/fastest_rc_python/blob/master/reverse_complement_tests.py#L91
    """
    if isinstance(seq, str):
        seq = seq.encode()
    rev_comp = seq.translate(tab_extended)[::-1]
    if return_str:
        return rev_comp.decode()
    else:
        return rev_comp



def generate_simulation_id(params: dict) -> str:
    """Generate a unique simulation ID based on input parameters.

    Args:
        params(dict): Dictionary of simulation parameters

    Returns:
        str: 8-character hex string unique to these parameters
    """
    # Sort the parameters to ensure consistent hashing
    param_str = json.dumps(params, sort_keys=True)
    # Generate hash
    hash_obj = hashlib.sha256(param_str.encode())
    # Return first 8 characters of hash
    return hash_obj.hexdigest()[:8]


def simulate_data_rust(
    contig_length_range: tuple[int, int],
    spacer_length_range: tuple[int, int],
    n_mismatch_range: tuple[int, int],
    sample_size_contigs: int,
    sample_size_spacers: int,
    insertion_range: tuple[int, int],
    n_insertion_range: tuple[int, int] = (0, 0),
    n_deletion_range: tuple[int, int] = (0, 0),
    contigs: Optional[Union[dict[str, str], str]] = None,
    spacers: Optional[Union[dict[str, str], str]] = None,
    prop_rc: float = 0.5,
    debug: bool = False,
    threads: Optional[int] = None,
    verify: bool = False,
    results_dir: Optional[str] = None,
    id_prefix: Optional[str] = None,
    # New parameters for base composition and distribution
    contig_distribution: Optional[Any] = None,
    spacer_distribution: Optional[Any] = None,
    base_composition: Optional[Any] = None,  # Deprecated: use contig/spacer specific params
    gc_content: Optional[float] = None,  # Deprecated: use contig/spacer specific params
    a_frac: Optional[float] = None,  # Deprecated: use contig/spacer specific params
    t_frac: Optional[float] = None,  # Deprecated: use contig/spacer specific params
    c_frac: Optional[float] = None,  # Deprecated: use contig/spacer specific params
    g_frac: Optional[float] = None,  # Deprecated: use contig/spacer specific params
    # New separate parameters
    contig_gc_content: Optional[float] = None,
    spacer_gc_content: Optional[float] = None,
    contig_a_frac: Optional[float] = None,
    contig_t_frac: Optional[float] = None,
    contig_c_frac: Optional[float] = None,
    contig_g_frac: Optional[float] = None,
    spacer_a_frac: Optional[float] = None,
    spacer_t_frac: Optional[float] = None,
    spacer_c_frac: Optional[float] = None,
    spacer_g_frac: Optional[float] = None,
    data_subdir: Optional[str] = None,
) -> tuple:
    """Simulate CRISPR spacer data using Rust implementation.

    Args:
        contig_length_range(tuple[int, int]): Range of contig lengths
        spacer_length_range(tuple[int, int]): Range of spacer lengths
        n_mismatch_range(tuple[int, int]): Range of mismatches
        sample_size_contigs(int): Number of contigs to sample
        sample_size_spacers(int): Number of spacers to sample
        insertion_range(tuple[int, int]): Range for insertions
        n_insertion_range(tuple[int, int]): Range for number of insertions
        n_deletion_range(tuple[int, int]): Range for number of deletions
        contigs(dict[str, str]): Pre-existing contigs
        spacers(dict[str, str]): Pre-existing spacers
        prop_rc(float): Proportion of reverse complement
        debug(bool): Enable debug mode
        threads(int): Number of threads
        verify(bool): Verify simulation
        results_dir(str): Results directory
        id_prefix(str): Simulation ID prefix
        contig_distribution: Distribution for contigs
        spacer_distribution: Distribution for spacers
        base_composition: Base composition (deprecated)
        gc_content(float): GC content (deprecated)
        a_frac(float): A fraction (deprecated)
        t_frac(float): T fraction (deprecated)
        c_frac(float): C fraction (deprecated)
        g_frac(float): G fraction (deprecated)
        contig_gc_content(float): GC content for contigs
        spacer_gc_content(float): GC content for spacers
        contig_a_frac(float): A fraction for contigs
        contig_t_frac(float): T fraction for contigs
        contig_c_frac(float): C fraction for contigs
        contig_g_frac(float): G fraction for contigs
        spacer_a_frac(float): A fraction for spacers
        spacer_t_frac(float): T fraction for spacers
        spacer_c_frac(float): C fraction for spacers
        spacer_g_frac(float): G fraction for spacers
        data_subdir(str): Data subdirectory

    Returns:
        tuple: (contigs_dict, spacers_dict, ground_truth_df)
    """
    # Use the Rust implementation for input file simulation
    if contigs is not None or spacers is not None:
        if debug:
            print("Using Rust simulate_from_input_files...")
            start_time = time.time()

        # Generate simulation ID if not provided
        if id_prefix is None:
            sim_params = {
                "contig_length_range": contig_length_range,
                "spacer_length_range": spacer_length_range,
                "n_mismatch_range": n_mismatch_range,
                "sample_size_contigs": sample_size_contigs,
                "sample_size_spacers": sample_size_spacers,
                "insertion_range": insertion_range,
                "prop_rc": prop_rc,
            }
            id_prefix = generate_simulation_id(sim_params)
            print(f"Using generated simulation ID: {id_prefix}")
        else:
            print(f"Using provided ID prefix: {id_prefix}")

        simulator = Simulator()

        # Handle base compositions for contigs and spacers separately
        contig_comp_obj = None
        spacer_comp_obj = None

        # Contig composition (use new params if provided, fallback to legacy params)
        if contig_gc_content is not None:
            contig_comp_obj = BaseComposition.from_gc_content(contig_gc_content)
        elif all(
            x is not None
            for x in [contig_a_frac, contig_t_frac, contig_c_frac, contig_g_frac]
        ):
            contig_comp_obj = BaseComposition.from_fractions(
                contig_a_frac, contig_t_frac, contig_c_frac, contig_g_frac
            )
        elif gc_content is not None:  # Fallback to legacy param
            contig_comp_obj = BaseComposition.from_gc_content(gc_content)
        elif all(
            x is not None for x in [a_frac, t_frac, c_frac, g_frac]
        ):  # Fallback to legacy param
            contig_comp_obj = BaseComposition.from_fractions(
                a_frac, t_frac, c_frac, g_frac
            )
        elif base_composition is not None:  # Fallback to legacy param
            contig_comp_obj = base_composition

        # Spacer composition (use new params if provided, fallback to legacy params)
        if spacer_gc_content is not None:
            spacer_comp_obj = BaseComposition.from_gc_content(spacer_gc_content)
        elif all(
            x is not None
            for x in [spacer_a_frac, spacer_t_frac, spacer_c_frac, spacer_g_frac]
        ):
            spacer_comp_obj = BaseComposition.from_fractions(
                spacer_a_frac, spacer_t_frac, spacer_c_frac, spacer_g_frac
            )
        elif gc_content is not None:  # Fallback to legacy param
            spacer_comp_obj = BaseComposition.from_gc_content(gc_content)
        elif all(
            x is not None for x in [a_frac, t_frac, c_frac, g_frac]
        ):  # Fallback to legacy param
            spacer_comp_obj = BaseComposition.from_fractions(
                a_frac, t_frac, c_frac, g_frac
            )
        elif base_composition is not None:  # Fallback to legacy param
            spacer_comp_obj = base_composition

        # Handle distribution types
        contig_dist_obj = None
        if contig_distribution is not None:
            contig_dist_obj = DistributionType(contig_distribution)

        spacer_dist_obj = None
        if spacer_distribution is not None:
            spacer_dist_obj = DistributionType(spacer_distribution)

        contigs_dict, spacers_dict, ground_truth = simulator.simulate_from_input_files(
            tuple(contig_length_range),
            tuple(spacer_length_range),
            tuple(n_mismatch_range),
            sample_size_contigs,
            sample_size_spacers,
            tuple(insertion_range),
            prop_rc,
            results_dir,
            id_prefix,
            contigs,  # Pass file path if string, None otherwise
            spacers,  # Pass file path if string, None otherwise
            data_subdir,
            debug,
            contig_dist_obj,
            spacer_dist_obj,
            contig_comp_obj,
            spacer_comp_obj,
            threads if threads is not None else 1,  # Add threads parameter
        )

        if debug:
            end_time = time.time()
            print(f"Simulation time when using Rust: {end_time - start_time} seconds")

        # Convert the ground truth data types before creating DataFrame
        processed_ground_truth = [
            [
                row[0],  # spacer_id (already has prefix from Rust)
                row[1],  # contig_id (already has prefix from Rust)
                int(row[2]),  # start (int)
                int(row[3]),  # end (int)
                row[4].lower() == "true",  # strand (bool)
                int(row[5]),  # mismatches (int)
            ]
            for row in ground_truth
        ]

        ground_truth_df = pl.DataFrame(
            processed_ground_truth,
            schema={
                "spacer_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "start": pl.UInt32,
                "end": pl.UInt32,
                "strand": pl.Boolean,
                "mismatches": pl.UInt32,
            },
            orient="row",
        )

        if verify:
            verrify_df = verify_simulated_data(
                contigs_dict, spacers_dict, ground_truth_df
            )
            if (
                verrify_df.filter(
                    pl.col("alignment_test") > pl.col("mismatches")
                ).height
                > 0
            ):
                print(
                    verrify_df.filter(pl.col("alignment_test") > pl.col("mismatches"))
                )
                raise ValueError("Simulated data is not correct")
            else:
                print("Simulated data is correct")

        return contigs_dict, spacers_dict, ground_truth_df

    if threads is None:
        import multiprocessing

        threads = multiprocessing.cpu_count()

    if debug:
        print(f"Using {threads} threads for Rust simulation...")
        start_time = time.time()

    # Generate simulation ID if not provided
    if id_prefix is None:
        sim_params = {
            "contig_length_range": contig_length_range,
            "spacer_length_range": spacer_length_range,
            "n_mismatch_range": n_mismatch_range,
            "sample_size_contigs": sample_size_contigs,
            "sample_size_spacers": sample_size_spacers,
            "insertion_range": insertion_range,
            "n_insertion_range": n_insertion_range,
            "n_deletion_range": n_deletion_range,
            "prop_rc": prop_rc,
        }
        id_prefix = generate_simulation_id(sim_params)
        print(f"Using generated simulation ID: {id_prefix}")
    else:
        print(f"Using provided ID prefix: {id_prefix}")

    simulator = Simulator()

    # Handle base compositions for contigs and spacers separately
    contig_comp_obj = None
    spacer_comp_obj = None

    # Contig composition (use new params if provided, fallback to legacy params)
    if contig_gc_content is not None:
        contig_comp_obj = BaseComposition.from_gc_content(contig_gc_content)
    elif all(
        x is not None
        for x in [contig_a_frac, contig_t_frac, contig_c_frac, contig_g_frac]
    ):
        contig_comp_obj = BaseComposition.from_fractions(
            contig_a_frac, contig_t_frac, contig_c_frac, contig_g_frac
        )
    elif gc_content is not None:  # Fallback to legacy param
        contig_comp_obj = BaseComposition.from_gc_content(gc_content)
    elif all(
        x is not None for x in [a_frac, t_frac, c_frac, g_frac]
    ):  # Fallback to legacy param
        contig_comp_obj = BaseComposition.from_fractions(a_frac, t_frac, c_frac, g_frac)
    elif base_composition is not None:  # Fallback to legacy param
        contig_comp_obj = base_composition

    # Spacer composition (use new params if provided, fallback to legacy params)
    if spacer_gc_content is not None:
        spacer_comp_obj = BaseComposition.from_gc_content(spacer_gc_content)
    elif all(
        x is not None
        for x in [spacer_a_frac, spacer_t_frac, spacer_c_frac, spacer_g_frac]
    ):
        spacer_comp_obj = BaseComposition.from_fractions(
            spacer_a_frac, spacer_t_frac, spacer_c_frac, spacer_g_frac
        )
    elif gc_content is not None:  # Fallback to legacy param
        spacer_comp_obj = BaseComposition.from_gc_content(gc_content)
    elif all(
        x is not None for x in [a_frac, t_frac, c_frac, g_frac]
    ):  # Fallback to legacy param
        spacer_comp_obj = BaseComposition.from_fractions(a_frac, t_frac, c_frac, g_frac)
    elif base_composition is not None:  # Fallback to legacy param
        spacer_comp_obj = base_composition

    # Handle distribution types
    contig_dist_obj = None
    if contig_distribution is not None:
        contig_dist_obj = DistributionType(contig_distribution)

    spacer_dist_obj = None
    if spacer_distribution is not None:
        spacer_dist_obj = DistributionType(spacer_distribution)

    contigs, spacers, ground_truth = simulator.simulate_data(
        tuple(contig_length_range),
        tuple(spacer_length_range),
        tuple(n_mismatch_range),
        sample_size_contigs,
        sample_size_spacers,
        tuple(insertion_range),
        tuple(n_insertion_range),
        tuple(n_deletion_range),
        prop_rc,
        threads,
        debug,
        results_dir,
        id_prefix,
        contig_dist_obj,
        spacer_dist_obj,
        contig_comp_obj,
        spacer_comp_obj,
        data_subdir,
    )

    if debug:
        end_time = time.time()
        print(f"Simulation time when using Rust: {end_time - start_time} seconds")

    # Convert the ground truth data types before creating DataFrame
    processed_ground_truth = [
        [
            row[0],  # spacer_id (already has prefix from Rust)
            row[1],  # contig_id (already has prefix from Rust)
            int(row[2]),  # start (int)
            int(row[3]),  # end (int)
            row[4].lower() == "true",  # strand (bool)
            int(row[5]),  # mismatches (int)
        ]
        for row in ground_truth
    ]
    # breakpoint()

    ground_truth_df = pl.DataFrame(
        processed_ground_truth,
        schema={
            "spacer_id": pl.Utf8,
            "contig_id": pl.Utf8,
            "start": pl.UInt32,
            "end": pl.UInt32,
            "strand": pl.Boolean,
            "mismatches": pl.UInt32,
        },
        orient="row",
    )
    if verify:
        verrify_df = verify_simulated_data(contigs, spacers, ground_truth_df)
        if (
            verrify_df.filter(pl.col("alignment_test") > pl.col("mismatches")).height
            > 0
        ):
            print(verrify_df.filter(pl.col("alignment_test") > pl.col("mismatches")))
            raise ValueError("Simulated data is not correct")
        else:
            print("Simulated data is correct")
    return contigs, spacers, ground_truth_df


def verify_simulated_data(
    contigs: Union[dict[str, str], str],
    spacers: Union[dict[str, str], str],
    ground_truth: Union[pl.DataFrame, str],
    return_fraction: bool = False,
    return_bool: bool = False,
    return_positive_only: bool = False,
) -> Any:
    """Verify the simulated data by checking the ground truth against the contigs and spacers.

    Args:
        contigs(dict[str, str] | str): Contig sequences or path to file
        spacers(dict[str, str] | str): Spacer sequences or path to file
        ground_truth(pl.DataFrame | str): Ground truth dataframe or path to TSV
        return_fraction(bool): Return fraction instead of count
        return_bool(bool): Return boolean result
        return_positive_only(bool): Return only positive results

    Returns:
        Verification results
    """
    if not isinstance(ground_truth, pl.DataFrame):
        if isinstance(ground_truth, str):
            ground_truth = pl.read_csv(ground_truth, separator="\t")
        else:
            raise ValueError("ground_truth must be a pandas DataFrame or a string")
    # load seqs
    if isinstance(contigs, str):
        ground_truth = populate_pldf_withseqs_needletail(
            ground_truth,
            seqfile=contigs,
            seqcol="contig_seq",
            idcol="contig_id",
            trim_to_region=False,
            reverse_by_strand_col=False,
        )
    else:
        contig_df = pl.DataFrame(
            data={
                "contig_id": list(contigs.keys()),
                "contig_seq": list(contigs.values()),
            }
        )
        ground_truth = ground_truth.join(contig_df, on="contig_id", how="left")
    if isinstance(spacers, str):
        ground_truth = populate_pldf_withseqs_needletail(
            ground_truth,
            seqfile=spacers,
            seqcol="spacer_seq",
            idcol="spacer_id",
            trim_to_region=False,
            reverse_by_strand_col=False,
        )
    else:
        spacer_df = pl.DataFrame(
            data={
                "spacer_id": list(spacers.keys()),
                "spacer_seq": list(spacers.values()),
            }
        )
        ground_truth = ground_truth.join(spacer_df, on="spacer_id", how="left")

    # check if the ground truth is correct
    verify_df = test_alignment_polars(ground_truth, ignore_region_strands=False)
    if return_bool:
        return (
            verify_df.filter(pl.col("alignment_test") != pl.col("mismatches")).height
            == ground_truth.height
        )
    if return_fraction:
        error_df = verify_df.filter(pl.col("alignment_test") != pl.col("mismatches"))
        return error_df.height / ground_truth.height
    elif return_positive_only:
        return verify_df.filter(pl.col("alignment_test") == pl.col("mismatches"))
    return verify_df


def write_fasta(sequences: dict[str, str], filename: str) -> None:
    """Write sequences to a FASTA file.

    Args:
        sequences(dict[str, str]): Dictionary of sequence_id -> sequence
        filename(str): Output filename
    """
    with open(filename, "w") as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n{seq}\n")


def read_fasta(filename: str) -> dict[str, str]:
    """Read sequences from a FASTA file.

    Args:
        filename(str): Path to the FASTA file

    Returns:
        dict[str, str]: Dictionary of sequence_id -> sequence
    """
    sequences = {}
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip()
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line.strip()
    return sequences


def run_tool(tool: dict, results_dir: str, debug: bool = False) -> dict:
    """
    Run a tool either via its bash script (with hyperfine) or directly (debug mode).

    Args:
        tool(dict): Tool configuration dictionary
        results_dir(str): Results directory
        debug(bool): If True, run command directly without hyperfine wrapper for better error messages

    Returns:
        dict: Timing data dictionary (empty dict in debug mode)
    """
    if debug:
        # Debug mode: run the actual command directly, not the hyperfine wrapper
        print(f"\n{'=' * 60}")
        print(f"DEBUG: Running {tool['name']} directly")
        print(f"{'=' * 60}")

        # Get the mamba environment if specified
        mamba_env = tool.get("mamba_env", None)

        # Build the command
        if mamba_env:
            # Activate mamba env and run command
            cmd = f'eval "$(micromamba shell hook --shell bash)" && micromamba activate {mamba_env} && {" ".join(tool["command"])}'
        else:
            cmd = " ".join(tool["command"])

        print(f"Command: {cmd}\n")

        # Run the command and capture output
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, executable="/bin/bash"
        )

        # Print stdout and stderr
        if result.stdout:
            print("STDOUT:")
            print(result.stdout)
        if result.stderr:
            print("STDERR:")
            print(result.stderr)

        # Check return code
        if result.returncode != 0:
            print(f"\n❌ Command failed with exit code {result.returncode}")
            raise RuntimeError(
                f"Tool {tool['name']} failed with exit code {result.returncode}"
            )
        else:
            print("\n✓ Command succeeded")

        print(f"{'=' * 60}\n")
        return {}  # No timing data in debug mode
    else:
        # Normal mode: run via bash script with hyperfine
        print(
            f"Running {tool['script_name']} in {results_dir}/bash_scripts/{tool['script_name']}"
        )
        subprocess.run(f"{results_dir}/bash_scripts/{tool['script_name']}", shell=True)

        with open(f"{results_dir}/raw_outputs/{tool['script_name']}.json", "r") as f:
            data = json.load(f)
        return data


def create_bash_script(tool: dict, results_dir: str, max_runs: int = 1, warmups: int = 0, hyperfine: bool = True) -> None:
    """Create a bash script for running a tool with optional hyperfine benchmarking.

    Args:
        tool(dict): Tool configuration dictionary
        results_dir(str): Results directory path
        max_runs(int): Maximum number of runs for hyperfine
        warmups(int): Number of warmup runs for hyperfine
        hyperfine(bool): Whether to use hyperfine for benchmarking
    """
    script_name = tool["script_name"]
    mamba_env = tool.get("mamba_env", None)
    if hyperfine:
        hyperfine_command = [
            "hyperfine",
            "--warmup",
            str(warmups),
            "--max-runs",
            str(max_runs),
            "--prepare",
            f'"rm -rf  {results_dir}/raw_outputs/tmp* || true"',
            "--output",
            f"{results_dir}/raw_outputs/hyperfine_output_{tool['script_name']}.txt",
            "--export-json",
            f"{results_dir}/raw_outputs/{tool['script_name']}.json",
        ]

    with open(f"{results_dir}/bash_scripts/{tool['script_name']}", "w") as f:
        f.write("#!/bin/bash\n")
        # Always initialize mamba shell
        f.write('eval "$(micromamba shell hook --shell bash)"\n')
        if mamba_env:
            f.write(f"micromamba activate {mamba_env}\n")
        if hyperfine:
            f.write(f"{' '.join(hyperfine_command)} ")
            f.write(f"'{' '.join(tool['command'])}'")
        else:
            f.write(f"{' '.join(tool['command'])}\n")

        os.chmod(f"{results_dir}/bash_scripts/{script_name}", 0o755)


def clean_everything(results_dir: str) -> None:
    """Clean all output directories in the results directory.

    Args:
        results_dir(str): Path to the results directory
    """
    os.system(f"rm -rf {results_dir}/raw_outputs/")
    os.system(f"rm -rf {results_dir}/bash_scripts/")
    os.system(f"rm -rf {results_dir}/results/")
    os.system(f"rm -rf {results_dir}/simulated_data/")

def clean_before_rerun(tool_name: str, results_dir: str) -> None:
    """Clean output files for a specific tool before rerunning.

    Args:
        tool_name(str): Name of the tool
        results_dir(str): Path to the results directory
    """
    os.system(f"rm -rf {results_dir}/raw_outputs/{tool_name}*tsv")
    os.system(f"rm -rf {results_dir}/raw_outputs/{tool_name}*sam")
    # os.system(f"rm -rf {results_dir}/bash_scripts/{tool_name}*") # generating scripts moved to generate_scripts.py
    os.system(f"rm -rf {results_dir}/raw_outputs/tmp*")


def get_aln_len_from_cigar(cigar: str) -> int:
    """Calculate reference span from CIGAR string.
    Only M, D, N, =, X consume reference positions.
    I, S, H do not consume reference positions.

    Args:
        cigar(str): CIGAR string

    Returns:
        int: Length of alignment on reference
    """
    # Match operations that consume reference: M, D, N, =, X
    ref_consuming = re.findall(r"(\d+)[MDN=X]", cigar)
    return sum(int(num) for num in ref_consuming)


def fix_cigar_for_pysamtools(cigar: str) -> str:
    """Fix CIGAR string to ensure all operations have numeric lengths for pysamtools.

    Args:
        cigar(str): CIGAR string that may have missing lengths

    Returns:
        str: Fixed CIGAR string with all operations having lengths
    """
    # if any character is not preceded by a digit, prepend 1 to it
    # handle both standard CIGAR operators (A-Z) and extended operators (=,X)
    t = re.sub(r"([A-Z|=|X])(?!\d)([A-Z|=|X])", r"\g<1>1\g<2>", cigar, count=0)
    if re.search(r"[A-Z|=|X][A-Z|=|X]", t):  # Check for consecutive operators
        return fix_cigar_for_pysamtools(t)
    else:
        return t


def get_mismatches_from_cigar(cigar: str) -> int:
    """Extract the number of mismatches from CIGAR string (X operations).

    Args:
        cigar(str): CIGAR string

    Returns:
        int: Number of mismatches

    Note:
        Cigar string should be in extended format (version 1.4 cigar) with X for mismatches
    """
    return sum(int(num[:-1]) for num in re.findall(r"\d+X", cigar) or [0])


def parse_samVext(sam_file: str, max_mismatches: int = 5) -> pl.DataFrame:

    """Parse SAM file using the version 1.4 cigar strings, where = for matches, X for mismatches, I for insertions, D for deletions.

    Args:
        sam_file (str): Path to the SAM file.
        max_mismatches (int): Maximum mismatches allowed.

    Returns:
        pl.DataFrame: Parsed results.
    """
    results = []
    with open(sam_file, "r") as f:
        for line in f:
            # print(line)
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if fields[2] == "*":  # unmapped
                continue
            if len(fields) >= 3:
                spacer_id = str(fields[0])
                flag = int(fields[1])
                is_reverse = bool(flag & 16)
                strand = is_reverse  # Changed from '-' if is_reverse else '+'
                contig_id = str(fields[2])
                start = int(fields[3]) - 1  # SAM is 1-based
                cigar = fields[5]
                seq_length = get_aln_len_from_cigar(cigar)
                mismatches = -1  # default value
                if len(fields) >= 12:  # has tags
                    # get which tag is NM:i:``
                    nm_tag = [x for x in fields[11:] if x.startswith("NM:i:")]
                    if len(nm_tag) > 0:
                        mismatches = int(nm_tag[0].split(":")[2])

                if (
                    mismatches == -1
                ):  # if we don't find NM:i: we try to get it from the cigar string.
                    mismatches = get_mismatches_from_cigar(cigar)

                if mismatches > max_mismatches:
                    continue
                end = start + seq_length
                results.append(
                    [
                        str(spacer_id),
                        str(contig_id),
                        str(strand),
                        int(start),
                        int(end),
                        int(mismatches),
                    ]
                )

    if len(results) == 0:
        return pl.DataFrame(
            schema={
                "spacer_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "strand": pl.Boolean,
                "start": pl.UInt32,
                "end": pl.UInt32,
                "mismatches": pl.UInt32,
            },
            orient="row",
        )  # empty dataframe, so we can still join with other results

    # Create lists for each column
    spacer_ids, contig_ids, strands, starts, ends, mismatches = zip(*results)

    df = pl.DataFrame(
        {
            "spacer_id": spacer_ids,
            "contig_id": contig_ids,
            "strand": strands,
            "start": starts,
            "end": ends,
            "mismatches": mismatches,
        }
    )
    return df.unique()


def order_columns_to_match(df1_to_order: pl.DataFrame, df2_to_match: pl.DataFrame) -> pl.DataFrame:
    """Select columns in df1 to match the order in df2.

    Args:
        df1_to_order: Dataframe to reorder.
        df2_to_match: Reference dataframe.

    Returns:
        pl.DataFrame: Reordered dataframe.
    """
    return df1_to_order.select(df2_to_match.columns)


def cast_cols_to_match(df1_to_cast: pl.DataFrame, df2_to_match: pl.DataFrame) -> pl.DataFrame:
    """Cast columns in df1 to match types in df2.

    Args:
        df1_to_cast: Dataframe to cast.
        df2_to_match: Reference dataframe.

    Returns:
        pl.DataFrame: Casted dataframe.
    """
    for col in df2_to_match.columns:
        df1_to_cast = df1_to_cast.with_columns(
            pl.col(col).cast(df2_to_match.schema[col])
        )
    return df1_to_cast


def vstack_easy(df1_to_stack: pl.DataFrame, df2_to_stack: pl.DataFrame) -> pl.DataFrame:
    """Vstack two dataframes after matching column types and order.

    Args:
        df1_to_stack: Base dataframe.
        df2_to_stack: Dataframe to append.

    Returns:
        pl.DataFrame: Combined dataframe.
    """
    df2_to_stack = cast_cols_to_match(df2_to_stack, df1_to_stack)
    df2_to_stack = order_columns_to_match(df2_to_stack, df1_to_stack)
    return df1_to_stack.vstack(df2_to_stack)


def parse_sassy_duckdb(
    sassy_file: str,
    max_mismatches: int = 5,
    spacer_lendf: Optional[pl.DataFrame] = None,
    max_gaps: int = 1,
    threads: int = 10,
    output_prefix: str = "sassy_parsed",
    output_dir: str = ".",
    memory_limit: str = "50GB",
    **kwargs,
) -> pl.LazyFrame:
    """
    Parses Sassy TSV using DuckDB with caching capabilities.
    If a valid Parquet file already exists for the given prefix/date, it is loaded directly.

    Args:
        sassy_file (str): Input path.
        max_mismatches (int): Mismatch threshold.
        spacer_lendf (pl.DataFrame): Spacer lengths.
        max_gaps (int): Gap threshold.
        threads (int): CPU cores.
        output_prefix (str): Cache filename prefix.
        output_dir (str): Cache directory.
        memory_limit (str): DuckDB RAM limit.

    Returns:
        pl.LazyFrame: Parsed results.
    """
    from pathlib import Path
    import datetime

    output_dir = Path(sassy_file).parent

    # 1. CONSTRUCT DETERMINISTIC PATH
    # Naming convention: {output_dir}/{prefix}_{YYYYMMDD}.parquet
    date_str = datetime.datetime.now().strftime("%Y%m%d")
    parquet_filename = f"{output_prefix}_{date_str}.parquet"
    parquet_path = os.path.join(output_dir, parquet_filename)

    # 2. CACHE CHECK (RESUME LOGIC)
    if os.path.exists(parquet_path):
        print(f"Checking existing cache at: {parquet_path}")
        try:
            # We try to scan the file. If the footer is missing (incomplete write),
            # this will throw an error immediately.
            # Using scan_parquet().schema is a very fast metadata-only check.
            existing_schema = pl.scan_parquet(parquet_path).collect_schema()

            # (Optional) Verify essential columns exist to ensure it's not an old version
            if "spacer_id" in existing_schema and "mismatches" in existing_schema:
                print(
                    ">> Valid cache found! Returning lazy frame (no materialization)."
                )
                # Return lazy frame - let read_results/DuckDB handle streaming
                return pl.scan_parquet(parquet_path)
            else:
                print(">> Cache exists but schema is incorrect. Reprocessing...")
        except Exception as e:
            print(
                f">> Cache file exists but appears corrupt/incomplete ({e}). Reprocessing..."
            )

    # =========================================================================
    # START PROCESSING (Only if Cache Miss)
    # =========================================================================

    # STRICT thread enforcement via environment variables
    # These are set BEFORE DuckDB connection to take effect
    # os.environ["OMP_NUM_THREADS"] = str(threads)  # OpenMP
    # os.environ["OPENBLAS_NUM_THREADS"] = str(threads)  # OpenBLAS
    # os.environ["MKL_NUM_THREADS"] = str(threads)  # Intel MKL
    # os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)  # macOS veclib
    # os.environ["NUMEXPR_NUM_THREADS"] = str(threads)  # NumExpr

    print(
        f"\n[ThreadControl] Setting environment variables for {threads} threads (Sassy)"
    )

    con = duckdb.connect(database=":memory:")

    # Hardware Tuning - Control query parallelism
    con.execute(f"SET threads TO {threads};")
    con.execute("SET preserve_insertion_order = false;")
    con.execute(f"SET memory_limit = '{memory_limit}';")

    # Logic Definitions
    cigar_check = "length(cigar) - length(replace(replace(replace(cigar, 'I', ''), 'D', ''), 'N', ''))"

    if max_gaps == 1:
        gap_logic = f"{cigar_check} = 0"
    else:
        gap_logic = f"""
            (CASE 
                WHEN (cigar LIKE '%I%' OR cigar LIKE '%D%' OR cigar LIKE '%N%') 
                THEN {cigar_check} <= {max_gaps}
                ELSE TRUE 
            END)
        """

    if spacer_lendf is not None:
        con.register("spacer_lookup", spacer_lendf)
        join_step = "INNER JOIN spacer_lookup s ON g.spacer_id = s.spacer_id"
        len_selection = "s.length AS spacer_length"
    else:
        join_step = ""
        len_selection = "0 AS spacer_length"

    try:
        print(f"Starting DuckDB Pipeline (Writing to: {parquet_path})...")

        query = f"""
            COPY (
                WITH 
                -- STEP 1: SCAN & CHEAP FILTER
                fast_filter AS (
                    SELECT 
                        spacer_id, contig_id, mismatches, strand, "start", "end", cigar
                    FROM read_csv(
                        '{sassy_file}', 
                        sep='\t', 
                        header=True, 
                        skip=1,
                        parallel=true,
                        auto_detect=false,
                        hive_partitioning=false,
                        columns={{
                            'spacer_id': 'VARCHAR', 'contig_id': 'VARCHAR', 
                            'mismatches': 'UINTEGER', 'strand': 'VARCHAR', 
                            'start': 'UINTEGER', 'end': 'UINTEGER',
                            'slice_str': 'VARCHAR', 'cigar': 'VARCHAR'
                        }}
                    )
                    WHERE mismatches <= {max_mismatches}
                ),

                -- STEP 2: EXPENSIVE FILTER
                gap_filter AS (
                    SELECT * FROM fast_filter
                    WHERE {gap_logic}
                )

                -- STEP 3: JOIN & SELECT
                SELECT 
                    g.spacer_id,
                    g.contig_id,
                    {len_selection},
                    (g.strand = '-') AS strand,
                    g."start",
                    g."end",
                    g.mismatches
                FROM gap_filter g
                {join_step}

            -- Write Config
            ) TO '{parquet_path}' (FORMAT 'PARQUET', COMPRESSION 'SNAPPY', ROW_GROUP_SIZE 100000);
        """

        # Execute the write
        con.execute(query)

        # Return lazy frame without materialization
        print("Returning lazy frame (no materialization).")
        if os.path.exists(parquet_path):
            # Return lazy frame - let read_results/DuckDB handle streaming
            results = pl.scan_parquet(parquet_path)
        else:
            print(
                "Warning: DuckDB finished but no file was created (empty result set?)."
            )
            results = pl.DataFrame(
                schema={
                    "spacer_id": pl.String,
                    "contig_id": pl.String,
                    "spacer_length": pl.UInt32,
                    "strand": pl.Boolean,
                    "start": pl.UInt32,
                    "end": pl.UInt32,
                    "mismatches": pl.UInt32,
                },
                orient="row",
            ).lazy()

    except Exception as e:
        print(f"Failed to process Sassy file: {e}")
        # If the write crashed halfway, the file might be corrupt. Clean it up?
        # Typically better to leave it for inspection or manual deletion,
        # but you can uncomment this if you prefer auto-cleanup:
        # if os.path.exists(parquet_path): os.remove(parquet_path)

        return pl.DataFrame(
            schema={
                "spacer_id": pl.String,
                "contig_id": pl.String,
                "spacer_length": pl.UInt32,
                "strand": pl.Boolean,
                "start": pl.UInt32,
                "end": pl.UInt32,
                "mismatches": pl.UInt32,
            },
            orient="row",
        ).lazy()

    return results


def parse_sassy(
    sassy_file: str, max_mismatches: int = 5, spacer_lendf: Optional[pl.DataFrame] = None, max_gaps: int = 2, **kwargs
) -> pl.DataFrame:
    """Parse Sassy TSV output format and return standardized coordinates.
    Sassy output format:
    pat_id	text_id	cost	strand	start	end	match_region	cigar

    Args:
        sassy_file (str): Input path.
        max_mismatches (int): Maximum mismatches allowed.
        spacer_lendf (pl.DataFrame): Length info for spacers.
        max_gaps (int): Maximum gaps allowed in CIGAR.

    Returns:
        pl.DataFrame: (optionally filtered) results.
    """

    try:
        results = pl.scan_csv(
            sassy_file,
            separator="\t",
            skip_lines=1,
            schema={
                "spacer_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "mismatches": pl.UInt32,
                "strand": pl.Utf8,
                "start": pl.UInt32,
                "end": pl.UInt32,
                "slice_str": pl.Utf8,
                "cigar": pl.Utf8,
            },
        ).drop("slice_str")

        if max_gaps == 1:
            results = results.filter(
                ~pl.col("cigar").str.contains_any(patterns=["I", "D", "N"])
            )
        else:
            results = results.filter(
                ~pl.col("cigar").str.count_matches(r"I|D|N") > max_gaps
            )

        # Filter by max_mismatches
        results = results.filter(pl.col("mismatches") <= max_mismatches)

    except Exception as e:
        print(f"Failed to read Sassy file {sassy_file}: {e}, returning empty dataframe")
        return pl.DataFrame(
            schema={
                "spacer_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "spacer_length": pl.UInt32,
                "strand": pl.Boolean,
                "start": pl.UInt32,
                "end": pl.UInt32,
                "mismatches": pl.UInt32,
            },
            orient="row",
        )

    # Convert strand from string to boolean
    results = results.with_columns(
        pl.when(pl.col("strand") == "-")
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias("strand")
    )

    results = results.collect()

    # Join with spacer lengths if provided
    if spacer_lendf is not None:
        results = spacer_lendf.join(results, on="spacer_id", how="inner")
        results = results.rename({"length": "spacer_length"})
    else:
        # If no spacer_lendf provided, create a dummy spacer_length column
        results = results.with_columns(pl.lit(0).alias("spacer_length"))

    # Select and order columns to match standard format
    results = results.select(
        "spacer_id",
        "contig_id",
        "spacer_length",
        "strand",
        "start",
        "end",
        "mismatches",
    ).unique()

    return results


def parse_blastn(blastn_file: str, max_mismatches: int = 5, spacer_lendf: Optional[pl.DataFrame] = None, **kwargs) -> pl.DataFrame:
    """Parse a custom BLAST or mmseqs output format
    (query,target,nident,alnlen,mismatch,qlen,gapopen,qstart,qend,tstart,tend,evalue,bits)
    and return standardized coordinates (1-based).
    BLAST format uses 1-based coordinates.
    Filter out rows with more than max_mismatches mismatches. (i.e. retain up to (including) max_mismatches mismatches)

    Args:
        blastn_file (str): Input path.
        max_mismatches (int): Mismatch limit.
        spacer_lendf (pl.DataFrame): Length info.

    Returns:
        pl.DataFrame: Standardized results.
    """
    try:
        results = pl.scan_csv(
            blastn_file,
            separator="\t",
            has_header=False,
            schema={
                "spacer_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "nident": pl.UInt32,
                "alnlen": pl.UInt32,
                "mismatch": pl.UInt32,
                "qlen": pl.UInt32,
                "gapopen": pl.UInt32,
                "qstart": pl.UInt32,
                "qend": pl.UInt32,
                "tstart": pl.UInt32,
                "tend": pl.UInt32,
                "evalue": pl.Float64,
                "bits": pl.Float64,
            },
        )
        results = results.filter((pl.col("mismatch")  <= max_mismatches) ) # do here to push down the filter
        # results = results.filter((pl.col("qlen") - (pl.col("nident")+3) <= max_mismatches) ) # commented out, will be done later after adding gaps too.
        results = results.collect()
        # results = results.limit(1212121).collect()

        if any(
            results["qend"] < results["qstart"]
        ):  # mmseqs reports reverse on the subject, not the query.
            print("qend < qstart detected, assuming tend and tstart are reversed")
            results = results.with_columns(
                pl.when(pl.col("qend") < pl.col("qstart"))
                .then(pl.struct(tstart=pl.col("tend"), tend=pl.col("tstart")))
                .otherwise(pl.struct(tstart=pl.col("tstart"), tend=pl.col("tend")))
                .struct.field(["tstart", "tend"])
            )
        # results = results.with_columns(
        #      pl.col("spacer_id").cast(pl.Utf8)
        # )
    except Exception as e:
        print(
            f"Failed to read BLAST file {blastn_file}: {e}, returning empty dataframe"
        )
        return pl.DataFrame(
            schema={
                "spacer_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "spacer_length": pl.UInt32,
                "strand": pl.Utf8,
                "start": pl.UInt32,
                "end": pl.UInt32,
                "mismatches": pl.UInt32,
            },
            orient="row",
        )

    results = results.with_columns(pl.col("spacer_id").cast(pl.Utf8))
    results = results.with_columns(
        # Determine strand based on coordinate ordering
        pl.when(pl.col("tstart") <= pl.col("tend"))
        .then(pl.lit(False))  # Changed from '+' to False
        .otherwise(pl.lit(True))  # Changed from '-' to True
        .alias("strand"),
        pl.min_horizontal(["tstart", "tend"]).alias("start_min") - 1,
        pl.max_horizontal(["tstart", "tend"]).alias("end_max"),
        # Calculate mismatches as edit distance (non-identical positions + gaps)
        (pl.col("qlen") - pl.col("nident")).alias(
            "mismatches"
        ),  # Need to consider setting this as max(qlen,alnlen), as mmseqs can report seems to allow gaps in the query to be counted in the aligned length.
    )
    # print(results.filter(pl.col("mismatches") != pl.col("mismatch")).height)

    # Filter out rows with more than max_mismatches mismatches
    results = results.with_columns(
        pl.col("mismatches") + pl.col("gapopen").alias("mismatches")
    )
    results = results.filter(max_mismatches >= pl.col("mismatches"))
    results = spacer_lendf.join(results, on="spacer_id", how="inner")
    results = results.rename({"length": "spacer_length"})

    results = results.select(
        "spacer_id",
        "contig_id",
        "spacer_length",
        "strand",
        pl.col("start_min").alias("start"),
        pl.col("end_max").alias("end"),
        "mismatches",
    ).unique()
    results = results.with_columns(pl.col("spacer_id").cast(pl.Utf8))
    return results


def parse_lexicmap(tsv_file: str, max_mismatches: int = 5, spacer_lendf: Optional[pl.DataFrame] = None, **kwargs) -> pl.DataFrame:
    """Parse LexicMap TSV output format and return standardized coordinates.

    LexicMap output columns (1-based positions):
    1.  query,    Query sequence ID.
    2.  qlen,     Query sequence length.
    3.  hits,     Number of subject genomes.
    4.  sgenome,  Subject genome ID.
    5.  sseqid,   Subject sequence ID.
    6.  qcovGnm,  Query coverage (percentage) per genome: $(aligned bases in the genome)/$qlen.
    7.  cls,      Nth HSP cluster in the genome. (just for improving readability)
                  It's useful to show if multiple adjacent HSPs are collinear.
    8.  hsp,      Nth HSP in the genome.         (just for improving readability)
    9.  qcovHSP   Query coverage (percentage) per HSP: $(aligned bases in a HSP)/$qlen.
    10. alenHSP,  Aligned length in the current HSP.
    11. pident,   Percentage of identical matches in the current HSP.
    12. gaps,     Gaps in the current HSP.
    13. qstart,   Start of alignment in query sequence.
    14. qend,     End of alignment in query sequence.
    15. sstart,   Start of alignment in subject sequence.
    16. send,     End of alignment in subject sequence.
    17. sstr,     Subject strand.
    18. slen,     Subject sequence length.
    19. evalue,   Expect value.
    20. bitscore, Bit score.
    21. cigar,    CIGAR string of the alignment.                      (optional with -a/--all)
    22. qseq,     Aligned part of query sequence.                      (optional with -a/--all)
    23. sseq,     Aligned part of subject sequence.                    (optional with -a/--all)
    24. align,    Alignment text ("|" and " ") between qseq and sseq. (optional with -a/--all)

    Args:
        tsv_file (str): Input path.
        max_mismatches (int): Mismatch threshold.
        spacer_lendf (pl.DataFrame): Length lookup.

    Returns:
        pl.DataFrame: Results in standard format.
    """
    try:
        results = pl.read_csv(
            tsv_file,
            separator="\t",
            has_header=True,
            infer_schema_length=100000,
        )
        results = results.rename(
            {
                "query": "spacer_id",
                "sseqid": "contig_id",
                "sstr": "strand",
            }
        )
        if results.height == 0:
            raise ValueError(f"No results found in {tsv_file}")
    except Exception as e:
        print(
            f"Failed to create read from file {tsv_file}: {e}, returning empty dataframe"
        )
        return pl.DataFrame(
            schema={
                "spacer_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "spacer_length": pl.UInt32,
                "strand": pl.Utf8,
                "start": pl.UInt32,
                "end": pl.UInt32,
                "mismatches": pl.UInt32,
            },
            orient="row",
        )

    # results = results.filter(pl.col("alenHSP") >= 17, pl.col("gaps") == 0) # we don't really need this here, but keeping commented to remember it was done previously. Potentially, with a max mismatch >3 this could have over restrict the max mismatch arguments values.
    results = results.with_columns(pl.col("spacer_id").cast(pl.Utf8))
    results = results.with_columns(
        pl.col("align").str.count_matches("|", literal=True).alias("matches")
    )

    results = spacer_lendf.join(results, on="spacer_id", how="inner")
    results = results.with_columns(
        (pl.col("length") - pl.col("matches")).alias("mismatches"),
        (pl.col("sstart") - 1).alias("start"),  # 1-based
        (pl.col("send")).alias("end"),  # 1-based
    )
    results = results.filter(pl.col("mismatches") <= max_mismatches)
    results = results.rename({"length": "spacer_length"})
    results = results.select(
        "spacer_id",
        "contig_id",
        "spacer_length",
        "strand",
        "start",
        "end",
        "mismatches",
    ).unique()
    results = results.with_columns(pl.col("spacer_id").cast(pl.Utf8))
    return results


def parse_samVn_with_lens_polars(sam_file: str, spacer_lendf: pl.DataFrame, max_mismatches: int = 5) -> pl.DataFrame:
    """Parse SAM file using spacer lengths for qlen and numerics in CIGAR for alignment length
    uses polars to speed up things

    Args:
        sam_file (str): Path to SAM.
        spacer_lendf (pl.DataFrame): Spacer lengths.
        max_mismatches (int): Filter threshold.

    Returns:
        pl.DataFrame: Parsed alignments.
    """
    # First read mandatory columns
    mandatory_sam_cols = [
        "spacer_id",
        "flag",
        "contig_id",
        "start",
        "MAPQ",
        "cigar",
        "rnext",
        "pnext",
        "tlen",
        "seq",
        "qual",
    ]

    try:
        # Read first 11 mandatory columns
        sam_df = pl.read_csv(
            sam_file,
            separator="\t",
            comment_prefix="@",
            has_header=False,
            new_columns=mandatory_sam_cols,
            columns=list(range(11)),
            truncate_ragged_lines=False,
            infer_schema_length=1000,
        )
        # add index column
        sam_df = sam_df.with_row_index(name="index")
        sam_df = sam_df.drop(["rnext", "pnext", "MAPQ", "qual"])
        # get the index of the unmapped lines
        unmmaped_lines_idx = sam_df.filter(pl.col("contig_id") == "*").get_column(
            "index"
        )
        if unmmaped_lines_idx.len() > 0:
            print(
                f"Warning: {unmmaped_lines_idx.len()} unmapped lines found in {sam_file}"
            )
            sam_df = sam_df.filter(pl.col("contig_id") != "*")

        # Now read optional tags as a single column
        tags_df = pl.read_csv(
            sam_file,
            separator="%",  # risky...
            comment_prefix="@",
            has_header=False,
            new_columns=["tags"],
            truncate_ragged_lines=False,
        )
        tags_df = tags_df.with_row_index(name="index")
        if unmmaped_lines_idx.len() > 0:
            tags_df = tags_df.filter(~pl.col("index").is_in(unmmaped_lines_idx))

        # tags_df = tags_df.with_columns(
        #      pl.col("tags").str.split( "\t").list.gather_every(n=1,offset=11).alias("tags")
        # )
        # tags_df = tags_df.with_columns(
        #      pl.col("tags").list.len().alias("tag_lengths")
        # )
        nm_tag = "nm:i:" if tags_df["tags"][1].count("NM:i:") == 0 else "NM:i:"
        tags_df = tags_df.with_columns(
            pl.col("tags").str.extract(f"{nm_tag}(\d+)").cast(pl.Int32).alias("nm")
        )
        # Join with main dataframe
        sam_df = sam_df.with_columns(tags_df["nm"].alias("nm"))
        sam_df = sam_df.drop(["index"])
    except Exception as e:
        print(f"Failed to read SAM file {sam_file}: {e}, returning empty dataframe")
        return pl.DataFrame(
            schema={
                "spacer_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "spacer_length": pl.UInt32,
                "strand": pl.Utf8,
                "start": pl.UInt32,
                "end": pl.UInt32,
                "mismatches": pl.UInt32,
            },
            orient="row",
        )

    # Continue with existing processing
    sam_df = sam_df.with_columns((pl.col("start") - 1).alias("start"))  # sam is 1-based
    sam_df = sam_df.with_columns(
        [
            pl.when(pl.col("flag") == 16)
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias("strand")
        ]
    )
    sam_df = sam_df.drop("flag")

    sam_df = sam_df.with_columns(
        [
            pl.col("cigar")
            .map_elements(
                lambda x: sum(int(num) for num in re.findall(r"(\d+)", x)),
                return_dtype=pl.UInt32,
            )
            .alias("aligned_length"),
            pl.col("cigar")
            .map_elements(
                lambda x: sum(int(num) for num in re.findall(r"(\d+)S", x)),
                return_dtype=pl.UInt32,
            )
            .alias("soft_clipped_length"),
        ]
    )

    # Join with spacer lengths
    sam_df = spacer_lendf.join(sam_df, on="spacer_id", how="inner")

    # Calculate mismatches based on length and matches
    sam_df = sam_df.with_columns(
        [
            (
                (pl.col("length") - pl.col("aligned_length"))
                + pl.col("soft_clipped_length")
                + pl.col("nm")
            ).alias("mismatches")
        ]
    )

    sam_df = sam_df.with_columns(
        (pl.col("start") + pl.col("aligned_length")).alias("end")
    )
    sam_df = sam_df.rename({"length": "spacer_length"})
    sam_df = sam_df.filter(pl.col("mismatches") <= max_mismatches)
    sam_df = sam_df.filter(pl.col("nm") <= max_mismatches)
    sam_df = sam_df.select(
        "spacer_id",
        "contig_id",
        "spacer_length",
        "strand",
        "start",
        "end",
        "mismatches",
    ).unique()

    return sam_df


def parse_sam(
    sam_file: str,
    spacer_lendf: pl.DataFrame,
    max_mismatches: int = 5,
    threads: int = 4,
    ref_file: Optional[str] = None,
    gaps_as_mismatches: bool = False,
    **kwargs,
) -> pl.DataFrame:
    """Parse SAM file using pysam and spacer lengths to compute mismatches
    RENAMED FROM parse_samVn_with_lens_pysam

    Args:
        sam_file (str): Input path.
        spacer_lendf (pl.DataFrame): Length lookup.
        max_mismatches (int): Filter limit.
        threads (int): Worker threads.
        ref_file (str): Reference FASTA.
        gaps_as_mismatches (bool): NM calculation logic.

    Returns:
        pl.DataFrame: Parsed results.
    """

    try:
        sam_file = verify_sam_file(sam_file, ref_file)
    except Exception as e:
        print(f"Failed to verify SAM file {sam_file}: {e}, returning empty dataframe")
        return pl.DataFrame(
            schema={
                "spacer_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "spacer_length": pl.UInt32,
                "strand": pl.Boolean,
                "start": pl.UInt32,
                "end": pl.UInt32,
                "mismatches": pl.UInt32,
            },
            orient="row",
        )

    results = []
    # Create a dictionary for quick spacer length lookups - TODO: maybe there's a faster native polars way to do this.
    spacer_lens = dict(zip(spacer_lendf["spacer_id"], spacer_lendf["length"]))

    with pysam.AlignmentFile(sam_file, "r", check_sq=False) as samfile:
        for read in tqdm(samfile, desc="Parsing SAM file"):
            if read.is_unmapped:
                continue

            # spacer_len = read.query_length # 213768.22it/s
            # read.cigarstring = fix_cigar_for_pysamtools(read.cigarstring) # if missing integars in the cigar string...
            spacer_len = spacer_lens.get(
                read.query_name
            )  #  199485it/s, with the "fix_cigar_for_pysamtools" ~99k.

            # spacer_len = read.infer_read_length() #
            # ## Get hard clipped length
            # hard_clipped = sum(length for op, length in read.cigartuples
            #                   if op == 5)  # H operation
            # # if hard_clipped > 0:
            # #      print(f"Hard clipped {hard_clipped} for {read.query_name}")
            # #      # break
            # if spacer_len == 0: # some alignments get 0 length, maybe because some tools do not output the sequence for secondary/ambiguous alignments. I suspect pysam get's the query length from the nchar of the seq.
            #      # print (f"spacer_len == 0 for {read.query_name}, extracting alignment info from cigar")
            # # Get spacer length from reference table instead of query sequence
            #      spacer_len = spacer_lens.get(read.query_name)
            # else:
            #      spacer_len = spacer_len + hard_clipped
            #      # print(f"{spacer_len}")
            #      # break
            if spacer_len is None:
                print(f"Warning: spacer {read.query_name} not found in reference table")
                break
            # Get alignment infto_string()
            # strand = '-' if read.is_reverse else '+'
            start = read.reference_start  # pysam uses 0-based coordinates (i.e. it knows the sam file is 1-based, and corrects for it)

            # Get alignment length excluding soft clips
            # query_alignment_length = read.query_alignment_length
            query_alignment_length = len(
                read.get_reference_positions(full_length=False)
            )
            # aligned_length = sum(length for op, length in read.cigartuples
            #                    if op in [0, 7, 8])  # M, =, X operations

            # # # Get soft clipped length
            soft_clipped = sum(
                length for op, length in read.cigartuples if op == 4
            )  # "S" operation

            # debug
            # if soft_clipped > 0:
                #  print(f"Soft clipped {soft_clipped} for {read.query_name}")
                 # break
            # Get edit distance from NM tag  (counts non-identical positions including gaps)
            nm = (read.get_tag("NM") if read.has_tag("NM") else 0) + soft_clipped

            # Calculate mismatches based on length and matches
            mismatches = (spacer_len - query_alignment_length) + nm
            if mismatches < 0:
                print("mismatches < 0")
                break
            if mismatches > max_mismatches or nm > max_mismatches:
                continue

            results.append(
                [
                    str(read.query_name),
                    str(read.reference_name),
                    int(spacer_len),
                    read.is_reverse,
                    int(start),
                    int(read.reference_end),  # read.reference_end
                    int(mismatches),
                ]
            )

    return pl.DataFrame(
        results,
        schema={
            "spacer_id": pl.Utf8,
            "contig_id": pl.Utf8,
            "spacer_length": pl.UInt32,
            "strand": pl.Boolean,
            "start": pl.UInt32,
            "end": pl.UInt32,
            "mismatches": pl.UInt32,
        },
        orient="row",
    ).unique()


def plot_matrix(
    matrix: pl.DataFrame,
    title: str,
    filename: str,
    hor_axis: str = "pairs in tool",
    vert_axis: str = "pairs not in tool",
    return_chart: bool = True,
) -> Any:
    """Create a heatmap matrix plot using Altair.

    Args:
        matrix (pl.DataFrame): Input matrix data.
        title (str): Chart title.
        filename (str): Output file path base.
        hor_axis (str): Horizontal axis title.
        vert_axis (str): Vertical axis title.
        return_chart (bool): If True, returns the altair object.

    Returns:
        Optional[alt.Chart]: The generated chart.
    """
    import altair as alt

    # Rename the tool1 column to 'in_tool'
    matrix = matrix.rename({"tool1": "in_tool"})

    # Calculate total counts for each tool to use for sorting
    tool_cols = [col for col in matrix.columns if col != "in_tool"]
    matrix = matrix.with_columns(pl.sum_horizontal(tool_cols).alias("total"))

    # Sort the matrix by total counts (descending)
    matrix = matrix.sort("total", descending=True)
    sorted_tools = matrix["in_tool"].to_list()

    # Drop the total column after sorting
    matrix = matrix.drop("total")

    # Reorder the columns to match the sorted rows
    ordered_cols = ["in_tool"] + sorted_tools
    matrix = matrix.select([col for col in ordered_cols if col in matrix.columns])

    # Convert the matrix to long form for Altair
    melted = matrix.melt(
        id_vars=["in_tool"],
        value_vars=[col for col in matrix.columns if col != "in_tool"],
        variable_name="not_in_tool",
        value_name="count",
    )

    # Create the heatmap with Altair
    chart = (
        alt.Chart(melted)
        .mark_rect()
        .encode(
            x=alt.X(
                "not_in_tool:N",
                title=vert_axis,
                sort=sorted_tools,
                axis=alt.Axis(labelAngle=-45, labelFontSize=14, titleFontSize=16),
            ),
            y=alt.Y(
                "in_tool:N",
                title=hor_axis,
                sort=sorted_tools,
                axis=alt.Axis(labelFontSize=14, titleFontSize=16),
            ),
            color=alt.Color(
                "count:Q",
                scale=alt.Scale(scheme="blues"),
                legend=alt.Legend(title="Count", labelFontSize=14, titleFontSize=16),
            ),
            tooltip=["in_tool", "not_in_tool", "count"],
        )
        .properties(title=title, width=600, height=600)
    )

    # Add text labels on the cells
    text = chart.mark_text(baseline="middle").encode(
        text=alt.Text("count:Q", format="d"),
        color=alt.condition(
            alt.datum.count > melted["count"].mean(),
            alt.value("white"),
            alt.value("black"),
        ),
    )

    # Combine the heatmap and text
    final_chart = chart + text

    # Save both HTML and SVG versions
    final_chart.save(filename + ".html")

    # Configure for SVG export with higher DPI
    svg_chart = final_chart.configure_view(
        strokeWidth=0.1,  # Remove border
    ).properties(
        width=600,  # Double size for higher resolution
        height=600,
    )

    # Save as SVG with higher DPI
    svg_chart.save(filename + ".svg", scale_factor=1.0)

    return final_chart if return_chart else None


def parse_tab(
    output_file: str, max_mismatches: int = 5, spacer_lendf: Optional[pl.DataFrame] = None, **kwargs
) -> pl.DataFrame:
    """expecting tab file with columns: spacer_id, contig_id, start, end, strand - because Antonio insists.

    Args:
        output_file (str): Input path.
        max_mismatches (int): Ignored here as mismatches are set to 0.
        spacer_lendf (pl.DataFrame): Length info.

    Returns:
        pl.DataFrame: Cleaned results.
    """
    try:
        results = pl.read_csv(
            output_file,
            separator="\t",
            has_header=False,
            new_columns=["contig_id", "spacer_id", "start", "end", "strand"],
            infer_schema_length=100000,
        )
    except Exception as e:
        print(f"Failed to read TAB file {output_file}: {e}, returning empty dataframe")
        return pl.DataFrame(
            schema={
                "spacer_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "spacer_length": pl.UInt32,
                "strand": pl.Utf8,
                "start": pl.UInt32,
                "end": pl.UInt32,
                "mismatches": pl.UInt32,
            },
            orient="row",
        )

    results = results.with_columns(
        pl.lit(0).alias(
            "mismatches"
        )  # we assume this function is only used for spacer_containment, so we don't have mismatches
    )
    results = results.with_columns(pl.col("spacer_id").cast(pl.Utf8))
    results = spacer_lendf.join(results, on="spacer_id", how="inner")
    results = results.rename({"length": "spacer_length"})
    results = results.with_columns(
        pl.when(pl.col("strand") == "-")
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias("strand")  # Convert string strand to boolean
    )
    return results.select(
        "spacer_id",
        "contig_id",
        "spacer_length",
        "strand",
        "start",
        "end",
        "mismatches",
    ).unique()


def parse_hyperfine_output(json_file: str) -> Optional[pl.DataFrame]:
    """Parse hyperfine JSON output into a Polars DataFrame.

    Args:
        json_file (str): Path to the hyperfine results file.

    Returns:
        Optional[pl.DataFrame]: Dataframe with timing info or None on failure.
    """
    try:
        with open(json_file, "r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"Failed to read hyperfine results for {json_file}: {e}")
        return None

    # Extract the first result (assuming one script-ran per file)
    result = data["results"][0]

    # Extract relevant timing information
    time_info = pl.DataFrame(
        {
            "mean_time": result["mean"],
            "stddev_time": result["stddev"]
            if result["stddev"] is not None
            else 0.0,  # need to debug this
            "median_time": result["median"],
            "user_time": result["user"],
            "system_time": result["system"],
            "min_time": result["min"],
            "max_time": result["max"],
            "command": result["command"],
            "tool": json_file.split("/")[-1].replace(".sh.json", ""),
        }
    )
    # breakpoint()
    return time_info


def read_hyperfine_results(tools: dict, results_dir: str) -> pl.DataFrame:
    """Read and aggregate hyperfine results for all tools.

    Args:
        tools (dict): Tool configurations.
        results_dir (str): Working directory.

    Returns:
        pl.DataFrame: Aggregated timing metrics.
    """
    results = pl.DataFrame(
        schema={
            "mean_time": pl.Float64,
            "stddev_time": pl.Float64,
            "median_time": pl.Float64,
            "user_time": pl.Float64,
            "system_time": pl.Float64,
            "min_time": pl.Float64,
            "max_time": pl.Float64,
            "command": pl.Utf8,
            "tool": pl.Utf8,
        },
        orient="row",
    )
    for tool in tools.values():
        try:
            tool_results = parse_hyperfine_output(
                f"{results_dir}/raw_outputs/{tool['script_name']}.json"
            )
            results = results.vstack(tool_results)
        except Exception as e:
            print(f"Failed to read hyperfine results for {tool['name']}: {e}")
    return results


def run_tools(tools: dict, results_dir: str, debug: bool = False) -> None:
    """
    Run multiple tools.

    Args:
        tools(dict): Dictionary of tool configurations
        results_dir(str): Results directory
        debug(bool): If True, run commands directly without hyperfine for better error messages
    """
    for tool in tools.values():
        try:
            clean_before_rerun(tool["name"], results_dir)
            tool["time_info"] = run_tool(tool, results_dir, debug=debug)
            if not debug:
                print(f"\nSuccessfully ran tool: {tool['name']}")
        except Exception as e:
            print(f"{tool['name']} failed: {e}")


def get_seq_from_fastx(
    seqfile: str,
    seq_ids: Union[str, list[str], pl.DataFrame],
    return_dict: bool = False,
    return_df: bool = True,
    idcol: str = "contig_id",
    seqcol: str = "contig_seq",
    output_file: Optional[str] = None,
) -> Any:
    """Extract specific sequences from a FASTX file using pyfastx or seqkit.

    Args:
        seqfile (str): Input path.
        seq_ids: Identifiers to extract.
        return_dict (bool): Return as mapping.
        return_df (bool): Return as Polars DF.
        idcol (str): ID column name.
        seqcol (str): Sequence column name.
        output_file (str): If provided, writes to file instead of returning objects.

    Returns:
        Any: List, Dict, DataFrame, or file path.
    """
    if isinstance(seq_ids, str):
        seq_ids = [seq_ids]
    if isinstance(seq_ids, pl.DataFrame):
        seq_ids = seq_ids["contig_id"].to_list()
    # make a tmp file in memory with the seq_ids
    with tempfile.NamedTemporaryFile(delete=False, mode="w") as tmp_file:
        tmp_file.write("\n".join(seq_ids))
        tmp_file_path = tmp_file.name
    if output_file is None:
        seqs = (
            subprocess.run(
                f"pyfastx extract {seqfile} -l '{tmp_file_path}' | seqkit seq -w 0",
                shell=True,
                capture_output=True,
            )
            .stdout.decode("utf-8")
            .replace(">", "")
            .splitlines()
        )
        tmp_file.close()
        os.remove(tmp_file_path)
        if not return_df and not return_dict:
            return seqs
        found_ids = [x.removeprefix(">") for x in seqs[::2]]
        seq_dicts = dict(zip(found_ids, seqs[1::2]))
        if return_dict:
            return seq_dicts
        elif return_df:
            return pl.DataFrame({idcol: found_ids, seqcol: seqs[1::2]})
    else:
        import pyfastx as pfx

        with open(output_file, "w") as f:
            fa = pfx.Fasta(seqfile)
            # load all sequence names into a set object
            all_names = set(fa.keys())
            for seqname in seq_ids:
                if seqname in all_names:
                    seq = fa[seqname].seq
                    f.write(f">{seqname}\n{seq}\n")
        return output_file


def soft_fetch_fastx(fa: pfx.Fasta, id: str, start: Optional[int] = None, end: Optional[int] = None, strand: Optional[str] = None) -> Optional[str]:
    """Safely fetch a sequence or sub-sequence from a pyfastx Fasta object.

    Args:
        fa: Pyfastx Fasta object.
        id (str): Sequence ID.
        start (int): Start coordinate.
        end (int): End coordinate.
        strand (str): '+' or '-'.

    Returns:
        Optional[str]: Sequence string or None if not found.
    """
    if start is not None and end is not None and strand is not None:
        try:
            return fa.fetch(id, (start, end), strand)
        except Exception:
            return None
    else:
        try:
            return fa[id].seq
        except Exception:
            return None  # if the id is not found, return None


def populate_pldf_from_fastx(
    seqfile: str,
    pldf: pl.DataFrame,
    trim_to_region: bool = False,
    idcol: str = "contig_id",
    start_col: str = "start",
    end_col: str = "end",
    strand_col: str = "strand",
    seqcol: str = "contig_seq",
    drop_missing: bool = True,
    **kwargs,
) -> pl.DataFrame:
    """Get sequences from a fasta file using pyfastx.
    sequences will be extracted from the start-1 to the end position.
    If strand_col is provided, the extracted sequences will be reverse complemented if the strand is boolean True.
    For some reason, this is slower than the populate_pldf_withseqs function.
    """
    # og_cols = pldf.columns
    fa = pfx.Fasta(seqfile)
    cols = start_col is not None and end_col is not None and strand_col is not None
    if (not cols) or (trim_to_region):
        print(
            "No start, end, or strand column provided, or trim_to_region is True, returning the entire sequences"
        )
        tmpdf = pldf.select([idcol]).unique()
        tmpdf = tmpdf.with_columns(
            pl.lit(None).alias(start_col),
            pl.lit(None).alias(end_col),
            pl.lit(None).alias(strand_col),
        )
        tmpdf = tmpdf.with_columns(
            pl.col(idcol).map_elements(
                lambda x: soft_fetch_fastx(fa, x, None, None, None)
            )
        )
        return tmpdf
    # temporary dataframe with the seq_ids, start, end, and strand taht we will export as a bed file (0 based coordinates, so we need to subtract 1 from the start and end)
    tmpdf = pldf.select([idcol, start_col, end_col, strand_col]).unique()
    tmpdf = tmpdf.with_columns(
        pl.when(pl.col(strand_col))
        .then(pl.lit("-"))
        .otherwise(pl.lit("+"))
        .alias("strand4bed"),
        (pl.col(start_col) - 1).alias("start_0based"),
    )
    # the following is not the best way to fetch sequences... Need to think of a better way to do this...
    seqs = [
        None
    ] * tmpdf.height  # initialize the list of sequences the same size as the tmpdf
    for ix in tqdm(
        range(tmpdf.height), desc="Fetching sequences", total=tmpdf.height
    ):  # 9998/9998 [00:24<00:00, 410.95it/s, 00:15<00:00, 645.96it/s when list is initialized with None
        row = tmpdf[ix]
        seq = soft_fetch_fastx(
            fa,
            row[idcol].item(),
            row["start_0based"].item(),
            row[end_col].item(),
            row["strand4bed"].item(),
        )
        seqs[ix] = seq
    tmpdf = tmpdf.with_columns(pl.Series(seqs).alias(seqcol))
    # # benchmarked this against the following:
    # time_start = time.time()
    # tmpdf = tmpdf.with_columns(
    # pl.struct(pl.col(idcol), pl.col("start_0based"), pl.col(end_col), pl.col("strand4bed")).map_elements(lambda x: soft_fetch_fastx(fa, x[idcol], x["start_0based"], x[end_col], x["strand4bed"]),strategy="threading", return_dtype=pl.Utf8).alias(seqcol)
    #      )
    # time_end = time.time() # 17.89 seconds (no threading) 16.66 seconds (threading)
    # print(f"Time taken: {time_end - time_start:.2f} seconds")

    tmpdf = tmpdf.drop(["strand4bed", "start_0based"])

    if drop_missing:
        tmpdf = tmpdf.filter(~pl.col(seqcol).is_in([None, "", "nan"]))
    pldf = pldf.join(tmpdf, on=[idcol, start_col, end_col, strand_col], how="left")
    return pldf


def populate_pldf_withseqs(
    pldf: pl.DataFrame,
    seqfile: str,
    chunk_size: int = 2000,
    drop_missing: bool = True,
    trim_to_region: bool = False,
    reverse_by_strand_col: bool = False,
    idcol: str = "contig_id",
    seqcol: str = "contig_seq",
    strand_col: str = "strand",
    **kwargs,
) -> pl.DataFrame:
    """Populate a Polars DataFrame with sequences from a FASTX file using chunked extraction.

    Args:
        pldf: Target dataframe.
        seqfile: Sequence file path.
        chunk_size: Processing chunk size.
        drop_missing: Drop rows where sequence isn't found.
        trim_to_region: Trim sequence to coordinates.
        reverse_by_strand_col: Reverse complement if strand is True.
        idcol: ID column name.
        seqcol: Output sequence column name.
        strand_col: Strand column name.

    Returns:
        pl.DataFrame: Populated dataframe.
    """
    # Note! if a "strand" column is present, and reverse_by_strand_col is True, it will be used to get the reverse complement of the sequence.
    # If drop_missing is True, it will drop any contig_ids that are not present in the sequence file. If no seqids were found, it will return the input dataframe.
    # NOTE: trim_to_region and reverse_by_strand_col are bugged when using slice, even though it's probably much faster.
    # get the unique contig_ids and strands
    stranded = strand_col in pldf.columns and reverse_by_strand_col
    # if stranded:
    #      nr_contigids = pldf[[idcol,strand_col]].unique()
    #      # print("rev true")
    # else:
    nr_contigids = pl.DataFrame({idcol: pldf[idcol].unique()})
    # print("rev false")
    # drop nulls nas or empty strings
    nr_contigids = nr_contigids.filter(~pl.col(idcol).is_in([None, "", "nan"]))
    print(f"Total of unique {idcol} ids to fetch: {nr_contigids.height}")
    nr_contigs_list = nr_contigids[idcol].unique().to_list()
    minipldf = pl.DataFrame(
        {idcol: nr_contigs_list, seqcol: [None] * len(nr_contigs_list)}
    )
    # Process chunks and update DataFrame
    for i in tqdm(
        range(0, len(nr_contigs_list), chunk_size),
        desc="Populating df with sequences",
        total=len(nr_contigs_list) // chunk_size,
    ):
        chunk = nr_contigs_list[i : i + chunk_size]
        # Get sequences for current chunk
        chunk_seqs = get_seq_from_fastx(
            seqfile=seqfile, seq_ids=chunk, return_df=True, idcol=idcol, seqcol=seqcol
        )

        # Update minipldf with new sequences using join
        minipldf = (
            minipldf.join(chunk_seqs, on=idcol, how="left")
            .with_columns(
                # Coalesce keeps existing values where not null, uses new values only where null
                pl.coalesce([pl.col(f"{seqcol}"), pl.col(f"{seqcol}_right")]).alias(
                    seqcol
                )
            )
            .drop(f"{seqcol}_right")
        )

    print(
        f"Estimated size of temporary dataframe: {minipldf.estimated_size('mb'):.2f} MB"
    )
    # count the number of nulls in the seqs_list
    nulls = [x for x in minipldf[seqcol].to_list() if x is None]
    print(f"Number of unfound ids: {len(nulls)}")
    if minipldf.height == 0:
        print("No sequences found, returning input as is")
        return pldf

    if drop_missing:
        print("Dropping unfound ids")
        pldf = pldf.filter(pl.col(idcol).is_in(minipldf[idcol].unique().to_list()))

    pldf = pldf.join(minipldf, on=idcol, how="left")

    if trim_to_region and "start" in pldf.columns and "end" in pldf.columns:
        print("Trimming to region")
        # Benchmarking with slice
        # time_start = time.time()
        # pldf = pldf.with_columns(
        # pl.col(seqcol).str.slice(length= pl.col("end")-pl.col("start"),offset=pl.col("start")-1).alias(seqcol)
        # )
        # time_end = time.time()
        # print(f"Time taken: {time_end - time_start:.2f} seconds") # Time taken: 0.15 seconds ### UPDATE: THIS IS BUGGED, REVERTING TO map_elements.
        # Benchmarking with map_elements
        # time_start = time.time()
        pldf = pldf.with_columns(
            pl.struct(pl.col(seqcol), pl.col("start"), pl.col("end"))
            .map_elements(
                lambda x: str(x[seqcol][x["start"] : x["end"]]), return_dtype=pl.Utf8
            )
            .alias(seqcol)
        )
        # time_end = time.time()
        # print(f"Time taken: {time_end - time_start:.2f} seconds") # Time taken: 0.06 seconds

    if stranded:
        print("Reversing stranded sequences")
        # get the strand from the nr_contigids
        pldf = pldf.with_columns(
            pl.when(pl.col("strand"))
            .then(pl.col(seqcol).map_elements(reverse_complement, return_dtype=pl.Utf8))
            .otherwise(pl.col(seqcol))
            .alias(seqcol)
        )

    return pldf


def populate_pldf_withseqs_needletail(
    pldf: pl.DataFrame,
    seqfile: str,
    chunk_size: int = 20000000,
    trim_to_region: bool = True,
    reverse_by_strand_col: bool = True,
    idcol: str = "contig_id",
    seqcol: str = "contig_seq",
    start_col: str = "start",
    end_col: str = "end",
    strand_col: str = "strand",
) -> pl.DataFrame:
    """Populate Polars DF using needletail for fast parsing and processing.

    Args:
        pldf: Target dataframe.
        seqfile: Sequence file.
        chunk_size: Max records per processing batch.
        trim_to_region: Trim sequence.
        reverse_by_strand_col: Rev-comp if strand is True.
        idcol: ID col name.
        seqcol: Output sequence col.
        start_col: Start coordinate.
        end_col: End coordinate.
        strand_col: Strand boolean col.

    Returns:
        pl.DataFrame: Updated dataframe.
    """
    merge_cols = [idcol]
    if reverse_by_strand_col:
        merge_cols.append(strand_col)
    if trim_to_region:
        merge_cols.extend([start_col, end_col])

    print(f"Initial pldf shape: {pldf.shape}")
    minipldf = pldf.select(merge_cols).unique()
    print(f"Unique entries in minipldf: {minipldf.shape}")

    minipldf = minipldf.filter(~pl.col(idcol).is_in([None, "", "nan"]))
    print(f"After filtering nulls: {minipldf.shape}")

    minipldf = minipldf.with_columns(pl.lit(None).alias(seqcol))

    seqs = []
    seq_ids = []

    # Get actual sequence count from file
    if seqfile.endswith(".gz"):
        seq_count = int(
            subprocess.run(
                f"zgrep -c '>'  {seqfile} ", shell=True, capture_output=True, text=True
            ).stdout.strip()
        )
    else:
        seq_count = int(
            subprocess.run(
                f"grep -F '>'  {seqfile} -c ",
                shell=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        )
    # seq_count = 0
    # for _ in parse_fastx_file(seqfile):
    #      seq_count += 1
    print(f"Actual number of sequences in file: {seq_count}")

    # Reset file iterator
    index = 0
    for record in parse_fastx_file(seqfile):
        seqs.append(record.seq)
        seq_ids.append(record.id)
        index += 1

        # Process chunk when we hit chunk_size or end of file
        if len(seqs) >= chunk_size or index == seq_count:
            print(f"\nProcessing chunk {index}/{seq_count}")
            print(f"Number of sequences in chunk: {len(seqs)}")

            chunk_seqs = pl.DataFrame({idcol: seq_ids, seqcol: seqs})

            chunk_seqs = chunk_seqs.join(
                minipldf.select(merge_cols), on=idcol, how="inner"
            )  # this join get's the info columns (start, end, strand) if needed, only for the entires in this chunk that are in the minipldf.

            if trim_to_region:
                print("Trimming sequences")
                # print(chunk_seqs.columns)
                # chunk_seqs = chunk_seqs.with_columns(len_toextract=pl.col(end_col)-pl.col(start_col) + 1)
                chunk_seqs = chunk_seqs.with_columns(
                #     pl.col(seqcol).str.slice(
                #         offset=pl.col(start_col) -1,
                #         length="len_toextract"                        
                #     ).alias(seqcol)
                # )
                            # ,  # is polars 0-based coordinated?
                # print (chunk_seqs[seqcol].head(10))
                    pl.struct(pl.col(seqcol), pl.col(start_col), pl.col(end_col))
                    .map_elements(
                        lambda x: str(x[seqcol][x[start_col] : x[end_col]])
                        if x[seqcol] is not None
                        else None,
                        return_dtype=pl.Utf8,
                    )
                    .alias(seqcol)
                )

            if reverse_by_strand_col:
                print("Reversing sequences")
                # print(chunk_seqs.columns)
                chunk_seqs = chunk_seqs.with_columns(
                    pl.when(pl.col(strand_col) == "-")
                    .then(
                        pl.col(seqcol).map_elements(
                            lambda x: reverse_complement(x) if x is not None else None,
                            return_dtype=pl.Utf8,
                        )
                    )
                    .otherwise(pl.col(seqcol))
                    .alias(seqcol)
                )

            print("Joining with nascent df")
            minipldf = minipldf.join(chunk_seqs, on=merge_cols, how="left")
            minipldf = minipldf.with_columns(
                pl.coalesce([pl.col(seqcol), pl.col(f"{seqcol}_right")]).alias(seqcol)
            ).drop(f"{seqcol}_right")

            print(f"Null count in seqcol after chunk: {minipldf[seqcol].null_count()}")

            seqs = []
            seq_ids = []
            # get count for remaining nulls, if zero, break - should be useful when fetching just a few sequences from a large file, at least if the needed seqs are closer to the start of the input fasta.
            if minipldf[seqcol].null_count() == 0:
                break

    print("\nFinal merge with original df")
    pldf = pldf.join(minipldf, on=merge_cols, how="left")
    print(f"Final null count in seqcol: {pldf[seqcol].null_count()}")

    return pldf


def test_alignment_from_faidx(spacers_file: str, contigs_file: str, alignment: pl.DataFrame, **kwargs) -> bool:
    """Test if an alignment matches between spacer and contig files using pyfastx.
    alignment is a single row extracted from a pl.df results  (spacer_id, contig_id, spacer_length, start, end, strand, mismatches, tool)

    Args:
        spacers_file: Path to spacer FASTA.
        contigs_file: Path to contig FASTA.
        alignment: Single-row dataframe.

    Returns:
        bool: True if alignment verification passes.
    """

    contig_seq = get_seq_from_fastx(contigs_file, alignment["contig_id"][0])
    # Get the sequences
    spacer_seq = get_seq_from_fastx(spacers_file, alignment["spacer_id"][0])
    # Extract the aligned region from contig
    aligned_region = contig_seq[alignment["start"].item() : alignment["end"].item()]

    # Compare sequences accounting for strand
    if alignment["strand"]:  # True means reverse strand
        aligned_region = reverse_complement(aligned_region)

    alignment = ps.nw_stats_scan(
        spacer_seq, aligned_region, open=10, extend=10, matrix=ps.nuc44
    )
    mismatches = len(spacer_seq) - (alignment.matches)
    if mismatches != alignment["mismatches"]:
        # print(f"Mismatch between pyfastx and alignment: {mismatches} != {alignment['mismatches']}")
        return False
    return True


def prettify_alignment(
    spacer_seq: str,
    contig_seq: str,
    strand: bool = False,
    start: Optional[int] = None,
    end: Optional[int] = None,
    return_ref_region: bool = False,
    return_query_region: bool = False,
    return_comp_str: bool = False,
    gap_cost: int = 10,
    extend_cost: int = 5,
    cost_matrix: Any = ps.nuc44,
    **kwargs,
) -> str:
    """given a spacer and contig sequence and coordinates, return | for a match and . for a mismatch.

    Args:
        spacer_seq: Query sequence.
        contig_seq: Reference sequence.
        strand: Rev-comp if True.
        start: Start coord.
        end: End coord.
        return_ref_region: Return traceback ref.
        return_query_region: Return traceback query.
        return_comp_str: Return comparison string (| .).
        gap_cost: Penalty.
        extend_cost: Penalty.
        cost_matrix: Scoring matrix.

    Returns:
        str: Alignment string representation.
    """
    if start is not None and end is not None:
        aligned_region = contig_seq[start:end]
    else:
        aligned_region = contig_seq

    if strand:  # True means reverse strand
        aligned_region = reverse_complement(aligned_region)
    # alignment = ps.nw_stats_scan(spacer_seq,aligned_region,open=10,extend=10,matrix=ps.nuc44)
    # alignment_string = ""
    # for i in range(len(spacer_seq)):
    #      if spacer_seq[i] == aligned_region[i]:
    #          alignment_string += "|"
    #      else:
    #          alignment_string += "."
    # return f"{spacer_seq}\n{alignment_string}\n{aligned_region}"
    ali = ps.nw_trace(
        spacer_seq,
        aligned_region,
        open=gap_cost,
        extend=extend_cost,
        matrix=cost_matrix,
    ).traceback
    if return_ref_region:
        return ali.ref
    if return_query_region:
        return ali.query
    if return_comp_str:
        return ali.comp
    return f"{ali.query}\n{ali.comp}\n{ali.ref}"

def test_alignment(
    spacer_seq: str,
    contig_seq: str,
    strand: bool = False,
    start: Optional[int] = None,
    end: Optional[int] = None,
    gap_cost: int = 10,
    extend_cost: int = 5,
    gaps_as_mismatch: bool = True,
    **kwargs,
)-> int:
    """Test if an alignment matches between a spacer and contig using Needleman-Wunsch algorithm.
    Input is a spacer and contig sequence and coordinates, and optionally strand, start, end, gap_cost, extend_cost, and gaps_as_mismatch.
    Returns the number of mismatches.

    Args:
        spacer_seq (str): Query sequence.
        contig_seq (str): Subject sequence.
        strand (bool): Reverse complement flag.
        start (int): Start index.
        end (int): End index.
        gap_cost (int): Opening penalty.
        extend_cost (int): Extension penalty.
        gaps_as_mismatch(bool): If True, count gaps as mismatches (similar to edit distance, but not quite).
                         If False, count only substitutions (hamming-like distance).

    Returns:
        int: Mismatch count.
    """
    if start is not None and end is not None:
        aligned_region = contig_seq[start:end]
    else:
        aligned_region = contig_seq
    if strand:  # True means reverse strand
        aligned_region = reverse_complement(aligned_region)

    # Always use nw_trace to get the alignment string for accurate counting
    alignment = ps.nw_trace(
        spacer_seq,
        aligned_region,
        open=gap_cost,
        extend=extend_cost,
        matrix=ps.nuc44,
    )

    # In parasail's comp string:
    # '|' = match
    # '.' = mismatch (substitution)
    # ' ' (space) = gap/indel

    comp_string = alignment.traceback.comp
    matches = comp_string.count("|")

    if gaps_as_mismatch:
        # Edit distance: spacer length minus matched positions
        # This counts all differences including substitutions and indels
        mismatches = len(spacer_seq) - matches
    else:
        # Hamming distance: count only substitutions, exclude gaps
        # This counts mismatched positions but ignores indels
        mismatches = comp_string.count(".")

    return mismatches


def calculate_hamming_distance(
    spacer_seq: str, contig_seq: str, strand: bool = False, start=None, end=None
):
    """
    Calculate hamming distance (substitutions only, no indels) using ungapped alignment.

    Args:
        spacer_seq(str): Query sequence
        contig_seq(str): Target sequence (may be pre-trimmed to region)
        strand(bool): If True, reverse complement contig_seq
        start(int): Region start coordinate (only used for slicing if sequence is full-length)
        end(int): Region end coordinate (only used for slicing if sequence is full-length)

    Returns:
        Hamming distance (number of substitutions only)
    
    Note:
        If sequences are of different lengths, returns the length difference as an invalid alignment.
    """
    if start is not None and end is not None:
        # Check if sequence is already trimmed to region
        expected_length = end - start
        if len(contig_seq) == expected_length:
            # Already trimmed, use as-is
            aligned_region = contig_seq
        else:
            # Full sequence, need to trim
            aligned_region = contig_seq[start:end]
    else:
        aligned_region = contig_seq

    if strand:
        aligned_region = reverse_complement(aligned_region)

    # If lengths don't match, can't do proper hamming distance, trying with padding (the shorter) from each side as "good enough effort"
    if len(spacer_seq) != len(aligned_region):
        print("Sequences of different lengths, cannot compute hamming distance accurately.")
        if len(spacer_seq) < len(aligned_region):
            # Pad spacer_seq (RIGHT)
            spacer_seq_r = spacer_seq.ljust(len(aligned_region), "@")
            ham_rght = sum(1 for a, b in zip(spacer_seq_r, aligned_region) if a != b)
            # Pad spacer_seq (LEFT)
            spacer_seq_l = spacer_seq.rjust(len(aligned_region), "@")
            ham_lft = sum(1 for a, b in zip(spacer_seq_l, aligned_region) if a != b)
        else:
            # Pad aligned_region (RIGHT)
            aligned_region_r = aligned_region.ljust(len(spacer_seq), "@")
            ham_rght = sum(1 for a, b in zip(spacer_seq, aligned_region_r) if a != b)
            # Pad aligned_region (LEFT)
            aligned_region_l = aligned_region.rjust(len(spacer_seq), "@")
            ham_lft = sum(1 for a, b in zip(spacer_seq, aligned_region_l) if a != b)
            # Return the minimum of the two padding attempts
            return min(ham_rght, ham_lft)
    # if they're the same length, proceed normally
    # Simple character-by-character comparison (true hamming distance)
    hamming = sum(1 for a, b in zip(spacer_seq, aligned_region) if a != b)
    return hamming


def calculate_edit_distance(
    spacer_seq: str,
    contig_seq: str,
    strand: bool = False,
    start: int = None,
    end: Optional[int] = None,
    gap_open: int = 10,
    gap_extend: int = 5,
    matrix: ps.Matrix = ps.nuc44
):
    """
    Calculate edit distance (allowing indels) using Needleman-Wunsch gapped alignment.

    Args:
        spacer_seq(str): Query sequence
        contig_seq(str): Target sequence (may be pre-trimmed to region)
        strand(bool): If True, reverse complement contig_seq
        start(int): Region start coordinate (only used for slicing if sequence is full-length)
        end(int): Region end coordinate (only used for slicing if sequence is full-length)
        gap_open(int): Gap opening penalty
        gap_extend(int): Gap extension penalty

    Returns:
        True edit distance: count of substitutions + insertions + deletions
    """
    if start is not None and end is not None:
        # Check if sequence is already trimmed to region
        expected_length = end - start
        if len(contig_seq) == expected_length:
            # Already trimmed, use as-is
            aligned_region = contig_seq
        else:
            # Full sequence, need to trim
            aligned_region = contig_seq[start:end]
    else:
        aligned_region = contig_seq

    if strand:
        aligned_region = reverse_complement(aligned_region)

    # Use Needleman-Wunsch global alignment (allows indels)
    # NW ensures we align the full sequences, giving true edit distance
    alignment = ps.nw_trace(
        spacer_seq,
        aligned_region,
        open=gap_open,
        extend=gap_extend,
        matrix=matrix,
    )

    comp_string = alignment.traceback.comp

    # True edit distance: count all non-match operations
    # '|' = match (0 cost)
    # '.' = substitution (1 cost)
    # ' ' = gap/indel (1 cost per gap)
    substitutions = comp_string.count(".")
    indels = comp_string.count(" ")

    edit_distance = substitutions + indels
    return edit_distance


def test_row(row, ignore_region_strands=False):
    if ignore_region_strands:
        return test_alignment(
            row["spacer_seq"],
            row["contig_seq"],
            None,
            None,
            None,
            gaps_as_mismatch=True,
        )
    else:
        return test_alignment(
            row["spacer_seq"],
            row["contig_seq"],
            row["strand"],
            row["start"],
            row["end"],
            gaps_as_mismatch=True,
        )


def test_alignment_polars(results: pl.DataFrame, ignore_region_strands=False, **kwargs):
    """Test if an alignment matches between a spacer and contig.
    Input is a pl.df with columns: spacer_id, contig_id, spacer_length, start, end, strand, mismatches, tool, spacer_seq, contig_seq
    """
    if "spacer_seq" not in results.columns or "contig_seq" not in results.columns:
        raise ValueError("spacer_seq and contig_seq columns are required")

    results = results.with_columns(
        pl.struct(
            pl.col("spacer_seq"),
            pl.col("contig_seq"),
            pl.col("strand"),
            pl.col("start"),
            pl.col("end"),
        )
        .map_elements(
            lambda x: test_row(x, ignore_region_strands=ignore_region_strands),
            return_dtype=pl.Int32(),
        )
        .alias("alignment_test")
    )

    return results


def recalculate_mismatches_streaming(
    parquet_path: str,
    spacers_file: str,
    contigs_file: str,
    output_parquet: str,
    max_mismatches: int = 3,
    batch_size: int = 200000,
    threads: int = 12,
    memory_limit: str = "100GB",
    ignore_region_strands: bool = True,
):
    """
    Recalculate mismatches using streaming approach with DuckDB.
    
    This function:
    1. Uses DuckDB to extract unique regions (streaming, no full load)
    2. Batches sequence population and alignment validation
    3. Joins recalculated mismatches back to original data
    4. Filters based on recalculated mismatches
    
    Args:
        parquet_path: Path to parquet file with tool results
        spacers_file: Path to spacers FASTA file
        contigs_file: Path to contigs FASTA file
        output_parquet: Path to save filtered results
        max_mismatches: Maximum mismatches to keep after recalculation
        batch_size: Number of unique regions to process per batch
        threads: Number of threads for DuckDB
        memory_limit: Memory limit for DuckDB
        ignore_region_strands: Whether to ignore strand info during alignment
        
    Returns:
        None (writes to output_parquet)
    """
    import os
    
    # Check if output already exists and is valid
    if os.path.exists(output_parquet):
        try:
            with open(output_parquet, 'rb') as f:
                f.seek(-4, 2)  # Seek to last 4 bytes
                footer = f.read()
                if footer == b'PAR1':
                    print(f"\n✓ Valid parquet file already exists: {output_parquet}")
                    print("  Skipping recalculation. Delete the file to recompute.")
                    return None
                else:
                    print("\n⚠ Output file exists but appears corrupted (missing PAR1 footer)")
                    print("  Rerunning recalculation...")
        except Exception as e:
            print(f"\n⚠ Error checking output file: {e}")
            print("  Rerunning recalculation...")
    
    print("\n[Streaming Mismatch Recalculation]")
    print(f"  Input: {parquet_path}")
    print(f"  Output: {output_parquet}")
    print(f"  Batch size: {batch_size:,}")
    print(f"  Memory limit: {memory_limit}")
    
    # Create working directory for intermediate files next to the output
    temp_dir = os.path.join(os.path.dirname(output_parquet), "mismatch_recalc_temp")
    os.makedirs(temp_dir, exist_ok=True)
    # temp_dir = os.path.abspath(tempfile.mkdtemp(prefix="mismatch_recalc_"))
    # temp_dir = tempfile.mkdtemp(prefix="mismatch_recalc_")
    unique_regions_path = os.path.join(temp_dir, "unique_regions.parquet")
    validated_regions_path = os.path.join(temp_dir, "validated_regions.parquet")
    
    try:
        # Step 1: Extract unique regions using DuckDB (streaming, no memory load)
        # Check if unique_regions already exists and is valid
        region_count = None
        if os.path.exists(unique_regions_path):
            try:
                # Check if file is large enough to contain PAR1 footer
                file_size = os.path.getsize(unique_regions_path)
                if file_size >= 4:
                    with open(unique_regions_path, 'rb') as f:
                        f.seek(-4, 2)  # Seek to last 4 bytes
                        footer = f.read()
                        if footer == b'PAR1':
                            print("\n[Step 1/4] ✓ Valid unique_regions.parquet found, reusing...")
                            # Get count for progress tracking
                            con = duckdb.connect(database=":memory:")
                            con.execute(f"SET threads TO {threads};")
                            con.execute(f"SET memory_limit = '{memory_limit}';")
                            region_count = con.execute(f"SELECT COUNT(*) FROM read_parquet('{unique_regions_path}')").fetchone()[0]
                            print(f"  Found {region_count:,} unique regions")
                            con.close()
                        else:
                            print("\n[Step 1/4] ⚠ unique_regions.parquet exists but appears corrupted (invalid footer)")
                            print("  Re-extracting unique regions...")
                else:
                    print(f"\n[Step 1/4] ⚠ unique_regions.parquet exists but is too small ({file_size} bytes)")
                    print("  Re-extracting unique regions...")
            except Exception as e:
                print(f"\n[Step 1/4] ⚠ Error checking unique_regions.parquet: {e}")
                print("  Re-extracting unique regions...")
        
        # Only extract if we didn't find a valid file
        if region_count is None:
            print("\n[Step 1/4] Extracting unique regions via DuckDB...")
            con = duckdb.connect(database=":memory:")
            con.execute(f"SET threads TO {threads};")
            con.execute(f"SET memory_limit = '{memory_limit}';")
            
            # Extract unique regions and write to parquet
            con.execute(f"""
                COPY (
                    SELECT DISTINCT 
                        spacer_id, contig_id, strand, "start", "end"
                    FROM read_parquet('{parquet_path}')
                ) TO '{unique_regions_path}' (FORMAT PARQUET, COMPRESSION SNAPPY)
            """)
            
            # Get count for progress tracking
            region_count = con.execute(f"SELECT COUNT(*) FROM read_parquet('{unique_regions_path}')").fetchone()[0]
            print(f"  Found {region_count:,} unique regions")
            con.close()
        
        # Step 2: Populate sequences in batches
        print("\n[Step 2/4] Populating sequences in batches...")
        
        # Read unique regions in lazy mode
        unique_regions = pl.scan_parquet(unique_regions_path)
        
        # Process in batches to avoid OOM
        processed_batches = []
        total_batches = (region_count + batch_size - 1) // batch_size
        
        for batch_num in range(total_batches):
            offset = batch_num * batch_size
            print(f"  Processing batch {batch_num + 1}/{total_batches} (offset={offset:,})...")
            
            # Read batch
            batch = unique_regions.slice(offset, batch_size).collect()
            
            # Populate spacer sequences
            batch = populate_pldf_withseqs_needletail(
                seqfile=spacers_file,
                pldf=batch,
                chunk_size=500000,
                reverse_by_strand_col=False,
                trim_to_region=False,
                idcol="spacer_id",
                seqcol="spacer_seq"
            )
            
            # Populate contig sequences
            batch = populate_pldf_withseqs_needletail(
                seqfile=contigs_file,
                pldf=batch,
                chunk_size=100000,
                reverse_by_strand_col=True,
                trim_to_region=True,
                idcol="contig_id",
                start_col="start",
                end_col="end",
                strand_col="strand",
                seqcol="contig_seq"
            )
            
            # Recalculate mismatches using parasail
            batch = test_alignment_polars(
                batch,
                ignore_region_strands=ignore_region_strands
            )
            
            # Keep only needed columns
            batch = batch.select([
                "spacer_id", "contig_id", "strand", "start", "end",
                "alignment_test", "spacer_seq", "contig_seq"
            ])
            
            processed_batches.append(batch)
            
            # Optional: write intermediate results to avoid memory buildup
            if len(processed_batches) >= 10:  # Write every 10 batches
                print("    Writing intermediate results...")
                if os.path.exists(validated_regions_path):
                    # Append to existing file
                    pl.concat(processed_batches).write_parquet(
                        validated_regions_path, 
                        compression="snappy",
                        row_group_size=100000
                    )
                else:
                    # Create new file
                    pl.concat(processed_batches).write_parquet(
                        validated_regions_path,
                        compression="snappy"
                    )
                processed_batches = []
        
        # Write any remaining batches
        if processed_batches:
            print("\n  Writing final batch...")
            if os.path.exists(validated_regions_path):
                # This is tricky - need to append properly
                existing = pl.scan_parquet(validated_regions_path)
                new_data = pl.concat(processed_batches)
                pl.concat([existing.collect(), new_data]).write_parquet(
                    validated_regions_path,
                    compression="snappy"
                )
            else:
                pl.concat(processed_batches).write_parquet(
                    validated_regions_path,
                    compression="snappy"
                )
        
        print("\n[Step 3/4] Joining recalculated mismatches back to original data...")
        
        # Step 3: Join using DuckDB and filter
        con = duckdb.connect(database=":memory:")
        con.execute(f"SET threads TO {threads};")
        con.execute(f"SET memory_limit = '{memory_limit}';")
        
        # Join and filter in one DuckDB query (streaming)
        con.execute(f"""
            COPY (
                SELECT 
                    t.*,
                    v.alignment_test,
                    v.spacer_seq,
                    v.contig_seq
                FROM read_parquet('{parquet_path}') t
                LEFT JOIN read_parquet('{validated_regions_path}') v
                    ON t.spacer_id = v.spacer_id 
                    AND t.contig_id = v.contig_id
                    AND t.strand = v.strand
                    AND t."start" = v."start"
                    AND t."end" = v."end"
                WHERE v.alignment_test <= {max_mismatches}
            ) TO '{output_parquet}' (FORMAT PARQUET, COMPRESSION SNAPPY)
        """)
        
        final_count = con.execute(f"SELECT COUNT(*) FROM read_parquet('{output_parquet}')").fetchone()[0]
        print(f"  Filtered to {final_count:,} alignments (≤{max_mismatches} mismatches)")
        
        con.close()
        
        print("\n[Step 4/4] Cleanup...")
        
    finally:
        # Cleanup temp files
        import shutil
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
            print(f"  Removed temp directory: {temp_dir}")
    
    print(f"\n[Complete] Results written to: {output_parquet}")
    return None


def read_results(
    tools: dict,
    max_mismatches: int = 5,
    spacer_lendf: Optional[pl.DataFrame] = None,
    ref_file: Optional[str] = None,
    threads: int = 10,
    use_duckdb: bool = True,
    output_parquet: Optional[str] = None,
    memory_limit: str = str(min(50000, psutil.virtual_memory().available // (1024 * 1024 * 2))) + "MB", # set a default memory limit to 50gb or half of all available RAM (if less than 50gb)
):
    """
    Read and combine results from multiple tools with memory-efficient processing.

    Args:
        tools(dict): Dictionary of tool configurations
        max_mismatches(int): Maximum number of mismatches to allow
        spacer_lendf(polars.DataFrame): Optional spacer length dataframe
        ref_file(str): Optional reference file path
        threads(int): Number of threads for processing
        use_duckdb(bool): If True, use DuckDB with streaming; if False, use Polars lazy evaluation
        output_parquet(str): If provided, save results to this Parquet file path instead of returning in-memory DF

    Returns:
        Polars DataFrame with combined results from all tools (or None if output_parquet is specified)
    """

    if use_duckdb:
        return _read_results_duckdb(
            tools,
            max_mismatches,
            spacer_lendf,
            ref_file,
            threads,
            output_parquet,
            memory_limit=memory_limit,
        )
    else:
        return _read_results_polars_lazy(
            tools, max_mismatches, spacer_lendf, ref_file, threads, output_parquet
        )


def _read_results_duckdb(
    tools: dict,
    max_mismatches: int = 5,
    spacer_lendf: Optional[pl.DataFrame] = None,
    ref_file: Optional[str] = None,
    threads: int = 10,
    output_parquet: Optional[str] = None,
    memory_limit: str = "50GB",
):
    """
    Read results using DuckDB with native Python bindings for memory efficiency.
    Results are streamed to disk (Parquet) if output_parquet is specified.

    Thread behavior:
    - The 'threads' parameter controls query parallelism (worker threads for CPU-bound work).
    - DuckDB may spawn additional background threads for I/O operations.
    - This function enforces strict thread limiting via environment variables and OS controls.
    - Additional background I/O threads may still be created but should be minimal.
    """
    # STRICT thread enforcement via environment variables
    # These are set BEFORE DuckDB connection to take effect
    # os.environ["OMP_NUM_THREADS"] = str(threads)  # OpenMP (used by DuckDB and Polars)
    # os.environ["OPENBLAS_NUM_THREADS"] = str(threads)  # OpenBLAS
    # os.environ["MKL_NUM_THREADS"] = str(threads)  # Intel MKL
    # os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)  # macOS veclib
    # os.environ["NUMEXPR_NUM_THREADS"] = str(threads)  # NumExpr

    print(f"\n[ThreadControl] Setting environment variables for {threads} threads")
    print(f"  OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}")
    print(f"  OPENBLAS_NUM_THREADS={os.environ.get('OPENBLAS_NUM_THREADS')}")

    # Create a temporary directory for intermediate results if needed
    temp_dir = tempfile.mkdtemp(prefix="spacer_bench_duckdb_")
    parquet_files = []

    # Initialize DuckDB connection AFTER setting environment variables
    con = duckdb.connect(database=":memory:")

    # Set hardware tuning parameters - Control query parallelism
    con.execute(f"SET threads TO {threads};")
    con.execute("SET preserve_insertion_order = false;")
    con.execute(f"SET memory_limit = '{memory_limit}';")

    # Log actual thread setting
    thread_setting = con.execute("SELECT current_setting('threads')").fetchone()[0]
    print(f"  DuckDB query parallelism: {thread_setting} threads")
    print(f"  DuckDB memory limit: {memory_limit}")

    try:
        for tool in tools.values():
            try:
                print(f"\nReading results for {tool['name']}...")
                parse_function = tool.get("parse_function")

                # Check if file exists and has content
                if os.path.exists(tool["output_file"]):
                    file_size = os.path.getsize(tool["output_file"])
                    print(f"File size: {file_size / 1024 / 1024:.2f} MB")
                    if file_size == 0:
                        print(f"File {tool['output_file']} is empty, skipping")
                        continue
                else:
                    print(f"File {tool['output_file']} does not exist, skipping")
                    continue

                # Parse results using the tool-specific parser
                # (parsers may return Polars dataframes or lazy frames)
                parsed_results = (
                    parse_function(tool["output_file"], max_mismatches=max_mismatches)
                    if spacer_lendf is None
                    else parse_function(
                        tool["output_file"],
                        max_mismatches=max_mismatches,
                        spacer_lendf=spacer_lendf,
                        threads=threads,
                        ref_file=ref_file,
                    )
                )

                # Convert to lazy frame if it's eager (for consistent handling)
                if not hasattr(parsed_results, "collect"):
                    # It's an eager DataFrame, convert to lazy
                    parsed_results = parsed_results.lazy()

                # Add tool name column (works on lazy frames)
                parsed_results = parsed_results.with_columns(
                    pl.lit(tool["name"]).alias("tool")
                )

                # Write lazy frame to parquet without materializing
                # sink_parquet handles lazy frames efficiently
                tool_parquet = os.path.join(temp_dir, f"{tool['name']}.parquet")
                parsed_results.sink_parquet(tool_parquet)
                parquet_files.append(tool_parquet)
                print(f"Saved {tool['name']} results to {tool_parquet}")

            except Exception as e:
                print(f"Failed to read results for {tool['name']}: {e}")
                import traceback

                traceback.print_exc()
                continue

        # If no results were collected, return empty dataframe
        if not parquet_files:
            print("No results collected from any tool")
            return pl.DataFrame(
                schema={
                    "spacer_id": pl.Utf8,
                    "contig_id": pl.Utf8,
                    "spacer_length": pl.UInt32,
                    "strand": pl.Boolean,
                    "start": pl.UInt32,
                    "end": pl.UInt32,
                    "mismatches": pl.UInt32,
                    "tool": pl.Utf8,
                },
                orient="row",
            )

        # Combine all parquet files using DuckDB's union operation
        print(f"\nCombining results from {len(parquet_files)} tools using DuckDB...")

        # Build union by concatenating all parquet files
        union_query = " UNION ALL ".join(
            [f"SELECT * FROM read_parquet('{pf}')" for pf in parquet_files]
        )

        # If output_parquet is specified, write directly to file using DuckDB's native COPY
        if output_parquet:
            print(
                f"Writing combined results to {output_parquet} using DuckDB (no materialization)..."
            )
            # Use DuckDB's COPY command to write directly without materializing in Polars
            copy_query = f"""
                COPY (
                    {union_query}
                ) TO '{output_parquet}' (FORMAT 'PARQUET', COMPRESSION 'SNAPPY', ROW_GROUP_SIZE 100000)
            """
            con.execute(copy_query)
            print(f"Results saved to {output_parquet}")
            result_df = None
        else:
            # Only load to Polars if not saving to Parquet
            print("Loading results to Polars...")
            union_relation = con.execute(union_query).pl()
            print("Results ready")
            result_df = union_relation

        return result_df

    finally:
        # Clean up temporary files
        import shutil

        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
            print(f"Cleaned up temporary directory: {temp_dir}")


def _read_results_polars_lazy(
    tools : dict,
    max_mismatches: int = 5,
    spacer_lendf: Optional[pl.DataFrame] = None,
    ref_file: Optional[str] = None,
    threads: int = 4,
    output_parquet: Optional[str] = None,
):
    """
    Read results using Polars lazy evaluation for memory efficiency.
    This alternative uses Polars' lazy evaluation to defer computation.
    """
    lazy_frames = []

    try:
        for tool in tools.values():
            try:
                print(f"\nReading results for {tool['name']}...")
                parse_function = tool.get("parse_function")

                # Check if file exists and has content
                if os.path.exists(tool["output_file"]):
                    file_size = os.path.getsize(tool["output_file"])
                    print(f"File size: {file_size / 1024 / 1024:.2f} MB")
                    if file_size == 0:
                        print(f"File {tool['output_file']} is empty, skipping")
                        continue
                else:
                    print(f"File {tool['output_file']} does not exist, skipping")
                    continue

                # Parse results using the tool-specific parser
                parsed_results = (
                    parse_function(tool["output_file"], max_mismatches=max_mismatches)
                    if spacer_lendf is None
                    else parse_function(
                        tool["output_file"],
                        max_mismatches=max_mismatches,
                        spacer_lendf=spacer_lendf,
                        threads=threads,
                        ref_file=ref_file,
                    )
                )

                # Add tool name column and convert to lazy frame
                parsed_results = parsed_results.with_columns(
                    pl.lit(tool["name"]).alias("tool")
                ).lazy()

                lazy_frames.append(parsed_results)
                print(f"Loaded {tool['name']} results (lazy evaluation)")

            except Exception as e:
                print(f"Failed to read results for {tool['name']}: {e}")
                import traceback

                traceback.print_exc()
                continue

        # If no results were collected, return empty dataframe
        if not lazy_frames:
            print("No results collected from any tool")
            return pl.DataFrame(
                schema={
                    "spacer_id": pl.Utf8,
                    "contig_id": pl.Utf8,
                    "spacer_length": pl.UInt32,
                    "strand": pl.Boolean,
                    "start": pl.UInt32,
                    "end": pl.UInt32,
                    "mismatches": pl.UInt32,
                    "tool": pl.Utf8,
                },
                orient="row",
            )

        # Combine using lazy concatenation
        print(
            f"\nCombining results from {len(lazy_frames)} tools using Polars lazy evaluation..."
        )
        combined_lazy = pl.concat(lazy_frames, how="vertical_relaxed")

        # If output_parquet is specified, write directly to file
        if output_parquet:
            print(f"Writing combined results to {output_parquet}...")
            combined_lazy.sink_parquet(output_parquet)
            print(f"Results saved to {output_parquet}")
            return None
        else:
            # Collect and return
            print("Collecting results...")
            result_df = combined_lazy.collect()
            return result_df

    except Exception as e:
        print(f"Error in Polars lazy evaluation: {e}")
        import traceback

        traceback.print_exc()
        raise


def verify_sam_file(sam_file: str, ref_file: Optional[str] = None):
    """Verify if a SAM file is valid and add SQ lines if needed (from reference file if provided or from other sam files in the same directory)
    also potentially replaces whitespaces with tabs in the first line"""

    # Read the first few lines to analyze the header
    with open(sam_file, "r") as f:
        lines = f.readlines()

    if len(lines) < 2:
        print("SAM file has less than 2 lines")
        return None

    first_line = lines[0].strip()
    second_line = lines[1].strip() if len(lines) > 1 else ""

    print("first line of sam file", first_line)
    print("second line of sam file", second_line)

    # Check if this is a mummer file - always rewrite HD line for mummer
    is_mummer_file = "mummer" in sam_file.lower()

    # Check if we have a proper HD header
    has_hd_header = first_line.startswith("@HD")

    # Check if HD line has spaces instead of tabs (malformed)
    hd_needs_fixing = False
    if has_hd_header and " " in first_line and "\t" not in first_line:
        hd_needs_fixing = True
        print("HD line has spaces instead of tabs, will fix")

    # Check if we have SQ lines
    sq_lines = []
    pg_lines = []
    alignment_start = 0

    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith("@SQ"):
            sq_lines.append(line)
        elif line.startswith("@PG"):
            pg_lines.append(line)
        elif not line.startswith("@"):
            alignment_start = i
            break

    has_sq_lines = len(sq_lines) > 0
    has_pg_lines = len(pg_lines) > 0

    print(f"Found {len(sq_lines)} SQ lines, {len(pg_lines)} PG lines")

    # Determine if we need to fix the file
    needs_fixing = (
        hd_needs_fixing or not has_sq_lines or not has_pg_lines or is_mummer_file
    )

    if needs_fixing:
        print("SAM file needs fixing")

        # Create new file with proper header
        with open(sam_file + ".new", "w") as new_sam_file:
            # Handle HD line
            if is_mummer_file:
                # Always rewrite HD line for mummer files
                new_sam_file.write("@HD\tVN:1.6\tSO:unsorted\n")
                print("Rewriting HD line for mummer file")
            elif has_hd_header:
                if hd_needs_fixing:
                    # Replace spaces with tabs in HD line
                    fixed_hd = first_line.replace(" VN:", "\tVN:").replace(
                        " SO:", "\tSO:"
                    )
                    new_sam_file.write(fixed_hd + "\n")
                else:
                    new_sam_file.write(first_line + "\n")
            else:
                # Add default HD line if missing
                new_sam_file.write("@HD\tVN:1.6\tSO:unsorted\n")

            # Add SQ lines
            if has_sq_lines:
                for sq_line in sq_lines:
                    new_sam_file.write(sq_line + "\n")
            elif ref_file is not None:
                print("Reference file provided, using it to add SQ lines")
                try:
                    tmp = pysam.FastaFile(ref_file)
                    for ref_name, ref_length in zip(tmp.references, tmp.lengths):
                        new_sam_file.write(f"@SQ\tSN:{ref_name}\tLN:{ref_length}\n")
                except Exception as e:
                    print(f"Error reading reference file: {e}")
                    return None
            else:
                print(
                    "No reference file provided, trying to find SQ lines from other SAM files"
                )
                present_sam_files = glob.glob(
                    os.path.join(os.path.dirname(sam_file), "*.sam")
                )
                present_sam_files = [
                    f for f in present_sam_files if f != os.path.normpath(sam_file)
                ]
                # Drop mummer4 and minimap2 sam files as they don't print SQ lines when batching
                present_sam_files = [
                    f
                    for f in present_sam_files
                    if "mummer" not in f and "minimap" not in f
                ]
                present_sam_files = [
                    f for f in present_sam_files if os.path.getsize(f) > 1
                ]

                if len(present_sam_files) == 0:
                    print("No non-empty sam files found, returning None")
                    return None

                # Find a SAM file with SQ lines
                sq_found = False
                for sam_file_candidate in present_sam_files:
                    with open(sam_file_candidate, "r") as f:
                        for line in f:
                            if line.startswith("@SQ"):
                                new_sam_file.write(line)
                                sq_found = True
                            elif line.startswith("@PG") or (
                                not line.startswith("@") and line.strip()
                            ):
                                break
                    if sq_found:
                        break

                if not sq_found:
                    print("No SQ lines found in any SAM files")
                    return None

            # Add PG lines
            if has_pg_lines:
                for pg_line in pg_lines:
                    # Fix PG line if it has spaces instead of tabs
                    if " " in pg_line and "\t" not in pg_line:
                        fixed_pg = (
                            pg_line.replace(" ID:", "\tID:")
                            .replace(" PN:", "\tPN:")
                            .replace(" VN:", "\tVN:")
                            .replace(" CL:", "\tCL:")
                        )
                        new_sam_file.write(fixed_pg + "\n")
                    else:
                        new_sam_file.write(pg_line + "\n")
            else:
                # Add default PG line if missing
                new_sam_file.write("@PG\tID:unknown\tPN:unknown\tVN:unknown\n")

            # Add alignment lines
            for i in range(alignment_start, len(lines)):
                line = lines[i]
                if not line.startswith("@"):  # Keep only the alignments
                    new_sam_file.write(line)

        # Replace original file with fixed version
        os.rename(sam_file + ".new", sam_file)
        print("SAM file has been fixed")

    else:
        print("SAM file looks good, no changes needed")

    return sam_file


def create_comparison_matrix(tools_results: pl.DataFrame, n_mismatches: int = 3):
    # Filter for specific number of mismatches
    tmp = tools_results.filter(pl.col("mismatches") == n_mismatches)

    # Get unique tools
    tools = tmp.select("tool").unique()
    tools_list = tools.get_column("tool").to_list()

    # Create empty matrix dataframe
    matrix_data = []

    for tool1 in tools_list:
        row_data = []
        # Get unique pairs for tool1
        tool1_pairs = (
            tmp.filter(pl.col("tool") == tool1)
            .select(["contig_id", "spacer_id"])
            .unique()
        )

        for tool2 in tools_list:
            if tool1 == tool2:
                # Diagonal will show total number of pairs for the tool
                row_data.append(tool1_pairs.height)
            else:
                # Get unique pairs for tool2
                tool2_pairs = (
                    tmp.filter(pl.col("tool") == tool2)
                    .select(["contig_id", "spacer_id"])
                    .unique()
                )

                # Count pairs in tool1 that are not in tool2
                diff_count = tool1_pairs.join(
                    tool2_pairs, on=["contig_id", "spacer_id"], how="anti"
                ).height
                row_data.append(diff_count)

        matrix_data.append(row_data)

    # Create matrix dataframe
    matrix_df = pl.DataFrame(
        matrix_data,
        schema=tools_list,
    ).with_columns(pl.Series(name="tool", values=tools_list))

    return matrix_df


def get_parse_function(func_name: str):
    return globals()[func_name]


class DebugArgs:
    def __init__(self):
        self.n_mismatch_range = [0, 0]
        self.contig_length_range = [1200, 12000]
        self.spacer_length_range = [19, 60]
        self.sample_size_contigs = 154
        self.sample_size_spacers = 7
        self.insertion_range = [5, 10]
        self.threads = 1
        self.contigs = None
        self.spacers = None
        self.prop_rc = 0.5  # Add default value for prop_rc
        self.max_mismatches = 5
        self.max_runs = 1


def validate_intervals_with_polars_bio(
    planned_intervals: pl.DataFrame, tool_results: pl.DataFrame, max_mismatches: int = 5
):
    """
    Validate tool results against planned intervals using polars-bio interval operations.

    Args:
        planned_intervals(polars.DataFrame): DataFrame with planned spacer insertions (spacer_id, contig_id, start, end, strand, mismatches)
        tool_results(polars.DataFrame): DataFrame with tool results (spacer_id, contig_id, start, end, strand, mismatches)
        max_mismatches(int): Maximum allowed mismatches for validation

    Returns:
        DataFrame with validation results
    """
    # Filter tool results by max mismatches
    tool_results_filtered = tool_results.filter(pl.col("mismatches") <= max_mismatches)

    # Add chrom column for polars-bio compatibility (use contig_id as chrom)
    planned_df = planned_intervals.with_columns(
        [
            pl.col("start").cast(pl.UInt32),
            pl.col("end").cast(pl.UInt32),
            pl.col("strand").cast(pl.Utf8),
        ]
    )

    tool_df = tool_results_filtered.with_columns(
        [
            pl.col("start").cast(pl.UInt32),
            pl.col("end").cast(pl.UInt32),
            pl.col("strand").cast(pl.Utf8),
        ]
    )

    # Check if dataframes are empty
    if planned_df.height == 0:
        print("Warning: planned_df is empty")
        return pl.DataFrame(
            schema={
                "spacer_id_1": pl.Utf8,
                "contig_id_1": pl.Utf8,
                "start_1": pl.UInt32,
                "end_1": pl.UInt32,
                "strand_1": pl.Utf8,
                "spacer_id_2": pl.Utf8,
                "contig_id_2": pl.Utf8,
                "start_2": pl.UInt32,
                "end_2": pl.UInt32,
                "strand_2": pl.Utf8,
                "overlap_type": pl.Utf8,
            }
        )
    if tool_df.height == 0:
        print("Warning: tool_df is empty")
        return pl.DataFrame(
            schema={
                "spacer_id_1": pl.Utf8,
                "contig_id_1": pl.Utf8,
                "start_1": pl.UInt32,
                "end_1": pl.UInt32,
                "strand_1": pl.Utf8,
                "spacer_id_2": pl.Utf8,
                "contig_id_2": pl.Utf8,
                "start_2": pl.UInt32,
                "end_2": pl.UInt32,
                "strand_2": pl.Utf8,
                "overlap_type": pl.Utf8,
            },
            orient="row",
        )

    # Use polars-bio overlap operation
    try:
        overlap_result = pb.overlap(
            planned_df,
            tool_df,
            cols1=["contig_id", "start", "end"],
            cols2=["contig_id", "start", "end"],
        )
    except Exception as e:
        print(f"Error in pb.overlap: {e}")
        print(f"planned_df shape: {planned_df.shape}")
        print(f"tool_df shape: {tool_df.shape}")
        raise

    # Calculate validation metrics with strand consideration
    validation_results = overlap_result.with_columns(
        [
            # First check if strands match (both forward or both reverse)
            pl.when(pl.col("strand_1") == pl.col("strand_2"))
            .then(
                # If strands match, check for coordinate overlaps
                # Use symmetric overlap detection logic
                pl.when(
                    # Exact match: intervals are identical
                    (pl.col("start_1") == pl.col("start_2"))
                    & (pl.col("end_1") == pl.col("end_2"))
                )
                .then(pl.lit("exact_match"))
                .when(
                    # Partial overlap: intervals overlap but are not identical
                    (pl.col("start_1") < pl.col("end_2"))
                    & (pl.col("start_2") < pl.col("end_1"))
                )
                .then(pl.lit("partial_overlap"))
                .otherwise(pl.lit("no_overlap"))
            )
            .otherwise(pl.lit("strand_mismatch"))  # Different strands = no valid match
            .alias("overlap_type")
        ]
    ).collect()

    return validation_results


def validate_intervals_fallback(planned_intervals: pl.DataFrame, tool_results: pl.DataFrame, max_mismatches: int = 5):
    """
    Fallback interval validation without polars-bio.
    """
    # Filter tool results by max mismatches
    tool_results_filtered = tool_results.filter(pl.col("mismatches") <= max_mismatches)

    # Simple overlap detection using polars operations
    joined = planned_intervals.join(
        tool_results_filtered,
        on=["spacer_id", "contig_id", "strand"],
        how="inner",
        suffix="_tool",
    )

    # Check for overlaps
    overlap_conditions = (
        (
            (pl.col("start") <= pl.col("start_tool"))
            & (pl.col("end") >= pl.col("end_tool"))
        )
        | (
            (pl.col("start_tool") <= pl.col("start"))
            & (pl.col("end_tool") >= pl.col("end"))
        )
        | (
            (pl.col("start") < pl.col("end_tool"))
            & (pl.col("end") > pl.col("start_tool"))
        )
    )

    validation_results = joined.with_columns(
        [
            pl.when(overlap_conditions)
            .then(pl.lit("overlap"))
            .otherwise(pl.lit("no_overlap"))
            .alias("overlap_type")
        ]
    )

    return validation_results


def load_multiple_fractions_duckdb(
    fraction_parquet_files, threads=4, memory_limit="50GB"
):
    """
    Load and combine results from multiple fraction parquet files using DuckDB.
    This enables analysis of stratified subsamples without loading all data into memory.

    Args:
        fraction_parquet_files(dict): Dictionary mapping fraction values to parquet file paths
            e.g., {0.001: 'path/to/alignments_fraction_0.001.parquet', ...}
        threads(int): Number of threads for DuckDB to use
        memory_limit(str): Memory limit for DuckDB operations (e.g., "50GB")

    Returns
    duckdb.DuckDBPyConnection
        A DuckDB connection with a registered view "fractions" containing all loaded data
    """
    # STRICT thread enforcement via environment variables
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)

    print(
        f"\n[DuckDB] Loading {len(fraction_parquet_files)} fractions with {threads} threads"
    )

    con = duckdb.connect(database=":memory:")
    con.execute(f"SET threads TO {threads};")
    con.execute("SET preserve_insertion_order = false;")
    con.execute(f"SET memory_limit = '{memory_limit}';")

    # Build union query for all fractions
    union_queries = []
    for frac, pf in sorted(fraction_parquet_files.items()):
        if os.path.exists(pf) and os.path.getsize(pf) > 0:
            union_queries.append(
                f"SELECT *, {frac} as fraction FROM read_parquet('{pf}')"
            )
            print(f"  Fraction {frac}: {os.path.getsize(pf) / 1024 / 1024:.1f} MB")
        else:
            print(f"  Fraction {frac}: NOT FOUND or empty, skipping")

    if not union_queries:
        raise ValueError("No valid parquet files found")

    # Create view for all fractions
    union_query = " UNION ALL ".join(union_queries)
    con.execute(f"CREATE VIEW fractions AS {union_query}")
    print(f"[DuckDB] Created view 'fractions' with {len(union_queries)} fractions")

    return con


def query_fractions_duckdb(con, query_template, output_parquet=None):
    """
    Execute a query on the DuckDB fractions view and optionally save to parquet.

    Args:
        con(duckdb.DuckDBPyConnection): DuckDB connection from load_multiple_fractions_duckdb()
        query_template(str): SQL query, use 'fractions' as the table name.
            e.g., "SELECT tool, COUNT(*) as count FROM fractions GROUP BY tool"
        output_parquet(str): Path to save results as parquet. If None, returns Polars DataFrame.

    Returns
    polars.DataFrame or None
        If output_parquet is None, returns results as Polars DataFrame.
        If output_parquet is specified, saves to disk and returns None.
    """
    print("\n[DuckDB] Executing query...")

    if output_parquet:
        copy_query = f"""
            COPY (
                {query_template}
            ) TO '{output_parquet}' (FORMAT 'PARQUET', COMPRESSION 'SNAPPY', ROW_GROUP_SIZE 100000)
        """
        con.execute(copy_query)
        print(f"[DuckDB] Results saved to {output_parquet}")
        return None
    else:
        result = con.execute(query_template).pl()
        return result


def create_spacer_counts_with_tools_duckdb(
    parquet_path: str,
    tools_list: list,
    mismatches: int = 3,
    exact_or_max: str = "exact",
    output_parquet: str = None,
    threads: int = 12,
    memory_limit: str = "100GB",
):
    """
    DuckDB version of create_spacer_counts_with_tools.
    Calculates occurrence counts and tool detection fractions without loading full dataset.
    
    Args:
        parquet_path: Path to parquet file with recalculated mismatches
        tools_list: List of tool names
        mismatches: Number of mismatches (exact or max depending on exact_or_max)
        exact_or_max: "exact" for == mismatches, "max" for <= mismatches
        output_parquet: Optional path to save results
        threads: Number of threads for DuckDB
        memory_limit: Memory limit for DuckDB
        
    Returns:
        Polars DataFrame or None (if output_parquet specified)
    """
    con = duckdb.connect(database=":memory:")
    con.execute(f"SET threads TO {threads};")
    con.execute(f"SET memory_limit = '{memory_limit}';")
    
    # Build mismatch filter
    if exact_or_max == "max":
        mismatch_filter = f"alignment_test <= {mismatches}"
    else:
        mismatch_filter = f"alignment_test = {mismatches}"
    
    # Create tools list SQL
    tools_sql = ", ".join([f"'{tool}'" for tool in tools_list])
    
    query = f"""
        WITH spacer_counts AS (
            -- Get total occurrences per spacer
            SELECT 
                spacer_id,
                COUNT(DISTINCT contig_id) as n_occurrences
            FROM read_parquet('{parquet_path}')
            WHERE {mismatch_filter}
            GROUP BY spacer_id
        ),
        tool_matches AS (
            -- Get matches per tool and spacer
            SELECT 
                spacer_id,
                tool,
                COUNT(DISTINCT contig_id) as tool_matches
            FROM read_parquet('{parquet_path}')
            WHERE {mismatch_filter}
            GROUP BY spacer_id, tool
        ),
        all_combinations AS (
            -- Cross join spacers with all tools
            SELECT 
                sc.spacer_id,
                sc.n_occurrences,
                t.tool
            FROM spacer_counts sc
            CROSS JOIN (SELECT UNNEST([{tools_sql}]) as tool) t
        )
        SELECT 
            ac.spacer_id,
            ac.n_occurrences,
            ac.tool,
            COALESCE(tm.tool_matches, 0) as tool_matches,
            COALESCE(tm.tool_matches, 0)::DOUBLE / ac.n_occurrences as fraction
        FROM all_combinations ac
        LEFT JOIN tool_matches tm
            ON ac.spacer_id = tm.spacer_id 
            AND ac.tool = tm.tool
    """
    
    # Pivot to get tools as columns
    tool_columns = ', '.join([f"MAX(CASE WHEN tool = '{tool}' THEN fraction ELSE 0 END) as \"{tool}\"" for tool in tools_list])
    pivot_query = f"""
        WITH base_data AS ({query})
        SELECT 
            spacer_id,
            n_occurrences,
            {tool_columns}
        FROM base_data
        GROUP BY spacer_id, n_occurrences
    """
    
    if output_parquet:
        con.execute(f"""
            COPY ({pivot_query}) 
            TO '{output_parquet}' (FORMAT PARQUET, COMPRESSION SNAPPY)
        """)
        con.close()
        return None
    else:
        result = con.execute(pivot_query).pl()
        con.close()
        return result


def get_summary_stats_duckdb(
    parquet_path: str,
    threads: int = 12,
    memory_limit: str = "100GB",
):
    """
    Get summary statistics from parquet file using DuckDB.
    
    Returns summary by tool without loading full dataset into memory.
    """
    con = duckdb.connect(database=":memory:")
    con.execute(f"SET threads TO {threads};")
    con.execute(f"SET memory_limit = '{memory_limit}';")
    
    query = """
        SELECT 
            tool,
            AVG(alignment_test) as mean_mismatches,
            COUNT(DISTINCT spacer_id) as n_spacers,
            COUNT(DISTINCT contig_id) as n_contigs,
            COUNT(*) as total_alignments
        FROM read_parquet(?)
        GROUP BY tool
        ORDER BY tool
    """
    
    result = con.execute(query, [parquet_path]).pl()
    con.close()
    return result


def create_tool_comparison_matrix_duckdb(
    parquet_path: str,
    tools_list: list,
    n_mismatches: int,
    output_csv: str = None,
    threads: int = 12,
    memory_limit: str = "100GB",
):
    """
    Create tool comparison matrix using DuckDB.
    Cell(i,j) = number of unique pairs in tool i but not in tool j.
    
    Returns Polars DataFrame with matrix.
    """
    import numpy as np
    
    con = duckdb.connect(database=":memory:")
    con.execute(f"SET threads TO {threads};")
    con.execute(f"SET memory_limit = '{memory_limit}';")
    
    # Create empty matrix
    matrix_data = np.zeros((len(tools_list), len(tools_list)), dtype=int)
    
    # Get unique pairs for each tool
    for i, tool_x in enumerate(tools_list):
        for j, tool_y in enumerate(tools_list):
            if tool_x == tool_y:
                continue
            
            # Count pairs in x but not in y using DuckDB
            query = f"""
                WITH tool_x_pairs AS (
                    SELECT DISTINCT contig_id, spacer_id, strand, "start", "end"
                    FROM read_parquet('{parquet_path}')
                    WHERE tool = '{tool_x}' AND alignment_test = {n_mismatches}
                ),
                tool_y_pairs AS (
                    SELECT DISTINCT contig_id, spacer_id, strand, "start", "end"
                    FROM read_parquet('{parquet_path}')
                    WHERE tool = '{tool_y}' AND alignment_test = {n_mismatches}
                )
                SELECT COUNT(*) as count
                FROM tool_x_pairs
                WHERE NOT EXISTS (
                    SELECT 1 FROM tool_y_pairs
                    WHERE tool_x_pairs.contig_id = tool_y_pairs.contig_id
                        AND tool_x_pairs.spacer_id = tool_y_pairs.spacer_id
                        AND tool_x_pairs.strand = tool_y_pairs.strand
                        AND tool_x_pairs."start" = tool_y_pairs."start"
                        AND tool_x_pairs."end" = tool_y_pairs."end"
                )
            """
            
            count = con.execute(query).fetchone()[0]
            matrix_data[i, j] = count
    
    con.close()
    
    # Convert to DataFrame
    matrix = pl.DataFrame(matrix_data, schema=tools_list)
    matrix = matrix.with_columns(pl.Series(name="tool1", values=tools_list, dtype=pl.Utf8))
    
    if output_csv:
        matrix.write_csv(output_csv, separator='\t')
    
    return matrix


def load_fraction_results_lazy(fraction_parquet_files, add_fraction_col=True):
    """
    Load multiple fraction parquet files as lazy Polars frames for memory-efficient processing.

    Args:
        fraction_parquet_files(dict): Dictionary mapping fraction values to parquet file paths
        add_fraction_col(bool): If True, add fraction column to each frame

    Returns
    polars.LazyFrame
        Concatenated lazy frames from all fractions
    """
    lazy_frames = []

    for frac, pf in sorted(fraction_parquet_files.items()):
        if os.path.exists(pf) and os.path.getsize(pf) > 0:
            lf = pl.scan_parquet(pf)
            if add_fraction_col:
                lf = lf.with_columns(pl.lit(frac).alias("fraction"))
            lazy_frames.append(lf)
            print(f"  Fraction {frac}: loaded (lazy)")
        else:
            print(f"  Fraction {frac}: NOT FOUND or empty, skipping")

    if not lazy_frames:
        raise ValueError("No valid parquet files found")

    combined = pl.concat(lazy_frames)
    print(f"[Polars] Concatenated {len(lazy_frames)} fractions as lazy frame")

    return combined


def analyze_sassy_edit_distance_false_positives(
    sassy_file,
    max_edit_distance=5,
    spacer_lendf=None,
    ref_file=None,
    threads=4,
    threads_alignment=1,
):
    """
    Analyze false positive rate using sassy output (which uses edit distance only).
    Compare actual observed false positives to expected values based on dataset statistics.

    Args:
        sassy_file(str): Path to sassy TSV output
        max_edit_distance(int): Maximum edit distance to consider
        spacer_lendf(polars.DataFrame): DataFrame with 'spacer_id' and 'length' columns
        ref_file(str): Path to reference contig fasta file
        threads(int): Threads for DuckDB operations
        threads_alignment(int): Threads for parasail alignment calculations

    Returns
    dict
        Analysis results including false positive rates, observed vs expected, etc.
    """
    # STRICT thread enforcement
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)

    print("\n[EditDistance] Analyzing false positives in sassy output (edit distance)")

    # Parse sassy output with edit distance filter
    sassy_results = parse_sassy(
        sassy_file,
        max_mismatches=max_edit_distance,
        spacer_lendf=spacer_lendf,
        threads=threads,
        ref_file=ref_file,
    )

    print(
        f"Loaded {sassy_results.height:,} sassy alignments with edit distance ≤ {max_edit_distance}"
    )

    # Load spacer and contig sequences for recalculation with hamming distance
    print("Recalculating alignments with hamming distance for comparison...")

    # Convert to eager if needed
    if hasattr(sassy_results, "collect"):
        sassy_results = sassy_results.collect()

    # Get unique regions for verification
    unique_regions = sassy_results.select(
        ["spacer_id", "contig_id", "strand", "start", "end"]
    ).unique()

    print(f"Verifying {unique_regions.height:,} unique regions...")

    # Load sequences
    unique_regions = populate_pldf_withseqs_needletail(
        seqfile=spacer_lendf.to_dict(as_series=False)
        if isinstance(spacer_lendf, pl.DataFrame)
        else spacer_lendf,
        pldf=unique_regions,
        idcol="spacer_id",
        seqcol="spacer_seq",
    )

    # This would require more detailed implementation based on your specific needs
    # For now, return basic structure
    return {
        "sassy_edit_results": sassy_results,
        "unique_regions_analyzed": unique_regions.height,
        "message": "Edit distance false positive analysis - detailed implementation in progress",
    }


def compare_tool_results_to_ground_truth(
    tool_results, ground_truth, max_hamming=3, max_edit=5
):
    """
    Compare tool alignment results to ground truth using both hamming and edit distance.

    For simulated data, we have ground truth and can classify each tool result as:
    - TRUE_POSITIVE_HAMMING: In ground truth AND hamming distance <= max_hamming
    - TRUE_POSITIVE_EDIT: In ground truth AND edit distance <= max_edit
    - FALSE_POSITIVE_HAMMING: NOT in ground truth BUT hamming distance <= max_hamming
    - FALSE_POSITIVE_EDIT: NOT in ground truth BUT edit distance <= max_edit
    - MISSED (False Negative): In ground truth BUT distance > threshold

    Args:
        tool_results(polars.DataFrame): DataFrame with tool results (spacer_id, contig_id, start, end, strand, mismatches)
        ground_truth(polars.DataFrame): DataFrame with ground truth (spacer_id, contig_id, start, end, strand, mismatches)
        max_hamming(int): Maximum hamming distance threshold
        max_edit(int): Maximum edit distance threshold

    Returns:
        DataFrame with results categorized by type
    """
    # Create position-based keys for ground truth (accounting for possible multiple insertions)
    gt_set = set()
    for row in ground_truth.iter_rows(named=True):
        key = (row["spacer_id"], row["contig_id"], row["strand"])
        gt_set.add(key)

    # Load sequences if not already present
    if (
        "spacer_seq" not in tool_results.columns
        or "contig_seq" not in tool_results.columns
    ):
        raise ValueError("tool_results must have spacer_seq and contig_seq columns")

    # Recalculate distances and categorize
    categorized = []
    for row in tool_results.iter_rows(named=True):
        # Recalculate hamming and edit distance
        hamming = calculate_hamming_distance(
            row["spacer_seq"],
            row["contig_seq"],
            strand=row["strand"],
            start=row["start"],
            end=row["end"],
        )

        edit = calculate_edit_distance(
            row["spacer_seq"],
            row["contig_seq"],
            strand=row["strand"],
            start=row["start"],
            end=row["end"],
        )

        # Check if in ground truth (positional)
        key = (row["spacer_id"], row["contig_id"], row["strand"])
        in_ground_truth = key in gt_set

        # Categorize
        if in_ground_truth:
            hamming_cat = (
                "TRUE_POSITIVE_HAMMING" if hamming <= max_hamming else "MISSED_HAMMING"
            )
            edit_cat = "TRUE_POSITIVE_EDIT" if edit <= max_edit else "MISSED_EDIT"
        else:
            hamming_cat = (
                "FALSE_POSITIVE_HAMMING" if hamming <= max_hamming else "NOT_ALIGNED"
            )
            edit_cat = "FALSE_POSITIVE_EDIT" if edit <= max_edit else "NOT_ALIGNED"

        categorized.append(
            {
                "spacer_id": row["spacer_id"],
                "contig_id": row["contig_id"],
                "start": row["start"],
                "end": row["end"],
                "strand": row["strand"],
                "hamming_distance": hamming,
                "edit_distance": edit,
                "in_ground_truth": in_ground_truth,
                "hamming_category": hamming_cat,
                "edit_category": edit_cat,
            }
        )

    return pl.DataFrame(categorized)
