# some utility functions for the benchmark.
import glob

import hashlib
import json
import os
import random
import re
import subprocess
import tempfile
import ipdb
import parasail as ps
import polars as pl
import pyfastx as pfx
import pysam
from needletail import parse_fastx_file  # , NeedletailError, normalize_seq
from tqdm import tqdm
from collections import Counter
import math
import numpy as np
import time
import polars_bio as pb
import duckdb
import uuid

# so printing to stdout doesn't break, line wrap, or truncate.
pl.Config.set_tbl_rows(123123)
pl.Config.set_tbl_cols(123123) # should be large enough
pl.Config.set_fmt_str_lengths(2100) 
pl.Config.set_tbl_width_chars(2100) 

def read_fasta_needletail(fasta_file):
    seqs = []
    seq_ids = []
    for record in parse_fastx_file(fasta_file):
        seqs.append(record.seq)
        seq_ids.append(record.id)
    return seq_ids, seqs

def apply_mismatches(sequence, n_mismatches):
    mismatch_positions = random.sample(range(len(sequence)), n_mismatches)
    sequence_list = list(sequence)
    for pos in mismatch_positions:
        original_base = sequence_list[pos]
        new_base = random.choice([b for b in "ATCG" if b != original_base])
        sequence_list[pos] = new_base
    return "".join(sequence_list)

def calculate_shannon_entropy(s: str) -> float:
    if not s:
        return 0.0  #here in case of NULLs so that applting over ovject of size n still returns a number.  TODO: think if 0.0 is the best place holder value...

    # Count character frequencies
    char_counts = Counter(s)
    total_chars = len(s)

    entropy = 0.0
    for count in char_counts.values():
        probability = count / total_chars
        # entropy -= probability * math.log10(probability) # TODO: check which base other people use.
        entropy -= probability * math.log2(probability)
    return entropy

def count_kmers_df_explicit(df: pl.DataFrame, seq_col: str = "seq",id_col: str = "seqid", k: int = 3, relative: bool = False) -> pl.DataFrame:
    """Calculate ALL k-mers counts for all sequences in a DataFrame. all possible k-mers are counted, not just the complete ones."""
    # Split sequences into characters
    import itertools
    all_kmers = [''.join(p) for p in itertools.product('ATCGN', repeat=k)]
    count_df = df.with_columns(
        pl.col(seq_col).str.extract_many(all_kmers, overlapping=True,ascii_case_insensitive=True).alias('kmers')
    ).group_by(id_col).agg(
        pl.col('kmers').explode().value_counts(normalize=relative).alias(f'kmer_{k}_relative' if relative else f'kmer_{k}_counts'),
    )
    return count_df



def count_kmers_df(df: pl.DataFrame, seq_col: str = "seq",id_col: str = "seqid", k: int = 3, relative: bool = False) -> pl.DataFrame:
    """Calculate k-mer counts for all sequences in a DataFrame"""
    # Split sequences into characters
    split_chars_expr = pl.col(seq_col).str.split('').alias('chars')
    
    # Create k-mers by shifting and concatenating # TODO: look if this can be down in one step or if there is some sliding window function.
    create_kmers_expr = pl.concat_str(
        [pl.col('chars').shift(-i).over(id_col) for i in range(k)]
    ).alias('substrings')
    
    # Filter for complete k-mers only
    filter_complete_kmers_expr = pl.col('substrings').str.len_chars() == k
    
    # Aggregate expressions
    agg_exprs = [
        pl.first(seq_col),  # Keep the original sequence
        pl.col('substrings').value_counts(normalize=relative).alias(f'kmer_{k}_relative' if relative else f'kmer_{k}_counts'),
        pl.exclude(seq_col, 'chars', 'substrings').first()  # Keep all other original columns
    ]
    
    return (
        df
        .with_columns(split_chars_expr)
        .explode('chars')
        .with_columns(create_kmers_expr)
        .filter(filter_complete_kmers_expr)
        .group_by(id_col, maintain_order=True)
        .agg(*agg_exprs)
    )

def filter_repetitive_kmers(df: pl.DataFrame, seq_col: str = "seq",id_col: str = "seqid", k: int = 5, max_count: int = 4) -> pl.DataFrame:
    """Filter sequences that have any k-mer appearing more than max_count times"""
    # First get k-mer counts
    df_with_kmers = count_kmers_df(df, seq_col,id_col, k, relative=False)
    
    # Filter for sequences without highly repetitive k-mers
    filter_repetitive_expr = ~pl.col('kmer_counts').list.eval(
        pl.element().struct.field('count') > max_count
    ).list.any()
    
    return df_with_kmers.filter(filter_repetitive_expr)

# lcc_mult and lcc_simp are sourced from the biopython library - we don't need to depend on it here just for these two functions.
def lcc_mult(seq, wsize):
    """Calculate Local Composition Complexity (LCC) values over sliding window.
    sourced from: https://github.com/biopython/biopython/blob/e451db211bdd855a5d0f1f6bba18985ffee12696/Bio/SeqUtils/lcc.py#L13

    Returns a list of floats, the LCC values for a sliding window over
    the sequence.

    seq - an unambiguous DNA sequence (a string or Seq object)
    wsize - window size, integer

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
def lcc_simp(seq):
    """Calculate Local Composition Complexity (LCC) for a sequence.
    sourced from: https://github.com/biopython/biopython/blob/e451db211bdd855a5d0f1f6bba18985ffee12696/Bio/SeqUtils/lcc.py#L120

    seq - an unambiguous DNA sequence (a string or Seq object)

    Returns the Local Composition Complexity (LCC) value for the entire
    sequence (as a float).

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
    return -(term_a + term_c + term_t + term_g)


def reverse_complement(sequence):
    """reverse complement non+ambiguous bases."""
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N",
        "U": "A",
        "R": "Y",
        "Y": "R",
        "W": "W",
        "S": "S",
        "M": "K",
        "K": "M",
        "B": "V",
        "V": "B",
        "D": "H",
        "H": "D",
    }
    return "".join(complement[base] for base in reversed(sequence))


from rust_simulator import Simulator, BaseComposition, DistributionType
# import shlex
import subprocess


def generate_simulation_id(params):
    """Generate a unique simulation ID based on input parameters.

    Args:
        params (dict): Dictionary of simulation parameters

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
    contig_length_range,
    spacer_length_range,
    n_mismatch_range,
    sample_size_contigs,
    sample_size_spacers,
    insertion_range,
    n_insertion_range=(0, 0),
    n_deletion_range=(0, 0),
    contigs=None,
    spacers=None,
    prop_rc=0.5,
    debug=False,
    threads=None,
    verify=False,
    results_dir=None,
    id_prefix=None,
    # New parameters for base composition and distribution
    contig_distribution=None,
    spacer_distribution=None,
    base_composition=None,
    gc_content=None,
    a_frac=None,
    t_frac=None,
    c_frac=None,
    g_frac=None,
    data_subdir=None,
):
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
        )
        
        if verify:
            verrify_df = verify_simulated_data(contigs_dict, spacers_dict, ground_truth_df)
            if (
                verrify_df.filter(pl.col("alignment_test") > pl.col("mismatches")).height
                > 0
            ):
                print(verrify_df.filter(pl.col("alignment_test") > pl.col("mismatches")))
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
    
    # Handle base composition
    composition_obj = None
    if base_composition is not None:
        composition_obj = base_composition
    elif gc_content is not None:
        composition_obj = BaseComposition.from_gc_content(gc_content)
    elif all(x is not None for x in [a_frac, t_frac, c_frac, g_frac]):
        composition_obj = BaseComposition.from_fractions(a_frac, t_frac, c_frac, g_frac)
    
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
        composition_obj,
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
        orient="row"
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
    contigs,
    spacers,
    ground_truth,
    return_fraction=False,
    return_bool=False,
    return_positive_only=False,
):
    """Verify the simulated data by checking the ground truth against the contigs and spacers."""
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


def write_fasta(sequences, filename):
    with open(filename, "w") as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n{seq}\n")


def read_fasta(filename):
    sequences = {}
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip()
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line.strip()
    return sequences


def run_tool(tool, results_dir, debug=False):
    """
    Run a tool either via its bash script (with hyperfine) or directly (debug mode).
    
    Args:
        tool: Tool configuration dictionary
        results_dir: Results directory
        debug: If True, run command directly without hyperfine wrapper for better error messages
    
    Returns:
        Timing data dictionary (empty dict in debug mode)
    """
    if debug:
        # Debug mode: run the actual command directly, not the hyperfine wrapper
        print(f"\n{'='*60}")
        print(f"DEBUG: Running {tool['name']} directly")
        print(f"{'='*60}")
        
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
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            executable="/bin/bash"
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
            raise RuntimeError(f"Tool {tool['name']} failed with exit code {result.returncode}")
        else:
            print(f"\n✓ Command succeeded")
        
        print(f"{'='*60}\n")
        return {}  # No timing data in debug mode
    else:
        # Normal mode: run via bash script with hyperfine
        print(f"Running {tool['script_name']} in {results_dir}/bash_scripts/{tool['script_name']}")
        subprocess.run(f"{results_dir}/bash_scripts/{tool['script_name']}", shell=True)
        
        with open(f"{results_dir}/raw_outputs/{tool['script_name']}.json", "r") as f:
            data = json.load(f)
        return data


def create_bash_script(tool, results_dir, max_runs=1, warmups=0, hyperfine=True):
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


def clean_everything(results_dir):
    os.system(f"rm -rf {results_dir}/raw_outputs/")
    os.system(f"rm -rf {results_dir}/bash_scripts/")
    os.system(f"rm -rf {results_dir}/results/")
    os.system(f"rm -rf {results_dir}/simulated_data/")


def clean_before_rerun(tool_name, results_dir):
    os.system(f"rm -rf {results_dir}/raw_outputs/{tool_name}*tsv")
    os.system(f"rm -rf {results_dir}/raw_outputs/{tool_name}*sam")
    # os.system(f"rm -rf {results_dir}/bash_scripts/{tool_name}*") # generating scripts moved to generate_scripts.py
    os.system(f"rm -rf {results_dir}/raw_outputs/tmp*")


def get_aln_len_from_cigar(cigar):
    """Calculate reference span from CIGAR string.
    Only M, D, N, =, X consume reference positions.
    I, S, H do not consume reference positions."""
    # Match operations that consume reference: M, D, N, =, X
    ref_consuming = re.findall(r"(\d+)[MDN=X]", cigar)
    return sum(int(num) for num in ref_consuming)


def fix_cigar_for_pysamtools(cigar):
    # if any character is not preceded by a digit, append to it 1
    # handle both standard CIGAR operators (A-Z) and extended operators (=,X)
    t = re.sub(r"([A-Z|=|X])(?!\d)([A-Z|=|X])", r"\g<1>1\g<2>", cigar, count=0)
    if re.search(r"[A-Z|=|X][A-Z|=|X]", t):  # Check for consecutive operators
        return fix_cigar_for_pysamtools(t)
    else:
        return t


def get_mismatches_from_cigar(cigar):
    return sum(int(num[:-1]) for num in re.findall(r"\d+X", cigar) or [0])


def parse_samVext(sam_file, max_mismatches=5):
    """Parse SAM file using the version 1.4 cigar strings, where = for matches, X for mismatches, I for insertions, D for deletions."""
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
            }
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


def order_columns_to_match(df1_to_order, df2_to_match):
    return df1_to_order.select(df2_to_match.columns)


def cast_cols_to_match(df1_to_cast, df2_to_match):
    for col in df2_to_match.columns:
        df1_to_cast = df1_to_cast.with_columns(
            pl.col(col).cast(df2_to_match.schema[col])
        )
    return df1_to_cast


def vstack_easy(df1_to_stack, df2_to_stack):
    df2_to_stack = cast_cols_to_match(df2_to_stack, df1_to_stack)
    df2_to_stack = order_columns_to_match(df2_to_stack, df1_to_stack)
    return df1_to_stack.vstack(df2_to_stack)


from pathlib import Path
def parse_sassy(
    sassy_file, 
    max_mismatches=5, 
    spacer_lendf=None, 
    max_gaps=2, 
    threads=10, 
    output_prefix="sassy_parsed",
    output_dir=".",
    **kwargs
):
    """
    Parses Sassy TSV using DuckDB with caching capabilities.
    If a valid Parquet file already exists for the given prefix/date, it is loaded directly.
    """
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
            existing_schema = pl.scan_parquet(parquet_path).schema
            
            # (Optional) Verify essential columns exist to ensure it's not an old version
            if "spacer_id" in existing_schema and "mismatches" in existing_schema:
                print(">> Valid cache found! Skipping DuckDB processing.")
                return pl.read_parquet(parquet_path)
            else:
                print(">> Cache exists but schema is incorrect. Reprocessing...")
        except Exception as e:
            print(f">> Cache file exists but appears corrupt/incomplete ({e}). Reprocessing...")
    
    # =========================================================================
    # START PROCESSING (Only if Cache Miss)
    # =========================================================================
    
    con = duckdb.connect(database=':memory:')
    
    # Hardware Tuning
    con.execute(f"SET threads TO {threads};")
    con.execute("SET preserve_insertion_order = false;") 
    con.execute("SET memory_limit = '50GB';")

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
        con.register('spacer_lookup', spacer_lendf)
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
        
        # Final Read
        print("Reading new Parquet file into Polars...")
        if os.path.exists(parquet_path):
            results = pl.read_parquet(parquet_path)
        else:
            print("Warning: DuckDB finished but no file was created (empty result set?)")
            results = pl.DataFrame(schema={
                "spacer_id": pl.String, "contig_id": pl.String, "spacer_length": pl.UInt32,
                "strand": pl.Boolean, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32,
            })

    except Exception as e:
        print(f"Failed to process Sassy file: {e}")
        # If the write crashed halfway, the file might be corrupt. Clean it up?
        # Typically better to leave it for inspection or manual deletion, 
        # but you can uncomment this if you prefer auto-cleanup:
        # if os.path.exists(parquet_path): os.remove(parquet_path)
        
        return pl.DataFrame(schema={
            "spacer_id": pl.String, "contig_id": pl.String, "spacer_length": pl.UInt32,
            "strand": pl.Boolean, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32,
        })

    return results

# def parse_sassy(sassy_file, max_mismatches=5, spacer_lendf=None, max_gaps=2, threads=10, **kwargs):
#     #checking size
#     # os.path.(sassy_file).
#     con = duckdb.connect(database=':memory:')
    
#     # 1. Hardware Tuning
#     con.execute(f"SET threads TO {threads};")
#     con.execute("SET preserve_insertion_order = false;") # Crucial for parallel writing
#     con.execute("SET memory_limit = '50GB';")
    
#     # Use a single temp file instead of a directory
#     temp_file = f"temp_sassy_{uuid.uuid4()}.parquet"

#     # 2. Logic Definitions
#     cigar_check = "length(cigar) - length(replace(replace(replace(cigar, 'I', ''), 'D', ''), 'N', ''))"
    
#     if max_gaps == 1:
#         gap_logic = f"{cigar_check} = 0" 
#     else:
#         gap_logic = f"""
#             (CASE 
#                 WHEN (cigar LIKE '%I%' OR cigar LIKE '%D%' OR cigar LIKE '%N%') 
#                 THEN {cigar_check} <= {max_gaps}
#                 ELSE TRUE 
#             END)
#         """

#     if spacer_lendf is not None:
#         con.register('spacer_lookup', spacer_lendf)
#         join_step = "INNER JOIN spacer_lookup s ON g.spacer_id = s.spacer_id"
#         len_selection = "s.length AS spacer_length"
#     else:
#         join_step = ""
#         len_selection = "0 AS spacer_length"

#     try:
#         print(f"parsing sassy (duckdb) (Threads: {threads})...")
        
#         query = f"""
#             COPY (
#                 WITH 
#                 -- STEP 1: SCAN & CHEAP FILTER
#                 -- Pushed down to the CSV reader for maximum I/O speed
#                 fast_filter AS (
#                     SELECT 
#                         spacer_id, contig_id, mismatches, strand, "start", "end", cigar
#                     FROM read_csv(
#                         '{sassy_file}', 
#                         sep='\t', 
#                         header=True, 
#                         skip=1,
#                         parallel=true,
#                         auto_detect=false,
#                         hive_partitioning=false,
#                         columns={{
#                             'spacer_id': 'VARCHAR', 'contig_id': 'VARCHAR', 
#                             'mismatches': 'UINTEGER', 'strand': 'VARCHAR', 
#                             'start': 'UINTEGER', 'end': 'UINTEGER',
#                             'slice_str': 'VARCHAR', 'cigar': 'VARCHAR'
#                         }}
#                     )
#                     WHERE mismatches <= {max_mismatches}
#                 ),

#                 -- STEP 2: EXPENSIVE FILTER
#                 -- Runs only on survivors of Step 1
#                 gap_filter AS (
#                     SELECT * FROM fast_filter
#                     WHERE {gap_logic}
#                 )

#                 -- STEP 3: JOIN & SELECT
#                 SELECT 
#                     g.spacer_id,
#                     g.contig_id,
#                     {len_selection},
#                     (g.strand = '-') AS strand,
#                     g."start",
#                     g."end",
#                     g.mismatches
#                 FROM gap_filter g
#                 {join_step}

#             -- Write to a SINGLE Parquet file. 
#             -- DuckDB handles parallel row-group writing automatically.
#             ) TO '{temp_file}' (FORMAT 'PARQUET', COMPRESSION 'SNAPPY', ROW_GROUP_SIZE 100000);
#         """
        
#         con.execute(query)
        
#         print("Reading Parquet back into Polars...")
#         if os.path.exists(temp_file):
#             # scan_parquet is generally faster/lighter than read_parquet for large files
#             results = pl.scan_parquet(temp_file).collect()
#         else:
#             print("Warning: No results generated.")
#             results = pl.DataFrame(schema={
#                 "spacer_id": pl.String, "contig_id": pl.String, "spacer_length": pl.UInt32,
#                 "strand": pl.Boolean, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32,
#             })

#     except Exception as e:
#         print(f"Failed to process Sassy file: {e}")
#         return pl.DataFrame(schema={
#             "spacer_id": pl.String, "contig_id": pl.String, "spacer_length": pl.UInt32,
#             "strand": pl.Boolean, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32,
#         })
#     finally:
#         if os.path.exists(temp_file):
#             os.remove(temp_file)

#     return results

def parse_sassy_polars(sassy_file, max_mismatches=5, spacer_lendf=None, max_gaps=2, **kwargs):
    """Parse Sassy TSV output format and return standardized coordinates.
    Sassy output format:
    pat_id	text_id	cost	strand	start	end	match_region	cigar
    """

    try:
        results = pl.scan_csv(
            sassy_file,
            separator="\t",
            skip_lines=1,
            schema={
                "spacer_id":pl.Utf8,
                "contig_id":pl.Utf8,
                "mismatches":pl.UInt32,
                "strand":pl.Utf8,
                "start":pl.UInt32,
                "end":pl.UInt32,
                "slice_str":pl.Utf8,
                "cigar":pl.Utf8,
            },
        ).drop("slice_str")

        if max_gaps == 1:
                results = results.filter(~pl.col("cigar").str.contains_any(patterns=["I","D","N"]))
        else:
                results = results.filter(~pl.col("cigar").str.count_matches(r"I|D|N")>max_gaps)

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
            }
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


def parse_blastn_custom(blastn_file, max_mismatches=5, spacer_lendf=None, **kwargs):
    """Parse a custom BLAST or mmseqs output format
    (query,target,nident,alnlen,mismatch,qlen,gapopen,qstart,qend,tstart,tend,evalue,bits)
    and return standardized coordinates (1-based).
    BLAST format uses 1-based coordinates.
    Filter out rows with more than max_mismatches mismatches. (i.e. retain up to (including) max_mismatches mismatches)"""
    try:
        results = pl.read_csv(
            blastn_file,
            separator="\t",
            has_header=False,
            infer_schema_length=100000,
            # n_threads=kwargs.get("threads", 1),
            new_columns=[
                "spacer_id",
                "contig_id",
                "nident",
                "alnlen",
                "mismatch",
                "qlen",
                "gapopen",
                "qstart",
                "qend",
                "tstart",
                "tend",
                "evalue",
                "bits",
            ],
        )
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
        #     pl.col("spacer_id").cast(pl.Utf8)
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
            }
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


def parse_lexicmap(tsv_file, max_mismatches=5, spacer_lendf=None, **kwargs):
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
    22. qseq,     Aligned part of query sequence.                     (optional with -a/--all)
    23. sseq,     Aligned part of subject sequence.                   (optional with -a/--all)
    24. align,    Alignment text ("|" and " ") between qseq and sseq. (optional with -a/--all)

    """
    try:
        results = pl.read_csv(
            tsv_file,
            separator="\t",
            has_header=True,
            infer_schema_length=100000,
        )
        results = results.rename({
            "query": "spacer_id",
            "sseqid": "contig_id",
            "sstr": "strand",
        })
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
            }
        )
        
    # results = results.filter(pl.col("alenHSP") >= 17, pl.col("gaps") == 0) # we don't really need this here, but keeping commented to remember it was done previously. Potentially, with a max mismatch >3 this could have over restrict the max mismatch arguments values.
    results = results.with_columns(pl.col("spacer_id").cast(pl.Utf8))
    results = results.with_columns(
        pl.col("align").str.count_matches("|", literal=True).alias("matches")
    )

    results = spacer_lendf.join(results, on="spacer_id", how="inner")
    results = results.with_columns(
        (pl.col("length") - pl.col("matches")).alias("mismatches"),
        (pl.col("sstart") - 1).alias("start"), # 1-based
        (pl.col("send")).alias("end"), # 1-based
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


def parse_samVn_with_lens_polars(
    sam_file, spacer_lendf, max_mismatches=5
):
    """Parse SAM file using spacer lengths for qlen and numerics in CIGAR for alignment length
    uses polars to speed up things"""
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
        #     pl.col("tags").str.split( "\t").list.gather_every(n=1,offset=11).alias("tags")
        # )
        # tags_df = tags_df.with_columns(
        #     pl.col("tags").list.len().alias("tag_lengths")
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
            }
        )

    # Continue with existing processing
    sam_df = sam_df.with_columns((pl.col("start") - 1).alias("start")) # sam is 1-based
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
    sam_file, spacer_lendf, max_mismatches=5, threads=4, ref_file=None, gaps_as_mismatches=False, **kwargs
):
    """Parse SAM file using pysam and spacer lengths to compute mismatches
    RENAMED FROM parse_samVn_with_lens_pysam """
    
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
            }
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
            #                  if op == 5)  # H operation
            # # if hard_clipped > 0:
            # #     print(f"Hard clipped {hard_clipped} for {read.query_name}")
            # #     # break
            # if spacer_len == 0: # some alignments get 0 length, maybe because some tools do not output the sequence for secondary/ambiguous alignments. I suspect pysam get's the query length from the nchar of the seq.
            #     # print (f"spacer_len == 0 for {read.query_name}, extracting alignment info from cigar")
            # # Get spacer length from reference table instead of query sequence
            #     spacer_len = spacer_lens.get(read.query_name)
            # else:
            #     spacer_len = spacer_len + hard_clipped
            #     # print(f"{spacer_len}")
            #     # break
            if spacer_len is None:
                print(
                    f"Warning: spacer {read.query_name} not found in reference table"
                )
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
            #                   if op in [0, 7, 8])  # M, =, X operations

            # # # Get soft clipped length
            soft_clipped = sum(
                length for op, length in read.cigartuples if op == 4
            )  # "S" operation


            #debug
            # if soft_clipped > 0:
            #     print(f"Soft clipped {soft_clipped} for {read.query_name}")
            #     break
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
        },orient="row"
    ).unique()


def plot_matrix(
    matrix,
    title,
    filename,
    hor_axis="pairs in tool",
    vert_axis="pairs not in tool",
    return_chart=True,
):
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

    # Convert to percentages if requested
    # if as_percent:
    #     # Get total matches for each tool
    #     tool_totals = matrix.group_by('in_tool').agg(
    #         pl.len().alias('total_matches')
    #     ).to_dict(as_series=False)
    #     tool_totals = dict(zip(tool_totals['in_tool'], tool_totals['total_matches']))

    #     # Convert counts to percentages
    #     for tool in sorted_tools:
    #         if tool != 'in_tool':
    #             matrix = matrix.with_columns(
    #                 (pl.col(tool) * 100 / tool_totals[tool]).alias(tool)
    #             )

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
    output_file, max_mismatches=5, spacer_lendf=None, **kwargs
):  # expecting tab file with columns: spacer_id, contig_id, start, end, strand - because Antonio insists.
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
            }
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


def parse_hyperfine_output(json_file):
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
            "stddev_time": result["stddev"] if result["stddev"] is not None else 0.0, # need to debug this
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


def read_hyperfine_results(tools, results_dir):
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
        }
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


def run_tools(tools, results_dir, debug=False):
    """
    Run multiple tools.
    
    Args:
        tools: Dictionary of tool configurations
        results_dir: Results directory
        debug: If True, run commands directly without hyperfine for better error messages
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
    seqfile,
    seq_ids,
    return_dict=False,
    return_df=True,
    idcol="contig_id",
    seqcol="contig_seq",
    output_file=None,
):
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


def soft_fetch_fastx(fa, id, start=None, end=None, strand=None):
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
    seqfile,
    pldf: pl.DataFrame,
    trim_to_region=False,
    idcol="contig_id",
    start_col="start",
    end_col="end",
    strand_col="strand",
    seqcol="contig_seq",
    drop_missing=True,
    **kwargs,
):
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
    #     )
    # time_end = time.time() # 17.89 seconds (no threading) 16.66 seconds (threading)
    # print(f"Time taken: {time_end - time_start:.2f} seconds")

    tmpdf = tmpdf.drop(["strand4bed", "start_0based"])

    if drop_missing:
        tmpdf = tmpdf.filter(~pl.col(seqcol).is_in([None, "", "nan"]))
    pldf = pldf.join(tmpdf, on=[idcol, start_col, end_col, strand_col], how="left")
    return pldf


# # # test
# one_mismatch = pl.read_parquet('results/real_data/results/one_mismatch.parquet')
# seqfile = 'results/real_data/results/matched_contigs.fna'
# pldf = one_mismatch[range(10000)]
# idcol = "contig_id"
# seqcol = "contig_seq"
# start_col = "start"
# end_col = "end"
# strand_col = "strand"
# test = populate_pldf_from_fastx (seqfile=seqfile, pldf=pldf, idcol=idcol, start_col=start_col, end_col=end_col, strand_col=strand_col, seqcol=seqcol)


def populate_pldf_withseqs(
    pldf: pl.DataFrame,
    seqfile,
    chunk_size=2000,
    drop_missing=True,
    trim_to_region=False,
    reverse_by_strand_col=False,
    idcol="contig_id",
    seqcol="contig_seq",
    strand_col="strand",
    **kwargs,
):
    # Note! if a "strand" column is present, and reverse_by_strand_col is True, it will be used to get the reverse complement of the sequence.
    # If drop_missing is True, it will drop any contig_ids that are not present in the sequence file. If no seqids were found, it will return the input dataframe.
    # NOTE: trim_to_region and reverse_by_strand_col are bugged when using slice, even though it's probably much faster.
    # get the unique contig_ids and strands
    stranded = strand_col in pldf.columns and reverse_by_strand_col
    # if stranded:
    #     nr_contigids = pldf[[idcol,strand_col]].unique()
    #     # print("rev true")
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
        # print(f"Time taken: {time_end - time_start:.2f} seconds") # Time taken: 0.15 seconds
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

    # debug
    # test = tools_results[range(100000)]
    # test2= populate_pldf_withseqs(seqfile=seqfile,chunk_size=12000, pldf=test, drop_missing=True, reverse_by_strand_col=False)
    # pldf=tools_results[range(10000)]
    # and add a few a null, empty string, or nan
    # pldf_null = pldf[range(4)].with_columns(
    #     pl.lit(None).alias("contig_id"),
    # )
    # pldf_empty = pldf[range(4)].with_columns(
    #     pl.lit("").alias("contig_id"),
    # )
    # pldf = pldf.vstack(pldf_null).vstack(pldf_empty)

    # pldf = pldf.vstack()


def populate_pldf_withseqs_needletail(
    pldf: pl.DataFrame,
    seqfile,
    chunk_size=20000000,
    trim_to_region=True,
    reverse_by_strand_col=True,
    idcol="contig_id",
    seqcol="contig_seq",
    start_col="start",
    end_col="end",
    strand_col="strand",
):
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
    #     seq_count += 1
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
                chunk_seqs = chunk_seqs.with_columns(
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
                    pl.when(pl.col(strand_col))
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


def test_alignment_from_faidx(spacers_file, contigs_file, alignment, **kwargs):
    """Test if an alignment matches between spacer and contig files using pyfastx.
    alignment is a single row extracted from a pl.df results  (spacer_id, contig_id, spacer_length, start, end, strand, mismatches, tool)
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
    spacer_seq,
    contig_seq,
    strand=False,
    start=None,
    end=None,
    return_ref_region=False,
    return_query_region=False,
    return_comp_str=False,
    gap_cost=10,
    extend_cost=5,
    cost_matrix=ps.nuc44,
    **kwargs,
):
    """given a spacer and contig sequence and coordinates, return | for a match and . for a mismatch."""
    if start is not None and end is not None:
        aligned_region = contig_seq[start:end]
    else:
        aligned_region = contig_seq

    if strand:  # True means reverse strand
        aligned_region = reverse_complement(aligned_region)
    # alignment = ps.nw_stats_scan(spacer_seq,aligned_region,open=10,extend=10,matrix=ps.nuc44)
    # alignment_string = ""
    # for i in range(len(spacer_seq)):
    #     if spacer_seq[i] == aligned_region[i]:
    #         alignment_string += "|"
    #     else:
    #         alignment_string += "."
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
    spacer_seq,
    contig_seq,
    strand=False,
    start=None,
    end=None,
    gap_cost=10,
    extend_cost=5,
    gaps_as_mismatch=True,
    **kwargs,
):
    """Test if an alignment matches between a spacer and contig using Needleman-Wunsch algorithm.
    Input is a spacer and contig sequence and coordinates, and optionally strand, start, end, gap_cost, extend_cost, and gaps_as_mismatch.
    Returns the number of mismatches.
    
    Args:
        gaps_as_mismatch: If True, count gaps as mismatches (edit distance).
                         If False, count only substitutions (hamming-like distance).
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


def read_results(
    tools,
    max_mismatches=5,
    spacer_lendf=None,
    ref_file=None,
    threads=4,
):
    results_df = pl.DataFrame(
        schema={
            "spacer_id": pl.Utf8,
            "contig_id": pl.Utf8,
            "spacer_length": pl.UInt32,
            "strand": pl.Boolean,
            "start": pl.UInt32,
            "end": pl.UInt32,
            "mismatches": pl.UInt32,
            "tool": pl.Utf8,
        }
    )
    for tool in tools.values():
        try:
            print(f"\nReading results for {tool['name']}...")
            parse_function = tool.get("parse_function")
            # if the file exists, get the size of it
            if os.path.exists(tool["output_file"]):
                file_size = os.path.getsize(tool["output_file"])
                print(f"File size: {file_size / 1024 / 1024:.2f} MB")
                if file_size == 0:
                    print(f"File {tool['output_file']} is empty, skipping")
                    continue
            else:
                print(f"File {tool['output_file']} does not exist, skipping")
                continue
            print(" output file ", tool["output_file"], "exists")
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
            parsed_results = parsed_results.with_columns(
                pl.lit(tool["name"]).alias("tool")
            )
            tool_df = pl.DataFrame(
                data=parsed_results,
                orient="row",
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
            )
            # print(parsed_results)
            results_df = results_df.vstack(tool_df)
        except Exception as e:
            print(f"Failed to read results for {tool['name']}: {e}")
            ipdb.set_trace()
    return results_df

def verify_sam_file(sam_file, ref_file=None): 
    """Verify if a SAM file is valid and add SQ lines if needed (from reference file if provided or from other sam files in the same directory)
    also potentially replaces whitespaces with tabs in the first line"""

    # Read the first few lines to analyze the header
    with open(sam_file, 'r') as f:
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
    needs_fixing = hd_needs_fixing or not has_sq_lines or not has_pg_lines or is_mummer_file
    
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
                    fixed_hd = first_line.replace(" VN:", "\tVN:").replace(" SO:", "\tSO:")
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
                print("No reference file provided, trying to find SQ lines from other SAM files")
                present_sam_files = glob.glob(os.path.join(os.path.dirname(sam_file), "*.sam"))
                present_sam_files = [f for f in present_sam_files if f != os.path.normpath(sam_file)]
                # Drop mummer4 and minimap2 sam files as they don't print SQ lines when batching
                present_sam_files = [f for f in present_sam_files if "mummer" not in f and "minimap" not in f]
                present_sam_files = [f for f in present_sam_files if os.path.getsize(f) > 1]
                
                if len(present_sam_files) == 0:
                    print("No non-empty sam files found, returning None")
                    return None
                
                # Find a SAM file with SQ lines
                sq_found = False
                for sam_file_candidate in present_sam_files:
                    with open(sam_file_candidate, 'r') as f:
                        for line in f:
                            if line.startswith("@SQ"):
                                new_sam_file.write(line)
                                sq_found = True
                            elif line.startswith("@PG") or (not line.startswith("@") and line.strip()):
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
                        fixed_pg = pg_line.replace(" ID:", "\tID:").replace(" PN:", "\tPN:").replace(" VN:", "\tVN:").replace(" CL:", "\tCL:")
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

def prefilter_sam_with_sambamba(input_sam, output_sam, max_mismatches=5, threads=1):
    """Prefilter SAM file using sambamba to remove unmapped reads with mismatches (NM) <= max_mismatches"""
    # Filter expression:
    # - not unmapped: keep only mapped reads
    # - [NM] <= max_mismatches: keep reads with max_mismatches or fewer mismatches
    # - cigar !~ /[ID]/: exclude reads with insertions or deletions in CIGAR string # UPDATE: never used.

    filter_expression = f'"not (unmapped) and [NM] <= {max_mismatches}"'
    if not os.path.exists(input_sam):
        print(f"File {input_sam} does not exist, skipping")
        return None
    if os.path.getsize(input_sam) == 0:
        print(f"File {input_sam} is empty, skipping")
        return None

    input_sam = verify_sam_file(input_sam)

    command = [
        "sambamba",
        "view",
        "-F",
        filter_expression,
        "--show-progress",
        "-t",
        str(threads),
        "-S",  # SAM input
        "-f",
        "sam",  # SAM output
        "-h",
        "-o",
        output_sam,
        input_sam,
    ]
    subprocess.run(" ".join(command), check=True, shell=True)
    return output_sam


def create_comparison_matrix(tools_results, n_mismatches):
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




def get_parse_function(func_name):
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



def validate_intervals_with_polars_bio(planned_intervals, tool_results, max_mismatches=5):
    """
    Validate tool results against planned intervals using polars-bio interval operations.
    
    Args:
        planned_intervals: DataFrame with planned spacer insertions (spacer_id, contig_id, start, end, strand, mismatches)
    tool_results: DataFrame with tool results (spacer_id, contig_id, start, end, strand, mismatches)
        max_mismatches: Maximum allowed mismatches for validation
    
    Returns:
        DataFrame with validation results
    """
    # Filter tool results by max mismatches
    tool_results_filtered = tool_results.filter(pl.col("mismatches") <= max_mismatches)
    
    # Add chrom column for polars-bio compatibility (use contig_id as chrom)
    planned_df = planned_intervals.with_columns([
        pl.col("start").cast(pl.UInt32),
        pl.col("end").cast(pl.UInt32),
        pl.col("strand").cast(pl.Utf8)
    ])
    
    tool_df = tool_results_filtered.with_columns([
        pl.col("start").cast(pl.UInt32),
        pl.col("end").cast(pl.UInt32),
        pl.col("strand").cast(pl.Utf8)
    ])
    
    # Check if dataframes are empty
    if planned_df.height == 0:
        print("Warning: planned_df is empty")
        return pl.DataFrame(schema={
            "spacer_id_1": pl.Utf8, "contig_id_1": pl.Utf8, "start_1": pl.UInt32, "end_1": pl.UInt32, "strand_1": pl.Utf8,
            "spacer_id_2": pl.Utf8, "contig_id_2": pl.Utf8, "start_2": pl.UInt32, "end_2": pl.UInt32, "strand_2": pl.Utf8,
            "overlap_type": pl.Utf8
        })
    if tool_df.height == 0:
        print("Warning: tool_df is empty")
        return pl.DataFrame(schema={
            "spacer_id_1": pl.Utf8, "contig_id_1": pl.Utf8, "start_1": pl.UInt32, "end_1": pl.UInt32, "strand_1": pl.Utf8,
            "spacer_id_2": pl.Utf8, "contig_id_2": pl.Utf8, "start_2": pl.UInt32, "end_2": pl.UInt32, "strand_2": pl.Utf8,
            "overlap_type": pl.Utf8
        })
    
    # Use polars-bio overlap operation
    try:
        overlap_result = pb.overlap(planned_df, tool_df,cols1=["contig_id","start", "end"],cols2=["contig_id","start", "end"])
    except Exception as e:
        print(f"Error in pb.overlap: {e}")
        print(f"planned_df shape: {planned_df.shape}")
        print(f"tool_df shape: {tool_df.shape}")
        raise
    
    # Calculate validation metrics with strand consideration
    validation_results = overlap_result.with_columns([
        # First check if strands match (both forward or both reverse)
        pl.when(pl.col("strand_1") == pl.col("strand_2"))
        .then(
            # If strands match, check for coordinate overlaps
            # Use symmetric overlap detection logic
            pl.when(
                # Exact match: intervals are identical
                (pl.col("start_1") == pl.col("start_2")) & (pl.col("end_1") == pl.col("end_2"))
            )
            .then(pl.lit("exact_match"))
            .when(
                # Partial overlap: intervals overlap but are not identical
                (pl.col("start_1") < pl.col("end_2")) & (pl.col("start_2") < pl.col("end_1"))
            )
            .then(pl.lit("partial_overlap"))
            .otherwise(pl.lit("no_overlap"))
        )
        .otherwise(pl.lit("strand_mismatch"))  # Different strands = no valid match
        .alias("overlap_type")
    ]).collect()
    
    return validation_results


def validate_intervals_fallback(planned_intervals, tool_results, max_mismatches=5):
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
        suffix="_tool"
    )
    
    # Check for overlaps
    overlap_conditions = (
        (pl.col("start") <= pl.col("start_tool")) & 
        (pl.col("end") >= pl.col("end_tool"))
    ) | (
        (pl.col("start_tool") <= pl.col("start")) & 
        (pl.col("end_tool") >= pl.col("end"))
    ) | (
        (pl.col("start") < pl.col("end_tool")) & 
        (pl.col("end") > pl.col("start_tool"))
    )
    
    validation_results = joined.with_columns([
        pl.when(overlap_conditions)
        .then(pl.lit("overlap"))
        .otherwise(pl.lit("no_overlap"))
        .alias("overlap_type")
    ])
    
    return validation_results


def detect_spurious_alignments(contigs, spacers, planned_intervals, max_mismatches=5, 
                              alignment_threshold=0.8):
    """
    Detect spurious (non-planned) alignments in contigs that match spacer sequences.
    
    Args:
        contigs: Dictionary of contig_id -> sequence
        spacers: Dictionary of spacer_id -> sequence  
        planned_intervals: DataFrame with planned insertions
        max_mismatches: Maximum mismatches to consider
        alignment_threshold: Minimum alignment score threshold
    
    Returns:
        DataFrame with detected spurious alignments
    """
    spurious_alignments = []
    
    # Create a set of planned positions for quick lookup
    planned_positions = set()
    for row in planned_intervals.iter_rows(named=True):
        key = (row['spacer_id'], row['contig_id'], row['start'], row['end'], row['strand'])
        planned_positions.add(key)
    
    # Check each spacer against each contig
    for spacer_id, spacer_seq in spacers.items():
        for contig_id, contig_seq in contigs.items():
            # Perform alignment using parasail
            result = ps.sg_trace_scan_sat(
                spacer_seq.encode('utf-8'),
                contig_seq.encode('utf-8'),
                1,  # gap open penalty (must be >= 0)
                1,  # gap extend penalty (must be >= 0)
                ps.nuc44  # substitution matrix
            )
            
            # Check if alignment meets criteria
            if result.score >= len(spacer_seq) * alignment_threshold:
                # Find all alignment positions
                for i in range(len(contig_seq) - len(spacer_seq) + 1):
                    region = contig_seq[i:i + len(spacer_seq)]
                    
                    # Calculate mismatches
                    mismatches = sum(1 for a, b in zip(spacer_seq, region) if a != b)
                    
                    if mismatches <= max_mismatches:
                        # Check if this position was planned
                        key = (spacer_id, contig_id, i, i + len(spacer_seq), False)
                        if key not in planned_positions:
                            spurious_alignments.append({
                                'spacer_id': spacer_id,
                                'contig_id': contig_id,
                                'start': i,
                                'end': i + len(spacer_seq),
                                'strand': False,
                                'mismatches': mismatches,
                                'type': 'spurious'
                            })
    
    return pl.DataFrame(spurious_alignments) if spurious_alignments else pl.DataFrame()


def _calculate_match_prob_with_base_composition(base_comp):
    """
    Calculate probability of a random match between two bases given base composition.
    
    Args:
        base_comp: BaseComposition object or dict with a_frac, t_frac, c_frac, g_frac
    
    Returns:
        Probability that two random bases match
    """
    if isinstance(base_comp, dict):
        a_frac = base_comp.get('a_frac', 0.25)
        t_frac = base_comp.get('t_frac', 0.25)
        c_frac = base_comp.get('c_frac', 0.25)
        g_frac = base_comp.get('g_frac', 0.25)
    else:
        # Assume BaseComposition object
        a_frac = base_comp.a_frac if hasattr(base_comp, 'a_frac') else 0.25
        t_frac = base_comp.t_frac if hasattr(base_comp, 't_frac') else 0.25
        c_frac = base_comp.c_frac if hasattr(base_comp, 'c_frac') else 0.25
        g_frac = base_comp.g_frac if hasattr(base_comp, 'g_frac') else 0.25
    
    # Probability of match = sum of P(both are A) + P(both are T) + ...
    match_prob = (a_frac ** 2 + t_frac ** 2 + c_frac ** 2 + g_frac ** 2)
    return match_prob

def _calculate_edit_distance_probability_hamming(spacer_len, max_edits, match_prob, total_db_size):
    """
    Calculate probability of finding a random string within Hamming distance.
    This gives the probability per position in the database.
    
    Args:
        spacer_len: Length of spacer sequence
        max_edits: Maximum edit distance allowed
        match_prob: Probability that two random bases match
        total_db_size: Total size of sequence database
    
    Returns:
        Probability of finding match at any given position
    """
    from scipy.special import comb
    
    mismatch_prob = 1 - match_prob
    total_prob = 0
    
    # Sum probabilities for 0 to max_edits mismatches
    for k in range(max_edits + 1):
        # Binomial probability: choose k positions to mismatch
        prob_k_mismatches = comb(spacer_len, k, exact=True) * \
                           (mismatch_prob ** k) * \
                           (match_prob ** (spacer_len - k))
        total_prob += prob_k_mismatches
    
    # Multiply by total_db_size to get expected number of matches
    return total_prob * total_db_size


def _calculate_blast_evalue(spacer_len, max_edits, match_prob, total_db_size,
                            match_reward=1, mismatch_penalty=-3,  # Standard BLASTN scoring
                            k_param=0.14, lambda_param=1.28):  # Lambda for +1/-3 scoring
    """
    Calculate expected spurious alignments using BLAST-style E-value approach.
    
    This follows BLAST's statistical theory more closely:
    1. Uses bit scores instead of raw scores
    2. Applies proper scaling factors
    3. Uses standard BLASTN scoring (+1/-3)
    
    Args:
        spacer_len: Length of spacer sequence
        max_edits: Maximum edit distance allowed
        match_prob: Probability that two random bases match (not used in BLAST calculation)
        total_db_size: Total size of sequence database in bases
        match_reward: Score for matching base (default: +1 like BLASTN)
        mismatch_penalty: Penalty for mismatching base (default: -3 like BLASTN)
        k_param: BLAST K parameter (default 0.14 for DNA)
        lambda_param: BLAST lambda parameter (default 1.28 for +1/-3 scoring)
    
    Returns:
        Expected number of spurious alignments (E-value)
    """
    # Calculate minimum score threshold
    min_matches = spacer_len - max_edits
    min_score = (min_matches * match_reward) + (max_edits * mismatch_penalty)
    
    # Convert to bit score (BLAST's statistical framework)
    # S_bit = (lambda * S - ln(K)) / ln(2)
    bit_score = (lambda_param * min_score - np.log(k_param)) / np.log(2)
    
    # Apply length-based corrections
    # These help account for edge effects and multiple testing
    H = 0.5772156649  # Euler's constant
    m_prime = spacer_len / (spacer_len - H)  # Effective query length
    n_prime = total_db_size / (spacer_len - H)  # Effective db size
    
    # Calculate E-value using bit score formulation
    # E = m'n' * 2^(-S_bit)
    e_value = m_prime * n_prime * np.power(2, -bit_score)
    
    return e_value


def calculate_expected_spurious_alignments(spacer_len, max_edits, total_db_size,
                                          base_composition=None, method='both', blast_evalue=False):
    """
    Calculate expected number of spurious alignments using different methods.
    
    Args:
        spacer_len: Length of spacer sequences
        max_edits: Maximum edit distance/mismatches allowed
        total_db_size: Total size of sequence database in bases
        base_composition: BaseComposition object or dict with base frequencies
        method: 'hamming', 'blast', or 'both'
    
    Returns:
        Dictionary with expected spurious alignment statistics for both Hamming and Edit distance
    """
    # Calculate match probability
    if base_composition is None:
        match_prob = 0.25  # Uniform distribution
    else:
        match_prob = _calculate_match_prob_with_base_composition(base_composition)
    
    results = {'blast_evalue': blast_evalue}
    
    if method in ['hamming', 'both'] and not blast_evalue:
        # Hamming distance (substitutions only) - use ungapped parameters
        expected_matches_hamming = _calculate_edit_distance_probability_hamming(
            spacer_len, max_edits, match_prob, total_db_size
        )
        results['expected_spurious_hamming'] = expected_matches_hamming
    
    if method in ['blast', 'both'] and not blast_evalue:
        # BLAST E-value approach
        expected_blast = _calculate_blast_evalue(
            spacer_len, max_edits, match_prob, total_db_size
        )
        results['expected_spurious_blast_evalue'] = expected_blast
    
    if blast_evalue:
        results['expected_spurious_blast_evalue'] = _calculate_blast_evalue(
            spacer_len, max_edits, match_prob, total_db_size
        )
    return results


def demonstrate_blast_based_spurious_estimation(contigs, spacers, max_edits_range=(2, 10), 
                                              base_composition=None, verbose=True):
    """
    Demonstrate BLAST-based spurious alignment estimation with educational output.
    
    This function shows how BLAST E-values work for different distance metrics
    and parameters, making it easier to understand the relationship between
    edit distance and spurious alignments.
    
    Args:
        contigs: Dictionary of contig_id -> sequence
        spacers: Dictionary of spacer_id -> sequence  
        max_edits_range: Tuple of (min, max) edit distance thresholds to test
        base_composition: Base composition dict or None for uniform
        verbose: Whether to print detailed output
    
    Returns:
        Dictionary with results for further analysis
    """
    if verbose:
        print("=== BLAST-based Spurious Alignment Demonstration ===")
        print("This shows how E-values change with different distance metrics")
        print()
    
    # Calculate total search space
    total_contig_length = sum(len(seq) for seq in contigs.values())
    total_spacer_length = sum(len(seq) for seq in spacers.values())
    avg_spacer_length = total_spacer_length / len(spacers)
    
    if verbose:
        print(f"Search space: {total_contig_length:,} bp in {len(contigs)} contigs")
        print(f"Query space: {len(spacers)} spacers, avg length {avg_spacer_length:.1f} bp")
        print()
    
    results = {
        'search_space': total_contig_length,
        'num_contigs': len(contigs),
        'num_spacers': len(spacers),
        'avg_spacer_length': avg_spacer_length,
        'hamming_results': [],
        'edit_results': [],
        'comparisons': []
    }
    
    # Test different edit distance thresholds
    for max_edits in range(max_edits_range[0], max_edits_range[1] + 1):
        # Hamming distance (substitutions only)
        expected_matches_hamming = _calculate_edit_distance_probability_hamming(
            avg_spacer_length, max_edits, 0.25, total_contig_length
        )
        
        # BLAST E-value calculation
        expected_matches_edit = _calculate_edit_distance_probability_with_indels(
            avg_spacer_length, max_edits, 0.25, total_contig_length,
            k_param=0.14, lambda_param=0.69  # These parameters only apply to BLAST E-value calculation
        )
        
        results['hamming_results'].append({
            'max_edits': max_edits,
            'mean_spurious': expected_matches_hamming,
            'std_spurious': expected_matches_hamming ** 0.5
        })
        
        results['edit_results'].append({
            'max_edits': max_edits,
            'mean_spurious': expected_matches_edit,
            'std_spurious': expected_matches_edit ** 0.5
        })
        
        # Calculate ratio
        ratio = (expected_matches_edit / expected_matches_hamming 
                if expected_matches_hamming > 0 else 0)
        
        results['comparisons'].append({
            'max_edits': max_edits,
            'hamming_evalue': expected_matches_hamming,
            'edit_evalue': expected_matches_edit,
            'ratio': ratio
        })
        
        if verbose:
            print(f"Max edits {max_edits:2d}: Hamming={expected_matches_hamming:.2e}, "
                  f"Edit={expected_matches_edit:.2e}, Ratio={ratio:.1f}x")
    
    if verbose:
        print()
        print("Key Insights:")
        print("• BLAST E-values are extremely conservative (very small numbers)")
        print("• Edit distance can be MORE restrictive than Hamming distance")
        print("• Gap penalties make indels expensive in BLAST scoring")
        print("• The relationship depends on edit distance to sequence length ratio")
        print()
    
    return results


def _calculate_score_range_for_edit_distance(length, max_edits, match_reward=2, mismatch_penalty=-3, 
                                           gap_open_penalty=-5, gap_extend_penalty=-2):
    """
    Calculate the range of possible raw scores for a given edit distance threshold.
    
    This simulates all possible combinations of substitutions, insertions, and deletions
    that result in edit distance <= max_edits.
    
    Args:
        length: Length of the query sequence
        max_edits: Maximum allowed edit distance
        match_reward: Score for matching characters
        mismatch_penalty: Score for mismatching base
        gap_open_penalty: Penalty for opening a gap
        gap_extend_penalty: Penalty for extending a gap
    
    Returns:
        Tuple of (min_score, max_score) for the given edit distance threshold
    """
    scores = []
    
    # Enumerate all valid combinations of substitutions, insertions, and deletions
    for total_ops in range(max_edits + 1):
        for subs in range(total_ops + 1):
            for ins in range(total_ops - subs + 1):
                dels = total_ops - subs - ins
                
                # Skip impossible combinations
                if dels > length:
                    continue
                
                # Calculate aligned length after operations
                aligned_len_source = length
                aligned_len_target = length - dels + ins
                
                if aligned_len_source <= 0 or aligned_len_target <= 0:
                    continue
                
                # Calculate raw score for this combination
                # Matches: positions that are neither substituted nor deleted
                matches = aligned_len_source - subs - dels
                score = (
                    matches * match_reward +           # matches
                    subs * mismatch_penalty +          # substitutions
                    (ins + dels) * gap_open_penalty +  # gap open penalty for each indel
                    (ins + dels) * gap_extend_penalty  # gap extend penalty
                )
                scores.append(score)
    
    return (min(scores), max(scores)) if scores else (0, 0)

def _calculate_edit_distance_probability_with_indels(spacer_len, max_edits, match_prob, 
                                                   total_db_size,  # Total database size in bases
                                                   match_reward=1, mismatch_penalty=-3,  # Changed to +1/-3 like BLASTN default
                                                   gap_open_penalty=-5, gap_extend_penalty=-2,
                                                   k_param=0.14, lambda_param=1.28):  # Lambda adjusted for +1/-3 scoring
    """
    Calculate expected spurious alignments using BLAST-style E-value approach for edit distance.
    
    This uses the Karlin-Altschul statistics to estimate the expected number of
    spurious alignments for edit distance-based matching (with indels).
    
    Args:
        spacer_len: Length of spacer
        max_edits: Maximum allowed edit distance
        match_prob: Probability that two random bases match
        total_db_size: Total size of sequence database in bases
        match_reward: Score for matching characters
        mismatch_penalty: Penalty for mismatching base
        gap_open_penalty: Penalty for opening a gap
        gap_extend_penalty: Penalty for extending a gap
        k_param: Karlin-Altschul K parameter for gapped alignments
        lambda_param: Karlin-Altschul lambda parameter for gapped alignments
    
    Returns:
        Expected number of spurious alignments for this spacer across the entire database
    """
    # Calculate the range of possible scores for this edit distance threshold
    min_score, max_score = _calculate_score_range_for_edit_distance(
        spacer_len, max_edits, match_reward, mismatch_penalty, 
        gap_open_penalty, gap_extend_penalty
    )
    
    # Convert to bit score (like BLAST does)
    # S' = (lambda * S - ln(K)) / ln(2)
    # where S is the raw score
    bit_score = (lambda_param * min_score - np.log(k_param)) / np.log(2)
    
    # Calculate effective lengths (simplified from BLAST's approach)
    # In BLAST, these account for edge effects and statistical adjustments
    effective_query_len = max(1, spacer_len - 11)  # BLAST-like adjustment
    effective_db_size = max(1, total_db_size - (11 * (total_db_size // spacer_len)))
    
    # Calculate E-value using bit score and effective lengths
    # E = mn * 2^(-S')
    # where m and n are effective lengths
    e_value = (effective_query_len * effective_db_size) * np.power(2, -bit_score)
    
    return e_value

def estimate_expected_spurious_alignments_fast(
    contigs, 
    spacers, 
    max_mismatches=5, 
    max_edit_distance=None,
    base_composition=None,
    use_edit_distance=False
):
    """
    Fast statistical estimate of expected spurious alignments using two approaches:
    1. Hamming distance probability calculation (theoretical number of sequences within Hamming distance)
    2. BLAST E-value calculation (statistically expected number of alignments by chance)
    
    Note: BLAST E-values are typically much lower than Hamming-based estimates because:
    - They account for gap penalties and scoring matrices
    - They are calibrated for biological relevance
    - They are intentionally conservative to minimize false positives
    
    Args:
        contigs: Dictionary of contig_id -> sequence (or just contig count and total length)
        spacers: Dictionary of spacer_id -> sequence
        max_mismatches: Maximum mismatches to consider (for Hamming distance)
        max_edit_distance: Maximum edit distance to consider (for edit distance with indels)
                          If None and use_edit_distance=True, uses max_mismatches
        base_composition: BaseComposition object or dict with a_frac, t_frac, c_frac, g_frac
                         If None, assumes uniform (0.25 each)
        use_edit_distance: If True, the 'mean_spurious' field uses edit distance
                          If False, the 'mean_spurious' field uses Hamming distance
                          (Both are always calculated and returned)
    
    Returns:
        Dictionary with both Hamming-based and BLAST-based estimates, including:
        - expected_spurious_hamming: Theoretical count of sequences within Hamming distance
        - blast_evalue: BLAST E-value (expected number of chance alignments)
        - std_spurious_hamming: Standard deviation of Hamming-based estimate
        - std_spurious_edit: Standard deviation of BLAST-based estimate
        - max_spurious_hamming: Upper bound for Hamming-based estimate (mean + 2*std)
        - min_spurious_hamming: Lower bound for Hamming-based estimate (mean - 2*std)
        - max_blast_evalue: Upper bound for BLAST-based estimate (mean + 2*std)
        - min_blast_evalue: Lower bound for BLAST-based estimate (mean - 2*std)
        - per_spacer_expected_hamming: List of Hamming-based estimates per spacer
        - per_spacer_blast_evalue: List of BLAST-based estimates per spacer
        - method: Description of the methods used
        - base_composition: Base composition used for the calculation
        - base_comp_note: Note on the base composition uniformity
        - distance_type: Type of distance used for the estimate (hamming or blast_evalue)
        - distance_note: Description of the distance metric used
    """
    # Handle base composition
    if base_composition is None:
        base_comp_dict = {'a_frac': 0.25, 't_frac': 0.25, 'c_frac': 0.25, 'g_frac': 0.25}
    elif isinstance(base_composition, dict):
        base_comp_dict = base_composition
    elif hasattr(base_composition, 'a_frac'):
        base_comp_dict = {
            'a_frac': base_composition.a_frac,
            't_frac': base_composition.t_frac,
            'c_frac': base_composition.c_frac,
            'g_frac': base_composition.g_frac
        }
    else:
        base_comp_dict = {'a_frac': 0.25, 't_frac': 0.25, 'c_frac': 0.25, 'g_frac': 0.25}
    
    match_prob = _calculate_match_prob_with_base_composition(base_comp_dict)
    
    # Calculate total search space
    if isinstance(contigs, dict):
        total_contig_length = sum(len(seq) for seq in contigs.values())
        num_contigs = len(contigs)
    else:
        # If just numbers provided
        total_contig_length = contigs
        num_contigs = 1
    
    # Calculate spacer statistics
    if isinstance(spacers, dict):
        spacer_lengths = [len(seq) for seq in spacers.values()]
        num_spacers = len(spacers)
    else:
        # If just numbers provided
        spacer_lengths = [spacers]
        num_spacers = 1
    
    # Determine max edits for both calculations
    max_edits = max_mismatches
    if max_edit_distance is not None:
        max_edits = max_edit_distance
    
    # For Hamming distance: Account for both strands
    total_db_size_hamming = total_contig_length * 2
    
    # For BLAST E-value: Use effective database size (with length adjustments)
    # Based on BLAST's statistical theory for local alignment
    avg_spacer_len = sum(spacer_lengths) / len(spacer_lengths)
    H = 0.5772156649  # Euler's constant
    effective_db_size = total_contig_length / (avg_spacer_len - H)
    
    # Calculate BOTH Hamming and Edit distance estimates
    expected_matches_hamming_per_spacer = []
    expected_matches_edit_per_spacer = []
    
    for spacer_len in spacer_lengths:
        # Hamming distance (substitutions only)
        expected_matches_hamming = _calculate_edit_distance_probability_hamming(
            spacer_len, max_edits, match_prob, total_db_size_hamming
        )
        expected_matches_hamming_per_spacer.append(expected_matches_hamming)
        
        # BLAST E-value calculation with corrected parameters
        # Using proper nucleotide BLAST parameters
        expected_matches_edit = _calculate_edit_distance_probability_with_indels(
            spacer_len, max_edits, match_prob, effective_db_size,
            match_reward=1, mismatch_penalty=-3,  # Standard BLASTN scoring
            k_param=0.14, lambda_param=1.28)  # Lambda adjusted for +1/-3 scoring
        expected_matches_edit_per_spacer.append(expected_matches_edit)
    
    # Total expected spurious alignments (sum over all spacers)
    total_expected_hamming = sum(expected_matches_hamming_per_spacer)
    total_expected_edit = sum(expected_matches_edit_per_spacer)
    
    # Standard deviations
    std_dev_hamming = np.sqrt(total_expected_hamming)
    std_dev_edit = np.sqrt(total_expected_edit)
    
    base_comp_note = 'uniform' if all(abs(v - 0.25) < 1e-6 for v in base_comp_dict.values()) else 'non-uniform'
    
    # Return both estimates
    result = {
        'expected_spurious_hamming': total_expected_hamming,
        'std_spurious_hamming': std_dev_hamming,
        'blast_evalue': total_expected_edit,  # Renamed from expected_spurious_edit
        'blast_evalue_std': std_dev_edit,     # Renamed from std_spurious_edit
        'max_spurious_hamming': total_expected_hamming + 2 * std_dev_hamming,
        'min_spurious_hamming': max(0, total_expected_hamming - 2 * std_dev_hamming),
        'max_blast_evalue': total_expected_edit + 2 * std_dev_edit,
        'min_blast_evalue': max(0, total_expected_edit - 2 * std_dev_edit),
        'per_spacer_expected_hamming': expected_matches_hamming_per_spacer,
        'per_spacer_blast_evalue': expected_matches_edit_per_spacer,  # Renamed
        'method': 'combined_estimates',
        'base_composition': base_comp_dict,
        'base_comp_note': base_comp_note,
        'note': f'Combined Hamming-based and BLAST E-value estimates with {base_comp_note} base composition'
    }
    
    # For backward compatibility, also include the single estimate based on use_edit_distance
    if use_edit_distance:
        result['mean_spurious'] = total_expected_edit
        result['std_spurious'] = std_dev_edit
        result['distance_type'] = 'blast_evalue'  # Renamed from edit_distance
        result['distance_note'] = f'BLAST E-value estimate (conservative)'
    else:
        result['mean_spurious'] = total_expected_hamming
        result['std_spurious'] = std_dev_hamming
        result['distance_type'] = 'hamming'
        result['distance_note'] = f'Hamming distance <= {max_mismatches} (theoretical count)'
    
    return result
