# THIS IS THE GitHub repository for continued devlopment, for the original repository see: [https://code.jgi.doe.gov/spacersdb/spacer_matching_bench](https://code.jgi.doe.gov/spacersdb/spacer_matching_bench)


# Benchmarking protospacer identification tools
Scripts, code, and data for the manuscript: "Tool Choice drastically Impacts CRISPR Spacer-Protospacer Detection"   
Raw data files available in zenodo: [https://zenodo.org/doi/10.5281/zenodo.15171878](https://zenodo.org/doi/10.5281/zenodo.15171878).  
This benchmark is a companion project to the [SpacerDB](https://spacers.jgi.doe.gov/) and the [SpacerExtractor](https://code.jgi.doe.gov/SRoux/spacerextractor) tools.  


## Overview
CLI Python script for simulating spacer insertions into contigs and benchmarking sequence aligners. Generates synthetic data with known ground truth and evaluates aligner performance metrics (precision, recall, F1) and runtime. Note that for proper in-depth analysis, the jupyter notebooks are preferable (less heuristics, more detailed, and runs additional tests to verify the different tool reported alignments (so more consistent approach)).

## How should I analyse my own samples? (tl;dr)
For analyzing your own samples, we recommend using [SpacerExtractor](https://code.jgi.doe.gov/SRoux/spacerextractor#extra-functions---map-spacers-to-potential-targets-at-scale) rather than this benchmarking repository. SpacerExtractor is designed for practical analysis, using only bowtie1 for alignments, and includes features like:
- Wrapping bowtie1 commands for indexing and alignment of contigs and spacers.
- Extraction of flanking sequences (10 bases up/downstream of matches) for PAM analysis
- Built-in SAM file parsing (and filtering) and alignment visualization.

### Still, when should I use code from this benchmarking repository for my own (non-benchmarking) analysis?
(Apart from running benchmarks, of course), You might find the code useful for if you are interested in:
1. Utilizing multiple search tools,
2. Double-checking the reported alignments
3. Matches with >3 mismatches, or with gaps (bowtie1, wrapped by SpacerExtractor, is limited to 3).
4. The performance/resource usage of a tool (for example, you are optimizing a tool for a specific range of datasets size and parameter ranges).

### Key Features
- Simulates contigs and spacers at given frequency and specifications.
- Hyperfine integration for runtime benchmarking.
- Configurable parameters: contig count (default: 1500), spacer count (default: 50), mismatches (default: 0-5), reverse complements (default: 0.5), length of contigs (default: 1000-32000bp), length of spacers (default: 20-60bp), number of threads to use (default: 1), and number of hyperfine warmups and max runs (default: 0 and 5).
- Parallel processing support (most simulations is done in the rust code (see [rust_simulator](https://code.jgi.doe.gov/spacersdb/spacer_matching_bench/-/tree/main/rust_simulator?ref_type=heads/)).
- Modular tool configuration system.

## Aligners and search-engines tested so far 
For a complete list, see the [tool_versions.txt](tool_configs/tool_versions.txt) and the json files in the `tool_configs` folder. Note, not all of these were tested on the real dataset, some are later additions or variants of the tools tested.
1. [Bowtie1](https://github.com/BenLangmead/bowtie)
2. [Bowtie2](https://github.com/BenLangmead/bowtie2)
3. [Minimap2](https://github.com/lh3/minimap2)
4. [BBMap-skimmer](https://sourceforge.net/projects/bbmap/)
5. [StrobeAlign](https://github.com/ksahlin/StrobeAlign)
6. [BLASTN-short](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
7. [MMseqs2](https://github.com/soedinglab/MMseqs2)
8. [Naive impl. in Rust](https://github.com/apcamargo/spacer-containment)
9. [nucmer](https://github.com/mummer4/mummer)

## Installation
The code is wrapped as a python package with a rust module. These require some dependencies. I recomend using [pixi](https://pixi.sh/) to manage the default environment used for the simulation, but you can use mamba or conda to install to create a basic environment too (see the [`benchy_env.sh`](./src/bench/utils/benchy_env.sh) script).
The basic/default environment includes Python <=3.9, polars, hyperfine, the simplistic implementation of substring containment in Rust ([spacer-containment](https://github.com/apcamargo/spacer-containment)), and other general dependencies used in the notebooks (like matplotlib, seaborn, altair, pyfastx, needletail, etc). IIRC all used tools should be listed in the [references bibtex](draft/main/references.bib).
After the base environment is created, you could use the `tool_env_maker.sh` script to create micromamba environments for the aligner/search-engines tested in the manuscript. The idea is to avoid conflicts between the different versions of the tools, and have one environment per tool.  
Note - all these different environments can take up a lot of disk space.    
Note2 - the `tool_env_maker.sh` script is not mandatory - you can create the environments manually, but keep in mind that the environment names are used in the tool configs json files (e.g. `bowtie1_env`), so you might need to change them in the json files.  
Note3 - This is designed (tested) to run on a Linux system with a bash shell.  
Note4 - the full `simulated_data` and `raw_outputs` directories are not uploaded to gitlab to conserve space. A static copy is available in the [zenodo](https://zenodo.org/record/15171878) deposit.   

tl;dr - with pixi: 
```bash
# clone the repo
git clone https://code.jgi.doe.gov/spacersdb/spacer_matching_bench.git
cd spacer_matching_bench
# get pixi if you don't have it
if ! command -v pixi &> /dev/null; then
    # linux and mac
    curl -fsSL https://pixi.sh/install.sh | sh 
    # windows
  # powershell -ExecutionPolicy ByPass -c "irm -useb https://pixi.sh/install.ps1 | iex"
fi
# Create the default enviroment
pixi install 
# install the default environment and the "bench" python package (plus the rust simulation module)
pixi run build-all
# create the tool specific mamba envs
bash ./src/bench/utils/tool_env_maker.sh
```

If you don't want to use pixi, you can use `benchy_env.sh` to create a base environment with the general dependencies.

## Quick Start / Usage

The tool provides a `full-run` command that executes the entire pipeline, or you can run individual steps:

### One-Command Pipeline

Run the complete benchmarking pipeline in a single command:

```bash
spacer_bencher full-run \
  -o results/my_benchmark \  # Output directory
  -nc 100 \                  # 100 contigs
  -ns 50 \                   # 50 spacers
  -lm 0 2 \                  # 0-2 mismatches
  -t 4                       # 4 threads
```

This internally runs all four pipeline steps (simulate, generate-scripts, execute-tools, compare-results).

### Individual Pipeline Steps

For more control, run each step separately:

#### Step 1: Generate Simulated Data

```bash
spacer_bencher simulate \
  -nc 100 \              # Number of contigs
  -ns 50 \               # Number of spacers
  -ir 1 1 \              # Insert each spacer 1 time
  -lm 0 0 \              # 0 mismatches
  -nir 0 0 \             # 0 insertion indels
  -ndr 0 0 \             # 0 deletion indels
  -o tests/validation    # Output directory
```

#### Step 2: Generate Tool Execution Scripts

```bash
spacer_bencher generate-scripts \
  -i tests/validation \  # Input directory with simulated data
  -t 4                   # Number of threads
```

#### Step 3: Execute Alignment Tools

```bash
spacer_bencher execute-tools \
  -i tests/validation    # Input directory (scripts already have settings)
```

#### Step 4: Compare Results

```bash
spacer_bencher compare-results \
  -i tests/validation \      # Input directory
  -mm 3 \                    # Max mismatches to consider
  -o results.tsv             # Output file
```

## Command Reference

### `full-run` - Complete Pipeline

Execute the entire benchmarking pipeline in a single command.

```bash
spacer_bencher full-run [OPTIONS]

# Example
spacer_bencher full-run -o results/benchmark -nc 100 -ns 50 -lm 0 2 -t 4

# Get detailed help
spacer_bencher full-run --help
```

#### Key Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `-o` | `--output-dir` | Output directory for all pipeline results | **required** |
| `-nc` | `--num-contigs` | Number of contigs to generate | 100 |
| `-ns` | `--num-spacers` | Number of spacers to generate | 50 |
| `-cl` | `--contig-length` | Contig length range (min max) | 2000 5000 |
| `-ls` | `--spacer-length` | Spacer length range (min max) | 20 60 |
| `-ir` | `--spacer-insertions` | Number of times to insert each spacer | 1 1 |
| `-nir` | `--indel-insertions` | Number of insertion mutations per spacer | 0 0 |
| `-ndr` | `--indel-deletions` | Number of deletion mutations per spacer | 0 0 |
| `-lm` | `--mismatch-range` | Substitution mismatches per spacer | 0 0 |
| `-prc` | `--reverse-complement` | Proportion of reverse complement (0-1) | 0.5 |
| `-t` | `--threads` | Number of threads | 4 |
| `-mr` | `--max-runs` | Benchmark runs (hyperfine) | 1 |
| `-w` | `--warmups` | Warmup runs (hyperfine) | 0 |
| `-st` | `--skip-tools` | Comma-separated tools to skip | "" |
| `-rt` | `--only-tools` | Only run these tools | None |
| `-mm` | `--max-mismatches` | Max mismatches for comparison | 5 |

#### Pipeline Steps

The `full-run` command internally executes:
1. **Simulate** - Generates synthetic data with ground truth
2. **Generate-scripts** - Creates bash scripts for each tool
3. **Execute-tools** - Runs all alignment tools
4. **Compare-results** - Validates and compares tool outputs

#### Output Structure

```
results/benchmark/
├── simulated_data/
│   ├── simulated_contigs.fa
│   ├── simulated_spacers.fa
│   └── planned_ground_truth.tsv
├── raw_outputs/
│   ├── bowtie1_output.sam
│   ├── minimap2_output.sam
│   └── ...
└── comparison_results.tsv
```

### `simulate` - Generate Simulated Data

Generate synthetic CRISPR spacer and contig sequences with configurable mutations.

```bash
spacer_bencher simulate [OPTIONS]

# Get detailed help
spacer_bencher simulate --help
```

#### Key Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `-nc` | `--num-contigs` | Number of contigs to generate | 100 |
| `-ns` | `--num-spacers` | Number of spacers to generate | 50 |
| `-cl` | `--contig-length` | Contig length range (min max) | 2000 5000 |
| `-ls` | `--spacer-length` | Spacer length range (min max) | 20 60 |
| `-ir` | `--spacer-insertions` | **Number of times to insert each spacer** | 1 1 |
| `-nir` | `--indel-insertions` | **Number of insertion mutations per spacer** | 0 0 |
| `-ndr` | `--indel-deletions` | **Number of deletion mutations per spacer** | 0 0 |
| `-lm` | `--mismatch-range` | Substitution mismatches per spacer | 0 0 |
| `-prc` | `--reverse-complement` | Proportion of reverse complement (0-1) | 0.5 |
| `-o` | `--output-dir` | Output directory | **required** |
| `-t` | `--threads` | Number of threads | 4 |

#### Understanding the Parameters

**CRITICAL: Insertion Terminology**

There are TWO different meanings of "insertion" in this tool:

1. **`--spacer-insertions` (-ir)**: How many TIMES each spacer is placed into contigs
   - Example: `-ir 1 2` means each spacer will be inserted 1 or 2 times
   - This creates the ground truth data

2. **`--indel-insertions` (-nir)**: How many INSERTION MUTATIONS (adding bases) to apply
   - Example: `-nir 0 2` means 0-2 insertion indel events within the spacer sequence
   - This creates sequence-level mutations (like biological indels)

3. **`--indel-deletions` (-ndr)**: How many DELETION MUTATIONS (removing bases) to apply
   - Example: `-ndr 0 1` means 0-1 deletion indel events within the spacer sequence

#### Example Use Cases

**Perfect matches (validation baseline)**:
```bash
spacer_bencher simulate \
  -nc 100 -ns 50 \
  -ir 1 1 \          # Insert each spacer exactly once
  -lm 0 0 \          # No mismatches
  -nir 0 0 \         # No insertion indels
  -ndr 0 0 \         # No deletion indels
  -o tests/perfect_matches
```

**multiple variantion data**:
```bash
spacer_bencher simulate \
  -nc 500 -ns 100 \
  -ir 1 3 \          # Each spacer inserted 1-3 times
  -lm 0 2 \          # 0-2 substitution mutations
  -nir 0 1 \         # 0-1 insertion indels
  -ndr 0 1 \         # 0-1 deletion indels
  -prc 0.5 \         # 50% reverse complement
  -o tests/realistic
```

**Custom GC content**:
```bash
spacer_bencher simulate \
  -nc 100 -ns 50 \
  --gc-content 60 \  # 60% GC content
  -o tests/high_gc
```

#### Output Files

The simulation creates a `simulated_data` subdirectory containing:

- **`simulated_contigs.fa`**: FASTA file with contig sequences
- **`simulated_spacers.fa`**: FASTA file with spacer sequences  
- **`planned_ground_truth.tsv`**: Ground truth annotations with columns:
  - `spacer_id`: Spacer identifier
  - `contig_id`: Contig identifier
  - `start`: Start position (0-based)
  - `end`: End position (0-based, exclusive)
  - `strand`: Boolean (true=reverse complement, false=forward)
  - `mismatches`: Number of substitution mismatches

### `generate-scripts` - Create Tool Execution Scripts

Generate bash scripts for running alignment tools.

```bash
spacer_bencher generate-scripts [OPTIONS]

# Example
spacer_bencher generate-scripts -i tests/validation -t 4 -mr 1

# Get detailed help
spacer_bencher generate-scripts --help
```

| Parameter | Long Option | Description | Default |
|-----------|-------------|-------------|---------|
| `-i` | `--input-dir` | Directory with simulated data | **required** |
| `-t` | `--threads` | Threads for tool execution | 4 |
| `-mr` | `--max-runs` | Benchmark runs (hyperfine) | 1 |
| `-w` | `--warmups` | Warmup runs (hyperfine) | 0 |
| `-st` | `--skip-tools` | Comma-separated tools to skip | "" |
| `-rt` | `--only-tools` | Only run these tools | None |

### `execute-tools` - Run Alignment Tools

Execute the generated bash scripts and collect results.

```bash
spacer_bencher execute-tools [OPTIONS]

# Example
spacer_bencher execute-tools -i tests/validation

# Get detailed help
spacer_bencher execute-tools --help
```

| Parameter | Long Option | Description | Default |
|-----------|-------------|-------------|---------|
| `-i` | `--input-dir` | Directory with generated scripts | **required** |
| `-st` | `--skip-tools` | Comma-separated tools to skip | "" |
| `-rt` | `--only-tools` | Only run these tools | None |
| `--debug` | | Run commands directly (no hyperfine) for better error messages | False |

**Note:** Thread counts and benchmark settings are already baked into the generated scripts from `generate-scripts`, so you don't need to specify them here.

**Debug Mode:** Use `--debug` flag to run tools directly without hyperfine wrapper, which provides full stdout/stderr output for troubleshooting failures:
```bash
spacer_bencher execute-tools -i tests/validation -rt lexicmap --debug
```

**Output Files** (saved to `raw_outputs/` subdirectory):
- `TOOL_NAME_output.sam` - Alignment results in SAM format
- `TOOL_NAME.json` - Benchmarking results (if max-runs > 1)

### `compare-results` - Validate and Compare Results

Compare tool results against ground truth and calculate performance metrics.

```bash
spacer_bencher compare-results [OPTIONS]

# Example
spacer_bencher compare-results -i tests/validation -mm 3 -o results.tsv

# Get detailed help
spacer_bencher compare-results --help
```

| Parameter | Long Option | Description | Default |
|-----------|-------------|-------------|---------|
| `-i` | `--input-dir` | Directory with tool outputs | **required** |
| `-mm` | `--max-mismatches` | Max mismatches to consider | 5 |
| `-o` | `--output-file` | Output file for results | stdout |
| `-t` | `--threads` | Threads for processing | 4 |

**Metrics Calculated:**
- **True Positives (TP)** - Correctly identified spacer locations
- **False Positives (FP)** - Incorrectly reported locations  
- **False Negatives (FN)** - Missed spacer locations
- **Precision** - TP / (TP + FP)
- **Recall** - TP / (TP + FN)
- **F1 Score** - Harmonic mean of precision and recall

## Complete Workflow Example

### Option 1: Single Command (Recommended for Quick Runs)

```bash
# Run the complete pipeline with one command
spacer_bencher full-run \
  -o tests/my_validation \
  -nc 100 -ns 50 \
  -lm 0 1 \
  -ir 1 1 \
  -nir 0 0 -ndr 0 0 \
  -t 4 \
  -mm 3
```

### Option 2: Step-by-Step (More Control)

```bash
# 1. Create test data with 0-1 mismatches
spacer_bencher simulate \
  -nc 100 -ns 50 \
  -lm 0 1 \
  -ir 1 1 \
  -nir 0 0 -ndr 0 0 \
  -t 4 \
  -o tests/my_validation

# 2. Generate execution scripts
spacer_bencher generate-scripts -i tests/my_validation -t 4

# 3. Run alignment tools (no parameters needed - uses generated scripts)
spacer_bencher execute-tools -i tests/my_validation

# 4. Compare and validate results
spacer_bencher compare-results -i tests/my_validation -mm 3 -o results.tsv
```

## Additional Features

### Logging and Debugging

Enable verbose logging with timestamps and detailed tracebacks:

```bash
spacer_bencher --verbose simulate -nc 100 -ns 50 -o tests/debug
```

### Tool Configuration

Tool configurations are stored as JSON files in `tool_configs/` directory. Each tool has:
- Command template with parameter placeholders
- Environment name (micromamba environment)
- Output file format
- Parser function name

See existing tool configs for examples when adding new tools.

### Documentation

- **[src/rust_simulator/README.md](src/rust_simulator/README.md)** - Detailed simulator documentation

## Citation

If you use this tool in your research, please cite:

<!-- ## License
Code required for running the benchmark is available under the [MIT License](LICENSE).   
Note, each individual tool tested by the benchmark is likely to have its own license and copyright.   
 -->
