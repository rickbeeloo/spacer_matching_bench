
# Benchmarking protospacer identification tools
Scripts, code, and data for the manuscript: "Tool Choice _drastically_ Impacts CRISPR Spacer-Protospacer Detection"   
Raw data files available in zenodo: [https://zenodo.org/doi/10.5281/zenodo.15171878](https://zenodo.org/doi/10.5281/zenodo.15171878).  
This benchmark is a companion project to the [SpacerDB](https://spacers.jgi.doe.gov/) and the [SpacerExtractor](https://code.jgi.doe.gov/SRoux/spacerextractor) tools.  

## Overview
CLI programs and scripts for:
1. simulating CRISPR spacer and contig sequences with configurable mutations,
2. running multiple sequence aligners/search-engines to identify spacer locations in contigs,
3. comparing the reported alignments to a ground truth/planned insertions,
4. calculating performance metrics using different criteria (edit/hamming, out of planned/observed) and collecting runtime and resource usage information (using hyperfine or from SLURM logs).

## Aligners and search-engines tested so far 
For a complete list, see the [tool_versions.csv](tool_configs/tool_versions.csv) and the json files in the `tool_configs` folder. Note, not all of these were tested on the "full" real dataset, some are later additions or variants of the tools tested.
1. [Bowtie1](https://github.com/BenLangmead/bowtie)
2. [Bowtie2](https://github.com/BenLangmead/bowtie2)
3. [Minimap2](https://github.com/lh3/minimap2)
4. [indelfree.sh](https://sourceforge.net/projects/bbmap/)
5. [StrobeAlign](https://github.com/ksahlin/StrobeAlign)
6. [BLASTN-short](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
7. [MMseqs2](https://github.com/soedinglab/MMseqs2)
8. [sassy](https://github.com/ragnargk/sassy)
9. [nucmer](https://github.com/mummer4/mummer)

## Installation
The code is wrapped as a python package with a rust module. These require some dependencies. I recomend using [pixi](https://pixi.sh/) to manage the default environment used for the simulation, but you can use mamba or conda to install to create a basic environment too (see the [`benchy_env.sh`](./src/bench/utils/benchy_env.sh) script).
The basic/default environment includes Python >=3.10, polars, hyperfine, sassy (note - requires nighly rust iirc) and some other general dependencies used in the notebooks (like matplotlib, seaborn, altair, pyfastx, needletail, etc). IIRC all used tools should be listed in the [references bibtex](draft/main/references.bib).  
After the base environment is created, you could use the `tool_env_maker.sh` script to create micromamba environments for the aligner/search-engines tested in the manuscript. The idea is to avoid conflicts between the different versions of the tools, and have one environment per tool.  
**Note** - all these different environments can take up a lot of disk space.    
**Note2** - the `tool_env_maker.sh` script is not mandatory - you can create the environments manually, but keep in mind that the environment names are used in the tool configs json files (e.g. `bowtie1_env`), so you might need to change them in the json files.  
**Note3** - This is designed (tested) to run on a Linux system with a bash shell.  
**Note4** - Some directories (like the full `simulated_data` and `raw_outputs` directories) are listed in the .gitignore as they're too big for github. A static copy is/will-be available in the [zenodo](https://zenodo.org/record/15171878) deposit (# TODO: REMINDER reupload before submitting revision).

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

## CLI commands and options:
`Usage: spacer_bencher [OPTIONS] COMMAND [ARGS]...`
All CLI subcommands have a `--help` flag for with specific information (copy-pasted below too).  
You can run a full pipeline with a single command, or run each step separately, altough some commands reuse previous command outputs/config.  
**NOTE** For analyses of large datasets, we recommend using the SLURM based workflow used in the notebooks (a la `notebooks/*.ipynb` folder). The CLI version uses hyperfine and runs locally - hyperfine is great for clearing I/O biases/caching effects. Also, the CLI version's method for comparing results is much more simplistic "quick and dirty" process, while in the notebooks we also validate the reported mismatches by recalculating them from the sequences + some spot checks. The CLI can support this but keep in mind.
**Note2** use the `--verbose` flag to enable detailed logging with timestamps and tracebacks of some spot checked alignemnts.  

#### Commands:
  1. `simulate`          Generate simulated CRISPR spacer and contig sequences.
```bash
  Options:
  -nc, --num-contigs INTEGER      Number of contigs to generate  [default: 1000]
  -ns, --num-spacers INTEGER      Number of spacers to generate  [default: 200]
  -cl, --contig-length <INTEGER INTEGER>...
                                  Range of contig lengths (min max)  [default: 2000, 15000]
  -ls, --spacer-length <INTEGER INTEGER>...
                                  Range of spacer lengths (min max)  [default: 20, 60]
  -lm, --mismatch-range <INTEGER INTEGER>...
                                  Range of substitution mismatches per spacer (min max)  [default: 0, 3]
  -ir, --spacer-insertions <INTEGER INTEGER>...
                                  Number of times to insert each spacer into contigs (min max) [default: 1, 4]
  -nir, --indel-insertions <INTEGER INTEGER>...
                                  Number of insertion mutations (indels) to add within each spacer (min max)  [default: 0, 0]
  -ndr, --indel-deletions <INTEGER INTEGER>...
                                  Number of deletion mutations (indels) to add within each spacer (min max)  [default: 0, 0]
  -prc, --reverse-complement FLOAT
                                  Proportion of spacers to reverse complement (0.0-1.0)  [default: 0.5]
  -t, --threads INTEGER           Number of threads for parallel processing [default: 4]
  -o, --output-dir PATH           Output directory for simulated data [required]
  -id, --id-prefix TEXT           Prefix for sequence IDs (default: auto-generated)
  --gc-content FLOAT              GC content percentage (0-100) - DEPRECATED: use --contig-gc-content and --spacer-gc-content. if this is set, it will override the spacer/contig specific ones.
  --contig-gc-content FLOAT       GC content percentage for contigs (0-100)
  --spacer-gc-content FLOAT       GC content percentage for spacers (0-100)
  --contig-distribution [uniform|normal|bell]
                                  Distribution for contig lengths  [default: uniform]
  --spacer-distribution [uniform|normal|bell]
                                  Distribution for spacer lengths  [default: uniform]
  --verify                        Verify simulation after generation
  --contigs PATH                  Path to existing contigs FASTA file (optional, for read-and-insert mode)
  --spacers PATH                  Path to existing spacers FASTA file (optional, for read-and-insert mode)
  ```
  2. `generate-scripts`  Generate execution scripts for alignment tools.
```bash
  Options:
  -i, --input-dir PATH           Input directory containing simulated data [required]
  -o, --output-dir PATH          Output directory for generated scripts (defaults to input-dir if specified)
  --contigs PATH                 Path to custom contigs file
  --spacers PATH                 Path to custom spacers file
  -t, --threads INTEGER          Number of threads for tool execution [default: 4]
  -mr, --max-runs INTEGER        Maximum number of benchmark runs (hyperfine) [default: 1]
  -w, --warmups INTEGER          Number of warmup runs (hyperfine) [default: 0]
  -st, --skip-tools TEXT         Comma-separated list of tools to skip [default: ""]
  -rt, --only-tools TEXT         Comma-separated list of tools to run (overrides skip-tools)
  -mm, --max-mismatches INTEGER  Maximum number of mismatches/edit-distance parameter for tools that support it (i.e. subs in indelfree.sh, k in sassy, -v in bowtie1) [default: 5]
  ```
  3. `execute-tools`     Execute alignment tools using generated bash scripts.
```bash
  Options:
  -i, --input-dir PATH    Input directory containing scripts and data [required]
  -st, --skip-tools TEXT  Comma-separated list of tools to skip [default: ""]
  -rt, --only-tools TEXT  Comma-separated list of tools to run
  --debug                 Run commands directly without hyperfine for better error messages
  --help                  Show this message and exit.
```
  4. `compare-results`   Compare and validate alignment tool results.
```bash
  Options:
  -i, --input-dir PATH           Input directory containing tool outputs and ground truth [required]
  -mm, --max-mismatches INTEGER  Maximum number of mismatches to consider
  -o, --output-file PATH         Output file for comparison results (default: stdout)
  -t, --threads INTEGER          Number of threads for processing
  --augment-ground-truth         Count verified non-planned alignments as true positives
  -st, --skip-tools TEXT         Comma-separated list of tools to skip
  -rt, --only-tools TEXT         Comma-separated list of tools to run
  --distance [hamming|edit]      Distance metric for validation: hamming (substitutions only) or edit (substitutions + indels). Default: hamming
  --gap-open-penalty INTEGER     Gap open penalty for alignment validation
  --gap-extend-penalty INTEGER   Gap extension penalty for alignment validation
  --help                         Show this message and exit.
```
  5. `full-run`          Runs 1-4 above.
```bash
  Options:
  all of the above...
```

#### Output Structure
for `-o blabla`, the output structure will be: 
```
blabla/
├── simulated_data/
│   ├── simulated_contigs.fa
│   ├── simulated_spacers.fa
│   └── planned_ground_truth.tsv
├── bash_scripts/
│   ├── <tool_name>.sh
│   └── ...
├── job_scripts/   # only if using SLURM workflow as in the notebooks...
│   ├── <tool_name>.sh
│   └── ...
├── slurm_logs/   # only if using SLURM workflow as in the notebooks...
│   ├── <tool_name>-jobid.out
│   ├── <tool_name>-jobid.err
│   └── ...
├── raw_outputs/
│   ├── bowtie1_output.sam
│   ├── minimap2_output.sam
│   └── ...
└── comparison_results.tsv
```

#### Understanding the Parameters

##### **CRITICAL: Insertion Terminology** (a mutation type vs number of appearances)

There are TWO different meanings of "insertion" in this project:
1. **`--spacer-insertions` (-ir)**: How many TIMES each spacer is placed into contigs
   - Example: `-ir 1 2` means each spacer will be inserted 1 or 2 times. This is partially synonymous with "copy number" / "abundance" / "occurrences" / "appearances".
   - This creates the ground truth data

2. **`--indel-insertions` (-nir)**: How many INSERTION MUTATIONS (adding bases) to apply
   - Example: `-nir 0 2` means 0-2 insertion indel events within the spacer sequence. This is partially synonymous with "indel events" or "insertion mutations".
   - This creates sequence-level mutations (like biological indels)

3. **`--indel-deletions` (-ndr)**: How many DELETION MUTATIONS (removing bases) to apply
   - Example: `-ndr 0 1` means 0-1 deletion indel events within the spacer sequence. This is partially synonymous with "indel events" or "deletion mutations".   
For simplicity, the indel options are not really used now and might cause some bugs. In the **current** version of the project, we are focusing on substitution mismatches in the simulation, after noticing high false positive with indels. This might be revisited in the future, and we include some edit distance support in the other commands/notebooks.

##### Semi-Synthetic options: mixing real spacers/contigs with simulated ones
You can provide your own contig and/or spacer FASTA files using the `--contigs` and `--spacers` options. If a number of contigs/spacers is also specified, that much will be subsampled from the provided files, while all other sequence feature arguments (length, gc, distributions) for that set (spacer/contigs/both) will be ignored.  
The motivation for this is to allow estimating (for a given size of datasets) the amount of non-planned (a.k.a "spurious") alignments that could be reported tools when using real biological spacers, tested for occurrence in simulated contigs. This is trynig to answer "how many potential alignments of this spacer set should be expectd to occur by chance alone when allowing for xyz thersholds/tool-selection etc". Note, these are "false-positives" in the sense we didn't explicitly planned their occurences in the contigs, this term "false-positive" doesn't imply these are incorrect alignments - they must still pass the same validation criteria.  


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

### Tool Configuration

Tool configurations are stored as JSON files in `tool_configs/` directory. Each tool has:
- Command template with parameter placeholders
- Environment name (micromamba environment)
- Output file format
- Parser function name

See existing tool configs for examples when adding new tools.

## Citation

If you use this tool in your research, please cite:

# Benchmark Reproduction  (notebooks)
In the preprint/paper we use the CLI tool or code/fucntions/etc to test relatively large datasets and various scnarios (fraction of real data, semi-synthetic, simulated, hard-limiting reporting etc). These are tracked using jupyter notebooks (and SLURM logs), in the `notebooks` folder. The general execution order:  
1. [fetch_real_data](notebooks/fetch_real_data.ipynb) - download and prepare the "real" contig set (IMG/VR4 filtered to HQ contigs), and the real spacer set (iPHoP spacers).
2. [spacer_inspection](notebooks/spacer_inspection.ipynb) - Examines the spacer set sequence characterisics (GC/nucleobase range and distribution, lengths, complexity, repetition...) and does some (minimaly) filtering for length, redundancy, and some complexity metrics. This notebook also demonstrates that our simulation module can recapitulate the spacer length and GC content distributions and (roughly) the complexity metrics, as the real spacers.
2.5 Investigating edit/hamming distance validation differences:  
   - [distance_metric_analysis](notebooks/distance_metric_analysis.ipynb) - investigates the differences between edit distance and hamming distance thesholds and validation methods, and their impact on potentially non-planned/spurious alignments.  
3. [Prepare_all_jobs](notebooks/Prepare_all_jobs.ipynb) - The main SLURM workflow logic; this preparesthe dirctory structre for the "real", "synthetic", and "semi-synthetic" datasets. Briefly, internally, this generates the subsets (`fraction_1/0.nnn`), the different spacer X contig simulated runs, and the "real-baseline" run (real spacers in simulated contigs matching the real contigs sequence features). Note that for this also splits the allowed `max_distance` provided to some tools; specifically, large datasets are either not tested with the compute heavy, exustaive tools (namely indelfree_bruteforce and sassy) or they are limited to 3 (edit/hamming) distance only.  
4. Performence_analysis notebooks:
   - [subsamples_analysis](notebooks/subsamples_analysis.ipynb) - analyses the results from the various fractions of the real dataset (real spacers in real contigs).
   - [semi_synthetic_analysis](notebooks/semi_synthetic_analysis.ipynb) - analyses the results from the semi-synthetic datasets (real spacers in simulated contigs).
   - [simulated_analysis](notebooks/simulated_analysis.ipynb) - analyses the results from the fully simulated datasets (simulated spacers in simulated contigs).
5. Resource usage analysis:
   - [resource_usage_analysis](notebooks/resource_usage_analysis.ipynb) - collects and analyses the hyperfine results from the various runs, generating plots and summary tables.
...  

n. [get_tool_versions](notebooks/get_tool_versions.ipynb)


## FAQs:

### How should I analyse my own samples? (tl;dr)
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
5. Pre-filtering spacers sequences using customisable criteria (length, base composition, complexity metrics) before calling a search/alignment tool.

### How to add a new tool to the benchmark?
Ideally (patient) - open an issue with detailed instructions on the tool you want added, ideally noting the exact parameters / CLI args to use with it, and we will consider it for a future version (which might take some time).  
Otherwise, or to get a feel if you should even reccommend the tool for addition, you canuse one of these two options: "Quick and Dirty" (only get performence metrics using the pre-generated/subsampled datasets), or "Full" (generate new scripts, run the tool, compare results, AND measure runtime, but probably using a small subset and an approximate alignment validation step).  

"Quick and Dirty" - create a config json file in the `tool_configs` folder as described above. Fetch the pre-existing datasets from zenodo (be patient, could be slow), and execute your tool on them to produce a SAM/similar output. Then process the output following the steps in the notebook of choise (i.e. if you used any of the fractions, use the `subsamples_analysis` notebook, if you used the "full" fraction_1, follow the ... and so on for the simulated too).  Note that the processed results for the other tools (alignments verified etc) are also uploaded to zenodo, and that the "read_results" function will ignore any missing SAM/outputs by default, so you can manually add that parquet to your tools' processed results).  
"Full" - After adding your tool's json config, use the `full-run` command as follows (note, this example is for a simulated, relatively small set, with 0 warmup runs and max 3 actual runs):
```bash
# 500 spacers, 5_000 contigs, HIGH INSERTIONS! 100-2500 insertions, 25-40 spacer length, 10000-150000 contig length, 0-5 mismatches
pixi run "spacer_bencher full-run --contig-distribution normal  --spacer-distribution normal --num-spacers 500 --spacer-gc-content 49 --contig-gc-content 46 --spacer-insertions 1 2500  --num-contigs 10000 --spacer-length 25 40  --contig-length 10000 750000 --mismatch-range 0 5 --output-dir simulated/ns_500_nc_10000/ --threads 8"
```



## Older versions (preprint)
The original repository was hosted on JGI's gitlab: [https://code.jgi.doe.gov/spacersdb/spacer_matching_bench](https://code.jgi.doe.gov/spacersdb/spacer_matching_bench).
That version is now obsolete, and superceded by this repository. I only mention it for discloure purposes...



<!-- ## License
Code required for running the benchmark is available under the [MIT License](LICENSE).   
Note, each individual tool tested by the benchmark is likely to have its own license and copyright.   
 -->