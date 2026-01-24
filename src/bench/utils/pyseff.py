import os
# import re
import shlex
import subprocess
import tempfile

import polars as pl

def parse_time_to_seconds(time_str):
    """Convert Slurm time format to seconds"""
    if not time_str or time_str == "":
        return 0
    try:
        # Handle days-hours:minutes:seconds format
        if "-" in time_str:
            days, time_part = time_str.split("-")
            hours, minutes, seconds = time_part.split(":")
            return (
                int(days) * 24 * 3600
                + int(hours) * 3600
                + int(minutes) * 60
                + int(seconds)
            )
        # Handle hours:minutes:seconds format
        else:
            parts = time_str.split(":")
            if len(parts) == 3:
                hours, minutes, seconds = parts
                return int(hours) * 3600 + int(minutes) * 60 + int(seconds)
            return 0
    except (ValueError, AttributeError):
        return float(0)


def parse_mem(mem_str):
    if not mem_str:
        return float(0)
    try:
        if mem_str.endswith("K"):
            return float(mem_str[:-1]) / 1024
        elif mem_str.endswith("G"):
            return float(mem_str[:-1]) * 1024
        elif mem_str.endswith("T"):
            return float(mem_str[:-1]) * 1024 * 1024
        else:
            return float(mem_str)
    except (ValueError, AttributeError):
        return float(0)


def format_seconds(seconds):
    days = seconds // (24 * 3600)
    remaining = seconds % (24 * 3600)
    hours = remaining // 3600
    remaining %= 3600
    minutes = remaining // 60
    seconds = remaining % 60
    if days > 0:
        return f"{days}-{hours:02d}:{minutes:02d}:{seconds:02d}"
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def pyseff(sacct_df=None,start_date = "2025-10-01", remove_cancelled=True,remove_failed=True, calculate_memory_efficiency=True,calculate_cpu_efficiency=True):
    """
    Main function to analyze job efficiency from sacct output, similar to Slurm's seff

    Uses specific fields from sacct:
    - JobID: for identifying jobs
    - JobName: name of the job
    - State: job state
    - ExitCode: exit status
    - AllocCPUS: allocated CPUs
    - Elapsed: total job time
    - TotalCPU: total CPU time (user + system)
    - MaxRSS: maximum resident set size
    - ReqMem: requested memory
    """
    if sacct_df is None:
        # Create temporary file for sacct output with specific fields
        tmp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
        sacct_cmd = (
            "sacct -u $USER --parsable "
            "--format=JobID,JobName,State,ExitCode,AllocCPUS,Elapsed,"
            "TotalCPU,MaxRSS,ReqMem "
            f"--starttime {start_date} "
            "> " + tmp_file.name
        )
        subprocess.run(sacct_cmd, shell=True, text=True, capture_output=True)
        sacct_df = pl.read_csv(tmp_file.name, separator="|")
        tmp_file.close()
    sacct_df = sacct_df.with_columns(
        pl.col("JobID").str.split(".").list.first().alias("BaseJobID")
    )
    # Remove cancelled jobs if specified
    if remove_cancelled:
        cancelled_jobs = (
            sacct_df.filter(pl.col("State").str.contains_any(patterns=["CANCELLED"]))[
                "BaseJobID"
            ]
            .unique()
            .to_list()
        )

        # Get base job IDs and filter out cancelled jobs
        sacct_df = sacct_df.filter(~pl.col("BaseJobID").is_in(cancelled_jobs))
    # Remove failed jobs if specified
    if remove_failed:
        failed_jobs = (
            sacct_df.filter(pl.col("State").str.contains_any(patterns=["FAILED"]))[
                "BaseJobID"
            ]
            .unique()
            .to_list()
        )

        # Get base job IDs and filter out failed jobs
        sacct_df = sacct_df.filter(~pl.col("BaseJobID").is_in(failed_jobs))

    if sacct_df.is_empty():
        return sacct_df
    
    # convert CPU time to seconds
    sacct_df = sacct_df.with_columns(
            [
                pl.col("TotalCPU")
                .map_elements(parse_time_to_seconds, return_dtype=pl.UInt32)
                .alias("TotalCPU_Seconds"),
                pl.col("Elapsed")
                .map_elements(parse_time_to_seconds, return_dtype=pl.UInt32)
                .alias("Elapsed_Seconds"),
            ]
        )
    
    # Calculate CPU efficiency
    if calculate_cpu_efficiency:
        sacct_df = sacct_df.with_columns(
            [
                (
                    pl.col("TotalCPU_Seconds")
                / (pl.col("Elapsed_Seconds") * pl.col("AllocCPUS"))
                * 100
            )
            .round(2)
            .alias("CPU_Efficiency")
            ]
        )

    # Always create MaxRSS_MB and ReqMem_MB columns
    sacct_df = sacct_df.with_columns(
        [
            pl.col("MaxRSS").map_elements(parse_mem, return_dtype=pl.Float64,skip_nulls=True)
            .alias("MaxRSS_MB"),
            pl.col("ReqMem")
            .map_elements(parse_mem, return_dtype=pl.Float64,skip_nulls=True)
            .alias("ReqMem_MB"),
        ]
    )

    # Calculate memory efficiency if requested
    if calculate_memory_efficiency:
        sacct_df = sacct_df.with_columns(
            (
                (pl.col("MaxRSS_MB") / pl.col("ReqMem_MB")) * 100
            )
            .round(2)
            .alias("Memory_Efficiency")
        )

    # Aggregate by BaseJobID
    agg_cols = [
        pl.col("JobName").first().alias("JobName"),
        pl.col("AllocCPUS").max().alias("AllocCPUS"),
        pl.col("State").unique().alias("State"),
        pl.col("ExitCode").unique().alias("ExitCode"),
        pl.col("MaxRSS").max().alias("MaxRSS"),
        pl.col("ReqMem").first().alias("ReqMem"),
        pl.col("Elapsed_Seconds").sum().alias("Elapsed_Seconds"),
        pl.col("TotalCPU_Seconds").sum().alias("TotalCPU_Seconds"),
        # Always include memory columns
        pl.col("MaxRSS_MB").max().alias("MaxRSS_MB"),
        pl.col("ReqMem_MB").first().alias("ReqMem_MB"),
    ]
    
    # Add CPU efficiency to aggregation if it was calculated
    if calculate_cpu_efficiency:
        agg_cols.append(pl.col("CPU_Efficiency").mean().alias("CPU_Efficiency"))
    
    # Add Memory efficiency to aggregation if it was calculated
    if calculate_memory_efficiency:
        agg_cols.append(pl.col("Memory_Efficiency").mean().alias("Memory_Efficiency"))
    
    sacct_df = sacct_df.group_by("BaseJobID").agg(agg_cols)

    # Format time columns back to human-readable format
    result_df = sacct_df.with_columns(
        [
            pl.col("Elapsed_Seconds")
            .map_elements(format_seconds, return_dtype=pl.Utf8)
            .alias("Elapsed"),
            pl.col("TotalCPU_Seconds")
            .map_elements(format_seconds, return_dtype=pl.Utf8)
            .alias("TotalCPU"),
        ]
    ) #.select(
        # [
    #         "BaseJobID",
    #         "JobName",
    #         "State",
    #         "ExitCode",
    #         "AllocCPUS",
    #         "Elapsed",
    #         "TotalCPU",
    #         # "CPU_Efficiency",
    #         "MaxRSS_MB",
    #         "ReqMem_MB",
    #         # "Memory_Efficiency",
    #     ]
    # )

    return result_df



def create_slurm_job_scripts(results_dir, slurm_opts=None, job_subdir="job_scripts",threads=16):
    """Create SLURM job scripts for each bash script under results_dir/bash_scripts.

    Args:
        results_dir (str): Path to the results directory which contains `bash_scripts` and `slurm_logs`.
        slurm_opts (dict or None): Optional SBATCH overrides. Keys are SBATCH flags without leading dashes, e.g. {'c':16, 'mem':'168G'}.
        job_subdir (str): Subdirectory under results_dir to write job scripts into.

    Returns:
        list: list of created job script paths.
    """
    bash_dir = os.path.join(results_dir, "bash_scripts")
    out_job_dir = os.path.join(results_dir, job_subdir)
    log_dir = os.path.join(results_dir, "slurm_logs")
    os.makedirs(out_job_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    equal_opts = ['mail-user', 'mail-type']
    # space_opts = [ 'A', 'q', 't']

    default_opts = {
        'mail-user': 'uneri@lbl.gov',
        'mail-type': 'FAIL,END,BEGIN',
        'A': 'grp-org-sc-metagen',
        'q': 'jgi_normal',
        'c': threads,
        'mem': '168G',
        't': '72:00:00',
    }
    if slurm_opts:
        default_opts.update(slurm_opts)

    created = []
    if not os.path.isdir(bash_dir):
        return created

    for entry in sorted(os.listdir(bash_dir)):
        if not entry.endswith('.sh'):
            continue
        script_path = os.path.join(bash_dir, entry)
        tool_name = os.path.splitext(entry)[0]
        job_name = f"{tool_name}"
        job_file = os.path.join(out_job_dir, f"{tool_name}.sh")

        # Build sbatch header
        with open(job_file, 'w') as jf:
            jf.write('#!/bin/bash\n')
            for k, v in default_opts.items():
                # Use single-dash for single-letter opts, double-dash otherwise
                if k in equal_opts:
                    prefix = '-' if len(str(k)) == 1 else '--'
                    jf.write(f"#SBATCH {prefix}{k}={v}\n")
                else:
                    prefix = '-' if len(str(k)) == 1 else '--'
                    jf.write(f"#SBATCH {prefix}{k} {v}\n")

            jf.write(f"#SBATCH -J \"{job_name}\"\n")
            jf.write(f"#SBATCH -o {os.path.join(log_dir, job_name)}-%A.out\n")
            jf.write(f"#SBATCH -e {os.path.join(log_dir, job_name)}-%A.err\n\n")

            # Run command: call the generated script with number of cpus
            # Use shlex.quote to safely embed paths
            jf.write(f"bash {shlex.quote(script_path)} \"$SLURM_CPUS_PER_TASK\"\n")

        os.chmod(job_file, 0o755)
        created.append(job_file)

    return created



def print_dataset_slurm_summary(results_dict, title):
    """Print a summary of script generation results."""
    print(f"{title}")
    for dataset_name, result in results_dict.items():
        print(f"{dataset_name}:")
        print(f"  max_distance: {result.get('max_distance', 'N/A')}")
        print(f"  Bash scripts: {len(result.get('tools', []))}")
        print(f"  SLURM jobs: {len(result.get('jobs', []))}")
        if result.get('tools'):
            print(f"  Tools: {', '.join(result['tools'][:5])}{'...' if len(result['tools']) > 5 else ''}")

def submit_slurm_jobs(results_dir, job_subdir="job_scripts", dry_run=False):
    """Submit all job scripts in results_dir/job_subdir via sbatch.

    Returns list of (script, sbatch_output) tuples for jobs submitted. If dry_run is True, only return the scripts.
    """
    job_dir = os.path.join(results_dir, job_subdir)
    submitted = []
    if not os.path.isdir(job_dir):
        return submitted

    for job in sorted(os.listdir(job_dir)):
        job_path = os.path.join(job_dir, job)
        if not job.endswith('.sh'):
            continue
        if dry_run:
            submitted.append((job_path, None))
            continue
        # Use subprocess to call sbatch and capture output
        try:
            proc = subprocess.run(['sbatch', job_path], capture_output=True, text=True, check=True)
            submitted.append((job_path, proc.stdout.strip()))
        except subprocess.CalledProcessError as e:
            submitted.append((job_path, e.stderr.strip() if e.stderr else str(e)))

    return submitted



def generate_scripts_for_dataset(
    dataset_dir,
    contigs_file,
    spacers_file,
    max_distance=5,
    threads=16,
    skip_tools=None,
    slurm_threads=16,
    slurm_opts=None,
    hyperfine=False
):
    """
    Generate bash and SLURM scripts for a dataset with specified max_distance.
    
    Args:
        dataset_dir: Path to the dataset directory
        contigs_file: Path to contigs FASTA file
        spacers_file: Path to spacers FASTA file
        max_distance: Maximum edit distance/mismatches for tools (default: 5)
        threads: Number of threads for bash scripts (default: 16)
        skip_tools: List of tool names to skip (default: None)
        slurm_threads: Number of threads for SLURM jobs (default: 16)
        slurm_opts: Dict of SLURM options to override (default: None)
        hyperfine: Whether to wrap commands in hyperfine (default: False)
    
    Returns:
        dict with keys 'tools' (list of tool names) and 'jobs' (list of job script paths)
    """
    import logging
    from pathlib import Path
    logger = logging.getLogger(__name__)
    from .tool_commands import load_tool_configs
    from .functions import create_bash_script
    dataset_dir = Path(dataset_dir)
    
    # Create necessary directories
    (dataset_dir / 'bash_scripts').mkdir(exist_ok=True, parents=True)
    (dataset_dir / 'raw_outputs').mkdir(exist_ok=True, parents=True)
    (dataset_dir / 'slurm_logs').mkdir(exist_ok=True, parents=True)
    (dataset_dir / 'job_scripts').mkdir(exist_ok=True, parents=True)
    
    logger.info(f"Generating scripts for {dataset_dir.name} (max_distance={max_distance})")
    
    # Load tool configs with max_distance (now called max_mismatches in the function)
    tools = load_tool_configs(
        results_dir=str(dataset_dir),
        threads=threads,
        contigs_file=contigs_file,
        spacers_file=spacers_file,
        max_mismatches=max_distance  # Parameter is called max_mismatches in the function
    )
    
    # Filter out skipped tools
    if skip_tools:
        tools = {k: v for k, v in tools.items() if k not in skip_tools}
        logger.info(f"Skipping tools: {skip_tools}")
    
    # Generate bash scripts
    created_tools = []
    for tool_name, tool in tools.items():
        try:
            create_bash_script(tool, str(dataset_dir), max_runs=1, warmups=0, hyperfine=hyperfine)
            created_tools.append(tool_name)
            logger.debug(f"  ✓ Created bash script for {tool_name}")
        except Exception as e:
            logger.error(f"  ✗ Failed to create script for {tool_name}: {e}")
    
    logger.info(f"Created {len(created_tools)} bash scripts")
    
    # Generate SLURM job scripts
    created_jobs = create_slurm_job_scripts(
        str(dataset_dir),
        slurm_opts=slurm_opts,
        threads=slurm_threads
    )
    logger.info(f"Created {len(created_jobs)} SLURM job scripts")
    
    return {
        'tools': created_tools,
        'jobs': created_jobs,
        'max_distance': max_distance
    }



if __name__ == "__main__":
    result_df = pyseff()
    print(result_df)
    result_df.write_csv("slurm_efficiency.tsv", separator="\t")


