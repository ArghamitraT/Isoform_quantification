#!/usr/bin/env python3
"""
submit_JKMS_slurm.py
=====================
Generates and submits SLURM jobs for the JOLI-Kallisto multi-sample
real-data runs. Mirrors AT_code/generate_bash.py's pattern: for each run,
write a separate "prg" file (conda activate + the pipeline invocation) and a
"slurm" file (the #SBATCH resource header that calls the prg file), then
submit both via sbatch.

Each prg file calls run_multisample_joli.sh (not main_multisample_joli.py
directly) so the pipeline's own per-read-type logic stays the single source
of truth: it auto-resolves eff_len_mode ("uniform" for long, "fld" for
short), em_include_single_tx (false for long, true for short), kmer size,
and index per sample, and skips kallisto bus/bustools automatically when the
TCC cache already exists. This script only generates the small per-run
override (SAMPLES array + THREADS + EXPERIMENT_COMMENT + RUN_DIR) that
run_multisample_joli.sh sources via JOLI_MS_OVERRIDE.

Each run's exprmnt_YYYY_MM_DD__HH_MM_SS output folder name is decided here,
at generation time (one shared base timestamp + a 1-second offset per run),
and passed down via RUN_DIR -> --run_dir. This keeps the plain, standard
folder-name format (no PID/job-id suffix) while guaranteeing the runs in a
batch never race on an independently-generated timestamp the way concurrent
SLURM jobs calling datetime.now() themselves could.

Inputs:
    CONFIG section below -- SLURM resources (hour/memory/nthred), the LR+SR
    sample pairs to run (runs), conda env name, and raw-read directories.

Outputs:
    <job_path>/override_JKms_run<N>_<hash>_<timestamp>.sh
        -- per-run SAMPLES/THREADS/EXPERIMENT_COMMENT override
    <job_path>/prg_JKms_run<N>_<hash>_<timestamp>.sh
        -- per-run command script (conda activate + pipeline invocation)
    <job_path>/slurm_JKms_run<N>_<hash>_<timestamp>.sh
        -- per-run sbatch resource header
    <job_path>/out_JKms_run<N>_<hash>_<timestamp>.%j
        -- SLURM stdout/stderr for the job
    Each run writes to its own pre-assigned files/results/exprmnt_<timestamp>/
    folder (standard exprmnt_YYYY_MM_DD__HH_MM_SS format, no suffix).

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
    conda activate NanoCount_5
    python scripts/submit_JKMS_slurm.py
"""
import os
import time
import random
from datetime import datetime, timedelta

# ============================================================
# CONFIG -- edit these variables before running; do not edit below
# ============================================================
hour      = 4                  # SLURM wall-clock time limit (hours)
memory    = 80                 # GB
nthred    = 4                  # number of CPU (also passed to the pipeline as THREADS)
conda_env = "Joli_kallisto"
mail_user = "atalukder@nygenome.org"
mail_type = "END,FAIL"

pipeline_script = "/gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto/scripts/run_multisample_joli.sh"
job_path        = "/gpfs/commons/home/atalukder/RNA_Splicing/files/cluster_job_submission_files/JK_multisample"
results_base    = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

lr_dir = "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long"
sr_dir = "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/short"

# Run definitions: (lr_sample, sr_sample)
runs = [
    ("flnc_31", "ill_51"),
    ("flnc_32", "ill_52"),
    ("flnc_51", "ill_31"),
    ("flnc_52", "ill_32"),
]

slurm_job_prefix = "JKms"
# ============================================================

trimester = time.strftime("_%Y_%m_%d__%H_%M_%S")
base_time = datetime.now()   # shared base for pre-assigning each run's exprmnt_* folder


def assign_run_dir(run_id):
    """
    Description:
        Pre-assign this run's exprmnt_YYYY_MM_DD__HH_MM_SS output folder at
        generation time, offsetting by run_id seconds from a shared base
        timestamp so runs in the same batch never collide.
    Args:
        run_id (int): 1-based run index.
    Returns:
        str: absolute path, e.g. ".../results/exprmnt_2026_07_19__18_05_01".
    """
    run_time = base_time + timedelta(seconds=run_id)
    timestamp = run_time.strftime("%Y_%m_%d__%H_%M_%S")
    return os.path.join(results_base, f"exprmnt_{timestamp}")


def create_job_dir(path):
    """
    Description:
        Create the directory that holds generated prg/slurm/override files, if missing.
    Args:
        path (str): directory to create.
    Returns:
        str: the same path, guaranteed to exist.
    """
    os.makedirs(path, exist_ok=True)
    return path


def get_file_name(kind, ext=True):
    """
    Description:
        Build a timestamped file name for a generated script.
    Args:
        kind (str): short label identifying the file's role, e.g. "prg_JKms_run1_123".
        ext (bool): append ".sh" when True.
    Returns:
        str: file name, e.g. "prg_JKms_run1_123_2026_07_19__17_30_00.sh".
    """
    file_name = f"{kind}{trimester}"
    if ext:
        file_name = f"{file_name}.sh"
    return file_name


def create_override_file(override_path, run_id, lr, sr, run_dir):
    """
    Description:
        Write the SAMPLES/THREADS/EXPERIMENT_COMMENT/RUN_DIR override that
        run_multisample_joli.sh sources via JOLI_MS_OVERRIDE. All per-sample
        derivation (eff_len_mode, em_include_single_tx, kmer size, index,
        cache-hit skip) is left to the pipeline script itself -- not
        duplicated here.
    Args:
        override_path (str): output path for the override script.
        run_id (int): 1-based run index (used in the experiment comment).
        lr (str): long-read sample stem, e.g. "flnc_31".
        sr (str): short-read sample stem, e.g. "ill_51".
        run_dir (str): pre-assigned exprmnt_* output folder for this run.
    Returns:
        str: the override_path written.
    """
    content = (
        f'EXPERIMENT_COMMENT="JK MS real — Run{run_id}: {lr}(LR) + {sr}(SR)"\n'
        f'DATA_TYPE="real"\n'
        f'REPLICATE_PAIRS=()\n'
        f'GT_MAP=()\n'
        f'THREADS={nthred}\n'
        f'RUN_DIR="{run_dir}"\n'
        f'SAMPLES=(\n'
        f'    "{lr}  long   {lr_dir}  {lr}.fastq"\n'
        f'    "{sr}  short  {sr_dir}  {sr}_R1.fastq  {sr}_R2.fastq"\n'
        f')\n'
    )
    with open(override_path, "w") as f:
        f.write(content)
    return override_path


def create_prg_file(prg_file_path, override_path):
    """
    Description:
        Write the command script: activates the conda env, points the
        pipeline at this run's override file, and invokes it.
    Args:
        prg_file_path (str): output path for the command script.
        override_path (str): override file this run's pipeline call should source.
    Returns:
        str: the prg_file_path written.
    """
    header = (
        "#!/bin/bash\n"
        "set -e\n"
        "source ~/.bashrc\n"
        f"conda activate {conda_env}\n"
        f'export JOLI_MS_OVERRIDE="{override_path}"\n'
        f'bash "{pipeline_script}"\n'
    )
    with open(prg_file_path, "w") as f:
        f.write(header)
    return prg_file_path


def create_slurm_file(prg_file_path, job_name, slurm_file_path):
    """
    Description:
        Write the #SBATCH resource header that invokes the prg file.
    Args:
        prg_file_path (str): command script this job should run.
        job_name (str): SLURM job name.
        slurm_file_path (str): output path for the sbatch script.
    Returns:
        str: the slurm_file_path written.
    """
    header = (
        "#!/bin/bash\n"
        "##NECESSARY JOB SPECIFICATIONS\n"
        f"#SBATCH --job-name={job_name}\n"
        f"#SBATCH --time={hour}:00:00\n"
        f"#SBATCH --mem={memory}G\n"
        f"#SBATCH --cpus-per-task={nthred}\n"
        f"#SBATCH --mail-type={mail_type}\n"
        f"#SBATCH --mail-user={mail_user}\n"
        f"#SBATCH --output={job_path}/out_{job_name}.%j\n"
        f'bash "{prg_file_path}"\n'
    )
    with open(slurm_file_path, "w") as f:
        f.write(header)
    return slurm_file_path


def gen_and_submit():
    """
    Description:
        For each (lr, sr) pair in `runs`, generate override/prg/slurm files
        and submit the job via sbatch.
    """
    create_job_dir(job_path)

    for run_id, (lr, sr) in enumerate(runs, start=1):
        hash_obj = random.getrandbits(25)
        kind = f"{slurm_job_prefix}_run{run_id}_{hash_obj}"

        override_path   = os.path.join(job_path, get_file_name(kind=f"override_{kind}"))
        prg_file_path    = os.path.join(job_path, get_file_name(kind=f"prg_{kind}"))
        slurm_file_path  = os.path.join(job_path, get_file_name(kind=f"slurm_{kind}"))
        run_dir          = assign_run_dir(run_id)

        create_override_file(override_path, run_id, lr, sr, run_dir)
        create_prg_file(prg_file_path, override_path)
        create_slurm_file(prg_file_path, job_name=kind, slurm_file_path=slurm_file_path)

        os.system(f"chmod u+x {prg_file_path}")
        os.system(f"chmod u+x {slurm_file_path}")

        print(f"Submitting Run {run_id}: {lr}(LR) + {sr}(SR) -> {run_dir}")
        os.system(f"sbatch {slurm_file_path}")


def main():
    gen_and_submit()


if __name__ == "__main__":
    main()
