"""
generate_slurm_jobs.py
======================
SLURM job generator for the JOLI-Kallisto pipeline.

Follows the same structure as AT_code/generate_bash.py:
  1. Generates a "prg file" — a bash script that activates conda and runs the pipeline.
  2. Generates a "SLURM file" — a script with #SBATCH headers that calls the prg file.
  3. Saves both to the cluster job submission directory.
  4. Submits each via `sbatch`.

Run this script to submit jobs:
  conda activate Joli_kallisto
  cd JOLI_Kallisto/
  python generate_slurm_jobs.py

Inputs:
  - JOBS list in CONFIG: one dict per job to submit.
      Required key:
        "script"    (str) — bash script to run (filename relative to this script's
                            directory, or an absolute path).
      Optional per-job overrides (fall back to the defaults below if omitted):
        "job_name"  (str) — SLURM job label; auto-derived from script name if absent.
        "hour"      (int) — wall-clock limit in hours.
        "memory"    (int) — RAM in GB.
        "ncpus"     (int) — CPUs per task.
        "conda_env" (str) — conda environment to activate.

  - Default resource settings: HOUR, MEMORY, NCPUS, CONDA_ENV (applied when a job
    dict omits the corresponding key).

Outputs:
  - prg_<job_name>_<timestamp>.sh    — saved to JOB_DIR
  - slurm_<job_name>_<timestamp>.sh  — saved to JOB_DIR
  - SLURM job(s)                     — submitted via sbatch; prints job IDs
"""

import os
import time


# ============================================================
# CONFIG — edit everything here; do not touch logic below
# ============================================================

# --- Jobs to submit ---
# Each entry is a dict.  Required key: "script" (filename relative to this
# script's directory, or an absolute path).  All other keys are optional.
#
# Examples:
#   Submit both pipelines with default resources:
#     JOBS = [
#         {"script": "run_joli_kallisto.sh"},
#         {"script": "run_lr_kallisto.sh"},
#     ]
#
#   Submit only JOLI with a custom time limit:
#     JOBS = [
#         {"script": "run_joli_kallisto.sh", "hour": 48},
#     ]
#
#   Submit both, each with different memory:
#     JOBS = [
#         {"script": "run_joli_kallisto.sh", "memory": 150},
#         {"script": "run_lr_kallisto.sh",   "memory": 80},
#     ]
JOBS = [
    {"script": "run_joli_kallisto.sh"},
    {"script": "run_lr_kallisto.sh"},
]

# --- Default SLURM resources (used when a job dict omits the key) ---
HOUR      = 10          # wall-clock limit in hours
MEMORY    = 100         # RAM in GB
NCPUS     = 32          # CPUs per task
CONDA_ENV = "Joli_kallisto"

# --- Email notifications ---
MAIL_USER = "atalukder@nygenome.org"
MAIL_TYPE = "END,FAIL"

# --- Directories ---
# Where prg_*.sh and slurm_*.sh are saved
JOB_DIR = "/gpfs/commons/home/atalukder/RNA_Splicing/files/cluster_job_submission_files"

# Where SLURM stdout/stderr logs go (%j is replaced by SLURM job ID)
SLURM_LOG_DIR = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

# ============================================================
# END CONFIG
# ============================================================


THIS_DIR  = os.path.dirname(os.path.realpath(__file__))
TIMESTAMP = time.strftime("%Y_%m_%d__%H_%M_%S")


def resolve_script_path(script: str) -> str:
    """
    Return the absolute path to a pipeline script.

    If `script` is already absolute, return it as-is.
    Otherwise, treat it as relative to this script's directory.

    Args:
        script (str): Script filename or absolute path.

    Returns:
        str: Absolute path to the script.
    """
    if os.path.isabs(script):
        return script
    return os.path.join(THIS_DIR, script)


def job_name_from_script(script: str) -> str:
    """
    Derive a short SLURM job label from a script filename.

    Strips the directory and the .sh extension.
    Example: "run_joli_kallisto.sh" -> "run_joli_kallisto"

    Args:
        script (str): Script filename or path.

    Returns:
        str: Short label suitable for a SLURM job name.
    """
    return os.path.splitext(os.path.basename(script))[0]


def create_prg_file(pipeline_script_path: str, prg_file_path: str, conda_env: str) -> str:
    """
    Write a bash script that activates the conda environment and runs the pipeline.

    This is the actual work script. It mirrors the structure used in
    AT_code/generate_bash.py: source ~/.bashrc, conda activate, then run pipeline.

    Args:
        pipeline_script_path (str): Absolute path to the pipeline bash script.
        prg_file_path (str): Where to write the prg file.
        conda_env (str): Conda environment to activate before running the script.

    Returns:
        str: Path to the written prg file.
    """
    content = (
        "#!/bin/bash\n"
        "set -e\n"
        "cd $HOME\n"
        "source ~/.bashrc\n"
        f"conda activate {conda_env}\n"
        f"bash {pipeline_script_path}\n"
    )
    with open(prg_file_path, "w") as f:
        f.write(content)
    os.chmod(prg_file_path, 0o755)
    print(f"    prg  : {prg_file_path}")
    return prg_file_path


def create_slurm_file(
    prg_file_path: str,
    job_name: str,
    slurm_file_path: str,
    slurm_out_path: str,
    hour: int,
    memory: int,
    ncpus: int,
) -> str:
    """
    Write a SLURM submission script with #SBATCH headers.

    The script's only command is calling the prg file, keeping all pipeline
    logic separate. Mirrors AT_code/generate_bash.py::create_slurm_file().

    Args:
        prg_file_path (str): Absolute path to the prg file.
        job_name (str): SLURM job name (shown in squeue).
        slurm_file_path (str): Where to write the SLURM script.
        slurm_out_path (str): Path for SLURM stdout/stderr (%j = job ID).
        hour (int): Wall-clock limit in hours.
        memory (int): RAM in GB.
        ncpus (int): CPUs per task.

    Returns:
        str: Path to the written SLURM file.
    """
    content = (
        "#!/bin/bash\n"
        "##ENVIRONMENT SETTINGS; REPLACE WITH CAUTION\n"
        "##NECESSARY JOB SPECIFICATIONS\n"
        f"#SBATCH --job-name={job_name}\n"
        f"#SBATCH --time={hour}:00:00\n"
        f"#SBATCH --mem={memory}G\n"
        f"#SBATCH --cpus-per-task={ncpus}\n"
        f"#SBATCH --mail-type={MAIL_TYPE}\n"
        f"#SBATCH --mail-user={MAIL_USER}\n"
        f"#SBATCH --output={slurm_out_path}\n"
        f"{prg_file_path}\n"
    )
    with open(slurm_file_path, "w") as f:
        f.write(content)
    os.chmod(slurm_file_path, 0o755)
    print(f"    slurm: {slurm_file_path}")
    return slurm_file_path


def submit_job(job: dict) -> None:
    """
    Generate and submit one SLURM job from a job dict.

    Reads the required "script" key and optional resource override keys,
    creates a prg file and a SLURM file, then calls sbatch.

    Args:
        job (dict): Job specification with keys:
            script    (str, required) — script filename or absolute path.
            job_name  (str, optional) — SLURM job label.
            hour      (int, optional) — wall-clock hours.
            memory    (int, optional) — RAM in GB.
            ncpus     (int, optional) — CPU count.
            conda_env (str, optional) — conda environment name.

    Returns:
        None
    """
    script_path = resolve_script_path(job["script"])

    if not os.path.isfile(script_path):
        print(f"  [ERROR] script not found: {script_path}")
        return

    # Per-job overrides, falling back to module-level defaults
    label     = job.get("job_name",  job_name_from_script(job["script"]))
    hour      = job.get("hour",      HOUR)
    memory    = job.get("memory",    MEMORY)
    ncpus     = job.get("ncpus",     NCPUS)
    conda_env = job.get("conda_env", CONDA_ENV)

    full_name       = f"{label}_{TIMESTAMP}"
    prg_file_path   = os.path.join(JOB_DIR, f"prg_{full_name}.sh")
    slurm_file_path = os.path.join(JOB_DIR, f"slurm_{full_name}.sh")
    slurm_out_path  = os.path.join(SLURM_LOG_DIR, f"slurm_{full_name}.%j.out")

    print(f"\n  Job: {label}")
    print(f"    script   : {script_path}")
    print(f"    resources: {ncpus} CPUs  {memory}G RAM  {hour}h")
    print(f"    conda env: {conda_env}")

    create_prg_file(
        pipeline_script_path=script_path,
        prg_file_path=prg_file_path,
        conda_env=conda_env,
    )
    create_slurm_file(
        prg_file_path=prg_file_path,
        job_name=full_name,
        slurm_file_path=slurm_file_path,
        slurm_out_path=slurm_out_path,
        hour=hour,
        memory=memory,
        ncpus=ncpus,
    )

    ret = os.system(f"sbatch {slurm_file_path}")
    if ret == 0:
        print(f"    [OK] submitted")
    else:
        print(f"    [ERROR] sbatch exit code: {ret}")


def main():
    """
    Entry point: iterate over the JOBS list and submit each as a separate SLURM job.

    Returns:
        None
    """
    os.makedirs(JOB_DIR, exist_ok=True)

    print("=" * 60)
    print(f"SLURM job generator — {TIMESTAMP}")
    print(f"Jobs to submit: {len(JOBS)}")
    print("=" * 60)

    if not JOBS:
        print("JOBS list is empty — nothing to submit.")
        return

    for job in JOBS:
        submit_job(job)

    print("\n" + "=" * 60)
    print("All jobs submitted.")
    print("Monitor with: squeue -u $(whoami)")
    print("=" * 60)


if __name__ == "__main__":
    main()
