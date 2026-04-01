#!/usr/bin/env bash
# submit_multisample_joli.sh
# ==========================
# Thin SLURM wrapper for run_multisample_joli.sh.
# All pipeline logic and sample config lives in run_multisample_joli.sh.
#
# Usage:
#   bash submit_multisample_joli.sh        # submit to SLURM cluster
#   bash run_multisample_joli.sh           # run locally (no SLURM)
#   python main_multisample_joli.py ...    # run MAP EM only (bustools output exists)

# ============================================================
# CONFIG — SLURM resources only; edit pipeline config in run_multisample_joli.sh
# ============================================================
JOB_NAME="joli_multisample"
MEM="100G"
CPUS=32
TIME="24:00:00"
MAIL_USER="atalukder@nygenome.org"
MAIL_TYPE="END,FAIL"
CONDA_ENV="Joli_kallisto"

# Free-text description of this submission (for your records; echoed at submit time).
# The actual experiment comment that gets saved in results lives in run_multisample_joli.sh.
EXPERIMENT_COMMENT=""
# ============================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINE_SCRIPT="${SCRIPT_DIR}/run_multisample_joli.sh"

echo "Submitting JOLI multi-sample pipeline to SLURM..."
echo "  Pipeline script: ${PIPELINE_SCRIPT}"
echo "  Conda env:       ${CONDA_ENV}"
echo "  Resources:       ${CPUS} CPUs, ${MEM} RAM, ${TIME} time limit"
if [[ -n "${EXPERIMENT_COMMENT}" ]]; then
    echo "  Comment:         ${EXPERIMENT_COMMENT}"
fi

sbatch \
    --job-name="${JOB_NAME}" \
    --mem="${MEM}" \
    --cpus-per-task="${CPUS}" \
    --time="${TIME}" \
    --mail-user="${MAIL_USER}" \
    --mail-type="${MAIL_TYPE}" \
    --output="${SCRIPT_DIR}/slurm_%j.out" \
    --error="${SCRIPT_DIR}/slurm_%j.err" \
    --wrap="conda run -n ${CONDA_ENV} bash ${PIPELINE_SCRIPT}"

echo "Job submitted."
