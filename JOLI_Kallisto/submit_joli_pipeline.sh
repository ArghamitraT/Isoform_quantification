#!/usr/bin/env bash
# submit_joli_pipeline.sh
# =======================
# Thin SLURM wrapper for run_joli_kallisto.sh.
# All pipeline logic and sample config lives in run_joli_kallisto.sh.
#
# Usage:
#   bash submit_joli_pipeline.sh           # submit to SLURM cluster
#   bash run_joli_kallisto.sh              # run locally (standard)
#   python main_joli.py --sample_dir ...  # run JOLI EM only (breakpoints)

# ============================================================
# CONFIG — SLURM resources only; edit pipeline config in run_joli_kallisto.sh
# ============================================================
JOB_NAME="joli_pipeline"
MEM="100G"
CPUS=32
TIME="24:00:00"
MAIL_USER="atalukder@nygenome.org"
MAIL_TYPE="END,FAIL"
CONDA_ENV="NanoCount_5"
# ============================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINE_SCRIPT="${SCRIPT_DIR}/run_joli_kallisto.sh"

echo "Submitting JOLI-Kallisto pipeline to SLURM..."
echo "  Pipeline script: ${PIPELINE_SCRIPT}"
echo "  Conda env:       ${CONDA_ENV}"
echo "  Resources:       ${CPUS} CPUs, ${MEM} RAM, ${TIME} time limit"

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
