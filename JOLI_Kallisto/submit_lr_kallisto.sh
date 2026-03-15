#!/usr/bin/env bash
# =============================================================================
# submit_lr_kallisto.sh
#
# Description:
#   Thin SLURM submission wrapper for run_lr_kallisto.sh.
#   Creates the timestamped results directory first so that SLURM can write
#   the running.log there from the start, then calls sbatch with resource
#   flags pointing at run_lr_kallisto.sh.
#
#   The pipeline logic lives entirely in run_lr_kallisto.sh — this file only
#   sets SLURM resource flags, keeping the two concerns cleanly separated.
#
# Usage:
#   bash submit_lr_kallisto.sh
#
# Inputs:  None (resources configured in the CONFIG section below)
# Outputs: Submits run_lr_kallisto.sh as a SLURM job; prints the job ID.
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIG — SLURM resources; adjust as needed
# =============================================================================

JOB_NAME=lr_kallisto
MEM=100G
CPUS=32
TIME=24:00:00
MAIL_USER=atalukder@nygenome.org
MAIL_TYPE=END,FAIL

# Base output directory — must match OUTPUT_BASE inside run_lr_kallisto.sh
OUTPUT_BASE=/gpfs/commons/home/atalukder/RNA_Splicing/files/results

# Path to the pipeline script (relative to this file's location)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_SCRIPT="${SCRIPT_DIR}/run_lr_kallisto.sh"

# =============================================================================
# Pre-create the timestamped run directory so SLURM can log there immediately.
# The pipeline script will also call mkdir -p on the same path, which is safe.
# =============================================================================

TIMESTAMP=$(date +"%Y_%m_%d__%H_%M_%S")
RUN_DIR="${OUTPUT_BASE}/exprmnt_${TIMESTAMP}"
mkdir -p "${RUN_DIR}"

SLURM_OUT="${RUN_DIR}/running.log"

echo "Submitting job..."
echo "  Script    : ${PIPELINE_SCRIPT}"
echo "  Run dir   : ${RUN_DIR}"
echo "  SLURM log : ${SLURM_OUT}"

# =============================================================================
# Submit — #SBATCH flags are passed on the command line here, not in the script,
# so run_lr_kallisto.sh stays free of SLURM-specific headers and runs cleanly
# with plain `bash run_lr_kallisto.sh` locally.
# =============================================================================

sbatch \
    --job-name="${JOB_NAME}" \
    --mem="${MEM}" \
    --cpus-per-task="${CPUS}" \
    --time="${TIME}" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    --output="${SLURM_OUT}" \
    "${PIPELINE_SCRIPT}"

echo "Job submitted. Monitor with: squeue -u $(whoami)"
echo "Results will appear in: ${RUN_DIR}"
