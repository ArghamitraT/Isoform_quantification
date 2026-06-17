#!/usr/bin/env bash
# =============================================================================
# run.sh — Configure and launch the RNA-seq simulation pipeline.
#
# THIS is the file you edit. All pipeline options live here.
# main.py reads them as CLI arguments — do not edit main.py directly.
#
# Requirements:
#   Conda env  : lrgsp_simulation
#                  (numpy, pandas, scipy, pysam, biopython, nanocount)
#                  Create from: conda env create -f env/lrgsp_simulation.yml
#   Modules    : minimap2/2.29            — PacBio and ONT alignment
#                SAMtools/1.21            — BAM/FASTA processing
#                RSEM/1.3.3-foss-2022a   — Illumina simulation only
#   ONT only   : lrgsp_simulation_2 (clone with HTSeq)
#                  conda create -n lrgsp_simulation_2 --clone lrgsp_simulation
#                  conda activate lrgsp_simulation_2 && pip install HTSeq
#                  Switch CONDA_ENV below to lrgsp_simulation_2 if running ONT.
#
# Estimated runtime (16 threads, full read counts):
#   Phase 1 — prepare reference     : ~TBD
#   Phase 2 — generate abundances   : ~TBD
#   Phase 3a — Illumina (100M reads): ~TBD
#   Phase 3b — PacBio  (10M reads)  : ~TBD
#   Phase 3c — ONT     (30M reads)  : ~TBD (separate env — see above)
#   Total                           : ~TBD
#
# Run:
#   cd /gpfs/commons/home/atalukder/RNA_Splicing
#   bash code/Simulations/run.sh
# =============================================================================

# ============================================================
# CONFIG — edit everything here; do not touch the logic below
# ============================================================

# --- Environment ---
CONDA_ENV="lrgsp_simulation"

# --- HPC modules to load ---
MODULE_MINIMAP2="minimap2/2.29"
MODULE_SAMTOOLS="SAMtools/1.21"
MODULE_RSEM="RSEM/1.3.3-foss-2022a"

# --- Input files ---
GENOME="data/ReqToSimulateData/human.genome.fasta"
GTF="data/ReqToSimulateData/human.annotation.gtf"
RSEM_MODEL="data/ReqToSimulateData/Illumina150.RNA.model"
SQANTI_PREFIX=""          # leave empty to skip novel isoforms
N_NOVEL=0                 # number of novel isoforms to add from SQANTI
POLYA_LEN=100             # poly-A tail length (bp)

# --- Phase toggles (1 = run, 0 = skip) ---
RUN_PHASE1=1              # 0 = skip if ref_data/ already exists from a prior run
RUN_PHASE2=1              # 0 = skip if abundances.tsv already exists
RUN_ILLUMINA=1            # Illumina paired-end reads (needs RSEM module)
RUN_PACBIO=1              # PacBio CCS reads (needs minimap2 + SAMtools)
RUN_ONT=0                 # ONT reads (needs lrgsp_simulation_2 with HTSeq)

# --- Abundance parameters ---
ABUNDANCE_MODE="lognormal"  # lognormal | custom | uniform
CUSTOM_TPM_FILE=""          # TSV (transcript_id, TPM) — required if mode=custom
CANCER_GENES_FILE=""        # gene list to up-regulate — optional, mode=custom only
CANCER_FOLD_CHANGE=5.0      # fold-change multiplier for cancer genes
DROPOUT_FRACTION=0.3        # fraction of transcripts silenced (lognormal only)
ABUNDANCE_MU=2.0            # log-normal mean  (lognormal only)
ABUNDANCE_SIGMA=2.0         # log-normal sigma (lognormal only)

# --- Read counts (full run) ---
ILLUMINA_COUNT=100000000    # paired-end read pairs
PACBIO_COUNT=10000000       # reads
ONT_COUNT=30000000          # reads
ONT_TYPE="cDNA"             # cDNA or dRNA

# --- Run settings ---
THREADS=16
SEED=42
KEEP_ISOFORM_IDS=0          # 1 = preserve original isoform IDs in PacBio read names
NOISE_READS=0               # 1 = include unaligned/noise reads in ONT output

# --- Test mode ---
# Set TEST_MODE=1 to override read counts with tiny values for a quick sanity check.
TEST_MODE=1
TEST_ILLUMINA_COUNT=1000
TEST_PACBIO_COUNT=1000
TEST_ONT_COUNT=1000

# ============================================================
# Do not edit below this line
# ============================================================

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MAIN_PY="${PROJECT_ROOT}/code/Simulations/main.py"

echo "============================================================"
echo "  RNA-seq Simulation Pipeline — run.sh"
echo "  Project root : ${PROJECT_ROOT}"
echo "  Date         : $(date)"
echo "  Test mode    : ${TEST_MODE}"
echo "============================================================"

# ── Activate conda env first (must come before module loads) ─────────────────
# Modules (especially RSEM) inject their own Python into PATH and will override
# conda if activated after. By activating conda first and capturing its Python
# path explicitly, we guarantee the right interpreter is used regardless of
# what modules do to PATH afterwards.
echo "[run.sh] Activating conda env: ${CONDA_ENV}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
PYTHON="${CONDA_PREFIX}/bin/python"   # explicit path — survives module PATH changes
echo "[run.sh] Python: ${PYTHON} ($(${PYTHON} --version 2>&1))"

# ── Load HPC modules (after conda, so they don't override Python) ─────────────
echo "[run.sh] Loading modules..."
module load "${MODULE_MINIMAP2}"
module load "${MODULE_SAMTOOLS}"
[[ "${RUN_ILLUMINA}" -eq 1 ]] && module load "${MODULE_RSEM}"
echo "[run.sh] Modules loaded."

# ── Build argument list ───────────────────────────────────────────────────────
ARGS=(
    --genome          "${GENOME}"
    --gtf             "${GTF}"
    --rsem_model      "${RSEM_MODEL}"
    --n_novel         "${N_NOVEL}"
    --polya_len       "${POLYA_LEN}"
    --abundance_mode  "${ABUNDANCE_MODE}"
    --cancer_fold_change "${CANCER_FOLD_CHANGE}"
    --dropout_fraction   "${DROPOUT_FRACTION}"
    --abundance_mu       "${ABUNDANCE_MU}"
    --abundance_sigma    "${ABUNDANCE_SIGMA}"
    --illumina_count  "${ILLUMINA_COUNT}"
    --pacbio_count    "${PACBIO_COUNT}"
    --ont_count       "${ONT_COUNT}"
    --ont_type        "${ONT_TYPE}"
    --threads         "${THREADS}"
    --seed            "${SEED}"
    --test_illumina_count "${TEST_ILLUMINA_COUNT}"
    --test_pacbio_count   "${TEST_PACBIO_COUNT}"
    --test_ont_count      "${TEST_ONT_COUNT}"
)

# Optional string args (only pass if non-empty)
[[ -n "${SQANTI_PREFIX}"     ]] && ARGS+=(--sqanti_prefix    "${SQANTI_PREFIX}")
[[ -n "${CUSTOM_TPM_FILE}"   ]] && ARGS+=(--custom_tpm_file  "${CUSTOM_TPM_FILE}")
[[ -n "${CANCER_GENES_FILE}" ]] && ARGS+=(--cancer_genes_file "${CANCER_GENES_FILE}")

# Phase skip flags
[[ "${RUN_PHASE1}"    -eq 0 ]] && ARGS+=(--skip_phase1)
[[ "${RUN_PHASE2}"    -eq 0 ]] && ARGS+=(--skip_phase2)
[[ "${RUN_ILLUMINA}"  -eq 0 ]] && ARGS+=(--skip_illumina)
[[ "${RUN_PACBIO}"    -eq 0 ]] && ARGS+=(--skip_pacbio)
[[ "${RUN_ONT}"       -eq 0 ]] && ARGS+=(--skip_ont)

# Boolean flags
[[ "${KEEP_ISOFORM_IDS}" -eq 1 ]] && ARGS+=(--keep_isoform_ids)
[[ "${NOISE_READS}"      -eq 1 ]] && ARGS+=(--noise_reads)
[[ "${TEST_MODE}"        -eq 1 ]] && ARGS+=(--test_mode)

# ── Launch ────────────────────────────────────────────────────────────────────
echo "[run.sh] Launching pipeline..."
cd "${PROJECT_ROOT}"
"${PYTHON}" "${MAIN_PY}" "${ARGS[@]}"

echo "[run.sh] Done."
