#!/usr/bin/env bash
# =============================================================================
# run_lr_kallisto.sh
#
# Description:
#   One-click lr-kallisto pipeline for long-read (PacBio/ONT) or short-read
#   RNA-seq isoform quantification.
#
#   Steps (in order):
#     0. [optional] kallisto index  — build index from transcriptome FASTA
#     1. kallisto bus               — align reads → output.bus
#     2. bustools sort              — sort .bus file
#     3. bustools count             — generate TCC count matrix
#     4. kallisto quant-tcc         — EM quantification → abundance.tsv
#
#   All steps run for every sample in the SAMPLES array.
#   Set BUILD_INDEX=1 to rebuild the index before the sample loop;
#   set BUILD_INDEX=0 to skip and use the existing INDEX_FILE.
#
# Inputs (configured in CONFIG section below):
#   - SAMPLES array: name, reads directory, reads file(s)
#   - kallisto binary, bustools binary
#   - TRANSCRIPTOME_FASTA, INDEX_FILE, T2G_FILE
#   - READ_TYPE, KMER_SIZE, BUILD_INDEX flag
#
# Outputs:
#   /gpfs/.../files/results/exprmnt_{timestamp}/
#   ├── experiment_description.log   # config dump + sample list
#   ├── running.log                  # combined stdout/stderr for all steps
#   ├── runtime.txt                  # total wall-clock time
#   ├── code_snapshot/
#   │   └── run_lr_kallisto.sh       # exact copy of this script
#   ├── {sample_name_1}/
#   │   ├── output.bus, sorted.bus
#   │   ├── count.mtx, count.ec.txt, transcripts.txt
#   │   ├── abundance.tsv
#   │   └── run_info.json
#   └── {sample_name_2}/  ...
#
# Usage:
#   bash run_lr_kallisto.sh            # local run
#   sbatch run_lr_kallisto.sh          # SLURM run (headers ignored by bash)
#   bash submit_lr_kallisto.sh         # submit via thin SLURM wrapper
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIG — edit everything here; do not touch the pipeline logic below
# =============================================================================

# --- Binaries ---
KALLISTO=/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/kallisto/build/src/kallisto
BUSTOOLS=/gpfs/commons/home/atalukder/miniconda3/envs/mpaqt/bin/bustools

# --- Reference files ---
LR_KALLISTO_BASE=/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto

# Transcriptome FASTA used to build the index (only needed when BUILD_INDEX=1)
# Real data:
TRANSCRIPTOME_FASTA="${LR_KALLISTO_BASE}/transcriptome.fasta"
# Simulation (uncomment to switch):
# TRANSCRIPTOME_FASTA=/gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result/ref_data/human.transcripts.fasta

# Path for the kallisto index (built here if BUILD_INDEX=1, otherwise must already exist)
INDEX_FILE="${LR_KALLISTO_BASE}/new_index.idx"

# Transcript-to-gene map
T2G_FILE="${LR_KALLISTO_BASE}/t2g.txt"

# --- Index settings ---
# Set to 1 to build/rebuild the index from TRANSCRIPTOME_FASTA before running samples. read this link: https://kallisto.readthedocs.io/en/latest/index/index_generation.html#index-generation
# Set to 0 to skip index building and use the existing INDEX_FILE.
BUILD_INDEX=0

# k-mer size for index building: 63 for long-read (PacBio/ONT), 31 for short-read
KMER_SIZE=63

# --- Run settings ---
READ_TYPE=long          # "long" (PacBio/ONT) or "short" (paired Illumina)
THREADS=32
THRESHOLD=0.8           # kallisto bus --threshold

# --- Output base directory ---
OUTPUT_BASE=/gpfs/commons/home/atalukder/RNA_Splicing/files/results

# --- Sample array ---
# Format (long-read):  "sample_name  reads_dir  reads_file"
# Format (short-read): "sample_name  reads_dir  reads_file_R1  reads_file_R2"
# Each entry is a single quoted string; fields are space-separated.
# Comment out any sample you do not want to run.
SAMPLES=(
    # "sim2  /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/human_simulated_job_correct/  PacBio.simulated.fasta"
    # "ds_52_furtherDownsampled  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long/downsampled/ds_52_furtherDownsampled.fastq"
    # "flnc_01  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long/  flnc_01.fastq"
    # "flnc_02  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long/  flnc_02.fastq"
    # "flnc_03  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long/  flnc_03.fastq"
    # "flnc_31  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long/  flnc_31.fastq"
    # "flnc_32  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long/  flnc_32.fastq"
    "sim1  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/  ds_100_num1_aln_01_long.fasta"
    "sim2  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/  ds_100_num1_aln_21_long.fasta"
)

# =============================================================================
# END CONFIG
# =============================================================================

# =============================================================================
# SETUP — create timestamped run directory and log initial config
# =============================================================================

TIMESTAMP=$(date +"%Y_%m_%d__%H_%M_%S")
RUN_DIR="${OUTPUT_BASE}/exprmnt_${TIMESTAMP}"
mkdir -p "${RUN_DIR}"

LOG="${RUN_DIR}/running.log"
EXPERIMENT_LOG="${RUN_DIR}/experiment_description.log"

# First line of every running.log: which script produced this log
echo "Script: $(realpath "$0")" > "${LOG}"

# Start wall-clock timer
RUN_START=$(date +%s)

# Activate conda environment (no-op when already active or not available)
source ~/.bashrc 2>/dev/null || true
conda activate Joli_kallisto 2>/dev/null || true

# Write experiment description / config dump
{
    echo "===== lr-kallisto experiment ====="
    echo "Timestamp     : ${TIMESTAMP}"
    echo "Run directory : ${RUN_DIR}"
    echo ""
    echo "--- Config ---"
    echo "KALLISTO           : ${KALLISTO}"
    echo "BUSTOOLS           : ${BUSTOOLS}"
    echo "TRANSCRIPTOME_FASTA: ${TRANSCRIPTOME_FASTA}"
    echo "INDEX_FILE         : ${INDEX_FILE}"
    echo "T2G_FILE           : ${T2G_FILE}"
    echo "BUILD_INDEX        : ${BUILD_INDEX}"
    echo "KMER_SIZE          : ${KMER_SIZE}"
    echo "READ_TYPE          : ${READ_TYPE}"
    echo "THREADS            : ${THREADS}"
    echo "THRESHOLD          : ${THRESHOLD}"
    echo ""
    echo "--- Samples ---"
    for entry in "${SAMPLES[@]}"; do
        echo "  ${entry}"
    done
    echo ""
    echo "--- Environment ---"
    echo "HOSTNAME      : $(hostname)"
    echo "USER          : $(whoami)"
    echo "CONDA_ENV     : ${CONDA_DEFAULT_ENV:-not set}"
    echo "kallisto ver  : $(${KALLISTO} version 2>&1 || echo 'version check failed')"
    echo "bustools ver  : $(${BUSTOOLS} version 2>&1 || echo 'version check failed')"
    echo "======================================"
} > "${EXPERIMENT_LOG}"

echo "Run directory   : ${RUN_DIR}"
echo "Experiment log  : ${EXPERIMENT_LOG}"
echo "Running log     : ${LOG}"

# =============================================================================
# STEP 0 (optional): BUILD INDEX
# Runs only when BUILD_INDEX=1. Skipped when BUILD_INDEX=0.
# =============================================================================

if [ "${BUILD_INDEX}" = "1" ]; then
    echo ""
    echo "=== Step 0: kallisto index ===" | tee -a "${LOG}"
    echo "  Transcriptome: ${TRANSCRIPTOME_FASTA}"  | tee -a "${LOG}"
    echo "  Output index : ${INDEX_FILE}"            | tee -a "${LOG}"
    echo "  k-mer size   : ${KMER_SIZE}"             | tee -a "${LOG}"

    if [ ! -f "${TRANSCRIPTOME_FASTA}" ]; then
        echo "ERROR: TRANSCRIPTOME_FASTA not found: ${TRANSCRIPTOME_FASTA}" | tee -a "${LOG}"
        exit 1
    fi

    "${KALLISTO}" index \
        --index "${INDEX_FILE}" \
        -k "${KMER_SIZE}" \
        "${TRANSCRIPTOME_FASTA}" \
        >> "${LOG}" 2>&1

    echo "  [OK] Index built: ${INDEX_FILE}" | tee -a "${LOG}"
else
    echo ""
    echo "=== Step 0: skipping index build (BUILD_INDEX=0) — using existing ${INDEX_FILE} ===" | tee -a "${LOG}"
    if [ ! -f "${INDEX_FILE}" ]; then
        echo "ERROR: INDEX_FILE not found: ${INDEX_FILE}. Set BUILD_INDEX=1 to build it." | tee -a "${LOG}"
        exit 1
    fi
fi

# =============================================================================
# READ_TYPE FLAGS — set once, applied to every sample
# =============================================================================

if [ "${READ_TYPE}" = "long" ]; then
    # Long-read mode: --long --unmapped for kallisto bus;
    # -P PacBio for quant-tcc (PacBio error model)
    BUS_MODE_FLAGS="--long --unmapped"
    QUANT_MODE_FLAGS="--long -P PacBio"
else
    # Short-read / paired-end mode
    BUS_MODE_FLAGS="--paired"
    QUANT_MODE_FLAGS=""
fi

# =============================================================================
# HELPERS
# =============================================================================

# Strip known read-file extensions to get a clean sample stem for cache naming
get_stem() {
    local fname="$1"
    local stem="${fname}"
    for ext in .fastq.gz .fasta.gz .fastq .fasta .fa .fq; do
        if [[ "${stem}" == *"${ext}" ]]; then stem="${stem%${ext}}"; break; fi
    done
    echo "${stem}"
}

# Returns 0 (cache hit) if count.mtx + count.ec.txt + transcripts.txt all exist
cache_complete() {
    local cache_dir="$1"
    [[ -f "${cache_dir}/count.mtx"      ]] && \
    [[ -f "${cache_dir}/count.ec.txt"   ]] && \
    [[ -f "${cache_dir}/transcripts.txt" ]] && \
    return 0 || return 1
}

# =============================================================================
# PIPELINE LOOP — one iteration per sample
# =============================================================================

for entry in "${SAMPLES[@]}"; do
    # Parse fields from the space-delimited entry string
    read -r SAMPLE_NAME READS_DIR READS_FILE1 READS_FILE2_OPT <<< "${entry}" 2>/dev/null || true

    echo ""
    echo "======================================================" | tee -a "${LOG}"
    echo "Processing sample: ${SAMPLE_NAME}"                      | tee -a "${LOG}"
    echo "======================================================" | tee -a "${LOG}"

    # Per-sample result subfolder inside the timestamped run dir
    SAMPLE_OUT="${RUN_DIR}/${SAMPLE_NAME}"
    mkdir -p "${SAMPLE_OUT}"

    # Cache dir — persists across runs next to the reads file (same convention as JK)
    STEM=$(get_stem "${READS_FILE1}")
    CACHE_DIR="${READS_DIR}/kallisto_output/${STEM}"

    # Build the reads argument (single file for long-read; two files for paired short-read)
    READS_ARG="${READS_DIR}${READS_FILE1}"
    if [ -n "${READS_FILE2_OPT:-}" ]; then
        READS_ARG="${READS_ARG} ${READS_DIR}${READS_FILE2_OPT}"
    fi

    echo "  Cache dir   : ${CACHE_DIR}"  | tee -a "${LOG}"
    echo "  Result dir  : ${SAMPLE_OUT}" | tee -a "${LOG}"

    # -----------------------------------------------------------------
    # Steps 1-3: kallisto bus + bustools sort + bustools count
    # Skipped on cache hit — saves time when re-running with same reads
    # -----------------------------------------------------------------
    if cache_complete "${CACHE_DIR}"; then
        echo "  [CACHE HIT] Skipping bus/sort/count — using: ${CACHE_DIR}" | tee -a "${LOG}"
    else
        echo "  [CACHE MISS] Running bus/sort/count → ${CACHE_DIR}" | tee -a "${LOG}"
        mkdir -p "${CACHE_DIR}"

        # Step 1: kallisto bus
        echo ""                                                            >> "${LOG}"
        echo "=== [${SAMPLE_NAME}] Step 1: kallisto bus ==="              >> "${LOG}"
        "${KALLISTO}" bus \
            -x bulk \
            --threshold "${THRESHOLD}" \
            -t "${THREADS}" \
            ${BUS_MODE_FLAGS} \
            -i "${INDEX_FILE}" \
            ${READS_ARG} \
            -o "${CACHE_DIR}/" \
            >> "${LOG}" 2>&1
        echo "  [OK] kallisto bus done"                                   | tee -a "${LOG}"

        # Step 2: bustools sort
        echo ""                                                            >> "${LOG}"
        echo "=== [${SAMPLE_NAME}] Step 2: bustools sort ==="             >> "${LOG}"
        "${BUSTOOLS}" sort \
            -t "${THREADS}" \
            "${CACHE_DIR}/output.bus" \
            -o "${CACHE_DIR}/sorted.bus" \
            >> "${LOG}" 2>&1
        echo "  [OK] bustools sort done"                                  | tee -a "${LOG}"

        # Step 3: bustools count
        echo ""                                                            >> "${LOG}"
        echo "=== [${SAMPLE_NAME}] Step 3: bustools count ==="            >> "${LOG}"
        "${BUSTOOLS}" count \
            "${CACHE_DIR}/sorted.bus" \
            -t "${CACHE_DIR}/transcripts.txt" \
            -e "${CACHE_DIR}/matrix.ec" \
            -o "${CACHE_DIR}/count" \
            --cm -m \
            -g "${T2G_FILE}" \
            >> "${LOG}" 2>&1
        echo "  [OK] bustools count done. TCC files cached: ${CACHE_DIR}" | tee -a "${LOG}"
    fi

    # Copy TCC files to experiment result dir for reproducibility
    for f in count.mtx count.ec.txt transcripts.txt matrix.ec run_info.json; do
        [ -f "${CACHE_DIR}/${f}" ] && cp "${CACHE_DIR}/${f}" "${SAMPLE_OUT}/${f}"
    done

    # -----------------------------------------------------------------
    # Step 4: kallisto quant-tcc — EM quantification from TCC matrix
    # -----------------------------------------------------------------
    echo ""                                                        >> "${LOG}"
    echo "=== [${SAMPLE_NAME}] Step 4: kallisto quant-tcc ==="    >> "${LOG}"

    "${KALLISTO}" quant-tcc \
        -t "${THREADS}" \
        ${QUANT_MODE_FLAGS} \
        "${CACHE_DIR}/count.mtx" \
        -i "${INDEX_FILE}" \
        -e "${CACHE_DIR}/count.ec.txt" \
        -o "${SAMPLE_OUT}/" \
        >> "${LOG}" 2>&1

    echo "  [OK] kallisto quant-tcc done for ${SAMPLE_NAME}"     | tee -a "${LOG}"

    # -----------------------------------------------------------------
    # Step 5: mtx_to_tsv.py — convert matrix.abundance.tpm.mtx → abundance.tsv
    # Produces a human-readable TSV (transcript_id, tpm) for downstream analysis
    # -----------------------------------------------------------------
    echo ""                                                        >> "${LOG}"
    echo "=== [${SAMPLE_NAME}] mtx → abundance.tsv ==="           >> "${LOG}"

    python "$(dirname "$(realpath "$0")")/mtx_to_tsv.py" \
        "${SAMPLE_OUT}" \
        >> "${LOG}" 2>&1

    echo "  [OK] abundance.tsv written for ${SAMPLE_NAME}"       | tee -a "${LOG}"

    # Report non-zero transcript counts from abundance.tsv (tpm column 2)
    if [ -f "${SAMPLE_OUT}/abundance.tsv" ]; then
        TOTAL_TX=$(awk 'NR>1' "${SAMPLE_OUT}/abundance.tsv" | wc -l)
        NONZERO_TX=$(awk 'NR>1 && $2>0' "${SAMPLE_OUT}/abundance.tsv" | wc -l)
        echo "  Transcript counts for ${SAMPLE_NAME}:"          | tee -a "${LOG}"
        echo "    Total transcripts : ${TOTAL_TX}"               | tee -a "${LOG}"
        echo "    Non-zero (est_counts > 0) : ${NONZERO_TX}"     | tee -a "${LOG}"
    fi

    echo ""                                                        >> "${LOG}"
    echo "  Output written to: ${SAMPLE_OUT}/"                   | tee -a "${LOG}"

done  # end sample loop

# =============================================================================
# TEARDOWN — runtime, code snapshot
# =============================================================================

# Save total wall-clock runtime
RUN_END=$(date +%s)
ELAPSED=$(( RUN_END - RUN_START ))
echo "Total elapsed time: ${ELAPSED} seconds" > "${RUN_DIR}/runtime.txt"
echo ""
echo "Total runtime: ${ELAPSED} seconds"

# Save code snapshot — copy this script into the results folder
mkdir -p "${RUN_DIR}/code_snapshot"
cp "$(realpath "$0")" "${RUN_DIR}/code_snapshot/run_lr_kallisto.sh"
echo "Code snapshot saved to: ${RUN_DIR}/code_snapshot/"

echo ""
echo "===== Pipeline complete ====="
echo "Results: ${RUN_DIR}"
