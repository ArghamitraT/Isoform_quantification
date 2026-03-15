#!/usr/bin/env bash
# =============================================================================
# run_joli_kallisto.sh
# ====================
# JOLI-Kallisto full pipeline: pseudoalignment + JOLI EM quantification.
#
# Steps per sample:
#   [cache check] Does <reads_dir>/kallisto_output/<sample_stem>/ already have
#                 count.mtx + matrix.ec + transcripts.txt?
#      YES -> skip to JOLI EM (saves time on re-runs with different EM settings)
#      NO  -> Step 1: kallisto bus  -> output.bus
#             Step 2: bustools sort -> sorted.bus
#             Step 3: bustools count -> count.mtx, matrix.ec, transcripts.txt
#                     (saved to <reads_dir>/kallisto_output/<sample_stem>/)
#   Step 4: python main_joli.py (JOLI EM) -> abundance.tsv
#
# Cache directory convention:
#   <reads_dir>/kallisto_output/<sample_stem>/
#   e.g. reads/ds_52.fastq -> reads/kallisto_output/ds_52/count.mtx
#
# Run locally:   bash run_joli_kallisto.sh
# Run on SLURM:  bash submit_joli_pipeline.sh
# =============================================================================

set -euo pipefail

# ============================================================
# CONFIG — edit everything here; do not touch pipeline logic below
# ============================================================

# --- Tool paths ---
KALLISTO="/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/kallisto/build/src/kallisto"
BUSTOOLS="bustools"   # full path if not on PATH

# --- Reference files ---
TRANSCRIPTOME_FASTA="/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/transcriptome.fasta"
INDEX_FILE="/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/new_index.idx"
T2G_FILE="/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/t2g.txt"

# --- Index settings ---
# Set to 1 to build/rebuild the index from TRANSCRIPTOME_FASTA before running samples.
# Set to 0 to skip and use the existing INDEX_FILE.
BUILD_INDEX=0
# k-mer size: 63 for long-read (PacBio/ONT), 31 for short-read
KMER_SIZE=63

# --- Run settings ---
READ_TYPE="long"     # "long" (PacBio/ONT) or "short" (paired Illumina)
PLATFORM="PacBio"    # used when READ_TYPE=long: "PacBio" or "ONT"
THREADS=32
THRESHOLD=0.8        # kallisto bus alignment threshold

# --- Python / conda ---
CONDA_ENV="Joli_kallisto"
# Resolved below; override here if needed
PYTHON=""            # leave empty to auto-resolve from conda env

# --- Output base for experiment results ---
OUTPUT_BASE="/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

# --- JOLI EM settings ---
MAX_EM_ROUNDS=10000
MIN_ROUNDS=50
EFF_LEN_MODE="uniform"   # "uniform" (Phase 1) | "kallisto" (Phase 2+)
EM_TYPE="plain"           # "plain" | "MAP" | "VI"

# --- Samples ---
# Format (long-read):  "sample_name  reads_dir  reads_file"
# Format (short-read): "sample_name  reads_dir  R1_file  R2_file"
SAMPLES=(
    "ds_52_furtherDownsampled  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long/downsampled  ds_52_furtherDownsampled.fastq"
    # "sim2  /path/to/sim/reads  PacBio.simulated.fasta"
)

# ============================================================
# END CONFIG
# ============================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Fix GLIBCXX version mismatch: prefer conda's bundled libstdc++ over system one
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
fi

# Resolve Python interpreter from conda env if not set
if [[ -z "$PYTHON" ]]; then
    PYTHON="$(conda run -n "${CONDA_ENV}" which python 2>/dev/null || true)"
    if [[ -z "$PYTHON" ]]; then
        echo "ERROR: Could not find python in conda env '${CONDA_ENV}'. Set PYTHON in CONFIG."
        exit 1
    fi
fi

# ---- Create timestamped result directory ----
TIMESTAMP="exprmnt_$(date +%Y_%m_%d__%H_%M_%S)"
RUN_DIR="${OUTPUT_BASE}/${TIMESTAMP}"
mkdir -p "${RUN_DIR}"
LOG="${RUN_DIR}/running.log"

echo "============================================================" | tee -a "${LOG}"
echo "JOLI-Kallisto pipeline: ${TIMESTAMP}"                          | tee -a "${LOG}"
echo "Results: ${RUN_DIR}"                                           | tee -a "${LOG}"
echo "Python:  ${PYTHON}"                                            | tee -a "${LOG}"
echo "============================================================" | tee -a "${LOG}"

RUN_START=$(date +%s)

# ---- Save experiment description ----
cat > "${RUN_DIR}/experiment_description.log" <<EOF
JOLI-Kallisto Experiment
Date: $(date -Iseconds)
Script: ${BASH_SOURCE[0]}

=== CONFIG ===
KALLISTO:           ${KALLISTO}
BUSTOOLS:           ${BUSTOOLS}
TRANSCRIPTOME_FASTA:${TRANSCRIPTOME_FASTA}
INDEX_FILE:         ${INDEX_FILE}
T2G_FILE:           ${T2G_FILE}
BUILD_INDEX:        ${BUILD_INDEX}
KMER_SIZE:          ${KMER_SIZE}
READ_TYPE:          ${READ_TYPE}
PLATFORM:           ${PLATFORM}
THREADS:            ${THREADS}
THRESHOLD:          ${THRESHOLD}
MAX_EM_ROUNDS:  ${MAX_EM_ROUNDS}
MIN_ROUNDS:     ${MIN_ROUNDS}
EFF_LEN_MODE:   ${EFF_LEN_MODE}
EM_TYPE:        ${EM_TYPE}

=== SAMPLES ===
$(printf '%s\n' "${SAMPLES[@]}")
EOF

# ---- Code snapshot ----
SNAPSHOT_DIR="${RUN_DIR}/code_snapshot"
mkdir -p "${SNAPSHOT_DIR}"
for ext in py sh txt yml yaml; do
    find "${SCRIPT_DIR}" -maxdepth 1 -name "*.${ext}" -exec cp {} "${SNAPSHOT_DIR}/" \;
done
echo "Code snapshot saved to: ${SNAPSHOT_DIR}" | tee -a "${LOG}"

# ============================================================
# Helper: derive sample stem from reads filename (strip extensions)
# ============================================================
get_stem() {
    local fname="$1"
    local stem="${fname}"
    for ext in .fastq.gz .fasta.gz .fastq .fasta .fa .fq; do
        if [[ "${stem}" == *"${ext}" ]]; then
            stem="${stem%${ext}}"
            break
        fi
    done
    echo "${stem}"
}

# ============================================================
# Helper: check if bustools cache is complete
# Returns 0 (true) if cache hit, 1 (false) if cache miss
# ============================================================
cache_complete() {
    local cache_dir="$1"
    [[ -f "${cache_dir}/count.mtx"       ]] && \
    [[ -f "${cache_dir}/matrix.ec"       ]] && \
    [[ -f "${cache_dir}/transcripts.txt" ]] && \
    return 0 || return 1
}

# ============================================================
# Main per-sample processing
# ============================================================
ERRORS=()

# ============================================================
# Step 0 (optional): Build kallisto index
# Set BUILD_INDEX=1 in CONFIG to rebuild; 0 to skip.
# ============================================================
if [[ "${BUILD_INDEX}" == "1" ]]; then
    echo "Step 0: kallisto index" | tee -a "${LOG}"
    [[ ! -f "${TRANSCRIPTOME_FASTA}" ]] && echo "ERROR: TRANSCRIPTOME_FASTA not found: ${TRANSCRIPTOME_FASTA}" | tee -a "${LOG}" && exit 1
    "${KALLISTO}" index --index "${INDEX_FILE}" -k "${KMER_SIZE}" "${TRANSCRIPTOME_FASTA}" \
        2>&1 | tee -a "${LOG}"
    echo "[OK] Index built: ${INDEX_FILE}" | tee -a "${LOG}"
else
    echo "Step 0: skipping index build (BUILD_INDEX=0), using: ${INDEX_FILE}" | tee -a "${LOG}"
    [[ ! -f "${INDEX_FILE}" ]] && echo "ERROR: INDEX_FILE not found: ${INDEX_FILE}" | tee -a "${LOG}" && exit 1
fi

for SAMPLE_ENTRY in "${SAMPLES[@]}"; do
    # Parse sample entry fields
    read -r SAMPLE_NAME READS_DIR READS_FILE1 <<< "${SAMPLE_ENTRY}"
    READS_FILE2=""
    # If 4th field present: paired short-read R2
    if [[ $(echo "${SAMPLE_ENTRY}" | awk '{print NF}') -ge 4 ]]; then
        READS_FILE2=$(echo "${SAMPLE_ENTRY}" | awk '{print $4}')
    fi

    STEM=$(get_stem "${READS_FILE1}")
    CACHE_DIR="${READS_DIR}/kallisto_output/${STEM}"
    SAMPLE_RESULT_DIR="${RUN_DIR}/${SAMPLE_NAME}"
    mkdir -p "${SAMPLE_RESULT_DIR}"

    echo "" | tee -a "${LOG}"
    echo "------------------------------------------------------------" | tee -a "${LOG}"
    echo "Sample:     ${SAMPLE_NAME}"          | tee -a "${LOG}"
    echo "Reads:      ${READS_DIR}/${READS_FILE1}" | tee -a "${LOG}"
    echo "Cache dir:  ${CACHE_DIR}"            | tee -a "${LOG}"
    echo "Result dir: ${SAMPLE_RESULT_DIR}"    | tee -a "${LOG}"
    echo "------------------------------------------------------------" | tee -a "${LOG}"

    # ---- CACHE CHECK ----
    # This is the single decision point: bustools or skip?
    if cache_complete "${CACHE_DIR}"; then
        echo "[CACHE HIT] Skipping kallisto bus + bustools steps." | tee -a "${LOG}"
        echo "  Using existing TCC files in: ${CACHE_DIR}"        | tee -a "${LOG}"

    else
        echo "[CACHE MISS] Running kallisto bus + bustools steps." | tee -a "${LOG}"
        mkdir -p "${CACHE_DIR}"

        # -- Step 1: kallisto bus --
        echo "Step 1: kallisto bus" | tee -a "${LOG}"
        if [[ "${READ_TYPE}" == "long" ]]; then
            "${KALLISTO}" bus \
                -x bulk \
                -i "${INDEX_FILE}" \
                -o "${CACHE_DIR}" \
                -t "${THREADS}" \
                --long \
                --threshold "${THRESHOLD}" \
                --unmapped \
                "${READS_DIR}/${READS_FILE1}" \
                2>&1 | tee -a "${LOG}"
        else
            # Short-read paired
            "${KALLISTO}" bus \
                -x bulk \
                -i "${INDEX_FILE}" \
                -o "${CACHE_DIR}" \
                -t "${THREADS}" \
                "${READS_DIR}/${READS_FILE1}" \
                "${READS_DIR}/${READS_FILE2}" \
                2>&1 | tee -a "${LOG}"
        fi

        # Verify output.bus was created
        if [[ ! -f "${CACHE_DIR}/output.bus" ]]; then
            echo "ERROR: kallisto bus did not produce output.bus" | tee -a "${LOG}"
            ERRORS+=("${SAMPLE_NAME}")
            continue
        fi

        # -- Step 2: bustools sort --
        echo "Step 2: bustools sort" | tee -a "${LOG}"
        "${BUSTOOLS}" sort \
            -t "${THREADS}" \
            -o "${CACHE_DIR}/sorted.bus" \
            "${CACHE_DIR}/output.bus" \
            2>&1 | tee -a "${LOG}"

        # -- Step 3: bustools count --
        echo "Step 3: bustools count" | tee -a "${LOG}"
        "${BUSTOOLS}" count \
            "${CACHE_DIR}/sorted.bus" \
            -t "${CACHE_DIR}/transcripts.txt" \
            -e "${CACHE_DIR}/matrix.ec" \
            -o "${CACHE_DIR}/count" \
            --cm -m \
            -g "${T2G_FILE}" \
            2>&1 | tee -a "${LOG}"

        # Verify required TCC files exist
        if ! cache_complete "${CACHE_DIR}"; then
            echo "ERROR: bustools count did not produce required files in ${CACHE_DIR}" \
                | tee -a "${LOG}"
            ERRORS+=("${SAMPLE_NAME}")
            continue
        fi

        echo "Bustools pipeline complete. TCC files cached at: ${CACHE_DIR}" \
            | tee -a "${LOG}"
    fi

    # ---- Copy TCC intermediates to experiment result dir ----
    # Always copy (whether cache hit or miss) so every experiment folder has
    # a complete snapshot of the bustools output used for that run.
    echo "Copying TCC intermediate files to: ${SAMPLE_RESULT_DIR}" | tee -a "${LOG}"
    for f in output.bus sorted.bus matrix.ec transcripts.txt \
              count.mtx count.ec.txt count.barcodes.txt run_info.json flens.txt; do
        if [[ -f "${CACHE_DIR}/${f}" ]]; then
            cp "${CACHE_DIR}/${f}" "${SAMPLE_RESULT_DIR}/${f}"
        fi
    done

    # ---- Step 4: JOLI EM ----
    echo "Step 4: JOLI EM (main_joli.py)" | tee -a "${LOG}"
    "${PYTHON}" "${SCRIPT_DIR}/main_joli.py" \
        --sample_dir    "${CACHE_DIR}" \
        --output_dir    "${SAMPLE_RESULT_DIR}" \
        --eff_len_mode  "${EFF_LEN_MODE}" \
        --max_em_rounds "${MAX_EM_ROUNDS}" \
        --min_rounds    "${MIN_ROUNDS}" \
        --em_type       "${EM_TYPE}" \
        2>&1 | tee -a "${LOG}"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: JOLI EM failed for ${SAMPLE_NAME}" | tee -a "${LOG}"
        ERRORS+=("${SAMPLE_NAME}")
        continue
    fi

    echo "Sample ${SAMPLE_NAME} complete." | tee -a "${LOG}"
done

# ---- Save total runtime ----
RUN_END=$(date +%s)
ELAPSED=$(( RUN_END - RUN_START ))
echo "${ELAPSED} seconds" > "${RUN_DIR}/runtime.txt"

echo "" | tee -a "${LOG}"
echo "============================================================" | tee -a "${LOG}"
echo "Pipeline complete in ${ELAPSED}s" | tee -a "${LOG}"
echo "Results: ${RUN_DIR}"              | tee -a "${LOG}"

if [[ ${#ERRORS[@]} -gt 0 ]]; then
    echo "FAILED samples: ${ERRORS[*]}" | tee -a "${LOG}"
    exit 1
fi
echo "All samples succeeded."           | tee -a "${LOG}"
echo "============================================================" | tee -a "${LOG}"
