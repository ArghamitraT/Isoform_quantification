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
# k-mer sizes: long-read uses larger k (63); short-read uses standard k (31).
# Both are used automatically based on per-sample read_type (see SAMPLES below).
KMER_SIZE_LONG=63    # for long reads (PacBio/ONT)
KMER_SIZE_SHORT=31   # for short reads (paired Illumina)

# --- Run settings ---
PLATFORM="PacBio"    # long-read platform: "PacBio" or "ONT" (applies to all long-read samples)
THREADS=32
THRESHOLD=0.8        # kallisto bus alignment threshold (long-read only)

# --- Python / conda ---
CONDA_ENV="Joli_kallisto"
# Resolved below; override here if needed
PYTHON=""            # leave empty to auto-resolve from conda env

# --- Output base for experiment results ---
OUTPUT_BASE="/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

# --- JOLI EM settings ---
MAX_EM_ROUNDS=10000
MIN_ROUNDS=50
EFF_LEN_MODE="uniform"     # "uniform" (Phase 1) | "kallisto" (Phase 2+)
                           # NOTE: kallisto quant-tcc --long uses uniform weights internally
                           # (no fragment-length distribution available for long reads).
                           # JK must match this or correlation drops from ~0.99 to ~0.89.
EM_TYPE="plain"           # "plain" | "MAP" | "VI"
SAVE_SNAPSHOTS=true    # true  = save theta snapshots every SNAPSHOT_INTERVAL rounds
                        #         (needed for plot_convergence_animation.py)
SNAPSHOT_INTERVAL=5    # save a snapshot every N rounds (only used when SAVE_SNAPSHOTS=true)

# --- Experiment comment ---
# Free-text description saved to experiment_description.log for this run.
# Describe what you are testing, what changed, or what you expect.
EXPERIMENT_COMMENT="we are running JK with snapshot on;now records the iteration 0 before any calculation"

# --- Samples ---
# Each entry: "sample_name  read_type  reads_dir  file1  [file2]"
#   read_type = "long"  → single file, --long --threshold flags, KMER_SIZE_LONG
#   read_type = "short" → two files (R1 R2), standard paired-end, KMER_SIZE_SHORT
# Mixed long + short samples in the same array is supported.
SAMPLES=(
    # Long-read examples:
    # "toy       long  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long/downsampled  toy.fastq"
    # "flnc_01   long  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long/  flnc_01.fastq"
    # "flnc_02   long  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long/  flnc_02.fastq"
    # Short-read examples:
    # "sim_sr_s1  short  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/  ds_100_num1_aln_01_short_1.fq  ds_100_num1_aln_01_short_2.fq"
    # "sim_sr_s2  short  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/  ds_100_num1_aln_21_short_1.fq  ds_100_num1_aln_21_short_2.fq"
    # "sr_s1  short  /path/to/short_reads/  sample1_R1.fastq.gz  sample1_R2.fastq.gz"
    # "sr_s2  short  /path/to/short_reads/  sample2_R1.fastq.gz  sample2_R2.fastq.gz"
    "sim1  long  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/  ds_100_num1_aln_01_long.fasta"
    "sim2  long  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/  ds_100_num1_aln_21_long.fasta"
)

# ============================================================
# END CONFIG
# ============================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Resolve Python interpreter from conda env if not set
if [[ -z "$PYTHON" ]]; then
    PYTHON="$(conda run -n "${CONDA_ENV}" which python 2>/dev/null || true)"
    if [[ -z "$PYTHON" ]]; then
        echo "ERROR: Could not find python in conda env '${CONDA_ENV}'. Set PYTHON in CONFIG."
        exit 1
    fi
fi

# Fix GLIBCXX version mismatch: use the *conda env's* libstdc++, not the base env.
# Derive the env prefix from the resolved PYTHON path (e.g. .../envs/Joli_kallisto/bin/python
# → .../envs/Joli_kallisto/lib).
CONDA_ENV_LIB="$(dirname "$(dirname "${PYTHON}")")/lib"
if [[ -d "${CONDA_ENV_LIB}" ]]; then
    export LD_LIBRARY_PATH="${CONDA_ENV_LIB}:${LD_LIBRARY_PATH:-}"
fi

# ---- Create timestamped result directory ----
TIMESTAMP="exprmnt_$(date +%Y_%m_%d__%H_%M_%S)"
RUN_DIR="${OUTPUT_BASE}/${TIMESTAMP}"
mkdir -p "${RUN_DIR}"
LOG="${RUN_DIR}/running.log"

# First line of every running.log: which script produced this log
echo "Script: $(realpath "$0")" > "${LOG}"

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

=== EXPERIMENT COMMENT ===
${EXPERIMENT_COMMENT:-"(none)"}

=== CONFIG ===
KALLISTO:           ${KALLISTO}
BUSTOOLS:           ${BUSTOOLS}
TRANSCRIPTOME_FASTA:${TRANSCRIPTOME_FASTA}
INDEX_FILE:         ${INDEX_FILE}
T2G_FILE:           ${T2G_FILE}
BUILD_INDEX:        ${BUILD_INDEX}
KMER_SIZE_LONG:     ${KMER_SIZE_LONG}
KMER_SIZE_SHORT:    ${KMER_SIZE_SHORT}
PLATFORM:           ${PLATFORM}  (long-read only)
THREADS:            ${THREADS}
THRESHOLD:          ${THRESHOLD}  (long-read only)
MAX_EM_ROUNDS:      ${MAX_EM_ROUNDS}
MIN_ROUNDS:         ${MIN_ROUNDS}
EFF_LEN_MODE:       ${EFF_LEN_MODE}
EM_TYPE:            ${EM_TYPE}
SAVE_SNAPSHOTS:     ${SAVE_SNAPSHOTS}
SNAPSHOT_INTERVAL:  ${SNAPSHOT_INTERVAL}

=== SAMPLES ===
$(printf '%s\n' "${SAMPLES[@]}")
EOF

# ---- Code snapshot ----
SNAPSHOT_DIR="${RUN_DIR}/code_snapshot"
mkdir -p "${SNAPSHOT_DIR}"
for ext in py sh txt yml yaml; do
    find "${SCRIPT_DIR}/.." -name "*.${ext}" -exec cp {} "${SNAPSHOT_DIR}/" \;
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
    echo "Step 0: kallisto index (k=${KMER_SIZE_LONG}, long-read k-mer)" | tee -a "${LOG}"
    [[ ! -f "${TRANSCRIPTOME_FASTA}" ]] && echo "ERROR: TRANSCRIPTOME_FASTA not found: ${TRANSCRIPTOME_FASTA}" | tee -a "${LOG}" && exit 1
    "${KALLISTO}" index --index "${INDEX_FILE}" -k "${KMER_SIZE_LONG}" "${TRANSCRIPTOME_FASTA}" \
        2>&1 | tee -a "${LOG}"
    echo "[OK] Index built: ${INDEX_FILE}" | tee -a "${LOG}"
else
    echo "Step 0: skipping index build (BUILD_INDEX=0), using: ${INDEX_FILE}" | tee -a "${LOG}"
    [[ ! -f "${INDEX_FILE}" ]] && echo "ERROR: INDEX_FILE not found: ${INDEX_FILE}" | tee -a "${LOG}" && exit 1
fi

for SAMPLE_ENTRY in "${SAMPLES[@]}"; do
    # Parse fields: sample_name  read_type  reads_dir  file1  [file2]
    read -r SAMPLE_NAME SAMPLE_READ_TYPE READS_DIR READS_FILE1 READS_FILE2_MAYBE <<< "${SAMPLE_ENTRY}"
    READS_FILE2="${READS_FILE2_MAYBE:-}"

    # Validate read_type
    if [[ "${SAMPLE_READ_TYPE}" != "long" && "${SAMPLE_READ_TYPE}" != "short" ]]; then
        echo "ERROR: Unknown read_type '${SAMPLE_READ_TYPE}' for sample '${SAMPLE_NAME}'." | tee -a "${LOG}"
        echo "       Must be 'long' or 'short'." | tee -a "${LOG}"
        ERRORS+=("${SAMPLE_NAME}")
        continue
    fi

    # Short-read requires two files
    if [[ "${SAMPLE_READ_TYPE}" == "short" && -z "${READS_FILE2}" ]]; then
        echo "ERROR: Sample '${SAMPLE_NAME}' has read_type=short but only one file provided." | tee -a "${LOG}"
        echo "       Format: \"sample_name  short  reads_dir  R1_file  R2_file\"" | tee -a "${LOG}"
        ERRORS+=("${SAMPLE_NAME}")
        continue
    fi

    # Select k-mer size based on read type
    if [[ "${SAMPLE_READ_TYPE}" == "long" ]]; then
        KMER_SIZE="${KMER_SIZE_LONG}"
    else
        KMER_SIZE="${KMER_SIZE_SHORT}"
    fi

    STEM=$(get_stem "${READS_FILE1}")
    CACHE_DIR="${READS_DIR}/kallisto_output/${STEM}"
    SAMPLE_RESULT_DIR="${RUN_DIR}/${SAMPLE_NAME}"
    mkdir -p "${SAMPLE_RESULT_DIR}"

    echo "" | tee -a "${LOG}"
    echo "------------------------------------------------------------" | tee -a "${LOG}"
    echo "Sample:     ${SAMPLE_NAME}"                              | tee -a "${LOG}"
    echo "Read type:  ${SAMPLE_READ_TYPE}  (k=${KMER_SIZE})"      | tee -a "${LOG}"
    echo "Reads:      ${READS_DIR}/${READS_FILE1}"                 | tee -a "${LOG}"
    if [[ -n "${READS_FILE2}" ]]; then
        echo "            ${READS_DIR}/${READS_FILE2}"             | tee -a "${LOG}"
    fi
    echo "Cache dir:  ${CACHE_DIR}"                                | tee -a "${LOG}"
    echo "Result dir: ${SAMPLE_RESULT_DIR}"                        | tee -a "${LOG}"
    echo "------------------------------------------------------------" | tee -a "${LOG}"

    # ---- CACHE CHECK ----
    if cache_complete "${CACHE_DIR}"; then
        echo "[CACHE HIT] Skipping kallisto bus + bustools steps." | tee -a "${LOG}"
        echo "  Using existing TCC files in: ${CACHE_DIR}"        | tee -a "${LOG}"

    else
        echo "[CACHE MISS] Running kallisto bus + bustools steps." | tee -a "${LOG}"
        mkdir -p "${CACHE_DIR}"

        # -- Step 1: kallisto bus --
        echo "Step 1: kallisto bus (read_type=${SAMPLE_READ_TYPE})" | tee -a "${LOG}"
        if [[ "${SAMPLE_READ_TYPE}" == "long" ]]; then
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
    "${PYTHON}" "${SCRIPT_DIR}/../main_joli.py" \
        --sample_dir    "${CACHE_DIR}" \
        --output_dir    "${SAMPLE_RESULT_DIR}" \
        --eff_len_mode  "${EFF_LEN_MODE}" \
        --max_em_rounds "${MAX_EM_ROUNDS}" \
        --min_rounds    "${MIN_ROUNDS}" \
        --em_type            "${EM_TYPE}" \
        --save_snapshots     "${SAVE_SNAPSHOTS}" \
        --snapshot_interval  "${SNAPSHOT_INTERVAL}" \
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
