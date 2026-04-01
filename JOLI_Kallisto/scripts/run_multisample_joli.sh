#!/usr/bin/env bash
# =============================================================================
# run_multisample_joli.sh
# =======================
# JOLI-Kallisto multi-sample pipeline: pseudoalignment + joint MAP EM.
#
# Steps per sample:
#   [cache check] Does <reads_dir>/kallisto_output/<sample_stem>/ already have
#                 count.mtx + matrix.ec + transcripts.txt?
#      YES -> skip bustools (cache hit)
#      NO  -> Step 1: kallisto bus  -> output.bus
#             Step 2: bustools sort -> sorted.bus
#             Step 3: bustools count -> count.mtx, matrix.ec, transcripts.txt
#                     (cached at <reads_dir>/kallisto_output/<sample_stem>/)
#
# After all samples are processed:
#   Step 4: python main_multisample_joli.py --sample_dirs <cache_dir_1> ...
#           Runs joint MAP EM with shared Dirichlet prior across all samples.
#           Creates its own timestamped output folder under OUTPUT_BASE.
#
# Cache directory convention:
#   <reads_dir>/kallisto_output/<sample_stem>/
#   e.g. reads/flnc_01.fastq -> reads/kallisto_output/flnc_01/count.mtx
#
# Supported read types per sample (set in SAMPLES array, field 2):
#   long  -> single file, --long --threshold flags, k=KMER_SIZE_LONG (63)
#   short -> two files (R1 R2), standard paired-end, k=KMER_SIZE_SHORT (31)
# Mixed long + short samples in one run are fully supported.
#
# Run locally:   bash run_multisample_joli.sh
# Run on SLURM:  bash submit_multisample_joli.sh
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
# Set to 1 to build/rebuild the index from TRANSCRIPTOME_FASTA before running.
# Set to 0 to skip and use the existing INDEX_FILE.
BUILD_INDEX=0
# k-mer sizes: long-read uses larger k (63); short-read uses standard k (31).
# Both are used automatically based on per-sample read_type (see SAMPLES below).
KMER_SIZE_LONG=63    # for long reads (PacBio/ONT)
KMER_SIZE_SHORT=31   # for short reads (paired Illumina)

# --- Run settings ---
PLATFORM="PacBio"    # long-read platform: "PacBio" or "ONT" (used for all long-read samples)
THREADS=32
THRESHOLD=0.8        # kallisto bus alignment threshold (long-read only)

# --- Python / conda ---
CONDA_ENV="Joli_kallisto"
PYTHON=""            # leave empty to auto-resolve from conda env

# --- Output base for experiment results ---
OUTPUT_BASE="/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

# --- JOLI MAP EM settings (passed to main_multisample_joli.py) ---
EFF_LEN_MODE="uniform"     # "uniform" | "kallisto"
CONVERGENCE_MODE="joli"    # "joli" (recommended for MAP) | "kallisto"
MAX_EM_ROUNDS=10000
MIN_EM_ROUNDS=1
MAX_GD_ROUNDS=500
GD_LR=0.01
ALPHA_INITIAL=1.0
GD_CONVERGENCE_TOL=1e-6
GD_STEPS_PER_ROUND=10
MIN_READ_SUPPORT=0     # Fix A: min expected reads n[t] to apply Dirichlet prior
#                        #   0.0 = disabled (original behaviour)
#                        #   try 0.1 to reduce false-positive multi-mapping leakage
SAVE_SNAPSHOTS=true    # true  = save alpha+theta snapshots every SNAPSHOT_INTERVAL rounds
#                        #          to snapshots.pkl — used by plot_convergence_animation.py
#                        # false = disabled (default; saves disk space)
SNAPSHOT_INTERVAL=5    # save a snapshot every N rounds (only used when SAVE_SNAPSHOTS=true)
LOOP_MODE="em_wrapper"  # Training loop structure:
#                        #   "gd_wrapper" — GD outer loop, EM to convergence per round (default)
#                        #   "em_wrapper" — EM convergence drives outer loop, 1 EM step +
#                        #                 GD_STEPS_PER_ROUND Adam steps per iteration (AT_code)

# Free-text description of this experiment run — printed at startup and saved in
# experiment_description.log. Leave empty for no comment.
EXPERIMENT_COMMENT="em_wrapper; with snapshot ON; now records the iteration 0 before any calculation"

# --- Samples ---
# Each entry: "sample_name  read_type  reads_dir  file1  [file2]"
#   read_type = "long"  → single file, kallisto --long --threshold flags, KMER_SIZE_LONG
#   read_type = "short" → two files (R1 R2), standard paired-end flags, KMER_SIZE_SHORT
# Mixed long + short samples in the same array is supported.
# Minimum 2 samples required for multi-sample MAP EM.
SAMPLES=(
    # Long-read examples:
    # "flnc_01  long  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long/  flnc_01.fastq"
    # "flnc_02  long  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long/  flnc_02.fastq"
    # Short-read examples:
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

# Fix GLIBCXX version mismatch: use the conda env's libstdc++
CONDA_ENV_LIB="$(dirname "$(dirname "${PYTHON}")")/lib"
if [[ -d "${CONDA_ENV_LIB}" ]]; then
    export LD_LIBRARY_PATH="${CONDA_ENV_LIB}:${LD_LIBRARY_PATH:-}"
fi

# ---- Create a preprocessing log (bustools steps only) ----
# The final experiment folder is created by main_multisample_joli.py.
PREP_TIMESTAMP="$(date +%Y_%m_%d__%H_%M_%S)"
PREP_LOG="/tmp/joli_multisample_prep_${PREP_TIMESTAMP}.log"
echo "Script: $(realpath "$0")" > "${PREP_LOG}"
echo "Preprocessing started: ${PREP_TIMESTAMP}" | tee -a "${PREP_LOG}"
if [[ -n "${EXPERIMENT_COMMENT}" ]]; then
    echo "EXPERIMENT_COMMENT: ${EXPERIMENT_COMMENT}" | tee -a "${PREP_LOG}"
fi

RUN_START=$(date +%s)

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
# Validate: need at least 2 samples
# ============================================================
if [[ ${#SAMPLES[@]} -lt 2 ]]; then
    echo "ERROR: Multi-sample MAP EM requires at least 2 samples." | tee -a "${PREP_LOG}"
    echo "       Found ${#SAMPLES[@]} sample(s) in CONFIG SAMPLES array." | tee -a "${PREP_LOG}"
    exit 1
fi

# ============================================================
# Step 0 (optional): Build kallisto index
# ============================================================
if [[ "${BUILD_INDEX}" == "1" ]]; then
    echo "Step 0: kallisto index (k=${KMER_SIZE_LONG}, long-read k-mer)" | tee -a "${PREP_LOG}"
    [[ ! -f "${TRANSCRIPTOME_FASTA}" ]] && \
        echo "ERROR: TRANSCRIPTOME_FASTA not found: ${TRANSCRIPTOME_FASTA}" | tee -a "${PREP_LOG}" && exit 1
    "${KALLISTO}" index --index "${INDEX_FILE}" -k "${KMER_SIZE_LONG}" "${TRANSCRIPTOME_FASTA}" \
        2>&1 | tee -a "${PREP_LOG}"
    echo "[OK] Index built: ${INDEX_FILE}" | tee -a "${PREP_LOG}"
else
    echo "Step 0: skipping index build (BUILD_INDEX=0), using: ${INDEX_FILE}" | tee -a "${PREP_LOG}"
    [[ ! -f "${INDEX_FILE}" ]] && \
        echo "ERROR: INDEX_FILE not found: ${INDEX_FILE}" | tee -a "${PREP_LOG}" && exit 1
fi

# ============================================================
# Per-sample bustools preprocessing
# Collect cache dirs for passing to main_multisample_joli.py
# ============================================================
ERRORS=()
CACHE_DIRS=()    # will be passed as --sample_dirs to the Python script
SAMPLE_NAMES=()  # will be passed as --sample_names to the Python script

for SAMPLE_ENTRY in "${SAMPLES[@]}"; do
    # Parse fields: sample_name  read_type  reads_dir  file1  [file2]
    read -r SAMPLE_NAME SAMPLE_READ_TYPE READS_DIR READS_FILE1 READS_FILE2_MAYBE <<< "${SAMPLE_ENTRY}"
    READS_FILE2="${READS_FILE2_MAYBE:-}"

    # Validate read_type
    if [[ "${SAMPLE_READ_TYPE}" != "long" && "${SAMPLE_READ_TYPE}" != "short" ]]; then
        echo "ERROR: Unknown read_type '${SAMPLE_READ_TYPE}' for sample '${SAMPLE_NAME}'." \
            | tee -a "${PREP_LOG}"
        echo "       Must be 'long' or 'short'." | tee -a "${PREP_LOG}"
        ERRORS+=("${SAMPLE_NAME}")
        continue
    fi

    # Short-read requires two files
    if [[ "${SAMPLE_READ_TYPE}" == "short" && -z "${READS_FILE2}" ]]; then
        echo "ERROR: Sample '${SAMPLE_NAME}' has read_type=short but only one file provided." \
            | tee -a "${PREP_LOG}"
        echo "       Format: \"sample_name  short  reads_dir  R1_file  R2_file\"" \
            | tee -a "${PREP_LOG}"
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

    echo "" | tee -a "${PREP_LOG}"
    echo "------------------------------------------------------------" | tee -a "${PREP_LOG}"
    echo "Sample:    ${SAMPLE_NAME}"                          | tee -a "${PREP_LOG}"
    echo "Read type: ${SAMPLE_READ_TYPE}  (k=${KMER_SIZE})"  | tee -a "${PREP_LOG}"
    echo "Reads:     ${READS_DIR}/${READS_FILE1}"             | tee -a "${PREP_LOG}"
    if [[ -n "${READS_FILE2}" ]]; then
        echo "           ${READS_DIR}/${READS_FILE2}"         | tee -a "${PREP_LOG}"
    fi
    echo "Cache dir: ${CACHE_DIR}"                            | tee -a "${PREP_LOG}"
    echo "------------------------------------------------------------" | tee -a "${PREP_LOG}"

    if cache_complete "${CACHE_DIR}"; then
        echo "[CACHE HIT] Skipping kallisto bus + bustools steps." | tee -a "${PREP_LOG}"
        echo "  Using existing TCC files in: ${CACHE_DIR}"        | tee -a "${PREP_LOG}"
    else
        echo "[CACHE MISS] Running kallisto bus + bustools steps." | tee -a "${PREP_LOG}"
        mkdir -p "${CACHE_DIR}"

        # -- Step 1: kallisto bus --
        echo "Step 1: kallisto bus (read_type=${SAMPLE_READ_TYPE})" | tee -a "${PREP_LOG}"
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
                2>&1 | tee -a "${PREP_LOG}"
        else
            "${KALLISTO}" bus \
                -x bulk \
                -i "${INDEX_FILE}" \
                -o "${CACHE_DIR}" \
                -t "${THREADS}" \
                "${READS_DIR}/${READS_FILE1}" \
                "${READS_DIR}/${READS_FILE2}" \
                2>&1 | tee -a "${PREP_LOG}"
        fi

        if [[ ! -f "${CACHE_DIR}/output.bus" ]]; then
            echo "ERROR: kallisto bus did not produce output.bus for ${SAMPLE_NAME}" \
                | tee -a "${PREP_LOG}"
            ERRORS+=("${SAMPLE_NAME}")
            continue
        fi

        # -- Step 2: bustools sort --
        echo "Step 2: bustools sort" | tee -a "${PREP_LOG}"
        "${BUSTOOLS}" sort \
            -t "${THREADS}" \
            -o "${CACHE_DIR}/sorted.bus" \
            "${CACHE_DIR}/output.bus" \
            2>&1 | tee -a "${PREP_LOG}"

        # -- Step 3: bustools count --
        echo "Step 3: bustools count" | tee -a "${PREP_LOG}"
        "${BUSTOOLS}" count \
            "${CACHE_DIR}/sorted.bus" \
            -t "${CACHE_DIR}/transcripts.txt" \
            -e "${CACHE_DIR}/matrix.ec" \
            -o "${CACHE_DIR}/count" \
            --cm -m \
            -g "${T2G_FILE}" \
            2>&1 | tee -a "${PREP_LOG}"

        if ! cache_complete "${CACHE_DIR}"; then
            echo "ERROR: bustools count did not produce required files for ${SAMPLE_NAME}" \
                | tee -a "${PREP_LOG}"
            ERRORS+=("${SAMPLE_NAME}")
            continue
        fi

        echo "Bustools pipeline complete. TCC files cached at: ${CACHE_DIR}" \
            | tee -a "${PREP_LOG}"
    fi

    CACHE_DIRS+=("${CACHE_DIR}")
    SAMPLE_NAMES+=("${SAMPLE_NAME}")
done

# ---- Fail fast if any bustools step failed ----
if [[ ${#ERRORS[@]} -gt 0 ]]; then
    echo "ERROR: Bustools failed for: ${ERRORS[*]}" | tee -a "${PREP_LOG}"
    echo "Aborting before MAP EM step."             | tee -a "${PREP_LOG}"
    exit 1
fi

if [[ ${#CACHE_DIRS[@]} -lt 2 ]]; then
    echo "ERROR: Fewer than 2 samples passed bustools. Cannot run multi-sample MAP EM." \
        | tee -a "${PREP_LOG}"
    exit 1
fi

echo "" | tee -a "${PREP_LOG}"
echo "============================================================" | tee -a "${PREP_LOG}"
echo "All bustools preprocessing complete (${#CACHE_DIRS[@]} samples)." | tee -a "${PREP_LOG}"
echo "Sample dirs passed to MAP EM:"                                   | tee -a "${PREP_LOG}"
for d in "${CACHE_DIRS[@]}"; do echo "  $d" | tee -a "${PREP_LOG}"; done
echo "============================================================" | tee -a "${PREP_LOG}"

# ============================================================
# Step 4: Multi-sample MAP EM
# main_multisample_joli.py creates its own timestamped run dir
# under OUTPUT_BASE and writes all results there.
# ============================================================
echo "" | tee -a "${PREP_LOG}"
echo "Step 4: Multi-sample MAP EM (main_multisample_joli.py)" | tee -a "${PREP_LOG}"

"${PYTHON}" "${SCRIPT_DIR}/../main_multisample_joli.py" \
    --sample_dirs        "${CACHE_DIRS[@]}" \
    --sample_names       "${SAMPLE_NAMES[@]}" \
    --results_base       "${OUTPUT_BASE}" \
    --eff_len_mode       "${EFF_LEN_MODE}" \
    --convergence_mode   "${CONVERGENCE_MODE}" \
    --max_em_rounds      "${MAX_EM_ROUNDS}" \
    --min_em_rounds      "${MIN_EM_ROUNDS}" \
    --max_gd_rounds      "${MAX_GD_ROUNDS}" \
    --gd_lr              "${GD_LR}" \
    --alpha_initial      "${ALPHA_INITIAL}" \
    --gd_convergence_tol "${GD_CONVERGENCE_TOL}" \
    --gd_steps_per_round "${GD_STEPS_PER_ROUND}" \
    --min_read_support   "${MIN_READ_SUPPORT}" \
    --loop_mode           "${LOOP_MODE}" \
    --save_snapshots      "${SAVE_SNAPSHOTS}" \
    --snapshot_interval   "${SNAPSHOT_INTERVAL}" \
    --experiment_comment  "${EXPERIMENT_COMMENT}" \
    2>&1 | tee -a "${PREP_LOG}"

RUN_END=$(date +%s)
ELAPSED=$(( RUN_END - RUN_START ))

echo "" | tee -a "${PREP_LOG}"
echo "============================================================" | tee -a "${PREP_LOG}"
echo "run_multisample_joli.sh complete in ${ELAPSED}s"              | tee -a "${PREP_LOG}"
echo "Results created by main_multisample_joli.py under: ${OUTPUT_BASE}" | tee -a "${PREP_LOG}"
echo "Preprocessing log: ${PREP_LOG}"                               | tee -a "${PREP_LOG}"
echo "============================================================" | tee -a "${PREP_LOG}"
