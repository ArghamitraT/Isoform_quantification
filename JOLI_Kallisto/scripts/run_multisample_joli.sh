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
SHORT_INDEX_FILE="/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/new_index_sr.idx"
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
THREADS=4
THRESHOLD=0.8        # kallisto bus alignment threshold (long-read only)

# --- Python / conda ---
CONDA_ENV="Joli_kallisto"
PYTHON=""            # leave empty to auto-resolve from conda env

# --- Output base for experiment results ---
OUTPUT_BASE="/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

# --- JOLI MAP EM settings (passed to main_multisample_joli.py) ---
EFF_LEN_MODE="auto"        # "auto" = resolved per sample (long→uniform, short→fld) | "uniform" | "fld" | "kallisto"
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
EM_INCLUDE_SINGLE_TX=auto   # "auto" = resolved per sample (long→false, short→true)
#                           # true  = single-tx counts included in EM M-step every round
#                           # false = original kallisto behaviour (single-tx added post-convergence)
GD_INCLUDE_SINGLE_TX=true  # false (default) = single-tx excluded from GD theta
#                           #   GD gradient driven only by multi-tx (ambiguous) ECs
#                           # true  = old behaviour (single-tx baked into GD theta)
LOOP_MODE="gd_wrapper"  # Training loop structure:
#                        #   "gd_wrapper" — GD outer loop, EM to convergence per round (default)
#                        #   "em_wrapper" — EM convergence drives outer loop, 1 EM step +
#                        #                 GD_STEPS_PER_ROUND Adam steps per iteration (AT_code)

# Free-text description of this experiment run — saved in experiment_description.log.
# Overridden automatically when submitted via submit_JKMS_slurm.py.
EXPERIMENT_COMMENT="JK MS real — edit or use submit_JKMS_slurm.py to set per run"

# Exact exprmnt_* output directory to use, pre-assigned by a batch submitter
# (submit_JKMS_slurm.py sets this in the override file) so concurrent jobs
# never race on independently-generated timestamps. Leave empty for manual/
# standalone runs -- main_multisample_joli.py falls back to its own timestamp.
RUN_DIR=""

# --- Post-EM analysis ---
# DATA_TYPE="sim"   → run GT comparison after EM (simulation data with known ground truth)
# DATA_TYPE="real"  → run replicate comparison after EM (biological replicates, no GT)
DATA_TYPE="real"

# For sim: map sample_name to its ground truth TSV ("sample_name:path" pairs)
GT_MAP=(
    # Uncomment to switch back to simulation:
    # "sim2:/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample2_gt.tsv"
    # "sim_sr_s1:/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/ill_sample1_gt.tsv"
    # "sim1:/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample1_gt.tsv"
    # "sim_sr_s2:/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/ill_sample2_gt.tsv"
)

# For real: within-run replicate pairs — left empty because each run mixes LR and SR
# from different days (not biological replicates of each other).
# Cross-run comparison (e.g. Run1/flnc_31 vs Run2/flnc_32) is printed at job end.
REPLICATE_PAIRS=()

# --- Samples ---
# Each entry: "sample_name  read_type  reads_dir  file1  [file2]"
#   read_type = "long"  → single file, kallisto --long --threshold flags, KMER_SIZE_LONG
#   read_type = "short" → two files (R1 R2), standard paired-end flags, KMER_SIZE_SHORT
# When submitting via submit_all_ms_runs.sh, SAMPLES is set automatically per run.
# To run manually, uncomment one block below:
_LR="/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long"
_SR="/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/short"
SAMPLES=(
    # Run 1: flnc_31(LR) + ill_51(SR)
    "flnc_31  long   ${_LR}  flnc_31.fastq"
    "ill_51   short  ${_SR}  ill_51_R1.fastq  ill_51_R2.fastq"
    # Run 2: flnc_32(LR) + ill_52(SR)
    # "flnc_32  long   ${_LR}  flnc_32.fastq"
    # "ill_52   short  ${_SR}  ill_52_R1.fastq  ill_52_R2.fastq"
    # Run 3: flnc_51(LR) + ill_31(SR)
    # "flnc_51  long   ${_LR}  flnc_51.fastq"
    # "ill_31   short  ${_SR}  ill_31_R1.fastq  ill_31_R2.fastq"
    # Run 4: flnc_52(LR) + ill_32(SR)
    # "flnc_52  long   ${_LR}  flnc_52.fastq"
    # "ill_32   short  ${_SR}  ill_32_R1.fastq  ill_32_R2.fastq"
    # Simulation (uncomment to switch back):
    # "sim2       long   /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/  ds_100_num1_aln_21_long.fasta"
    # "sim_sr_s1  short  /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/  ds_100_num1_aln_01_short_1.fq  ds_100_num1_aln_01_short_2.fq"
)

# ============================================================
# END CONFIG
# ============================================================

# If launched via submit_all_ms_runs.sh, source the per-run override (sets SAMPLES + EXPERIMENT_COMMENT)
if [[ -n "${JOLI_MS_OVERRIDE:-}" ]] && [[ -f "${JOLI_MS_OVERRIDE}" ]]; then
    # shellcheck source=/dev/null
    source "${JOLI_MS_OVERRIDE}"
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
JOLI_KALLISTO_DIR="${SCRIPT_DIR}/.."

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
    # Args: cache_dir [require_flens=false]
    # Pass require_flens=true for short-read samples: flens.txt must exist for eff_len_mode=fld.
    local cache_dir="$1"
    local require_flens="${2:-false}"
    [[ -f "${cache_dir}/count.mtx"       ]] && \
    [[ -f "${cache_dir}/matrix.ec"       ]] && \
    [[ -f "${cache_dir}/transcripts.txt" ]] && \
    { [[ "${require_flens}" != "true" ]] || [[ -f "${cache_dir}/flens.txt" ]]; } && \
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
CACHE_DIRS=()      # will be passed as --sample_dirs to the Python script
SAMPLE_NAMES=()    # will be passed as --sample_names to the Python script
EFF_LEN_MODES_LIST=()   # per-sample eff_len_mode resolved from read_type
EM_SINGLE_TX_LIST=()    # per-sample em_include_single_tx resolved from read_type

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

    # Resolve per-sample settings from read_type
    if [[ "${SAMPLE_READ_TYPE}" == "long" ]]; then
        KMER_SIZE="${KMER_SIZE_LONG}"
        SAMPLE_INDEX="${INDEX_FILE}"
        SAMPLE_EFF_LEN_MODE="uniform"   # long reads: uniform effective lengths
        SAMPLE_EM_SINGLE_TX="false"     # long reads: exclude single-tx from EM loop
    else
        KMER_SIZE="${KMER_SIZE_SHORT}"
        SAMPLE_INDEX="${SHORT_INDEX_FILE}"
        SAMPLE_EFF_LEN_MODE="fld"       # short reads: FLD-based effective lengths
        SAMPLE_EM_SINGLE_TX="true"      # short reads: include single-tx in EM loop
    fi

    STEM=$(get_stem "${READS_FILE1}")
    CACHE_DIR="${READS_DIR}/kallisto_output/${STEM}"

    echo "" | tee -a "${PREP_LOG}"
    echo "------------------------------------------------------------" | tee -a "${PREP_LOG}"
    echo "Sample:      ${SAMPLE_NAME}"                                  | tee -a "${PREP_LOG}"
    echo "Read type:   ${SAMPLE_READ_TYPE}  (k=${KMER_SIZE})"          | tee -a "${PREP_LOG}"
    echo "eff_len_mode:${SAMPLE_EFF_LEN_MODE}"                         | tee -a "${PREP_LOG}"
    echo "Reads:       ${READS_DIR}/${READS_FILE1}"                     | tee -a "${PREP_LOG}"
    if [[ -n "${READS_FILE2}" ]]; then
        echo "             ${READS_DIR}/${READS_FILE2}"                 | tee -a "${PREP_LOG}"
    fi
    echo "Cache dir:   ${CACHE_DIR}"                                    | tee -a "${PREP_LOG}"
    echo "------------------------------------------------------------" | tee -a "${PREP_LOG}"

    # Short reads require flens.txt for eff_len_mode=fld
    REQUIRE_FLENS="false"
    [[ "${SAMPLE_READ_TYPE}" == "short" ]] && REQUIRE_FLENS="true"

    if cache_complete "${CACHE_DIR}" "${REQUIRE_FLENS}"; then
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
                -i "${SAMPLE_INDEX}" \
                -o "${CACHE_DIR}" \
                -t "${THREADS}" \
                --long \
                --threshold "${THRESHOLD}" \
                --unmapped \
                "${READS_DIR}/${READS_FILE1}" \
                2>&1 | tee -a "${PREP_LOG}"
        else
            # --paired writes flens.txt (FLD histogram) needed for eff_len_mode="fld"
            "${KALLISTO}" bus \
                -x bulk \
                -i "${SAMPLE_INDEX}" \
                -o "${CACHE_DIR}" \
                -t "${THREADS}" \
                --paired \
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

        if ! cache_complete "${CACHE_DIR}" "${REQUIRE_FLENS}"; then
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
    EFF_LEN_MODES_LIST+=("${SAMPLE_EFF_LEN_MODE}")
    EM_SINGLE_TX_LIST+=("${SAMPLE_EM_SINGLE_TX}")
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
# We pass --run_dir_output to a temp file so we can read the path back.
# ============================================================
echo "" | tee -a "${PREP_LOG}"
echo "Step 4: Multi-sample MAP EM (main_multisample_joli.py)" | tee -a "${PREP_LOG}"

RUN_DIR_TMP="$(mktemp /tmp/joli_ms_rundir_XXXXXX.txt)"

# Build optional fallback args (omit when "auto" — per-sample lists handle it)
EXTRA_ARGS=()
[[ "${EFF_LEN_MODE}" != "auto" ]]        && EXTRA_ARGS+=(--eff_len_mode "${EFF_LEN_MODE}")
[[ "${EM_INCLUDE_SINGLE_TX}" != "auto" ]] && EXTRA_ARGS+=(--em_include_single_tx "${EM_INCLUDE_SINGLE_TX}")
[[ -n "${RUN_DIR}" ]]                     && EXTRA_ARGS+=(--run_dir "${RUN_DIR}")

"${PYTHON}" "${SCRIPT_DIR}/../main_multisample_joli.py" \
    --sample_dirs             "${CACHE_DIRS[@]}" \
    --sample_names            "${SAMPLE_NAMES[@]}" \
    --results_base            "${OUTPUT_BASE}" \
    --eff_len_modes           "${EFF_LEN_MODES_LIST[@]}" \
    --transcriptome_fasta     "${TRANSCRIPTOME_FASTA}" \
    --em_include_single_tx_list "${EM_SINGLE_TX_LIST[@]}" \
    --convergence_mode        "${CONVERGENCE_MODE}" \
    --max_em_rounds           "${MAX_EM_ROUNDS}" \
    --min_em_rounds           "${MIN_EM_ROUNDS}" \
    --max_gd_rounds           "${MAX_GD_ROUNDS}" \
    --gd_lr                   "${GD_LR}" \
    --alpha_initial           "${ALPHA_INITIAL}" \
    --gd_convergence_tol      "${GD_CONVERGENCE_TOL}" \
    --gd_steps_per_round      "${GD_STEPS_PER_ROUND}" \
    --min_read_support        "${MIN_READ_SUPPORT}" \
    --loop_mode               "${LOOP_MODE}" \
    --save_snapshots          "${SAVE_SNAPSHOTS}" \
    --snapshot_interval       "${SNAPSHOT_INTERVAL}" \
    --gd_include_single_tx    "${GD_INCLUDE_SINGLE_TX}" \
    --experiment_comment      "${EXPERIMENT_COMMENT}" \
    --run_dir_output          "${RUN_DIR_TMP}" \
    "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}" \
    2>&1 | tee -a "${PREP_LOG}"

# Read the run_dir that Python created
RUN_DIR=""
if [[ -f "${RUN_DIR_TMP}" ]]; then
    RUN_DIR="$(tr -d '[:space:]' < "${RUN_DIR_TMP}")"
    rm -f "${RUN_DIR_TMP}"
fi

# ============================================================
# Step 5: Post-EM analysis
# ============================================================
if [[ -n "${RUN_DIR}" ]]; then
    echo "" | tee -a "${PREP_LOG}"
    echo "Step 5: Post-EM analysis (DATA_TYPE=${DATA_TYPE})" | tee -a "${PREP_LOG}"

    if [[ "${DATA_TYPE}" == "sim" ]]; then
        # Build --gt_override args for each sample's GT file
        GT_OVERRIDE_ARGS=()
        for entry in "${GT_MAP[@]}"; do
            GT_OVERRIDE_ARGS+=("--gt_override" "${entry}")
        done

        "${PYTHON}" "${JOLI_KALLISTO_DIR}/analysis/run_gt_comparison.py" \
            "${RUN_DIR}" \
            "${GT_OVERRIDE_ARGS[@]+"${GT_OVERRIDE_ARGS[@]}"}" \
            2>&1 | tee -a "${PREP_LOG}" \
            || echo "WARNING: GT comparison failed (non-fatal)" | tee -a "${PREP_LOG}"

    elif [[ "${DATA_TYPE}" == "real" ]]; then
        if [[ ${#REPLICATE_PAIRS[@]} -eq 0 ]]; then
            echo "" | tee -a "${PREP_LOG}"
            echo "  [INFO] No within-run replicate pairs configured." | tee -a "${PREP_LOG}"
            echo "  Cross-run comparison must be run after all jobs complete. Example:" | tee -a "${PREP_LOG}"
            echo "    cd ${JOLI_KALLISTO_DIR}" | tee -a "${PREP_LOG}"
            echo "    conda activate NanoCount_5" | tee -a "${PREP_LOG}"
            echo "    python analysis/run_replicate_comparison.py \\" | tee -a "${PREP_LOG}"
            echo "        ${RUN_DIR}/<sample1>/abundance.tsv \\" | tee -a "${PREP_LOG}"
            echo "        <OTHER_RUN_DIR>/<sample2>/abundance.tsv \\" | tee -a "${PREP_LOG}"
            echo "        --names <sample1> <sample2>" | tee -a "${PREP_LOG}"
        else
            for pair in "${REPLICATE_PAIRS[@]}"; do
                s1="${pair%%:*}"
                s2="${pair##*:}"
                ab1="${RUN_DIR}/${s1}/abundance.tsv"
                ab2="${RUN_DIR}/${s2}/abundance.tsv"
                if [[ ! -f "${ab1}" || ! -f "${ab2}" ]]; then
                    echo "" | tee -a "${PREP_LOG}"
                    echo "  [SKIP] Analysis NOT done for ${s1} vs ${s2} — abundance.tsv missing for one or both samples." | tee -a "${PREP_LOG}"
                    echo "  Run analysis later with:" | tee -a "${PREP_LOG}"
                    echo "    cd ${JOLI_KALLISTO_DIR}" | tee -a "${PREP_LOG}"
                    echo "    conda activate NanoCount_5" | tee -a "${PREP_LOG}"
                    echo "    python analysis/run_replicate_comparison.py \\" | tee -a "${PREP_LOG}"
                    echo "        ${ab1} \\" | tee -a "${PREP_LOG}"
                    echo "        ${ab2} \\" | tee -a "${PREP_LOG}"
                    echo "        --names ${s1} ${s2}" | tee -a "${PREP_LOG}"
                    continue
                fi
                echo "  Replicate comparison: ${s1} vs ${s2}" | tee -a "${PREP_LOG}"
                "${PYTHON}" "${JOLI_KALLISTO_DIR}/analysis/run_replicate_comparison.py" \
                    "${ab1}" \
                    "${ab2}" \
                    --names "${s1}" "${s2}" \
                    2>&1 | tee -a "${PREP_LOG}" \
                    || echo "WARNING: Replicate comparison failed for ${s1}:${s2}" \
                        | tee -a "${PREP_LOG}"
            done
        fi
    fi
else
    echo "WARNING: Could not determine RUN_DIR from Python output; skipping Step 5 analysis." \
        | tee -a "${PREP_LOG}"
fi

RUN_END=$(date +%s)
ELAPSED=$(( RUN_END - RUN_START ))

echo "" | tee -a "${PREP_LOG}"
echo "============================================================" | tee -a "${PREP_LOG}"
echo "run_multisample_joli.sh complete in ${ELAPSED}s"              | tee -a "${PREP_LOG}"
if [[ -n "${RUN_DIR}" ]]; then
    echo "Results: ${RUN_DIR}"                                       | tee -a "${PREP_LOG}"
fi
echo "Preprocessing log: ${PREP_LOG}"                               | tee -a "${PREP_LOG}"
echo "============================================================" | tee -a "${PREP_LOG}"
