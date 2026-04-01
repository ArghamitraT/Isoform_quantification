# Simulation Pipeline — Next Steps

**Status as of 2026-03-30:** All src modules written and tested. Remaining: entry point + SLURM wrappers.

---

## What Is Done

| Step | File | Tests |
|------|------|-------|
| ✅ 1 | `src/utility.py`             | 3/3 pass |
| ✅ 2 | `src/prepare_reference.py`   | 6/6 pass |
| ✅ 3 | `src/generate_abundances.py` | 8/8 pass |
| ✅ 4 | `src/simulate_illumina.py`   | 5/5 pass |
| ✅ 5 | `src/simulate_pacbio.py`     | 6/6 pass |
| ✅ 6 | `src/simulate_ont.py`        | 7/7 pass |

---

## Step 7 — `main.py`

Entry point that ties all phases together via a single CONFIG section.

**Flow:**
1. Phase 1: call `prepare_reference.py` logic (or subprocess it)
2. Phase 2: call `generate_abundances.py` logic
3. Phase 3: call each simulator (Illumina / PacBio / ONT) — each as subprocess so they can also be run independently via their SLURM scripts
4. Use `utility.py` helpers: `create_run_dir`, `save_runtime`, `save_code_snapshot`
5. Write `experiment_description.log` with full CONFIG dump
6. First line of `running.log`: `Script: <path/to/main.py>`

**CONFIG block** (at top of `main()`, clearly delimited):
```python
# ============================================================
# CONFIG — edit these variables before running; do not edit below
# ============================================================
reference_genome   = "/path/to/GRCh38.fa"
reference_gtf      = "/path/to/gencode.vXX.annotation.gtf"
output_prefix      = "human_sim"
abundance_mode     = "lognormal"    # "lognormal" | "custom" | "uniform"
custom_tpm_file    = ""
cancer_genes_file  = ""
cancer_fold_change = 5.0
dropout_fraction   = 0.3
illumina_count     = 100_000_000
pb_count           = 10_000_000
ont_count          = 30_000_000
ont_type           = "cDNA"
threads            = 16
seed               = 42
noise_reads        = False
test_mode          = False          # True = small read counts for quick test
rsem_model_file    = ""             # required for Illumina; see simulate_illumina.py docstring
# ============================================================
```

**test_mode values** (override read counts when `test_mode = True`):
- illumina_count → 10_000
- pb_count → 1_000
- ont_count → 5_000

**Output folder:** `files/results/exprmnt_{timestamp}/` — created by `utility.create_run_dir()`

---

## Step 8 — SLURM Scripts

Three independent scripts so each technology can be re-run without re-running the others.

### `submit_illumina.sh`
- `#SBATCH --cpus-per-task 16`
- `#SBATCH --mem 64G`
- `module load RSEM/1.3.3-foss-2022a SAMtools/1.21`
- `conda activate lrgsp_simulation`
- CONFIG section at top: paths, read count, threads, seed, rsem_model_file
- Calls: `python code/Simulations/src/simulate_illumina.py ...`
- First log line: `echo "Script: $(realpath "$0")" > "${LOG}"`

### `submit_pacbio.sh`
- `#SBATCH --cpus-per-task 16`
- `#SBATCH --mem 64G`
- `module load minimap2/2.29 SAMtools/1.21`
- `conda activate lrgsp_simulation`
- CONFIG section at top: paths, read count, threads, seed, error rates
- Calls: `python code/Simulations/src/simulate_pacbio.py ...`

### `submit_ont.sh`
- `#SBATCH --cpus-per-task 16`
- `#SBATCH --mem 128G`
- `module load minimap2/2.29 SAMtools/1.21`
- `conda activate lrgsp_simulation_2`   ← needs HTSeq (see note below)
- CONFIG section at top: paths, read count, threads, seed, ont_type
- Calls: `python code/Simulations/src/simulate_ont.py ...`

---

## Outstanding Issue — HTSeq for ONT

NanoSim requires `HTSeq`, which is **not** in `lrgsp_simulation`.
Before running `submit_ont.sh`, create a new env:

```bash
conda create -n lrgsp_simulation_2 --clone lrgsp_simulation
conda activate lrgsp_simulation_2
pip install HTSeq
conda env export > code/Simulations/env/lrgsp_simulation_2.yml
```

`submit_ont.sh` must activate `lrgsp_simulation_2`, not `lrgsp_simulation`.

---

## Outstanding Issue — RSEM Model File for Illumina

`simulate_illumina.py` requires a `*.model` file from a real RSEM run.
If you don't have one, generate it once from real Illumina RNA-seq data:

```bash
module load RSEM/1.3.3-foss-2022a SAMtools/1.21
rsem-calculate-expression --paired-end \
    read_1.fastq read_2.fastq \
    <rsem_reference_prefix> sample_name
# Model file → sample_name.stat/sample_name.model
```

Then set `rsem_model_file` in `main.py` CONFIG and `submit_illumina.sh` CONFIG.

---

## After Steps 7 & 8 — End-to-End Test

Run a quick end-to-end test using `test_mode = True`:

```bash
# Phase 1 + 2 on login node
conda activate lrgsp_simulation
python code/Simulations/main.py   # with test_mode = True

# Phase 3 via SLURM (submit all three)
sbatch code/Simulations/submit_illumina.sh
sbatch code/Simulations/submit_pacbio.sh
sbatch code/Simulations/submit_ont.sh   # requires lrgsp_simulation_2
```

Check that `files/results/exprmnt_{timestamp}/` contains all expected outputs.
