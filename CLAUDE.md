# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Multi-sample RNA isoform quantification from long-read sequencing (PacBio/ONT). Two parallel workstreams:

1. **`AT_code/`** — Existing custom EM algorithm (VI or MAP) with Bayesian Dirichlet hyperparameter optimization via gradient descent.
2. **`JOLI_Kallisto/`** — Active development: merging the AT_code EM/VI approach with [lr-kallisto](https://github.com/pachterlab/kallisto) to make quantification faster. New code goes here.

---

## Environment Setup

```bash
conda activate NanoCount_5
```

- yml/txt: `AT_code/NanoCount_5.yml`, `AT_code/NanoCount_5.txt`
- Key packages: PyTorch 2.3.1, Pyro-PPL 1.9.1, NanoCount 1.0.0, pysam 0.22.1, NumPy, SciPy, pandas, Matplotlib

> Per global rules: new dependencies go in a new environment (e.g. `NanoCount_6`) with an incremented yml saved in `AT_code/`.

---

## AT_code — How to Run

### Architecture

```
BAM files
  → process_bam_files.py       # Parse alignments → pickle dicts
  → EM_VIorMAP_GD_vector.py    # Core EM (VI or MAP)
      E-step: read-to-isoform assignment weights
      M-step: update theta (isoform abundances)
      Dirichlet opt: update alpha via GD (PyTorch, DirichletOptimizer_vector.py)
  → TSV output files
  → result_analysis/           # Downstream stats, correlation, plots
```

### Key files

| File | Role |
|------|------|
| `AT_code/main_EM_VIorMAP_GD_vector.py` | Primary entry point |
| `AT_code/EM_VIorMAP_GD_vector.py` | Core EM logic — VI or MAP, multi-sample |
| `AT_code/DirichletOptimizer_vector.py` | PyTorch Dirichlet GD optimizer |
| `AT_code/process_bam_files.py` | BAM → pickle conversion pipeline |
| `AT_code/generate_bash.py` | Generates + submits SLURM batch jobs |
| `AT_code/result_analysis/util.py` | `ExperimentFileProcessor` — parses output filenames via regex |
| `AT_code/result_analysis/generate_result_stat.py` | Spearman/Pearson correlation vs. ground truth |
| `AT_code/result_analysis/ThetaCorrWGrndTruth.py` | Compare theta estimates to SIRV ground truth |
| `AT_code/result_analysis/plt_experiment_stats.py` | Plot EM convergence and GD loss trajectories |
| `AT_code/simulation_code/simulation.py` | Synthetic data generation for algorithm validation |

---

### 1. BAM → Pickle (pre-processing)

Converts BAM alignment files into pickle dicts used by the EM algorithm.

```bash
conda activate NanoCount_5
cd AT_code
python process_bam_files.py
```

- Edit `process_bam_files.py` to set input BAM paths and output pickle directory before running.
- Output: `.pkl` files in the specified output folder.

---

### 2. EM Run — local / interactive

Run the EM algorithm directly from the command line (no SLURM).

```bash
conda activate NanoCount_5
cd AT_code
python main_EM_VIorMAP_GD_vector.py \
    --data_folder <pkl_files_dir> \
    --output_path <output_path> \
    --sample1 <file1> [<file2> ...] \
    --sample2 <file1> [<file2> ...] \
    --alpha_initial 1.0 \
    --GD_lr 0.01 \
    --max_em_rounds 30 \
    --experiment_num <1|2|4|5> \
    --EM_type <VI|MAP> \
    --dirichlet_process <theta|expectation_log_theta> \
    --process_bam_required <0|1>
```

**Experiment types:**
| `--experiment_num` | Meaning |
|--------------------|---------|
| `1` | Single sample, no GD |
| `2` | Two samples merged, no GD |
| `4` | Multi-sample + Dirichlet GD |
| `5` | Merged multi-sample + Dirichlet GD |

**Key flags:**
| Flag | Values | Notes |
|------|--------|-------|
| `--EM_type` | `VI` / `MAP` | Inference method |
| `--dirichlet_process` | `theta` / `expectation_log_theta` | GD optimization target |
| `--process_bam_required` | `0` / `1` | `1` = run BAM→pickle first; `0` = pkl files already exist |
| `--load` | `0` / `1` | `1` = resume from a saved checkpoint |
| `--load_filename` | path | Path to `.pkl` checkpoint; required when `--load 1` |

**Quick trial with SIRV data:**
```bash
python main_EM_VIorMAP_GD_vector.py \
    --data_folder ../SIRV_data/aln_E0/ \
    --output_path /tmp/sirv_test/ \
    --sample1 <sirv_file> --sample2 NA \
    --experiment_num 1 --EM_type MAP \
    --alpha_initial 1.0 --GD_lr 0.01 --max_em_rounds 10 \
    --process_bam_required 0
```

---

### 3. EM Run — SLURM batch (multiple hyperparameter combinations)

`generate_bash.py` generates one SLURM job per combination of samples × alpha × GD_lr × EM_rounds and submits them all.

**Step 1 — edit the CONFIG block inside `AT_code/generate_bash.py`:**

| Variable | What to set |
|----------|-------------|
| `samples_file_names` | List of sample pairs/groups (LR first, SR second) |
| `alpha_val_arr` | List of `alpha_initial` values to sweep |
| `GDlr_val_arr` | List of GD learning rates to sweep |
| `EM_round_arr` | List of max EM round counts to sweep |
| `experiment_num` | `1`, `2`, `4`, or `5` |
| `EM_type` | `'VI'` or `'MAP'` |
| `dirichlet_process` | `'theta'` or `'expectation_log_theta'` |
| `simulation` | `1` for simulation data, `0` for real PacBio data |
| `process_bam_required` | `0` (pkl already exists) or `1` |
| `hour`, `memory`, `nthred` | SLURM resource limits |
| `slurm_file_name` | Short label for the job name in SLURM |
| `readme_comment` | Human-readable note logged in the results readme |

**Step 2 — run:**
```bash
conda activate NanoCount_5
cd AT_code
python generate_bash.py
```

Generated `.sh` scripts land in:
`/gpfs/commons/home/atalukder/RNA_Splicing/files/cluster_job_submission_files/`

They are also submitted to SLURM automatically by the script (`sbatch` is called inside `gen_combination()`).

**To resume from a checkpoint**, set `copy_needed = 1` and `from_where_to_copy = "exprmnt_<timestamp>"` before running.

---

### 4. Result Analysis

All analysis scripts are in `AT_code/result_analysis/`. Edit the config variables at the top of each script, then run directly.

#### 4a. Correlation vs. ground truth — `generate_result_stat.py`

Computes Spearman and Pearson correlation between EM estimates and ground truth across replicas.

```bash
conda activate NanoCount_5
cd AT_code/result_analysis
# Edit main_result_dir, experiment_file, simulation, experiment inside the script
python generate_result_stat.py
```

#### 4b. Theta vs. SIRV ground truth — `ThetaCorrWGrndTruth.py`

Compares theta estimates to SIRV spike-in ground truth values.

```bash
conda activate NanoCount_5
cd AT_code/result_analysis
# Edit main_result_dir, experiment_file, groundTruth_main_dir, groundTruth_file inside the script
python ThetaCorrWGrndTruth.py
```

Key config variables:
| Variable | Meaning |
|----------|---------|
| `main_result_dir` | Top-level results directory |
| `experiment_file` | Timestamped subfolder name (`exprmnt_...`) |
| `simulation` | `1` = simulation data, `0` = real data |
| `experiment` | Experiment type (`1`, `2`, `4`, `5`) |
| `groundTruth_main_dir` | Path to SIRV ground truth files |
| `groundTruth_file` | Dict mapping sample IDs to ground truth filenames |

#### 4c. EM convergence and GD loss plots — `plt_experiment_stats.py`

Plots alpha/convergence/loss trajectories from saved `.pkl` stats files.

```bash
conda activate NanoCount_5
cd AT_code/result_analysis
# Edit final_result_dir, main_dir, experiment_file, fig_generic_name inside the script
python plt_experiment_stats.py
```

Figures are saved to `<experiment_dir>/figures/`.

---

### Output filename convention (AT_code)

```
output_[PacIllu|Simulation]_VIGD_token_<hash>_sample<N>_file<idx>_ds<pct>num<id>aln<day><replica><len>_GDlr_<rate>_AlphaInitial_<alpha>_EMround_<rounds>_
```

Use `ExperimentFileProcessor` in `AT_code/result_analysis/util.py` to parse these — do not parse manually.

---

## JOLI_Kallisto — How to Run

### Architecture

```
FASTQ/FASTA reads
  → run_lr_kallisto.sh
      kallisto bus       # align reads → output.bus
      bustools sort      # sort .bus file
      bustools count     # count matrix (count.mtx, count.ec.txt)
      kallisto quant-tcc # EM quantification → abundance.tsv
  → per-sample output subfolder
```

### Key files

| File | Role |
|------|------|
| `JOLI_Kallisto/run_lr_kallisto.sh` | Full pipeline logic — local or SLURM compatible |
| `JOLI_Kallisto/submit_lr_kallisto.sh` | Thin SLURM submission wrapper |

### lr-kallisto reference files

Base directory: `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/`

| File | Path |
|------|------|
| kallisto binary | `.../kallisto/build/src/kallisto` (v0.51.1) |
| Index (new, used by pipeline) | `.../new_index.idx` |
| Transcript-to-gene map | `.../t2g.txt` |
| Transcriptome FASTA | `.../transcriptome.fasta` |
| Simulation index | `.../sim_index` |

---

### 1. Configure the pipeline

Edit the **CONFIG section** at the top of `JOLI_Kallisto/run_lr_kallisto.sh`:

| Variable | What to set |
|----------|-------------|
| `KALLISTO` | Full path to kallisto binary |
| `BUSTOOLS` | `bustools` (if on PATH) or full path |
| `INDEX_FILE` | Path to kallisto index |
| `T2G_FILE` | Path to transcript-to-gene map |
| `READ_TYPE` | `long` (PacBio/ONT) or `short` (paired Illumina) |
| `THREADS` | CPU threads (default: 32) |
| `THRESHOLD` | kallisto bus threshold (default: 0.8) |
| `OUTPUT_BASE` | Base results directory (default: `/gpfs/.../files/results`) |
| `SAMPLES` | Bash array of `"sample_name  reads_dir  reads_file"` entries |

**Sample array format:**
```bash
# Long-read (single file per sample)
SAMPLES=(
    "sim2  /path/to/sim/reads/  PacBio.simulated.fasta"
    "ds52  /path/to/real/reads/ ds_52.fastq"
)

# Short-read (paired: two files per sample)
SAMPLES=(
    "sample1  /path/to/reads/  R1.fastq  R2.fastq"
)
```

---

### 2. Run locally

```bash
bash JOLI_Kallisto/run_lr_kallisto.sh
```

No SLURM required. The script creates the timestamped output folder and runs all samples sequentially.

---

### 3. Run on SLURM

**Option A — via wrapper (recommended):**
```bash
bash JOLI_Kallisto/submit_lr_kallisto.sh
```

The wrapper (`submit_lr_kallisto.sh`) pre-creates the timestamped output folder and calls `sbatch` with resource flags. Edit the CONFIG block inside `submit_lr_kallisto.sh` to change SLURM resources:

| Variable | Default |
|----------|---------|
| `MEM` | `100G` |
| `CPUS` | `32` |
| `TIME` | `24:00:00` |
| `MAIL_USER` | `atalukder@nygenome.org` |

**Option B — direct sbatch:**
```bash
sbatch JOLI_Kallisto/run_lr_kallisto.sh
```

(SLURM `#SBATCH` headers inside `run_lr_kallisto.sh` are picked up automatically.)

---

### Per-sample output layout

```
/gpfs/.../files/results/exprmnt_{timestamp}/
├── experiment_description.log   # config dump + sample list
├── running.log                  # combined stdout/stderr for all steps
├── runtime.txt                  # total wall-clock time
├── code_snapshot/
│   └── run_lr_kallisto.sh       # exact copy of this script at run time
├── {sample_name_1}/
│   ├── output.bus               # raw bus file
│   ├── sorted.bus               # sorted bus file
│   ├── count.mtx                # TCC count matrix
│   ├── count.ec.txt             # equivalence classes
│   ├── transcripts.txt          # transcript list
│   ├── abundance.tsv            # final isoform quantification
│   └── run_info.json            # kallisto run metadata
└── {sample_name_2}/
    └── ...
```

---

## Data Paths

| Data | Path |
|------|------|
| Real PacBio pkl files | `/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln_pklfiles/` |
| Simulation pkl files | `/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/pklfiles/` |
| SIRV trial data | `SIRV_data/aln_E0/`, `SIRV_data/aln_E2/` |
| lr-kallisto base | `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/` |

---

## Shared Output Format

All results (AT_code and JOLI_Kallisto) save to:
```
/gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmnt_{timestamp}/
├── experiment_description.log   # config + what/why/expected outcome
├── running.log                  # full stdout/stderr of the run
├── best_checkpoint/             # model weights (AT_code: use --load 1 --load_filename to resume)
├── results_summary.txt          # accuracy tables, key metrics
├── runtime.txt                  # total wall-clock time
└── code_snapshot/               # selective copy of the code folder (see rules below)
```

### Code snapshot rules

`code_snapshot/` must mirror the full directory structure of the code folder, but **only copy files with these extensions**:

| Extension | Examples |
|-----------|---------|
| `.py` | all Python scripts |
| `.sh` | all shell/SLURM scripts |
| `.txt` | requirements, notes |
| `.yml` / `.yaml` | conda environment files |

**Do NOT copy** large or non-code files: `.pkl`, `.idx`, `.bam`, `.fastq`, `.fasta`, `.bus`, `.mtx`, `.log`, or any other data/binary assets. These are large and not needed for reproducibility.

**Implementation** — bash:
```bash
find code/ \( -name "*.py" -o -name "*.sh" -o -name "*.txt" -o -name "*.yml" -o -name "*.yaml" \) \
    | while read f; do
        dest="code_snapshot/${f#code/}"
        mkdir -p "$(dirname "$dest")"
        cp "$f" "$dest"
    done
```
Python (`save_code_snapshot()` in `utility.py`):
```python
SNAPSHOT_EXTS = {".py", ".sh", ".txt", ".yml", ".yaml"}
for src in Path("code/").rglob("*"):
    if src.suffix in SNAPSHOT_EXTS:
        dest = run_dir / "code_snapshot" / src.relative_to("code/")
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
```

---

## Utility Helpers

Per global rules, a `utility.py` should exist providing `create_run_dir()`, `save_runtime()`, and `save_code_snapshot()`. Add this to `AT_code/` (for AT_code runs) and `JOLI_Kallisto/` (for kallisto-merged runs) as development progresses.
