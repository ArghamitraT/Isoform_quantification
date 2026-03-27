# experiments.md

Detailed script descriptions, run instructions, and configuration reference for both workstreams.

---

## AT_code — Detailed Run Guide

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

### Step 1 — BAM → Pickle (pre-processing)

Converts BAM alignment files into pickle dicts used by the EM algorithm.

```bash
conda activate NanoCount_5
cd AT_code
python process_bam_files.py
```

- Edit `process_bam_files.py` to set input BAM paths and output pickle directory before running.
- Output: `.pkl` files in the specified output folder.

---

### Step 2a — EM Run: local / interactive

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

### Step 2b — EM Run: SLURM batch (hyperparameter sweep)

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

They are also submitted to SLURM automatically (`sbatch` is called inside `gen_combination()`).

**To resume from a checkpoint**, set `copy_needed = 1` and `from_where_to_copy = "exprmnt_<timestamp>"` before running.

---

### Step 3 — Result Analysis

All analysis scripts are in `AT_code/result_analysis/`. Edit the config variables at the top of each script, then run directly.

#### 3a. Correlation vs. ground truth — `generate_result_stat.py`

Computes Spearman and Pearson correlation between EM estimates and ground truth across replicas.

```bash
conda activate NanoCount_5
cd AT_code/result_analysis
# Edit main_result_dir, experiment_file, simulation, experiment inside the script
python generate_result_stat.py
```

#### 3b. Theta vs. SIRV ground truth — `ThetaCorrWGrndTruth.py`

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

#### 3c. EM convergence and GD loss plots — `plt_experiment_stats.py`

Plots alpha/convergence/loss trajectories from saved `.pkl` stats files.

```bash
conda activate NanoCount_5
cd AT_code/result_analysis
# Edit final_result_dir, main_dir, experiment_file, fig_generic_name inside the script
python plt_experiment_stats.py
```

Figures are saved to `<experiment_dir>/figures/`.

---

### AT_code Output Filename Convention

```
output_[PacIllu|Simulation]_VIGD_token_<hash>_sample<N>_file<idx>_ds<pct>num<id>aln<day><replica><len>_GDlr_<rate>_AlphaInitial_<alpha>_EMround_<rounds>_
```

Use `ExperimentFileProcessor` in `AT_code/result_analysis/util.py` to parse these — do not parse manually.

---

## JOLI_Kallisto — Detailed Run Guide

### Architecture

Two parallel pipelines — lr-kallisto (C++ EM baseline) and JOLI (Python MAP EM):

```
FASTQ/FASTA reads
  → scripts/run_lr_kallisto.sh         # lr-kallisto: C++ EM baseline
      kallisto bus                     #   align reads → output.bus
      bustools sort/count              #   count matrix (TCC format)
      kallisto quant-tcc               #   EM → abundance.tsv
  → scripts/run_joli_kallisto.sh       # JOLI single-sample: bustools + Python EM
      kallisto bus                     #   align reads → output.bus
      bustools sort/count              #   count matrix (TCC format)
      python main_joli.py              #   JOLI EM → abundance.tsv
  → main_multisample_joli.py           # JOLI multi-sample MAP EM (Phase 2)
      core/multi_sample_em.py          #   per-sample MAP EM + shared Dirichlet GD
      core/dirichlet_optimizer.py      #   Adam optimizer on shared alpha
```

### Folder Structure

```
JOLI_Kallisto/
├── core/                        # Library modules (not run directly)
│   ├── load_tcc.py              #   Parse bustools output → TCCData
│   ├── weights.py               #   Compute effective lengths and EC weights
│   ├── em_algorithm.py          #   JoliEM class — plain EM + MAP EM
│   ├── output_writer.py         #   Write abundance.tsv
│   ├── dirichlet_optimizer.py   #   Adam optimizer for shared Dirichlet alpha
│   └── multi_sample_em.py       #   MultiSampleJoliEM — outer GD + per-sample MAP EM
├── scripts/                     # Executable pipeline scripts
│   ├── run_joli_kallisto.sh     #   Full single-sample pipeline (bustools + JOLI EM)
│   ├── run_lr_kallisto.sh       #   Full lr-kallisto baseline pipeline
│   ├── submit_joli_pipeline.sh  #   SLURM wrapper for run_joli_kallisto.sh
│   └── submit_lr_kallisto.sh    #   SLURM wrapper for run_lr_kallisto.sh
├── analysis/                    # Analysis and comparison tools
│   ├── compare_abundance_files.py   # Compare JK vs LK abundance.tsv
│   └── run_batch_comparison.py      # Batch comparison across experiment folders
├── test/                        # Unit tests (run from JOLI_Kallisto/ root)
│   ├── test_JolitoKallisto.py   #   Tests for load_tcc, weights, output_writer
│   ├── test_em_algorithm.py     #   Tests for JoliEM (plain + MAP)
│   ├── test_dirichlet_optimizer.py  # Tests for DirichletOptimizer
│   ├── test_multi_sample_em.py  #   Tests for MultiSampleJoliEM
│   ├── test_main_multisample_joli.py # Tests for main_multisample_joli helpers
│   └── test_helpers.py          #   Shared make_sample_dir() and check() utilities
├── plans/                       # Design documents
├── main_joli.py                 # Entry point: single-sample JOLI EM (Python only)
├── main_multisample_joli.py     # Entry point: multi-sample MAP EM (Phase 2)
└── main_pipeline.py             # Full pipeline entry point
```

### Key Files

| File | Role |
|------|------|
| `main_joli.py` | Single-sample JOLI EM entry point (assumes bustools output exists) |
| `main_multisample_joli.py` | Multi-sample MAP EM entry point — Phase 2 |
| `scripts/run_joli_kallisto.sh` | Full single-sample pipeline: bustools → JOLI EM |
| `scripts/run_lr_kallisto.sh` | Full lr-kallisto baseline pipeline |
| `scripts/submit_joli_pipeline.sh` | SLURM wrapper for `run_joli_kallisto.sh` |
| `scripts/submit_lr_kallisto.sh` | SLURM wrapper for `run_lr_kallisto.sh` |
| `core/em_algorithm.py` | `JoliEM` class — plain EM and MAP EM with Dirichlet prior |
| `core/multi_sample_em.py` | `MultiSampleJoliEM` — outer GD loop + per-sample MAP EM |
| `core/dirichlet_optimizer.py` | Adam optimizer for shared Dirichlet `alpha` |
| `analysis/compare_abundance_files.py` | Compare JK vs LK abundance.tsv outputs |

### lr-kallisto Reference Files

Base directory: `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/`

| File | Path |
|------|------|
| kallisto binary | `.../kallisto/build/src/kallisto` (v0.51.1) |
| Index (new, used by pipeline) | `.../new_index.idx` |
| Transcript-to-gene map | `.../t2g.txt` |
| Transcriptome FASTA | `.../transcriptome.fasta` |
| Simulation index | `.../sim_index` |

---

### Pipeline A — lr-kallisto baseline (C++ EM)

#### Step 1 — Configure

Edit the **CONFIG section** at the top of `JOLI_Kallisto/scripts/run_lr_kallisto.sh`:

| Variable | What to set |
|----------|-------------|
| `KALLISTO` | Full path to kallisto binary |
| `BUSTOOLS` | `bustools` (if on PATH) or full path |
| `INDEX_FILE` | Path to kallisto index |
| `T2G_FILE` | Path to transcript-to-gene map |
| `READ_TYPE` | `long` (PacBio/ONT) or `short` (paired Illumina) |
| `THREADS` | CPU threads (default: 32) |
| `THRESHOLD` | kallisto bus threshold (default: 0.8) |
| `OUTPUT_BASE` | Base results directory |
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

#### Step 2 — Run locally

```bash
conda activate NanoCount_5
bash JOLI_Kallisto/scripts/run_lr_kallisto.sh
```

#### Step 3 — Run on SLURM

```bash
# Option A — via wrapper (recommended):
bash JOLI_Kallisto/scripts/submit_lr_kallisto.sh

# Option B — direct sbatch:
sbatch JOLI_Kallisto/scripts/run_lr_kallisto.sh
```

Edit the CONFIG block inside `submit_lr_kallisto.sh` for SLURM resources (`MEM`, `CPUS`, `TIME`, `MAIL_USER`).

---

### Pipeline B — JOLI single-sample (full pipeline: bustools + Python EM)

#### Step 1 — Configure

Edit the **CONFIG section** at the top of `JOLI_Kallisto/scripts/run_joli_kallisto.sh`:

| Variable | What to set |
|----------|-------------|
| `KALLISTO` | Full path to kallisto binary |
| `BUSTOOLS` | `bustools` (if on PATH) or full path |
| `INDEX_FILE` | Path to kallisto index |
| `T2G_FILE` | Path to transcript-to-gene map |
| `READ_TYPE` | `long` (PacBio/ONT) or `short` (paired Illumina) |
| `EFF_LEN_MODE` | `uniform` (all 1.0) or `kallisto` (length-based) |
| `CONVERGENCE_MODE` | `kallisto` (matches lr-kallisto exactly) or `joli` (faster, for MAP) |
| `EM_TYPE` | `plain` (Phase 1) — MAP/VI not yet supported in single-sample entry |
| `THREADS` | CPU threads (default: 32) |
| `OUTPUT_BASE` | Base results directory |
| `SAMPLES` | Same format as lr-kallisto above |

#### Step 2 — Run locally

```bash
conda activate NanoCount_5
bash JOLI_Kallisto/scripts/run_joli_kallisto.sh
```

#### Step 3 — Run on SLURM

```bash
# Option A — via wrapper (recommended):
bash JOLI_Kallisto/scripts/submit_joli_pipeline.sh

# Option B — direct sbatch:
sbatch JOLI_Kallisto/scripts/run_joli_kallisto.sh
```

---

### Pipeline C — JOLI single-sample (Python EM only, bustools output exists)

Use this when bustools TCC files already exist (e.g. from a previous run) and you only want to re-run the EM with different settings.

```bash
conda activate NanoCount_5
cd JOLI_Kallisto
python main_joli.py \
    --sample_dir  /path/to/bustools_output/ \
    --output_dir  /path/to/output/ \
    --eff_len_mode   uniform \
    --convergence_mode  kallisto \
    --em_type        plain \
    --max_em_rounds  10000 \
    --min_rounds     50
```

**Key flags:**

| Flag | Values | Notes |
|------|--------|-------|
| `--sample_dir` | path | Directory with `count.mtx`, `matrix.ec`, `transcripts.txt` |
| `--output_dir` | path | Where to write `abundance.tsv` (default: same as `--sample_dir`) |
| `--eff_len_mode` | `uniform` / `kallisto` | `uniform` = all eff_len=1.0; `kallisto` = length-based from `flens.txt` |
| `--convergence_mode` | `kallisto` / `joli` | `kallisto` = raw count threshold (matches lr-kallisto); `joli` = normalized theta (faster) |
| `--em_type` | `plain` | MAP and VI not yet wired to single-sample entry (use multi-sample for MAP) |
| `--max_em_rounds` | int | Max EM iterations (default: 10000) |
| `--min_rounds` | int | Min iterations before convergence check (default: 50) |

**Outputs:**
```
<output_dir>/
├── abundance.tsv   # per-transcript quantification (kallisto format)
└── runtime.txt     # wall-clock time
```

---

### Pipeline D — JOLI multi-sample MAP EM (Phase 2)

Shared Dirichlet prior `alpha` learned jointly across all samples via gradient descent. Per-sample `theta` estimated via MAP EM.

Two ways to run:
- **Full pipeline** (`scripts/run_multisample_joli.sh`) — runs bustools for each sample then MAP EM. Use when starting from raw reads.
- **MAP EM only** (`main_multisample_joli.py`) — skips bustools. Use when bustools TCC output already exists.

#### Option D1 — Full pipeline (reads → TCC → MAP EM)

Edit the **CONFIG section** at the top of `JOLI_Kallisto/scripts/run_multisample_joli.sh`:

| Variable | What to set |
|----------|-------------|
| `KALLISTO` / `BUSTOOLS` | Tool paths |
| `INDEX_FILE` / `T2G_FILE` | Reference files |
| `READ_TYPE` | `long` or `short` |
| `THREADS` | CPU threads |
| `OUTPUT_BASE` | Base results directory |
| `EFF_LEN_MODE` | `uniform` or `kallisto` |
| `CONVERGENCE_MODE` | `joli` (recommended) or `kallisto` |
| `MAX_EM_ROUNDS` / `MIN_EM_ROUNDS` | Inner EM iteration limits |
| `MAX_GD_ROUNDS` / `GD_LR` | Outer GD settings |
| `ALPHA_INITIAL` | Initial Dirichlet concentration |
| `GD_CONVERGENCE_TOL` / `GD_STEPS_PER_ROUND` | GD convergence settings |
| `SAMPLES` | Same array format as single-sample pipelines; minimum 2 entries |

```bash
conda activate NanoCount_5

# Run locally:
bash JOLI_Kallisto/scripts/run_multisample_joli.sh

# Run on SLURM:
bash JOLI_Kallisto/scripts/submit_multisample_joli.sh
```

Edit the CONFIG block inside `submit_multisample_joli.sh` for SLURM resources (`MEM`, `CPUS`, `TIME`, `MAIL_USER`).

#### Option D2 — MAP EM only (bustools output already exists)

Edit the **CONFIG section** at the top of `JOLI_Kallisto/main_multisample_joli.py`:

| Variable | What to set |
|----------|-------------|
| `SAMPLE_DIRS` | List of bustools output directories (Option A: explicit) |
| `SAMPLES_FOLDER` | Parent folder — auto-discovers all subdirs (Option B) |
| `RESULTS_BASE` | Base results directory |
| `READ_TYPE` | `long` or `short` |
| `EFF_LEN_MODE` | `uniform` or `kallisto` |
| `CONVERGENCE_MODE` | `joli` (recommended for MAP) or `kallisto` |
| `MAX_EM_ROUNDS` | Max inner EM rounds per sample per GD round (default: 10000) |
| `MIN_EM_ROUNDS` | Min inner EM rounds before convergence check (default: 50) |
| `MAX_GD_ROUNDS` | Max outer GD iterations (default: 500) |
| `GD_LR` | Adam learning rate for shared alpha (default: 0.01) |
| `ALPHA_INITIAL` | Initial Dirichlet concentration — uniform prior (default: 1.0) |
| `GD_CONVERGENCE_TOL` | Stop outer loop when `|loss_change| < tol` (default: 1e-6) |
| `GD_STEPS_PER_ROUND` | Adam steps per outer GD round (default: 10) |

**Two input options (set ONE, leave the other blank/empty):**
```python
# Option A — explicit list of bustools output directories
SAMPLE_DIRS = [
    "/path/to/kallisto_output/sample_1/",
    "/path/to/kallisto_output/sample_2/",
]
SAMPLES_FOLDER = ""

# Option B — auto-discover all subdirs in a parent folder
SAMPLE_DIRS    = []
SAMPLES_FOLDER = "/path/to/kallisto_output/"
```

> All samples must have been processed with the same transcriptome (same `transcripts.txt`).
> Run `scripts/run_joli_kallisto.sh` or `scripts/run_lr_kallisto.sh` first to produce bustools outputs.

```bash
conda activate NanoCount_5
cd JOLI_Kallisto

# Standalone (CONFIG provides sample dirs):
python main_multisample_joli.py

# Passing sample dirs from the command line:
python main_multisample_joli.py \
    --sample_dirs /path/to/sample1_tcc/ /path/to/sample2_tcc/ \
    --results_base /path/to/results/ \
    --convergence_mode joli \
    --max_gd_rounds 500
```

#### Outputs

```
<RESULTS_BASE>/exprmnt_{timestamp}/
├── experiment_description.log   # full CONFIG dump + sample list
├── running.log                  # combined stdout (teed from terminal)
├── runtime.txt                  # total wall-clock time
├── code_snapshot/               # all .py/.sh/.txt/.yml files (recursive)
├── alpha_final.npy              # final shared Dirichlet alpha (T,)
├── gd_loss_history.pkl          # list of per-round loss histories
└── {sample_name}/
    └── abundance.tsv            # per-transcript quantification for this sample
```

---

### Running Tests

All tests live in `JOLI_Kallisto/test/` and are run from the `JOLI_Kallisto/` root:

```bash
conda activate NanoCount_5
cd JOLI_Kallisto

python test/test_JolitoKallisto.py          # load_tcc, weights, output_writer
python test/test_em_algorithm.py            # JoliEM: plain EM + MAP EM
python test/test_dirichlet_optimizer.py     # DirichletOptimizer
python test/test_multi_sample_em.py         # MultiSampleJoliEM
python test/test_main_multisample_joli.py   # main_multisample_joli helpers
```

Each test prints `PASS` / `FAIL` per check and a summary line. No external test framework required.

---

### Per-sample Output Layout (single-sample pipelines)

```
/gpfs/.../files/results/exprmnt_{timestamp}/
├── experiment_description.log   # config dump + sample list
├── running.log                  # combined stdout/stderr for all steps
├── runtime.txt                  # total wall-clock time
├── code_snapshot/               # all .py/.sh/.txt/.yml files (recursive)
├── {sample_name_1}/
│   ├── output.bus               # raw bus file
│   ├── sorted.bus               # sorted bus file
│   ├── count.mtx                # TCC count matrix
│   ├── matrix.ec                # equivalence class definitions
│   ├── transcripts.txt          # transcript list
│   ├── abundance.tsv            # final isoform quantification
│   └── run_info.json            # kallisto run metadata
└── {sample_name_2}/
    └── ...
```

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

`code_snapshot/` mirrors the full directory structure but **only copies** these extensions:

| Extension | Examples |
|-----------|---------|
| `.py` | all Python scripts |
| `.sh` | all shell/SLURM scripts |
| `.txt` | requirements, notes |
| `.yml` / `.yaml` | conda environment files |

**Do NOT copy**: `.pkl`, `.idx`, `.bam`, `.fastq`, `.fasta`, `.bus`, `.mtx`, `.log`, or any data/binary assets.

**Implementation — bash:**
```bash
find code/ \( -name "*.py" -o -name "*.sh" -o -name "*.txt" -o -name "*.yml" -o -name "*.yaml" \) \
    | while read f; do
        dest="code_snapshot/${f#code/}"
        mkdir -p "$(dirname "$dest")"
        cp "$f" "$dest"
    done
```

**Implementation — Python** (`save_code_snapshot()` in `utility.py`):
```python
SNAPSHOT_EXTS = {".py", ".sh", ".txt", ".yml", ".yaml"}
for src in Path("code/").rglob("*"):
    if src.suffix in SNAPSHOT_EXTS:
        dest = run_dir / "code_snapshot" / src.relative_to("code/")
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
```
