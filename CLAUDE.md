# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.
For detailed run instructions, script descriptions, and config references see **[experiments.md](experiments.md)**.

---

## Project Overview

Multi-sample RNA isoform quantification from long-read sequencing (PacBio/ONT). Two parallel workstreams:

1. **`AT_code/`** — Custom EM algorithm (VI or MAP) with Bayesian Dirichlet hyperparameter optimization via gradient descent.
2. **`JOLI_Kallisto/`** — Active development: merging the AT_code EM/VI approach with [lr-kallisto](https://github.com/pachterlab/kallisto) for faster quantification.

---

## Environment Setup

```bash
conda activate NanoCount_5
```

- yml/txt: `AT_code/NanoCount_5.yml`, `AT_code/NanoCount_5.txt`
- Key packages: PyTorch 2.3.1, Pyro-PPL 1.9.1, NanoCount 1.0.0, pysam 0.22.1, NumPy, SciPy, pandas, Matplotlib
- New dependencies → new environment `NanoCount_6`, saved as `AT_code/NanoCount_6.yml`.

---

## AT_code — Summary

Pipeline: BAM → pickle (`process_bam_files.py`) → EM (`main_EM_VIorMAP_GD_vector.py`) → TSV outputs → analysis (`result_analysis/`).

| File | Role |
|------|------|
| `main_EM_VIorMAP_GD_vector.py` | Primary entry point |
| `EM_VIorMAP_GD_vector.py` | Core EM logic |
| `DirichletOptimizer_vector.py` | PyTorch Dirichlet GD optimizer |
| `process_bam_files.py` | BAM → pickle |
| `generate_bash.py` | SLURM batch job generator |
| `result_analysis/generate_result_stat.py` | Spearman/Pearson vs. ground truth |
| `result_analysis/ThetaCorrWGrndTruth.py` | Theta vs. SIRV ground truth |
| `result_analysis/plt_experiment_stats.py` | EM convergence / GD loss plots |
| `simulation_code/simulation.py` | Synthetic data generation |

Experiment types (`--experiment_num`): `1`=single sample, `2`=merged, `4`=multi-sample+GD, `5`=merged-multi+GD.
Output filenames parsed via `ExperimentFileProcessor` in `result_analysis/util.py`.

> See [experiments.md](experiments.md) for full CLI flags, config tables, and step-by-step instructions.

---

## JOLI_Kallisto — Summary

Four pipelines available — choose based on what you need:

| Pipeline | Entry point | When to use |
|----------|-------------|-------------|
| lr-kallisto baseline | `scripts/run_lr_kallisto.sh` | C++ EM, comparison baseline |
| JOLI full (single-sample) | `scripts/run_joli_kallisto.sh` | bustools + Python EM, end-to-end |
| JOLI EM only (single-sample) | `main_joli.py` | bustools output already exists |
| JOLI multi-sample full pipeline | `scripts/run_multisample_joli.sh` | Phase 2: bustools + joint MAP EM |
| JOLI multi-sample MAP EM only | `main_multisample_joli.py` | Phase 2: bustools output already exists |

### Folder structure

```
JOLI_Kallisto/
├── core/                        # Library modules (not run directly)
│   ├── load_tcc.py              #   bustools output → TCCData
│   ├── weights.py               #   effective lengths + EC weights
│   ├── em_algorithm.py          #   JoliEM: plain EM + MAP EM
│   ├── output_writer.py         #   write abundance.tsv
│   ├── dirichlet_optimizer.py   #   Adam optimizer for shared alpha
│   ├── multi_sample_em.py       #   MultiSampleJoliEM: GD + per-sample MAP EM
│   └── training_tracker.py      #   TrainingTracker: per-round metric collection
├── scripts/                     # Executable pipeline scripts
│   ├── run_joli_kallisto.sh         #   full single-sample pipeline
│   ├── run_lr_kallisto.sh           #   lr-kallisto baseline pipeline
│   ├── run_multisample_joli.sh      #   full multi-sample pipeline (bustools + MAP EM)
│   ├── submit_joli_pipeline.sh      #   SLURM wrapper for run_joli_kallisto.sh
│   ├── submit_lr_kallisto.sh        #   SLURM wrapper for run_lr_kallisto.sh
│   └── submit_multisample_joli.sh   #   SLURM wrapper for run_multisample_joli.sh
├── analysis/                    # Post-run analysis and plotting
│   ├── compare_abundance_files.py   # Compare JK vs LK abundance.tsv
│   ├── run_batch_comparison.py      # Batch comparison across experiment folders
│   └── plot_training.py             # Training diagnostic figures from training_stats.pkl
├── test/                        # Unit tests (run from JOLI_Kallisto/ root)
│   ├── test_JolitoKallisto.py
│   ├── test_em_algorithm.py
│   ├── test_dirichlet_optimizer.py
│   ├── test_multi_sample_em.py
│   ├── test_main_multisample_joli.py
│   ├── test_training_tracker.py
│   └── test_helpers.py
├── plans/                       # Design documents
├── main_joli.py                 # Single-sample JOLI EM entry point
├── main_multisample_joli.py     # Multi-sample MAP EM entry point (Phase 2)
└── main_pipeline.py             # Full pipeline entry point
```

### Key implementation notes

- `core/` modules are imported by entry points via `sys.path.insert(0, .../core)` — no `__init__.py` needed.
- Test files import from `../core` (not `..`); `test_main_multisample_joli.py` imports from `..` (root) for `main_multisample_joli`.
- `convergence_mode="kallisto"` matches lr-kallisto exactly (raw count threshold + zeroing on raw alpha); `"joli"` uses normalized theta (faster, for MAP).
- Multi-sample MAP EM: posterior mean M-step `theta = (n + alpha) / sum(n + alpha)`; shared `alpha` updated via Adam on `log_alpha`.
- `main_multisample_joli.py` accepts CLI args (`--sample_dirs`, `--results_base`, all EM/GD params) that override CONFIG; used by `run_multisample_joli.sh` to pass cache dirs dynamically.
- `TrainingTracker` records per-round metrics (GD loss, EM rounds, Spearman/Pearson inter-sample + theta-vs-alpha, alpha sum/entropy/change, nonzero transcripts); saved as `training_stats.pkl`; figures auto-generated to `figures/` via `plot_training.py`.
- `save_code_snapshot()` uses `rglob` — captures files from all subdirs including `core/`.

lr-kallisto base: `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/`

> See [EXPERIMENTS.md](../EXPERIMENTS.md) for full CONFIG variables, run commands, CLI flags, and output layouts.

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

All results save to `/gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmnt_{timestamp}/` containing:
`experiment_description.log`, `running.log`, `best_checkpoint/`, `results_summary.txt`, `runtime.txt`, `code_snapshot/`.

`code_snapshot/` copies only `.py`, `.sh`, `.txt`, `.yml`, `.yaml` files — no data/binary assets.

> See [experiments.md](experiments.md) for the full output tree and snapshot implementation.

---

## Utility Helpers

`utility.py` (add to `AT_code/` and `JOLI_Kallisto/` as development progresses) must provide:
- `create_run_dir()` — create timestamped output folder
- `save_runtime(run_dir, elapsed_seconds)` — write `runtime.txt`
- `save_code_snapshot(run_dir)` — copy code into `code_snapshot/`

---

## CRITICAL: Off-Limits Directories

**NEVER read from, write to, modify, create, or delete any files in:**

```
/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation
```

This directory belongs to collaborator Shree Raghavendra and is **read-only reference material**. Treat it as untouchable — no edits, no new files, no deletions, no moves, ever.

---

## Simulations (`code/Simulations/`)

A clean re-implementation of a 3-phase RNA-seq read simulation pipeline that generates synthetic reads for **Illumina**, **PacBio**, and **ONT** technologies. Based on the original pipeline in `data/Shree_stuff/Simulation/` (do not touch that folder).

> Full reference documentation: [`code/Simulations/REFERENCE.md`](Simulations/REFERENCE.md)
