# Project Handover — RNA Isoform Quantification (JOLI / JOLI_Kallisto)

**Purpose of this file:** a single entry point for the next developer (and for a fresh Claude
session) to understand where everything lives — code, data, results, environments, external
docs — and what state the project is in. All paths below were verified to exist on the NYGC
cluster (`pe2-login`) on **2026-10-27**. Where a path could not be verified it is marked
`[UNVERIFIED — from journal, not re-checked]`.

If you are a new Claude session reading this: read this file first, then
`/gpfs/commons/home/atalukder/RNA_Splicing/code/CLAUDE.md` (hard rules, folder structure) and
`/gpfs/commons/home/atalukder/RNA_Splicing/code/EXPERIMENTS.md` (CLI flags / run commands) before
touching anything.

---

## 1. One-paragraph project summary

The project quantifies RNA transcript-isoform abundance (θ) from long-read (PacBio/ONT) and
short-read (Illumina) RNA-seq, with two goals: (a) combine long- and short-read evidence, and
(b) share statistical strength **across multiple samples** (replicates / days / conditions) via a
Bayesian Dirichlet prior on the per-sample isoform-abundance vectors, optimized jointly with the
per-sample EM by gradient descent (MAP or VI). The current engineering effort (`JOLI_Kallisto/`)
re-implements this on top of `lr-kallisto`'s pseudoalignment/equivalence-class machinery instead
of `NanoCount`-style BAM alignment, for speed. The project has a ~2.5-year running research
journal (Jan 2024 → Jul 2026, pasted into this chat session in full) documenting the math,
literature review, and every experiment; this file is a condensed map + current-status summary,
**not** a replacement for that journal.

---

## 2. Repo layout (code)

Repo root: `/gpfs/commons/home/atalukder/RNA_Splicing/code`
GitHub: `https://github.com/ArghamitraT/Isoform_quantification` (branch `main`; the remote URL
in `git remote -v` has a personal access token embedded — do not paste that full URL anywhere
public).

```
/gpfs/commons/home/atalukder/RNA_Splicing/code/
├── CLAUDE.md                     # hard rules: env, folder conventions, off-limits dirs, output format
├── EXPERIMENTS.md                # full CLI flags / config tables / run commands for AT_code & JOLI_Kallisto
├── README.md
├── HANDOVER.md                   # this file
├── git_push.sh
├── download_Gdrive.sh
├── encode_2026_downloader.py     # newer (2026-07) ENCODE bulk-download workstream
├── encode_2026_experiments.csv
├── ENCODE_2026_README.md
│
├── AT_code/                      # ORIGINAL/legacy workstream: custom EM (MAP/VI) + Dirichlet GD on BAM-derived pickles
│   ├── main_EM_VIorMAP_GD_vector.py   # primary entry point (current/clean version)
│   ├── EM_VIorMAP_GD_vector.py        # core EM logic
│   ├── DirichletOptimizer_vector.py   # PyTorch Dirichlet GD optimizer
│   ├── process_bam_files.py           # BAM -> pickle
│   ├── generate_bash.py               # SLURM batch job generator
│   ├── NanoCount_5.yml / .txt         # conda env spec for this workstream
│   ├── result_analysis/               # Spearman/Pearson vs ground truth, EM/GD convergence plots
│   ├── simulation_code/simulation.py  # older synthetic-data generator (superseded by code/Simulations/)
│   ├── newly_developed/               # in-progress rewrite (JOLI.py, main.py) — precursor to JOLI_Kallisto
│   ├── extra_codes/                   # graveyard of old/experimental variants — do not build on these
│   ├── pkl_files/, results/, figures/ # local scratch data/output, mostly stale
│   └── .vscode/launch.json
│
├── JOLI_Kallisto/                # ACTIVE / CURRENT workstream — JOLI (custom MAP/VI EM) + lr-kallisto
│   ├── core/                     # library modules (imported via sys.path.insert, no __init__.py)
│   │   ├── load_tcc.py           #   bustools output -> TCCData
│   │   ├── weights.py            #   effective lengths + EC weights (incl. Python-side eff-len calc for SR)
│   │   ├── em_algorithm.py       #   JoliEM: plain EM + MAP EM; em_include_single_tx flag
│   │   ├── multi_sample_em.py    #   MultiSampleJoliEM: GD + per-sample MAP EM; gd_include_single_tx flag,
│   │   │                         #   loop_mode = "em_wrapper" | "gd_wrapper"
│   │   ├── dirichlet_optimizer.py#   Adam optimizer for shared alpha (log-alpha space)
│   │   ├── output_writer.py      #   write abundance.tsv
│   │   └── training_tracker.py   #   per-round metrics -> training_stats.pkl
│   ├── scripts/                  # executable pipelines
│   │   ├── run_lr_kallisto.sh        # baseline: lr-kallisto's own C++ EM
│   │   ├── run_joli_kallisto.sh      # single-sample: bustools + Python JOLI EM
│   │   ├── run_multisample_joli.sh   # multi-sample: bustools + joint MAP EM across samples
│   │   ├── submit_*_slurm / submit_JKMS_slurm.py
│   │   └── mtx_to_tsv.py, view_html.sh
│   ├── analysis/                 # post-run analysis/plots (run_gt_comparison.py, run_replicate_comparison.py,
│   │                              # run_eda.py, plot_training.py, plot_metric_curves.py, exp1_tpm_stealing.py, …)
│   ├── debug/                    # ad hoc debugging scripts (micro EM tests, LK eff-len checks)
│   ├── test/                     # unit tests (run from JOLI_Kallisto/ root)
│   ├── plans/                    # dated design docs — read these for design rationale on specific features
│   ├── main_joli.py              # single-sample JOLI EM entry point
│   ├── main_multisample_joli.py  # multi-sample MAP EM entry point
│   └── main_pipeline.py          # full pipeline entry point
│
├── Simulations/                  # clean re-implementation of the 3-phase (Illumina/PacBio/ONT) read simulator
│   ├── src/, tests/, commands/, plans/, REFERENCE.md, run.sh
│   └── (independent of data/Shree_stuff/Simulation — see §6 OFF-LIMITS)
│
├── NanoCount_codes/               # vendored/adapted copy of the original NanoCount tool (reference only)
├── SIRV_data/                     # pickled SIRV E0/E2 alignment + ref-length dicts (small, in-repo)
└── .claude/skills/view-html/      # project-local Claude skill
```

---

## 3. Conda environments

| Env name | Used by | Spec file | Notes |
|---|---|---|---|
| `NanoCount_5` | `AT_code/` | `AT_code/NanoCount_5.yml`, `AT_code/NanoCount_5.txt` | PyTorch 2.3.1, Pyro-PPL 1.9.1, NanoCount 1.0.0, pysam 0.22.1 |
| `Joli_kallisto` | `JOLI_Kallisto/` | **no yml committed in repo** — GAP, see below | activated as `conda activate Joli_kallisto` throughout the 2026 journal entries |
| `mpaqt` | Shree's simulation pipeline (`lrgasp-simulation`) | none in this repo | used in `run_sim.sh` for `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/AT_tried/simulation/` |
| `NanoCount_4`, `NanoCount_5_shree`, `NanoCount_trial` | one-off exploratory work (VS Code debugging, seqtk downsampling, bustools GLIBCXX fix) | none | historical, not needed for current work |

**Gap to fix:** `Joli_kallisto` env has no exported yml in the repo (violates the global CLAUDE.md
"store all environment ymls in `env/`" rule). Next step: `conda env export -n Joli_kallisto >
JOLI_Kallisto/env/Joli_kallisto.yml` and commit it.

---

## 4. Data locations (all verified to exist, 2026-10-27)

| Data | Path | Notes |
|---|---|---|
| Real PacBio (WTC-11) alignments, pre-pickled | `/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/` | subfolders: `transcriptome_aln/` (current/used), `transcriptome_aln_pklfiles/` (parsed read/ref-len dicts, incl. `day5_our_parsing` and `day5_NanoCount_parsing`), `spliced_STAR_aln/`, `unspliced_aln/` (old), `transc_sr_minimap/` |
| Simulated data (sim_real_data, RSEM/IsoSeq-based) | `/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/` **== identical to** `/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/` | both paths resolve to the same files (verified same contents/sizes) — use either; contains `ds_100_num1_aln_01_long.bam`, matching `.fasta`/`.fq` pairs, `dataStatlong/`, per-sample ground truth |
| SIRV pickled alignments (in-repo) | `/gpfs/commons/home/atalukder/RNA_Splicing/code/SIRV_data/` | `aln_E0_read_dicts.pkl`, `aln_E0_ref_len_dicts.pkl`, `aln_E2_read_dicts.pkl`, `aln_E2_ref_len_dicts.pkl` |
| SIRV raw BAMs (per journal) | `/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/SIRV/aln_E0.bam`, `aln_E2.bam` | `[UNVERIFIED — from journal, not re-checked]` |
| Collaborator Shree's SOTA tool installs + AT's adaptations | `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/` | subfolders: `lr-kallisto/`, `oarfish/`, `MPAQT/`, `salmon/`, `NanoCount/`, `IsoQuant/`, `h1_cell/`, `output_sim1/`, `generate_stat.py` |
| **lr-kallisto install + indices** (referenced throughout `JOLI_Kallisto/`) | `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/` | `kallisto/` (built binary at `kallisto/build/src/kallisto`), `bustools/`, index files `human_k-31.idx`, `human_k-63.idx`, `new_index.idx` (k=63, long-read), `new_index_k31.idx` (k=31, short-read), `human.t2g`, `cleaned_genomic.gtf`, `mtrx_to_csv.py` |
| Collaborator Shree's simulation pipeline source (reference only) | `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation/` | **🚫 OFF-LIMITS — see §6, do not read/write/modify** |
| AT's own adapted copies of SOTA/simulation tooling | `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/AT_tried/` | kallisto/bustools builds, `simulation/lrgasp-simulation/` (uses `mpaqt` env), `kallisto_stuff/run_pipeline_lr_kallisto.sh` |
| PacBio downsampled long-read fastq (used for toy/debug runs) | `/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq/PacBio/reads/long/downsampled/` | `[UNVERIFIED — from journal, not re-checked]` |
| **All experiment results** | `/gpfs/commons/home/atalukder/RNA_Splicing/files/results/` | 172 `exprmnt_YYYY_MM_DD__HH_MM_SS/` folders as of 2026-10-27 (verified); each contains `experiment_description.log`, `running.log`, `runtime.txt`, `code_snapshot/`, `figures/`, per-sample subfolders (e.g. `sim1/`, `sim2/`), and analysis outputs like `gt_comparison_*.txt` saved in-place per the CLAUDE.md "Analysis Output Convention" |
| SLURM job scripts submitted historically | `/gpfs/commons/home/atalukder/RNA_Splicing/files/cluster_job_submission_files/` | includes `JK_multisample/` subfolder and many `prg_main_<jobid>__<timestamp>.sh` scripts |
| Other `files/` subfolders | `/gpfs/commons/home/atalukder/RNA_Splicing/files/` | `bash_scripts/`, `environments/` (some historical env ymls may be here — check before recreating `Joli_kallisto` env from scratch), `presentations/`, `shareWJiayu/`, `shree_files/` |

---

## 5. External links (not re-verified — copied as given by the user/journal)

| What | Link |
|---|---|
| GitHub repo | `https://github.com/ArghamitraT/Isoform_quantification` |
| Overleaf (math writeup, edit) | `https://www.overleaf.com/2434624477tsdjrfmqymtc#477980` |
| Overleaf (math writeup, view-only) | `https://www.overleaf.com/read/vnrqsbyvfvxy#c010ad` |
| Overleaf project (referenced 2024-03-14, full MAP math derivation) | `https://www.overleaf.com/project/65b45b7b0649814843edee73` |
| Google Drive — main folder | `https://drive.google.com/drive/folders/1NEqOjM-do-T_odmlrJdiPyj8jeDffhK-` |
| Google Drive — all files | `https://drive.google.com/drive/folders/1UyRmBi_31EsnND7LGEjPoalNlJXiHDZS` |
| Google Drive — code backups | `https://drive.google.com/drive/folders/1UyEXXWl4snTqKUNvetedh_4Y5ur6BpDP` |
| Google Drive — data | `https://drive.google.com/drive/folders/1V16YHrr0lzJ0bERIq8ZAwwNwb8b9bXcB` |
| Google Drive — results | `https://drive.google.com/drive/folders/1kyRUbGa1teLagO8eyCOI5_woVygKrAjt` |
| Google Drive — conda env ymls | `https://drive.google.com/drive/folders/1RKbnfxaeSBfoRIYeaW9614qux-fsKvg9` |
| Result summary spreadsheet | `https://docs.google.com/spreadsheets/d/1YvvQWoJnev5OH_8skUx1JzW0BmroR3I6sC0nIvDORgM` |
| Shree Raghavendra's running journal | `https://docs.google.com/document/d/1VY7W-cthWH5iSjRKncdc4ZEyMAI9ZaPefyQ2pnkH0mI` |
| VI/EM math explainer doc (Hensman paper walkthrough) | `https://docs.google.com/document/d/1EP49BPa4FgvuKZONl_4kR-vkayI9Qfs_dJ1-POQKPRA` |
| Math-to-code implementation demo slides | `https://docs.google.com/presentation/d/1F2cUSfQLUgmVGqpf91Yu1W6OLFPK1B0g` |
| Hensman generative-model slides (key reference, cited repeatedly) | `https://docs.google.com/presentation/d/1NFozz5EM0oUFCC2V8O_yrrKEORrfw63x` |

**⚠️ Argha's own day-by-day research journal** (the ~3000-line text pasted into this chat session,
covering 2024-01-15 → 2026-07-18, with its own "Index" / "Imp Links" header) is the authoritative
narrative record of this project but **its Google Doc URL was not found in this session** — it was
only pasted as raw text. Ask Argha for the doc link and add it here. Until then, treat §7 below as
a lossy summary only.
https://docs.google.com/document/d/1EP49BPa4FgvuKZONl_4kR-vkayI9Qfs_dJ1-POQKPRA/edit?usp=sharing


---

## 6. 🚫 Off-limits directory (do not touch, per project CLAUDE.md)

```
/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation
```
Belongs to collaborator Shree Raghavendra Thavarekere — read-only reference material. No reads,
writes, edits, or deletions, ever. (The parallel, actively-developed reimplementation lives in
`code/Simulations/` instead.)

---

## 7. Condensed project narrative (full detail is in the journal, see §5)

### 7.1 Timeline of workstreams
1. **Spring–Summer 2024**: literature review (LIQA, NanoCount, ESPRESSO, FLAIR, Salmon, Kallisto,
   Scasa, XAEM, StringTie2, MPAQT); derived the generative model, EM, then full variational
   inference (mean-field, ELBO, CAVI) for isoform quantification following Hensman et al. 2015.
2. **Mid/late 2024**: `AT_code/` built — PyTorch-based EM with Dirichlet-prior MAP/VI, vectorized
   for speed (50min → 18min/iteration). Ran on SIRV spike-ins, real PacBio+Illumina ENCODE data,
   and custom simulations. Core finding across many experiments: multi-sample sharing (MS) helps
   recall/correlation vs single-sample, provided the Dirichlet `alpha` prior strength is tuned —
   too weak (small alpha) under-shares, too strong (huge alpha) forces replicate thetas to be
   identical (overshare), both hurting real performance.
3. **2025 (Recomb-seq paper)**: finalized MAP-vs-VI comparisons, real-data and simulated-data
   results tables (see journal "Recomb'25-seq results" section) establishing that MS (multi-sample
   MAP EM with gradient-descent-tuned Dirichlet alpha) beats single-sample on Spearman correlation
   consistently, at a modest false-positive cost.
4. **2026 — current phase**: rebuilding the pipeline on top of `lr-kallisto`'s C++
   pseudoalignment/EC machinery (`JOLI_Kallisto/`) instead of NanoCount/BAM parsing, nicknamed
   **JOLI** (JOLI+Kallisto). Goal: reproduce lr-kallisto's single-sample EM exactly in Python, then
   extend it with the same multi-sample Dirichlet-MAP machinery from `AT_code/` (**JK MS**).

### 7.2 Current architecture (JOLI_Kallisto)
- Pipeline: `kallisto bus` (pseudoalignment → equivalence classes) → `bustools sort/count` →
  `core/load_tcc.py` (parse TCC matrix) → `core/em_algorithm.py` (`JoliEM`, single-sample) or
  `core/multi_sample_em.py` (`MultiSampleJoliEM`, joint MAP EM + Adam-optimized shared Dirichlet
  `alpha` across samples) → `core/output_writer.py` (`abundance.tsv`).
- `JoliEM` reproduces lr-kallisto's own EM (`convergence_mode="kallisto"`) as a validated
  baseline; `convergence_mode="joli"` is the faster normalized-theta variant used for MAP/MS runs.
- Two independent boolean flags, added 2026-05-20 after a design bug was found:
  - `em_include_single_tx` — whether unambiguous (single-transcript EC) reads are folded into
    every EM M-step (affects which transcripts "win" ambiguous reads during iteration).
  - `gd_include_single_tx` — whether those same unambiguous counts are included in the theta
    fed to the Dirichlet-alpha gradient-descent step.
  - **Current best-known config for long reads**: `em_include_single_tx=False`,
    `gd_include_single_tx=True`, `loop_mode="gd_wrapper"` (or `em_wrapper`, both validated to match
    pre-fix JK MS results). Using `gd_include_single_tx=False` or running GD on non-converged
    multi-tx-only theta causes catastrophic Dirichlet-alpha explosion (observed up to α_max≈194M)
    — **do not reuse that combination without re-deriving the fix**.
  - **For short reads**, kallisto's own behavior includes single-tx counts in every EM round, so
    `em_include_single_tx` should auto-detect short vs long and default to `True` for short reads
    (`gd_wrapper` config above was tuned for long reads only — short-read tuning is less settled,
    see §7.3).
- `loop_mode`: `"em_wrapper"` (1 EM step → N GD steps → repeat) vs `"gd_wrapper"` (EM to full
  convergence → N GD steps → repeat). **`gd_wrapper` is the validated choice for metric-curve /
  snapshot analysis** — `em_wrapper` snapshots are taken mid-EM-convergence and are not reliable
  for tracking progress over GD rounds.
- Effective-length (`eff_len`) handling for short reads was buggy (relied on a separate
  `kallisto quant` call); fixed 2026-06-23 by computing eff_len directly in Python from the
  `kallisto bus --paired` fragment-length-distribution histogram (see `core/weights.py` and
  `plans/sr_efflen_python_plan_2026_06_23.md`).

### 7.3 Open problems / where to pick up (as of the latest journal entry, 2026-07-18)
1. **Evaluation-universe ambiguity**: three ways to compute GT-comparison metrics exist
   (`run_gt_comparison.py`): "All" (outer join incl. true negatives), "Active" (GT ∪ nonzero
   predictions, excludes TN), "GT non-zero only". Which is the fairest comparison across tools is
   **still an open question** — flagged explicitly in the journal ("Which one is more valid?
   Concord, all active or positive").
2. **False positives from the Dirichlet prior**: MS recovers substantially more true positives
   (lower FN) than single-sample, but at the cost of more FPs because the prior keeps
   zero-evidence transcripts' theta from decaying fully to zero. Investigated via "TPM stealing"
   experiments (`analysis/exp1_tpm_stealing.py`) — conclusion so far: **net benefit is positive**
   (MS wins on TP/FN/correlation; FP cost is small, ~0.4% of total TPM) but the question of whether
   stolen TPM concentrates on a few high-expression TPs (vs. spread thinly) was never finished.
3. **Short-read tuning is behind long-read tuning** — the `em_include_single_tx`/
   `gd_include_single_tx`/`loop_mode` combination validated for long reads has not been fully
   re-validated for short reads after the eff-len fix; last entries (`exprmnt_2026_06_23__*`) were
   mid-debugging this.
4. **Planned next phases** (explicitly listed in journal, "Future plan" / "TODO" sections, not yet
   started as of 2026-07-18):
   - CAVI (variational, not just MAP) implementation for the multi-sample model.
   - Linear/hierarchical model on top of per-day thetas (e.g. lasso to weight days/conditions).
   - Single-cell / multi-tissue extension.
   - More real datasets — matched long+short ENCODE experiments (a matrix query URL is saved in
     the journal under 07/15/2026), Megan Schertzer's cancer-cell dataset suggestion, lowly-
     expressed-isoform analysis tied to the "Perplexity as a Metric for Isoform Diversity" paper.
   - Differential isoform analysis.
   - Understanding what the learned Dirichlet `alpha` actually encodes biologically.
5. **Known collaborators/contacts** relevant to context: David Knowles (PI/advisor,
   `daknowles@nygenome.org`), Shree Raghavendra Thavarekere (built the SOTA-tool installs and
   original simulation pipeline), Megan Schertzer (`mschertzer@nygenome.org`, data/direction
   input), Liz Tseng (PacBio, WTC-11 + Kinnex data provider), Stella (ENCODE PacBio data contact),
   Gloria Kreitinger, Karin Isaev, Wayne, Tanya (kallisto/lr-kallisto internals discussions).

---

## 8. How to resume work (quick checklist for the next developer)

1. `ssh atalukder@pe2-login` (via Cisco AnyConnect VPN, `vpn.nygenome.org`), then request an
   interactive node (`srun -p pe2 --mem 128G --cpus-per-task=4 --pty -u bash`), per
   `CLAUDE.md`/journal "How to run" sections.
2. `cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto && conda activate Joli_kallisto`
   (if the env doesn't exist, check `/gpfs/commons/home/atalukder/RNA_Splicing/files/environments/`
   or ask Argha for an export — see §3 gap).
3. Read `JOLI_Kallisto/plans/` in date order to see design rationale for the single-tx flag split,
   the eff-len fix, and the multi-sample wrapper modes before changing any of that logic.
4. Look at the most recent `exprmnt_*` folders under
   `/gpfs/commons/home/atalukder/RNA_Splicing/files/results/` (sort by timestamp) and their
   `experiment_description.log` + `gt_comparison_*.txt` to see the latest validated numbers before
   re-running anything.
5. Per this repo's `CLAUDE.md`: never touch
   `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation`; always save new
   experiment outputs under `files/results/exprmnt_{timestamp}/` with the mandated file set
   (`experiment_description.log`, `running.log`, `runtime.txt`, `code_snapshot/`, `figures/`).
