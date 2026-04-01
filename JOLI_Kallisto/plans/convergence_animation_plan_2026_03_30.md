# Convergence Animation Plan
**Date:** 2026-03-30

## Goal

Visualise how isoform quantification converges toward ground truth across training
rounds, for three methods:
- **LK** — lr-kallisto (C++ EM, final theta only — no intermediate access)
- **JK** — JOLI single-sample plain EM (intermediate thetas via snapshot flag)
- **JK MS** — JOLI multi-sample MAP EM with shared Dirichlet prior (intermediate
  alpha + per-sample theta via snapshot flag)

Also compute:
1. GT vs GT inter-sample correlation (one-off comparison)
2. JK single-sample inter-sample Spearman across EM rounds (from snapshots)
3. JK MS inter-sample Spearman already tracked in `training_stats.pkl`

---

## Output: 2 GIFs (one per sample)

### `convergence_sim1.gif`
Four-panel page, sim1 perspective:

```
┌─────────────────────────────┬──────────────────────────────────────┐
│ Panel 1: GT sim1 (STATIC)   │ Panel 2: JK single-sample sim1 theta │
│ transcripts ranked by GT TPM│ same rank order, animated by round   │
├─────────────────────────────┼──────────────────────────────────────┤
│ Panel 3: JK MS alpha        │ Panel 4: JK MS theta sim1            │
│ shared alpha, animated      │ per-sample theta, animated           │
└─────────────────────────────┴──────────────────────────────────────┘
```

### `convergence_sim2.gif`
Same 4-panel layout for sim2 perspective.

---

## Transcript Universe Per Page

For each sample (sim1 or sim2):
1. **TP transcripts** — GT TPM > 0 AND predicted non-zero in any method.
   Sorted by GT TPM descending. Shown in the top portion of the strip.
2. **FP transcripts** — predicted non-zero in any method AND GT TPM = 0.
   Sorted by predicted TPM descending. Shown below a visible separator line.
3. **FN transcripts** — GT TPM > 0 but predicted zero — EXCLUDED (not shown).

The rank order (x-axis) is fixed to GT for that sample and does not change
across frames.

---

## Visual Encoding

- **Colormap:** `coolwarm` (blue = low, red = high)
- **Normalisation:** per-panel, across all frames — so round 0 colours are
  directly comparable to final-round colours
- **Panel 1 (GT):** static, shown in every frame; FP positions shown as grey
  (no GT value)
- **Frame title:** `"Round N"` (or `"Epoch N"`)
- **FPS:** configurable (default 2 fps = 0.5 s/frame)
- **Snapshot cadence:** every 5 rounds

---

## Item 1 — GT vs GT Correlation (standalone)

Simple one-off analysis in `plot_convergence_animation.py` (printed at startup):
- Load `PB_sample1_gt.tsv` and `PB_sample2_gt.tsv`
- Merge on `transcript_id`, fill missing with 0
- Compute and print Spearman + Pearson

---

## Item 2 — JK Single-Sample Inter-Sample Spearman Across Rounds

For each snapshot checkpoint (every 5 rounds):
- Load sim1 theta snapshot at round N
- Load sim2 theta snapshot at round N
- Align on shared transcript indices
- Compute Spearman(theta_sim1, theta_sim2)
- Plot: x = round, y = Spearman. Saved as `jk_inter_sample_spearman.png`
  in the JK MS experiment folder.

JK MS inter-sample Spearman is already recorded per round in
`training_stats.pkl` — just load and plot.

---

## Infrastructure Changes

### 1. `core/em_algorithm.py`

Add to `JoliEM.run()`:
```python
snapshot_interval: int = 0,   # 0 = disabled; save theta every N rounds
```

Inside the EM loop, every `snapshot_interval` rounds:
```python
if snapshot_interval > 0 and (round_num % snapshot_interval == 0):
    snapshots.append((round_num, theta.copy()))
```

Also append the final theta after the loop.

Add `snapshots` field to `EMResult`:
```python
@dataclass
class EMResult:
    alpha:     np.ndarray
    n_rounds:  int
    converged: bool
    snapshots: list = None   # list of (round_num, theta_array) or None
```

### 2. `core/training_tracker.py`

Add method:
```python
def record_snapshot(self, round_num: int, alpha: np.ndarray,
                    theta_list: list) -> None:
```

Appends `{"round": round_num, "alpha": alpha.copy(),
           "thetas": [t.copy() for t in theta_list]}` to `self.snapshots`.

Add `save_snapshots(path)` method — pickles the snapshot list.

Called from `_run_em_wrapper` every `snapshot_interval` rounds when
`save_snapshots=True`.

### 3. `core/multi_sample_em.py`

Add to `__init__`:
```python
save_snapshots:     bool = False,
snapshot_interval:  int  = 5,
```

In `_run_em_wrapper`, after `tracker.record(...)`:
```python
if self.save_snapshots and (round_num % self.snapshot_interval == 0):
    tracker.record_snapshot(round_num, self.dirichlet_opt.get_alpha(), thetas)
```

In `write_results`, after saving `training_stats.pkl`:
```python
if results.get("tracker") and self.save_snapshots:
    snap_path = os.path.join(output_dir, "snapshots.pkl")
    results["tracker"].save_snapshots(snap_path)
```

### 4. `main_joli.py`

Add to CONFIG:
```python
SAVE_SNAPSHOTS    = False   # set True to save per-round theta snapshots
SNAPSHOT_INTERVAL = 5       # save every N EM rounds (only when SAVE_SNAPSHOTS=True)
```

Pass to `JoliEM.run()`:
```python
snapshot_interval = SNAPSHOT_INTERVAL if SAVE_SNAPSHOTS else 0
```

Save snapshots to experiment folder:
```python
if SAVE_SNAPSHOTS and result.snapshots:
    snap_path = os.path.join(run_dir, sample_name, "theta_snapshots.pkl")
    with open(snap_path, "wb") as fh:
        pickle.dump(result.snapshots, fh)
```

### 5. `main_multisample_joli.py`

Add to CONFIG:
```python
SAVE_SNAPSHOTS    = False   # set True to save alpha + theta snapshots every N rounds
SNAPSHOT_INTERVAL = 5
```

Add `--save_snapshots` and `--snapshot_interval` CLI args.
Wire through to `MultiSampleJoliEM`.

### 6. `scripts/run_multisample_joli.sh`

Add to CONFIG section:
```bash
SAVE_SNAPSHOTS=false      # true = save alpha+theta snapshots every SNAPSHOT_INTERVAL rounds
SNAPSHOT_INTERVAL=5
```

Pass to Python:
```bash
--save_snapshots   "${SAVE_SNAPSHOTS}" \
--snapshot_interval "${SNAPSHOT_INTERVAL}" \
```

---

## New File: `analysis/plot_convergence_animation.py`

### CONFIG section
```python
# ============================================================
# CONFIG — edit before running
# ============================================================
JK_SINGLE_SIM1_DIR = "exprmnt_..."   # JK single-sample run for sim1
JK_SINGLE_SIM2_DIR = "exprmnt_..."   # JK single-sample run for sim2
JK_MS_DIR          = "exprmnt_..."   # JK MS run (contains snapshots.pkl)
GT_SIM1_PATH       = "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample1_gt.tsv"
GT_SIM2_PATH       = "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample2_gt.tsv"
RESULTS_BASE       = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"
FPS                = 2               # frames per second (2 = 0.5s/frame)
SAMPLE_NAMES       = ["sim1", "sim2"]
# ============================================================
```

### Key functions

```
load_gt(path)
    → pd.DataFrame[transcript_id, tpm_gt]

load_jk_single_snapshots(exp_dir, sample_name)
    → list of (round_num, theta: np.ndarray)  # indexed same as transcripts.txt

load_jk_ms_snapshots(exp_dir)
    → list of {"round", "alpha", "thetas": [theta_s1, theta_s2]}

build_transcript_index(jk_ms_dir, sample_name)
    → {transcript_id: array_index}  # from transcripts.txt

build_universe(gt_df, snapshot_thetas, transcript_index)
    → DataFrame with columns: transcript_id, gt_tpm, pred_tpm_final,
       is_tp (GT>0 and pred>0), is_fp (GT=0 and pred>0)
    → sort: TP by gt_tpm desc, then FP by pred_tpm desc
    → returns rank array and separator_index

build_frame(round_num, gt_strip, jk_single_theta, jk_ms_alpha,
            jk_ms_theta, universe, transcript_index, vmin, vmax)
    → matplotlib Figure with 4 panels (2×2), all sharing the same
      transcript rank axis

generate_gif(frames, output_path, fps)
    → saves animated GIF using PIL/imageio

main()
    1. Print GT vs GT correlation
    2. Load all snapshots
    3. Build universe for sim1 and sim2
    4. Compute per-panel vmin/vmax across all frames (for stable colormap)
    5. Render all frames for sim1 → save convergence_sim1.gif
    6. Render all frames for sim2 → save convergence_sim2.gif
    7. Plot JK single inter-sample Spearman curve → save jk_inter_sample_spearman.png
    8. Plot JK MS inter-sample Spearman curve (from training_stats.pkl) → save jkms_inter_sample_spearman.png
    All outputs saved inside JK_MS_DIR experiment folder.
```

### Panel rendering detail

Each panel is a 1D colour strip (horizontal bar, 1 row × N_transcripts columns):
- x-axis = transcript rank (fixed to GT order)
- colour = TPM/alpha value at that position
- TP region: coloured by value (coolwarm)
- FP region (below separator): coloured by predicted value; GT panel shows grey
- Colourbar on the right of each panel
- Panel title: method name + sample + round number

---

## Run instructions

```bash
# Step 1: re-run JK single-sample for sim1 + sim2 with SAVE_SNAPSHOTS=True
#         (set in main_joli.py CONFIG)
cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
conda activate NanoCount_5
python main_joli.py   # for sim1

# Step 2: re-run JK MS with SAVE_SNAPSHOTS=True in run_multisample_joli.sh

# Step 3: generate GIFs
python analysis/plot_convergence_animation.py
```

---

## File outputs (all saved inside JK_MS_DIR)

```
exprmnt_.../
├── snapshots.pkl                    # JK MS alpha + theta snapshots (every 5 rounds)
├── convergence_sim1.gif             # 4-panel animation, sim1 perspective
├── convergence_sim2.gif             # 4-panel animation, sim2 perspective
├── jk_inter_sample_spearman.png     # JK single-sample inter-sample Spearman curve
└── jkms_inter_sample_spearman.png   # JK MS inter-sample Spearman curve

exprmnt_jk_single_sim1/
└── sim1/theta_snapshots.pkl         # JK single-sample sim1 theta snapshots

exprmnt_jk_single_sim2/
└── sim2/theta_snapshots.pkl         # JK single-sample sim2 theta snapshots
```
