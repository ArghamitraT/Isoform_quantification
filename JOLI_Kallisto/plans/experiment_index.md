# JK Multi-Sample Experiment Index

Spearman columns show **all / active / GT-nonzero** universes (3-universe breakdown from `run_gt_comparison.py`).  
"active" = GT ∪ nonzero-pred; "GT-nonzero" = GT non-zero only.  
Alpha scale = max value of `alpha_final.npy` — healthy runs stay in the hundreds; exploding runs reach millions+.

---

## 2026-05-20 — Single-TX flag ablation

**Question:** Does including single-tx reads in the EM M-step or GD theta help or hurt?

| Folder | loop | em_stx | gd_stx | Spearman all/active/gt (sim1) | Spearman all/active/gt (sim2) | alpha max | Snapshots | Verdict |
|--------|------|--------|--------|-------------------------------|-------------------------------|-----------|-----------|---------|
| `exprmnt_2026_05_20__11_47_51` | em_wrapper | ✓ | ✗ | 0.086 / -0.430 / -0.106 | 0.072 / -0.405 / -0.106 | 194M | ✗ | ✗ RUNAWAY — em_stx=True destabilizes em_wrapper; feedback loop drives alpha to 194M |
| `exprmnt_2026_05_20__12_02_51` | gd_wrapper | ✓ | ✗ | 0.385 / -0.079 / 0.351 | 0.383 / -0.030 / 0.361 | 7M | ✗ | ✗ POOR — gd_wrapper more stable but alpha still explodes to 7M; active Spearman negative |
| `exprmnt_2026_05_20__12_07_44` | em_wrapper | ✗ | ✓ | 0.900 / 0.801 / 0.837 | 0.900 / 0.804 / 0.842 | 4M | ✗ | ✓ GOOD — old behavior restored; Spearman ~0.9 across both samples; alpha still large |
| `exprmnt_2026_05_20__13_04_38` | gd_wrapper | ✗ | ✓ | 0.903 / 0.833 / 0.814 | 0.912 / 0.830 / 0.817 | 17K | ✗ | ✓ BEST — old behavior + gd_wrapper; highest Spearman, alpha healthy (17K vs 4M); **no snapshots** (bug fixed 2026-05-21) |

**Key finding:** `em_include_single_tx=True` causes runaway feedback in both loop modes.  
`gd_include_single_tx=True` (old behavior) is stable. `gd_wrapper` keeps alpha ~17K vs em_wrapper's ~4M.

---

## Reference experiments (earlier runs)

| Folder | loop | em_stx | gd_stx | Notes | Snapshots |
|--------|------|--------|--------|-------|-----------|
| `exprmnt_2026_03_30__22_37_55` | gd_wrapper | ✗ | ✓ | Old reference MS run | ✗ |
| `exprmnt_2026_05_18__12_46_25` | — | — | — | **JK single-sample** plain EM (use as SS baseline) | ✓ sim1+sim2 |
| `exprmnt_2026_05_18__12_42_22` | — | — | — | lr-kallisto baseline (NOT a JK run — abundance.tsv has 2 cols only) | — |

---

## Quick lookup

| I want… | Use this folder |
|---------|-----------------|
| Best JK MS result (needs re-run for snapshots) | `exprmnt_2026_05_20__13_04_38` |
| JK SS baseline with snapshots | `exprmnt_2026_05_18__12_46_25` |
| Old MS reference (pre-May ablation) | `exprmnt_2026_03_30__22_37_55` |
| lr-kallisto baseline | `exprmnt_2026_05_18__12_42_22` |
| Next MS run (snapshots fixed, re-run needed) | run `scripts/run_multisample_joli.sh` |
