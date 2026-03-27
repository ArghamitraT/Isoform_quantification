# JOLI-Kallisto: Why JK Has More Isoforms Than Kallisto
**Date:** 2026-03-26
**Last updated:** 2026-03-26 (root cause confirmed from kallisto source)

---

## Current State

Two pipelines run on the same sample (`ds_52_furtherDownsampled`) with identical upstream
steps (kallisto bus → bustools sort → bustools count), producing identical TCC input:

| Pipeline | Experiment |
|----------|-----------|
| LK (lr-kallisto quant-tcc EM) | `exprmnt_2026_03_14__10_43_29` |
| JK (JOLI EM, uniform eff_lens, Fix A) | `exprmnt_2026_03_14__22_51_13` |

---

## What We Know

```
LK non-zero transcripts : 48,659
JK non-zero transcripts : 51,881
Non-zero in both        : 48,659
Non-zero ONLY in JK     : 3,222  ← mystery transcripts
Non-zero ONLY in LK     : 0
```

**Key facts about the 3,222 JK-unique transcripts:**
- All are **Cause B** (multi-tx EM only) — zero from Cause A or Cause C
- `matrix.ec` = `count.ec.txt` (29,851 ECs, identical) → no EC file mismatch
- `est_counts` range: 0.0007 → 12.17 (all fractional, from multi-tx EM)
- `tpm` range: 0.004 → 61.8

**Note on Jaccard = 0.223:** File format artifact.
- LK's `abundance.tsv`: 48,659 rows (non-zero only, filtered by `mtx_to_tsv.py`)
- JK's `abundance.tsv`: 218,255 rows (all transcripts including zeros)
- True Jaccard over non-zero sets = 48,659 / 51,881 = **0.938**

**From `check_lk_efflens.py`:**
- Both LK and JK use **uniform eff_lens = 1.0** (CV = 0.000002, ratio perfectly constant)
- `flens.txt` has real lengths (mean 3,591 bp) but **kallisto does NOT use them in --long mode**
- This explains why Fix B degraded correlation (0.993→0.770): JK got real lengths, LK stayed uniform

---

## ROOT CAUSE — CONFIRMED from `EMAlgorithm.h`

### Convergence criterion operates on different scales

**Kallisto long-read branch (`EMAlgorithm.h` line 302):**
```cpp
if (next_alpha[ec] > alpha_change_limit && ...)
//    ^^^^^^^^^^^ RAW EXPECTED COUNTS (unnormalized)
//    alpha_change_limit = 0.01  →  monitors transcripts with > 0.01 reads
```

**JK (`em_algorithm.py` line 289):**
```python
(theta_new > ALPHA_CHANGE_LIMIT) &
#  ^^^^^^^^^ NORMALIZED fraction (sums to 1)
#  ALPHA_CHANGE_LIMIT = 0.01  →  monitors transcripts with > 1% of total reads
#                                ≈ > ~1,969 reads  (with total_multi_reads ≈ 196,905)
```

**The threshold is ~196,905× too permissive in JK.**

All 3,222 mystery transcripts have est_counts 0.0007→12.17, so theta = 3.5e-9→6.2e-5.
All have `theta < 0.01`, so **JK never monitors them for convergence**. JK declares "done"
while these transcripts are still slowly drifting toward zero. Kallisto continues running
until transcripts with > 0.01 raw counts stabilize, at which point they have converged to 0
and are zeroed.

### Secondary difference: `finalRound` mechanism (kallisto lines 213-221)

When kallisto detects convergence, it does NOT stop immediately:
1. Sets `finalRound = true`
2. Applies zeroing (`alpha < 1e-8 → 0`) **immediately**
3. Runs **one more full EM round** with the zeroed alpha
4. Only then exits

This extra round lets remaining reads redistribute to surviving transcripts.
JK zeroes once post-loop and stops — no redistribution round.

---

## The Fix

### Fix 1 — Correct convergence criterion scale in `em_algorithm.py`

Change the convergence check from normalized theta to raw expected counts:

**Current (WRONG scale):**
```python
changed = int(np.sum(
    (theta_new > ALPHA_CHANGE_LIMIT) &                            # ← fraction threshold
    (np.abs(theta_new - theta) / np.maximum(theta_new, TOLERANCE) > ALPHA_CHANGE)
))
```

**Fixed (matches kallisto — raw count threshold):**
```python
alpha_new = theta_new * total_multi_reads    # convert to raw counts
alpha_old = theta     * total_multi_reads
changed = int(np.sum(
    (alpha_new > ALPHA_CHANGE_LIMIT) &                            # ← raw count threshold
    (np.abs(alpha_new - alpha_old) / np.maximum(alpha_new, TOLERANCE) > ALPHA_CHANGE)
))
```

With `total_multi_reads` available in `run()`, this is a small, targeted change.
The threshold remains `ALPHA_CHANGE_LIMIT = 1e-2` (matching kallisto exactly).

### Fix 2 — Add `finalRound` mechanism (one extra round after zeroing)

After convergence + zeroing, run one additional EM step to let reads redistribute:

```python
# After zeroing theta:
theta[theta < ALPHA_LIMIT / 10] = 0.0

# One extra round with zeroed theta (matches kallisto finalRound)
theta_final = self._em_step(theta)
total_final = theta_final.sum()
if total_final > 0:
    theta = theta_final / total_final
    theta[theta < ALPHA_LIMIT / 10] = 0.0  # re-zero after redistribution
```

---

## Implementation Plan

```
Step 1: Fix em_algorithm.py — convergence criterion scale (Fix 1)
        → change (theta_new > ALPHA_CHANGE_LIMIT) to (alpha_new > ALPHA_CHANGE_LIMIT)
        → requires passing total_multi_reads into the convergence check
        → estimated impact: most of the 3,222 → zero

Step 2: Fix em_algorithm.py — finalRound mechanism (Fix 2)
        → add one extra EM round after zeroing
        → estimated impact: minor, but completes the match to kallisto

Step 3: Rerun JK on ds_52_furtherDownsampled and compare
        → target: JK non-zero count ≈ LK non-zero count (both ~48,659)
        → target: Jaccard(nonzero) ≈ 1.0, Spearman > 0.993

Step 4: Once isoform count problem is resolved, revisit Fix B (flens.txt)
        → now that we know LK uses uniform, Fix B should NOT use raw flens.txt
        → kallisto uses uniform in --long mode → JK should too
        → Fix B as currently coded is WRONG (gives JK different weights than LK)
        → Either remove Fix B or implement it only for short-read mode
```

---

## Files

| File | Description |
|------|-------------|
| `JOLI_Kallisto/diagnose_unique_jk_transcripts.py` | Classifies 3,222 JK-unique transcripts. All Cause B. |
| `JOLI_Kallisto/check_lk_efflens.py` | Back-calculates eff_lens used by LK. Confirms uniform. |
| `JOLI_Kallisto/check_hypothesis1_duplicates.py` | (Not yet run — Hypothesis 1 ruled out by eff_len check) |
| `files/results/.../unique_jk_transcripts.tsv` | Detail TSV of 3,222 transcripts. |

---

## Additional Finding: Fix B is Conceptually Wrong for Long-Read Mode

`flens.txt` values (mean 3,591 bp) are real transcript lengths. Kallisto in `--long` mode
uses **uniform eff_lens = 1.0** (not flens.txt). Applying flens.txt in JK gives it
different weights than LK → different EM → lower correlation.

Fix B should be either removed or scoped only to short-read mode where kallisto DOES use
real eff_lens. For long-read mode, both JK and LK should use uniform weights.
