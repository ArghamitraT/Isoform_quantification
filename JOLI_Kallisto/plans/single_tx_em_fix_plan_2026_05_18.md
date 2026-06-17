# Plan: Include Single-tx Counts in the EM M-step

**Date:** 2026-05-18  
**File to change:** `core/em_algorithm.py` (only)

---

## Problem

Single-transcript ECs (reads that map unambiguously to one transcript) are currently
excluded from the EM loop and added to `alpha` only **after convergence**. This means
the M-step normalizes theta using only multi-tx read counts:

```
theta_i = n_i_multi / sum(n_multi)
```

A transcript A with 100 single-tx reads and transcript B with 30 single-tx reads look
**identical** to the EM when they share a multi-tx EC `{A, B}`. The 100 vs 30 asymmetry
is invisible to the optimizer — the E-step splits ambiguous reads without that evidence.

---

## Correct EM Formulation

Single-tx ECs have deterministic assignment (probability = 1.0), so their counts are
a **fixed addend** to the expected count of their transcript in every M-step. No change
to the E-step is needed.

**Correct M-step:**
```
theta_i = (n_i_multi + n_i_single) / (total_multi + total_single)
```

- `n_i_multi` = expected counts from multi-tx ECs (computed by E-step, changes each round)
- `n_i_single` = sum of raw counts from all single-tx ECs mapping to transcript i (fixed)
- `total_all_reads = total_multi + total_single` (fixed)

---

## Changes to `core/em_algorithm.py`

### 1. `_preprocess()` — pre-compute fixed per-transcript single-tx counts

After `self._total_multi_reads` is set, add:

```python
# Per-transcript total of single-tx read counts (constant — added every M-step)
self._single_tx_counts_per_tx = np.zeros(self.n_transcripts, dtype=np.float64)
if len(self._single_tx_ids) > 0:
    np.add.at(self._single_tx_counts_per_tx, self._single_tx_ids,
               self._single_ec_counts.astype(np.float64))
self._total_single_reads = float(self._single_tx_counts_per_tx.sum())
self._total_all_reads    = self._total_multi_reads + self._total_single_reads
```

### 2. `run()` — M-step: add single-tx counts every round

Inside the EM loop, change:
```python
n = self._em_step(theta)
numerator = n              # (or n + alpha_prior)
```
to:
```python
n = self._em_step(theta)
n_total = n + self._single_tx_counts_per_tx   # fixed single-tx contribution
numerator = n_total        # (or n_total + alpha_prior)
```

Fix A gating should check `n_total >= min_read_support` instead of `n >= min_read_support`.

### 3. `run()` — convergence check: use `_total_all_reads`

Kallisto-mode convergence scales theta to raw counts for comparison. Change:
```python
alpha_new = theta_new * total_multi_reads
alpha_old = theta     * total_multi_reads
```
to:
```python
alpha_new = theta_new * self._total_all_reads
alpha_old = theta     * self._total_all_reads
```

### 4. `run()` — final alpha: remove post-convergence single-tx addition

Change:
```python
alpha = theta * total_multi_reads
np.add.at(alpha, self._single_tx_ids, self._single_ec_counts)  # REMOVE
```
to:
```python
alpha = theta * self._total_all_reads   # single-tx already baked into theta
```

### 5. `em_step()` (public, used by MultiSampleJoliEM) — same fix

Apply the same `n_total = n + self._single_tx_counts_per_tx` addition and use
`self._total_all_reads` in the convergence check. No changes needed in
`multi_sample_em.py` — it calls `em_step()`, so the fix cascades automatically.

---

## What does NOT change

- `_em_step()` (private) — E-step is correct as-is; single-tx ECs don't need
  fractional assignment (they are always assigned with probability 1.0).
- `load_tcc.py`, `output_writer.py`, `multi_sample_em.py`, `dirichlet_optimizer.py`
  — no changes needed.
- `_active_tx_mask` — still correct; masks transcripts with zero reads.

---

## Expected Behavioral Change

- Transcripts with strong single-tx support get higher theta during EM → pull more
  ambiguous reads toward themselves correctly in the E-step.
- Transcripts appearing only in multi-tx ECs are unaffected (`_single_tx_counts_per_tx[i] = 0`).
- FP leakage should decrease: low-signal transcripts with near-zero single-tx counts
  compete against transcripts with real single-tx evidence and lose.
- JK SS and JK MS metrics (Spearman, Pearson, MAE, RMSE) expected to improve,
  especially for transcripts that are in both EC types.

---

## Tests to Update

`test/test_em_algorithm.py` — update expected values for:
- M-step numerator (now includes single-tx)
- Final `alpha` values (no longer a sum of two separate arrays)
- Convergence round counts (may change slightly)
