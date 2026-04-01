# False Positive Fix Plan — JK MS Leaked Transcripts
**Date:** 2026-03-29

---

## Problem Statement

JK MS (multi-sample MAP EM) produces **16,066 false positive** non-zero predictions for sim1
(transcripts predicted non-zero but absent from the ground truth). LK (plain EM) produces
only 6,092 false positives on the same data. The Dirichlet prior is the amplifying factor.

**Contingency table (sim1):**
```
                                    LK        JK MS
True  positives  (pred>0, GT>0) :  62,199    69,285
False positives  (pred>0, GT=0) :   6,092    16,066
False negatives  (pred=0, GT>0) :   8,956     1,870
True  negatives  (pred=0, GT=0) :       0   131,034
```

JK MS recovers ~7K transcripts LK misses (FN: 8,956 → 1,870) — this is a feature of the
Dirichlet prior. But it also introduces ~10K extra false positives relative to LK. The goal
is to keep the FN improvement while reducing the FP inflation.

---

## Root Cause Analysis

### What these false positives are

Transcripts that are:
- **Not expressed** in the simulation (absent from GT)
- **Present in the kallisto index** (218K reference, GT only covers 71-79K)
- **Active** — they appear in multi-transcript ECs alongside expressed transcripts,
  because they share k-mer structure with expressed isoforms

LK also has 6,092 FPs from the same mechanism — multi-mapping leakage is shared.
The Dirichlet prior makes it 2.6× worse in JK MS.

### Exact mechanism (tracing through `em_algorithm.py`)

For an active false-positive transcript `t` during JK MS training:

```
Round N:
  E-step:  weighted_theta[t] = theta[t] * weight[t,ec]  (small but > 0)
           → n[t] > 0  (fractional reads borrowed from co-occurring expressed transcripts)

  M-step:  numerator[t] = n[t] + alpha[t]
           → both terms > 0; alpha prevents convergence toward zero

  Mask:    _active_tx_mask[t] = True  (it IS in a multi-tx EC)
           → numerator NOT zeroed by mask

  Norm:    theta[t] ≈ 1/218K ≈ 4.6e-6

  Floor:   threshold = 1e-8 on normalized theta
           → 4.6e-6 >> 1e-8  → NOT zeroed
```

**These transcripts are never zeroed during training.** Alpha is not resurrecting
them — it is preventing them from converging toward zero in the first place.

In plain EM (LK, no alpha): fractional counts shrink each round as expressed
transcripts dominate → eventually falls below threshold → zeroed. Alpha
short-circuits this natural convergence.

### Secondary issue: `write_results` resurrection (bug)

`multi_sample_em.py:308-315` runs one additional MAP EM round during output writing.
For transcripts that were barely zeroed during training (theta crossed below 1e-8):
```
init_theta[t] = 0      (zeroed during training)
n[t]          = 0      (theta=0 → no contribution to E-step)
numerator[t]  = 0 + alpha[t] > 0   (prior resurrects)
_active_tx_mask[t] = True           (in multi-tx EC → mask doesn't help)
→ theta_new[t] > 0  → survives into output
```
This is a secondary contributor on top of root cause 1. Fix is separate but simpler.

---

## Why a Simple Threshold Fix Won't Work

Raising the zeroing threshold uniformly would reduce FP but would also re-introduce
the FN problem — exactly the transcripts that JK MS correctly recovers over LK.
The threshold cannot distinguish:
- Truly low-expressed transcript (e.g. GT TPM = 1.008, theta ≈ 1e-5): **keep**
- Spurious multi-mapping transcript (no GT expression, theta ≈ 4.6e-6): **zero**

Both have similar theta magnitude. A blanket threshold trades one error for the other.

---

## Proposed Fixes (ordered by priority)

### Fix A — Gate the prior on direct read support (recommended)
**File:** `core/em_algorithm.py`, M-step inside `run()`

**Idea:** Only apply alpha to transcripts where the E-step produced meaningful read
support. Transcripts with pure multi-mapping leakage (tiny n) get plain EM (no prior
boost) — plain EM naturally drives them to zero. Truly low-expressed transcripts with
real reads keep the prior benefit.

```python
# Current M-step:
numerator = n + alpha_prior

# Proposed M-step:
min_read_support = <threshold>   # e.g. 0.1 expected reads, to be tuned
if alpha_prior is not None:
    has_support = n >= min_read_support
    numerator = np.where(has_support, n + alpha_prior, n)
else:
    numerator = n
```

**Why this works:** Alpha amplifies real signal. Transcripts that only appear in
multi-tx ECs receive tiny fractional n (<<0.1 reads) — plain EM drives them toward
zero. Transcripts with real reads (single-tx ECs or dominant in multi-tx ECs) keep
the prior and are protected from over-zeroing.

**Tuning:** `min_read_support` needs to be chosen carefully. Start with 0.1 (fractional
reads) and validate against GT. Too high → reintroduces FN. Too low → doesn't help FP.

**Risk:** Low — the fix only removes the prior from transcripts with negligible read
support. The 69K true positives will all have n >> 0.1.

---

### Fix B — Tighten the active mask (alternative)
**File:** `core/em_algorithm.py`, `_preprocess()`

**Idea:** Instead of marking a transcript active if it appears in ANY EC, require
minimum read evidence — e.g., appears in a single-tx EC (unambiguous read) OR
appears in multi-tx ECs with total count ≥ threshold.

```python
# Current:
self._active_tx_mask[self._multi_flat_tx] = True

# Proposed: only mark active if total multi-tx EC count for this transcript >= min_count
tx_multi_counts = np.zeros(self.n_transcripts)
np.add.at(tx_multi_counts, self._multi_flat_tx,
           self._multi_ec_counts[self._multi_flat_ec_idx])
high_count_mask = tx_multi_counts >= min_multi_count   # e.g. 1.0
self._active_tx_mask |= high_count_mask
```

**Why this is less flexible than Fix A:** The mask is computed once at init (before
EM runs), so it can't adapt to how theta evolves. Fix A is per-round and more precise.

---

### Fix C — write_results resurrection (separate bug fix)
**File:** `core/multi_sample_em.py`, `write_results()`

**Idea:** Instead of running one MAP EM round with alpha_prior (which resurrects
zeroed transcripts), use the converged theta_s directly to build alpha for output_writer,
or run with `alpha_prior=None` (plain EM, no resurrection).

```python
# Option C1: skip the extra EM pass entirely — use converged theta_s directly
#            theta_s is already normalized; reconstruct alpha manually
total_multi_reads = em._multi_ec_counts.sum()
alpha_out = theta_s * total_multi_reads
# add single-tx counts
np.add.at(alpha_out, em._single_tx_ids, em._single_ec_counts)

# Option C2: run the extra pass with alpha_prior=None (plain EM, no resurrection)
result = em.run(
    max_em_rounds    = 1,
    min_rounds       = 1,
    convergence_mode = self.convergence_mode,
    alpha_prior      = None,          # ← no prior; prevent resurrection
    init_theta       = theta_s,
)
```

**Note:** Fix C alone does not address root cause 1 (transcripts that were never
zeroed during training). It only prevents the secondary resurrection of borderline
transcripts. Implement alongside Fix A, not instead of it.

---

## Implementation Order

1. **Run EDA first** (see `eda_plan_2026_03_29.md`, Priority 6 — leaked transcript
   analysis). Specifically need:
   - `fig_leaked_ec_sizes.png` — are FPs exclusively in large multi-tx ECs?
   - `fig_leaked_vs_lk.png` — does LK assign near-zero to the same transcripts?
   - The leaked TPM CDF — what threshold would eliminate most FPs without harming TPs?

2. **Implement Fix C** (write_results resurrection) — simple, low-risk, isolated.
   Quantify how many of the 16K FPs disappear. Remainder = root cause 1.

3. **Implement Fix A** (gate prior on read support) — main fix.
   Tune `min_read_support` against GT metrics. Target: FP ≈ LK's 6K while keeping
   FN ≈ JK MS's current 1,870 (not regressing to LK's 8,956).

4. **Rerun both experiments** and recompute contingency table + all three metric sets
   from `run_gt_comparison.py` to validate.

---

## Success Criteria

After fixes, on sim1:
```
                                    LK     JK MS (current)   JK MS (target)
True  positives  (pred>0, GT>0) : 62,199      69,285            ~69K  (keep)
False positives  (pred>0, GT=0) :  6,092      16,066            ~6-8K (reduce)
False negatives  (pred=0, GT>0) :  8,956       1,870            ~2-3K (keep low)
Active universe Spearman        :  0.???        0.???            > LK
GT-nonzero Spearman             :  0.824        0.814            > 0.824
```

The key constraint: **do not regress the FN improvement**. JK MS recovering 7K
transcripts that LK misses is the value of the multi-sample prior. The fix must
preserve that while eliminating the spurious multi-mapping FPs.
