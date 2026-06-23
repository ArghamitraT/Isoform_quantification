# Short-Read Effective Length: Python Replication Plan
**Date:** 2026-06-23

Goal: compute per-transcript effective lengths in Python (no `kallisto quant` call),
matching kallisto's output exactly, so JOLI EM never wastes a full quantification run.

---

## Part 1 — How kallisto computes effective lengths (source audit)

### Step A: Collect the empirical FLD (fragment length distribution)

**File:** `src/ProcessReads.cpp` lines 1177-1178

During pseudoalignment of paired-end reads, for every read pair that maps
**uniquely** to a single transcript, kallisto records the insert size `tl`:
```cpp
if (0 < tl && tl < flens.size()) {
    flens[tl]++;
}
```
`flens` is a `uint32_t` histogram of size `MAX_FRAG_LEN = 1000`.
`flens[i]` = number of paired-end reads with insert size `i`.

**Key fact: `kallisto bus` runs the same pseudoalignment and builds the same
`flens` histogram. It writes `flens.txt` to the output directory for short
paired-end reads (`main.cpp` lines 2511–2526).**

Format of `kallisto bus`'s `flens.txt`:
```
0 0 0 ... 42 158 1234 ... 0
```
One line, space-separated integers, 1000 values — `flens[0]` through `flens[999]`.

---

### Step B: Compute `mean_fl_trunc` — the length-conditional mean FLD

**File:** `src/MinCollector.cpp`, `compute_mean_frag_lens_trunc()`

```cpp
for (size_t i = 1; i < MAX_FRAG_LEN; ++i) {
    mass[i]   = flens[i] * i + mass[i-1];
    counts[i] = flens[i]     + counts[i-1];
    if (counts[i] > 0)
        mean_fl_trunc[i] = mass[i] / counts[i];
}
```

`mean_fl_trunc[i]` = E[fragment_length | fragment_length ≤ i]
= cumulative mean of the FLD up to length i.

For long transcripts (i >> mean), `mean_fl_trunc[i] ≈ global_mean`.
For short transcripts (i < mean), `mean_fl_trunc[i] << global_mean`.
This is the key trick: short transcripts cannot produce long fragments, so their
effective fragment length is smaller.

---

### Step C: Get per-transcript conditional means

**File:** `src/weights.cpp`, `get_frag_len_means()`

```cpp
for each transcript i of length L:
    if L >= MAX_FRAG_LEN:
        fl_mean[i] = mean_fl_trunc[MAX_FRAG_LEN - 1]  // global mean
    else:
        fl_mean[i] = mean_fl_trunc[L]                  // length-specific mean
```

---

### Step D: Compute effective length

**File:** `src/weights.cpp`, `calc_eff_lens(lengths, means)` (non-deprecated overload)

```cpp
eff_len[i] = lengths[i] - means[i] + 1.0
if eff_len[i] < 1.0:
    eff_len[i] = lengths[i]  // fallback: use raw transcript length
```

---

### Summary of what `kallisto quant` writes to `abundance.tsv`

The `eff_length` column in `abundance.tsv` = output of `calc_eff_lens()` above.
It is fully determined by the FLD histogram and the transcript lengths.
No EM result is needed to compute it.

---

## Part 2 — Python Implementation

### Where to add the code: `core/weights.py`

Add three functions mirroring the C++ exactly:

```python
def load_fld_histogram(flens_path: str) -> np.ndarray:
    """
    Load the FLD histogram written by kallisto bus for short paired-end reads.
    Format: one line, space-separated integers, 1000 values (flens[0]..flens[999]).
    Returns np.ndarray of uint32, shape (MAX_FRAG_LEN=1000,).
    """

def compute_mean_fl_trunc(fld: np.ndarray) -> np.ndarray:
    """
    Compute the length-conditional cumulative mean of the FLD.
    Mirrors MinCollector::compute_mean_frag_lens_trunc().

    mean_fl_trunc[i] = E[fragment_length | fragment_length <= i]
                     = cumulative_mass[i] / cumulative_counts[i]

    Returns np.ndarray of float64, shape (MAX_FRAG_LEN=1000,).
    """

def calc_eff_lens_from_fld(
    transcript_lengths: np.ndarray,
    mean_fl_trunc: np.ndarray,
) -> np.ndarray:
    """
    Compute per-transcript effective lengths from the truncated FLD means.
    Mirrors get_frag_len_means() + calc_eff_lens() in weights.cpp.

    For each transcript of length L:
        fl_mean = mean_fl_trunc[min(L, MAX_FRAG_LEN - 1)]
        eff_len = L - fl_mean + 1
        if eff_len < 1.0: eff_len = L   (fallback)

    Returns np.ndarray of float64, shape (n_transcripts,).
    """
```

Add a new mode `"fld"` (or rename existing `"kallisto"` path) in `compute_weights()`:
- Accept `fld_path` argument pointing to `flens.txt` from `kallisto bus`
- Call the three functions above
- Return `WeightData` as before

---

## Part 3 — FLD Source: `kallisto bus` already has it

**No `kallisto quant` call needed at all.**

`kallisto bus` writes `flens.txt` (the FLD histogram) to `CACHE_DIR` for
short paired-end reads. It is already present after Step 1 of the pipeline.

Pipeline change in `run_joli_kallisto.sh`:
- **Remove Step 3.5 entirely** (the `kallisto quant` run).
- Pass `fld_path="${CACHE_DIR}/flens.txt"` (the bus-produced histogram) to `main_joli.py`.
- `main_joli.py` calls `calc_eff_lens_from_fld()` → identical eff_lengths.

**Before (wasteful):**
```
kallisto bus  →  bustools sort/count  →  kallisto quant (Step 3.5, just for eff_lens)  →  JOLI EM
```

**After (clean):**
```
kallisto bus  →  bustools sort/count  →  JOLI EM (reads FLD from bus's flens.txt)
```

---

## Part 4 — Verification (must pass before merging)

After implementing the Python functions, run this check:

```python
# Load FLD histogram from kallisto bus output
fld = load_fld_histogram("CACHE_DIR/flens.txt")
mean_fl_trunc = compute_mean_fl_trunc(fld)

# Load transcript lengths (from transcripts.txt + index)
# Compute eff_lens in Python
eff_lens_python = calc_eff_lens_from_fld(transcript_lengths, mean_fl_trunc)

# Load kallisto quant's eff_lengths (ground truth)
import pandas as pd
abund = pd.read_csv("quant_tmp/abundance.tsv", sep="\t")
eff_lens_kallisto = abund["eff_length"].values

# Compare
import numpy as np
max_diff = np.abs(eff_lens_python - eff_lens_kallisto).max()
mean_diff = np.abs(eff_lens_python - eff_lens_kallisto).mean()
print(f"max_diff={max_diff:.6f}  mean_diff={mean_diff:.6f}")
# Expected: max_diff < 0.01, mean_diff ~ 0
```

Both should be effectively zero — the Python and C++ implementations should
give bit-for-bit identical results since both use integer accumulation on the
same `flens` array.

---

## Part 5 — Files to change

| File | Change |
|------|--------|
| `core/weights.py` | Add `load_fld_histogram`, `compute_mean_fl_trunc`, `calc_eff_lens_from_fld`; add `"fld"` branch in `compute_weights()` |
| `main_joli.py` | Pass `fld_path` instead of `flens_path`; update `eff_len_mode` logic |
| `scripts/run_joli_kallisto.sh` | Remove Step 3.5; pass `--fld_path ${CACHE_DIR}/flens.txt` |
| `test/test_em_algorithm.py` | Add test verifying eff_lens match kallisto output |

No changes needed to `core/load_tcc.py`, `core/em_algorithm.py`, or
`core/output_writer.py`.

---

## Open questions

1. Does `kallisto bus` always write `flens.txt` for short paired reads, or only
   sometimes? Verify by checking `CACHE_DIR` after a `kallisto bus` run. Look
   for `flens.txt` alongside `output.bus`.

2. Transcript lengths: `calc_eff_lens_from_fld` needs `transcript_lengths` in
   the same order as `transcripts.txt`. Confirm they come from the kallisto index
   (`index.target_lens_`) which matches `transcripts.txt` order.

3. The `flens.txt` written by `kallisto bus` has `MAX_FRAG_LEN=1000` entries.
   Confirm the file has exactly 1000 space-separated values before implementing
   the Python reader.
