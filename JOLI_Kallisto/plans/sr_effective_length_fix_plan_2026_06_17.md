# Short-Read Effective Length Fix Plan
**Date:** 2026-06-17

---

## Problem

When running `run_joli_kallisto.sh` on short-read (Illumina paired-end) data,
`_multi_flat_weights` in the EM are all 1.0 — meaning effective length
normalization is silently not happening.

### Root cause chain

1. `EFF_LEN_MODE="auto"` correctly resolves to `"kallisto"` for short reads.
2. Step 3.5 tries to generate `flens.txt` by running `kallisto quant-tcc`.
3. `quant-tcc` without `--fragment-file` or `-l`/`-s` does **no effective
   length normalization** (stated explicitly in its help text).
4. `flens.txt` is never written → fallback to `eff_len_mode="uniform"` → all
   eff_lens = 1.0 → all `ec_weights` = 1.0.

### Why `quant-tcc --long` is not the solution

Running `quant-tcc --long` does produce `flens.txt`, but the values are
long-read-specific effective lengths (small numbers like 9, 10, 41, 44),
not transcript lengths and not appropriate for short reads.

---

## Solution

Use `kallisto quant` (not `quant-tcc`) for short reads in Step 3.5.

From the kallisto manual:
> "For paired-end reads, the average fragment length can be directly
> estimated from the reads and the program will do so if -l is not used
> (this is the preferred run mode)."

`kallisto quant` on paired-end reads:
- Auto-estimates the fragment length distribution (FLD) from the reads.
- Outputs `abundance.tsv` with an `eff_length` column per transcript.
- Requires no external fragment length parameters.

We **ignore** `kallisto quant`'s EM output entirely. We only extract the
`eff_length` column from its `abundance.tsv` and write it as `flens.txt`
for JOLI EM.

---

## Implementation

### Change: Step 3.5 in `run_joli_kallisto.sh`

Replace the current short-read branch (which ran `quant-tcc` without FLD
info) with the following:

```bash
else
    # Short reads: run kallisto quant to auto-estimate FLD from paired-end
    # reads and extract per-transcript effective lengths from its abundance.tsv.
    # We discard kallisto's EM result — only eff_length column is used.
    QUANT_TMP="${CACHE_DIR}/quant_tmp"
    mkdir -p "${QUANT_TMP}"
    "${KALLISTO}" quant \
        -i "${INDEX_FILE}" \
        -o "${QUANT_TMP}" \
        -t "${THREADS}" \
        --plaintext \
        "${READS_DIR}/${READS_FILE1}" \
        "${READS_DIR}/${READS_FILE2}" \
        2>&1 | tee -a "${LOG}"

    # Extract eff_length column (col 3, skip header) → space-separated flens.txt
    # matching the format expected by load_flens() in core/load_tcc.py
    awk 'NR>1 {printf "%s ", $3}' "${QUANT_TMP}/abundance.tsv" \
        > "${CACHE_DIR}/flens.txt"
fi
```

### Why this works

- `kallisto quant` uses the same index as `kallisto bus`, so transcript order
  in its `abundance.tsv` matches `transcripts.txt` from bustools exactly.
- The `eff_length` column in `abundance.tsv` is the per-transcript effective
  length computed from the auto-estimated FLD — exactly what `load_flens()`
  needs (one value per transcript, space-separated).
- No fragment length parameters required from the user.
- `flens.txt` is cached in `CACHE_DIR`, so it is only computed once per sample.

---

## Files to change

| File | Change |
|------|--------|
| `scripts/run_joli_kallisto.sh` | Replace Step 3.5 short-read branch as above |

No changes needed to `core/weights.py`, `core/load_tcc.py`, or `main_joli.py`.

---

## Verification

After the fix, re-run with pdb breakpoint at `weights.py:204`:

```
(Pdb) eff_lens[:10]   # should show values like 856.0, 1204.0, not all 1.0
(Pdb) eff_lens.min()  # should be > 1.0
```

And `em_algorithm.py:269`:

```
(Pdb) self._multi_flat_weights[:10]   # should show 1/eff_len values, not all 1.0
```

---

## Open questions

- Does `kallisto quant`'s `abundance.tsv` for short reads list transcripts in
  the same order as the index / `transcripts.txt`? (Almost certainly yes, but
  verify on first run by cross-checking the first few `target_id` values
  against `transcripts.txt`.)
- Should the `quant_tmp/` directory be cleaned up after extracting `flens.txt`
  to save disk space? (It contains a full quantification run's output.)

---

## TODO — Next Steps

### TODO 1: Fix `run_lr_kallisto.sh` short-read effective length

**Problem:** `run_lr_kallisto.sh` uses `kallisto quant-tcc` for Step 4 for all
read types. For short reads, `quant-tcc` without a proper `--fld-file` gives
all eff_lens = 1.0 (same root cause as the JOLI fix above). The existing
`FLD_FILE_FLAG` logic is dead code — `flens.txt` is never produced by
`kallisto bus`.

**Fix:** For short reads, replace Step 4 `quant-tcc` with `kallisto quant`
directly. `kallisto quant` auto-estimates the FLD from paired-end reads and
computes correct effective lengths internally in a single call. Long reads keep
`quant-tcc --long -P PLATFORM` unchanged.

**Files to change:** `scripts/run_lr_kallisto.sh` — Step 4 short-read branch only.

---

### TODO 2: Eliminate the redundant `kallisto quant` run in `run_joli_kallisto.sh`

**Problem:** Step 3.5 runs `kallisto quant` on short reads solely to extract
eff_lengths, discarding its full EM result. This wastes time running a complete
quantification that we never use.

**Better approach:** Instead of running `kallisto quant` and throwing away its
abundance output, use `kallisto quant --fragment-file` to save only the FLD
histogram, then compute eff_lengths from the FLD directly (without running EM).
Alternatively, restructure so the `kallisto quant` result in `quant_tmp/` is
reused or at minimum the `quant_tmp/` directory is cleaned up post-extraction.

**Files to change:** `scripts/run_joli_kallisto.sh` Step 3.5 short-read branch;
possibly `core/weights.py` if FLD → eff_len computation is moved there.
