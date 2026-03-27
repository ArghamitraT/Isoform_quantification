# Dummy Bustools Dataset Plan

Small, hand-crafted dataset for debugging and unit-testing the JOLI EM pipeline
without needing real reads or a real kallisto/bustools run.

---

## Transcripts (10 total)

| Index | Name   | Length (bp) | Gene group |
|-------|--------|-------------|------------|
| 0     | TX_A1  | 1200        | Gene A (isoforms A1, A2, A3) |
| 1     | TX_A2  | 900         | Gene A |
| 2     | TX_A3  | 1500        | Gene A |
| 3     | TX_B1  | 800         | Gene B (isoforms B1, B2) |
| 4     | TX_B2  | 1100        | Gene B |
| 5     | TX_C1  | 700         | Gene C (isoforms C1, C2) |
| 6     | TX_C2  | 950         | Gene C |
| 7     | TX_D1  | 600         | Gene D (unique) |
| 8     | TX_E1  | 1300        | Gene E (unique) |
| 9     | TX_F1  | 850         | Gene F (unique) |

---

## Equivalence Classes (14 total)

### Ambiguous ECs (multiple transcripts)

| EC index | Transcripts (indices) | Transcript names      | # reads | Notes |
|----------|-----------------------|-----------------------|---------|-------|
| 0        | 0, 1                  | TX_A1, TX_A2          | 6       | A1/A2 ambiguous |
| 1        | 0, 1, 2               | TX_A1, TX_A2, TX_A3   | 4       | All Gene A ambiguous |
| 2        | 1, 2                  | TX_A2, TX_A3          | 3       | A2/A3 ambiguous |
| 3        | 3, 4                  | TX_B1, TX_B2          | 7       | All Gene B ambiguous |
| 4        | 5, 6                  | TX_C1, TX_C2          | 5       | All Gene C ambiguous |
| 5        | 0, 3                  | TX_A1, TX_B1          | 2       | Cross-gene ambiguity |
| 6        | 2, 4, 6               | TX_A3, TX_B2, TX_C2   | 2       | 3-way cross-gene |

### Unique ECs (single transcript — unambiguous)

| EC index | Transcript (index) | Transcript name | # reads | Notes |
|----------|--------------------|-----------------|---------|-------|
| 7        | 0                  | TX_A1           | 4       | unique to A1 |
| 8        | 2                  | TX_A3           | 3       | unique to A3 |
| 9        | 3                  | TX_B1           | 5       | unique to B1 |
| 10       | 4                  | TX_B2           | 3       | unique to B2 |
| 11       | 7                  | TX_D1           | 4       | Gene D only |
| 12       | 8                  | TX_E1           | 5       | Gene E only |
| 13       | 9                  | TX_F1           | 3       | Gene F only |

**Total reads: 56**

---

## Expected EM behavior

- **TX_A1** (EC 0, 1, 7 + cross): lots of ambiguous reads + 4 unique → should get
  non-trivial mass, pulled by unique evidence.
- **TX_A2** (EC 0, 1, 2): no unique reads → mass comes entirely from ambiguous ECs;
  EM must infer it from shared evidence.
- **TX_A3** (EC 1, 2, 8 + 3-way): moderate unique evidence (3 reads in EC 8).
- **TX_B1**, **TX_B2**: paired ambiguity (EC 3) + unique reads for each → should
  resolve reasonably.
- **TX_D1**, **TX_E1**, **TX_F1**: fully unambiguous → EM trivially assigns all reads;
  good sanity-check targets.

---

## Files to create under `.../kallisto_output/dummy/`

| File            | Format | Content |
|-----------------|--------|---------|
| `transcripts.txt` | one name per line | 10 transcript names |
| `matrix.ec`       | `<EC_idx>\t<comma-sep transcript indices>` | 14 ECs |
| `count.mtx`       | Matrix Market sparse | 1 barcode × 14 ECs, read counts |
| `count.ec.txt`    | `<EC_idx>\t<comma-sep transcript indices>` | same as matrix.ec (bustools count format) |
| `count.barcodes.txt` | single line | `"bulk"` |
| `flens.txt`       | one float per line | effective lengths (one per transcript, 10 lines) |
| `run_info.json`   | JSON | dummy run metadata |

---

## count.mtx format (Matrix Market)

```
%%MatrixMarket matrix coordinate real general
%
1 14 14
1 1 6.0
1 2 4.0
1 3 3.0
1 4 7.0
1 5 5.0
1 6 2.0
1 7 2.0
1 8 4.0
1 9 3.0
1 10 5.0
1 11 3.0
1 12 4.0
1 13 5.0
1 14 3.0
```

*(1-indexed: row=barcode, col=EC index+1, value=read count)*

---

---

## Pen-and-paper EM simulation (Round 0)

Using `eff_len_mode="uniform"` so all weights = 1.0.

### Internal flat arrays built by `_preprocess()`

Single-tx ECs (EC 7–13) are split out and never enter the loop.
Multi-tx ECs (EC 0–6) are expanded into flat parallel arrays:

```
flat_tx      = [0, 1,  0, 1, 2,  1, 2,  3, 4,  5, 6,  0, 3,  2, 4, 6]
flat_weights = [1, 1,  1, 1, 1,  1, 1,  1, 1,  1, 1,  1, 1,  1, 1, 1]
flat_ec_idx  = [0, 0,  1, 1, 1,  2, 2,  3, 3,  4, 4,  5, 5,  6, 6, 6]
multi_ec_counts = [6, 4, 3, 7, 5, 2, 2]   (one per multi-EC)
```

16 positions total, 7 multi-ECs, 7 single-ECs.

---

### Step 1 — Initialisation

```python
theta = np.full(10, 1.0 / 10)   # all = x = 0.1
```

---

### Step 2 — `theta_per_pos` and `weighted_theta`

```python
theta_per_pos  = theta[flat_tx]              # all x  (uniform theta, uniform flat_tx)
weighted_theta = theta_per_pos * flat_weights # all x  (weights = 1)
```

---

### Step 3 — `denominators = np.zeros(7)`

Seven zeros, one per multi-EC.

---

### Step 4 — `np.add.at(denominators, flat_ec_idx, weighted_theta)`

Accumulates `weighted_theta` (all x) grouped by `flat_ec_idx`:

| local EC | transcripts | tx count | denominator |
|----------|-------------|----------|-------------|
| 0 | A1, A2 | 2 | 2x |
| 1 | A1, A2, A3 | 3 | 3x |
| 2 | A2, A3 | 2 | 2x |
| 3 | B1, B2 | 2 | 2x |
| 4 | C1, C2 | 2 | 2x |
| 5 | A1, B1 | 2 | 2x |
| 6 | A3, B2, C2 | 3 | 3x |

`denominators = [2x, 3x, 2x, 2x, 2x, 2x, 3x]`

---

### Step 5 — Broadcast back to positions

```python
denom_per_pos = denominators[flat_ec_idx]
count_per_pos = multi_ec_counts[flat_ec_idx]
```

| pos | tx | EC | denom_per_pos | count_per_pos |
|-----|----|----|---------------|---------------|
| 0  | A1 | 0 | 2x | 6 |
| 1  | A2 | 0 | 2x | 6 |
| 2  | A1 | 1 | 3x | 4 |
| 3  | A2 | 1 | 3x | 4 |
| 4  | A3 | 1 | 3x | 4 |
| 5  | A2 | 2 | 2x | 3 |
| 6  | A3 | 2 | 2x | 3 |
| 7  | B1 | 3 | 2x | 7 |
| 8  | B2 | 3 | 2x | 7 |
| 9  | C1 | 4 | 2x | 5 |
| 10 | C2 | 4 | 2x | 5 |
| 11 | A1 | 5 | 2x | 2 |
| 12 | B1 | 5 | 2x | 2 |
| 13 | A3 | 6 | 3x | 2 |
| 14 | B2 | 6 | 3x | 2 |
| 15 | C2 | 6 | 3x | 2 |

---

### Step 6 — `contributions = count_per_pos * weighted_theta / denom_per_pos`

`weighted_theta` = all x → **x cancels** → contributions are independent of initial theta value:

| pos | tx | calc | contribution |
|-----|----|------|--------------|
| 0  | A1 | 6·x / 2x | **3** |
| 1  | A2 | 6·x / 2x | **3** |
| 2  | A1 | 4·x / 3x | **4/3** |
| 3  | A2 | 4·x / 3x | **4/3** |
| 4  | A3 | 4·x / 3x | **4/3** |
| 5  | A2 | 3·x / 2x | **3/2** |
| 6  | A3 | 3·x / 2x | **3/2** |
| 7  | B1 | 7·x / 2x | **7/2** |
| 8  | B2 | 7·x / 2x | **7/2** |
| 9  | C1 | 5·x / 2x | **5/2** |
| 10 | C2 | 5·x / 2x | **5/2** |
| 11 | A1 | 2·x / 2x | **1** |
| 12 | B1 | 2·x / 2x | **1** |
| 13 | A3 | 2·x / 3x | **2/3** |
| 14 | B2 | 2·x / 3x | **2/3** |
| 15 | C2 | 2·x / 3x | **2/3** |

---

### Step 7 — `np.add.at(n, flat_tx, contributions)`

Accumulate contributions per transcript:

| tx | positions summed | n[t] | ≈ |
|----|-----------------|------|---|
| A1 (0) | 3 + 4/3 + 1 | **16/3** | 5.33 |
| A2 (1) | 3 + 4/3 + 3/2 | **35/6** | 5.83 |
| A3 (2) | 4/3 + 3/2 + 2/3 | **7/2** | 3.50 |
| B1 (3) | 7/2 + 1 | **9/2** | 4.50 |
| B2 (4) | 7/2 + 2/3 | **25/6** | 4.17 |
| C1 (5) | 5/2 | **5/2** | 2.50 |
| C2 (6) | 5/2 + 2/3 | **19/6** | 3.17 |
| D1 (7) | (single-tx only) | **0** | — |
| E1 (8) | (single-tx only) | **0** | — |
| F1 (9) | (single-tx only) | **0** | — |

Sum = (32+35+21+27+25+15+19)/6 · (1/29·29) = **29** ✓  (= total multi-tx reads)

---

### Step 8 — M-step: `theta_new = n / 29`

| tx | n[t] | theta_new | ≈ |
|----|------|-----------|---|
| A1 | 16/3 | 16/87 | 0.1839 |
| A2 | 35/6 | 35/174 | 0.2011 |
| A3 | 7/2 | 7/58 | 0.1207 |
| B1 | 9/2 | 9/58 | 0.1552 |
| B2 | 25/6 | 25/174 | 0.1437 |
| C1 | 5/2 | 5/58 | 0.0862 |
| C2 | 19/6 | 19/174 | 0.1092 |
| D1 | 0 | 0 | — |
| E1 | 0 | 0 | — |
| F1 | 0 | 0 | — |

Sum = 174/174 = **1** ✓

---

### Step 9 — Convergence check

```python
changed = count of t where:
    theta_new[t] > 1e-2   AND
    |theta_new[t] - theta[t]| / theta_new[t] > 1e-2
```

theta_old = 0.1 for all:

| tx | theta_new | >1e-2? | \|Δ\|/theta_new | >1e-2? | changed? |
|----|-----------|--------|-----------------|--------|----------|
| A1 | 0.1839 | ✓ | 0.456 | ✓ | yes |
| A2 | 0.2011 | ✓ | 0.503 | ✓ | yes |
| A3 | 0.1207 | ✓ | 0.171 | ✓ | yes |
| B1 | 0.1552 | ✓ | 0.355 | ✓ | yes |
| B2 | 0.1437 | ✓ | 0.304 | ✓ | yes |
| C1 | 0.0862 | ✓ | 0.160 | ✓ | yes |
| C2 | 0.1092 | ✓ | 0.084 | ✓ | yes |
| D1 | 0 | ✗ | — | — | no |
| E1 | 0 | ✗ | — | — | no |
| F1 | 0 | ✗ | — | — | no |

**changed = 7** → not converged. Also `round_num=0 < min_rounds=50` so loop continues regardless.

---

### Key observations from Round 0

1. **x cancels** — the initial theta value (1/10) is irrelevant for round 0. Only EC
   structure (how many transcripts share each EC) determines the first update.

2. **A2 gets the highest theta (0.2011) despite zero unique reads** — pulled up by
   EC0 (6 reads, split 2 ways) and EC1 (4 reads, split 3 ways). EM correctly infers
   abundance from ambiguous evidence alone.

3. **D1, E1, F1 drop to zero in theta** — they only appear in single-tx ECs, so
   the multi-tx EM never sees them. Their counts (4, 5, 3) are added back
   post-convergence via Fix A:
   ```python
   np.add.at(alpha, self._single_tx_ids, self._single_ec_counts)
   ```

4. Round 1 starts with **non-uniform theta** → denominators will no longer be
   symmetric → ambiguous reads start redistributing toward transcripts with more
   unique evidence (e.g. A1 will pull ahead of A2 because A1 has 4 unique reads
   in EC7 that feed back into its alpha, whereas A2 has none).

---

## Notes / open questions

- `flens.txt` in the real data has one entry per transcript (not per EC).
  Use transcript lengths from the table above.
- `matrix.ec` uses 0-based EC indices; `count.mtx` is 1-based (Matrix Market standard).
- `count.ec.txt` appears to be identical in content to `matrix.ec` — confirm by
  checking `load_tcc.py` to see which file is actually read.
- Cross-gene ECs (EC 5, 6) are intentionally included to test that the EM doesn't
  assume all ECs are within a single gene.
