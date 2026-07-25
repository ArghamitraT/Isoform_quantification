# kallisto_ec_loader.py
import os
import numpy as np
import re

def read_count_mtx(mtx_path: str):
    """
    MatrixMarket reader for count.mtx.

    - Skips comment lines (starting with '%')
    - Skips blank lines
    - Skips malformed lines (e.g. '1 1255' lines with only 2 entries)
    """
    print(f"Reading MatrixMarket file: {mtx_path}")
    with open(mtx_path, "r") as f:
        # ---- read header: nrows ncols nnz ----
        nrows = ncols = nnz = None
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("%"):
                continue

            parts = s.split()
            if len(parts) != 3:
                print(f"Skipping malformed header-like line in {mtx_path}: {s!r}")
                continue

            nrows, ncols, nnz = map(int, parts)
            break

        if nrows is None:
            raise ValueError(f"Could not find a valid header line in {mtx_path}")

        # ---- read triplets ----
        rows = []
        cols = []
        vals = []

        for line in f:
            s = line.strip()
            if not s:
                continue
            parts = s.split()
            if len(parts) != 3:
                print(f"Skipping non-data line in {mtx_path}: {s!r}")
                continue
            r_str, c_str, v_str = parts
            rows.append(int(r_str))
            cols.append(int(c_str))
            vals.append(float(v_str))

    rows = np.array(rows, dtype=int)
    cols = np.array(cols, dtype=int)
    vals = np.array(vals, dtype=float)

    print(
        f"  nrows={nrows}, ncols={ncols}, nnz_header={nnz}, "
        f"nnz_parsed={len(vals)}"
    )
    return nrows, ncols, nnz, rows, cols, vals


def load_kallisto_ec_with_counts(kallisto_dir: str, max_ec_nonzero: int | None = None):
    """
    Load EC information from a Kallisto/BUS tools folder:

      - transcripts.txt
      - count.ec.txt
      - count.mtx

    Returns:
      all_read_names  : array of pseudo-read labels (EC_x)
      all_rnames      : array of transcript names per row
      compat_list     : 1 / |EC| for each row
      p_nm            : 1.0 for each row (no positional weighting)
      weights         : EC count (multiplicity) for each row
    """
    mtx_path = os.path.join(kallisto_dir, "count.mtx")
    ec_path  = os.path.join(kallisto_dir, "count.ec.txt")
    tx_path  = os.path.join(kallisto_dir, "transcripts.txt")

    # read count matrix
    nrows, ncols, nnz, rows, cols, vals = read_count_mtx(mtx_path)

    # BUS matrix is 1-based so convert to 0-based EC index
    ec_ids = cols - 1  # 0..(ncols-1)
    counts = vals  # multiplicities

    # Optionally limit to first N nonzero ECs
    order = np.arange(len(counts))
    if max_ec_nonzero is not None:
        order = order[:max_ec_nonzero]

    ec_ids = ec_ids[order]
    counts = counts[order]

    # transcripts
    with open(tx_path) as f:
        tx_names = [ln.strip() for ln in f if ln.strip()]

    # EC → list of transcript indices
    ec_to_txs: list[list[int]] = []
    with open(ec_path) as f:
        for line in f:
            s = line.strip()
            if not s:
                ec_to_txs.append([])
                continue

            # Handle both whitespace-separated and comma-separated formats
            # e.g. "1 2 3" or "1,2,3" or "1, 2, 3"
            tokens = re.split(r"[,\s]+", s)
            tokens = [t for t in tokens if t]  # drop empty strings
            try:
                tx_idx = [int(x) for x in tokens]
            except ValueError as e:
                print(f"[EC parse] Skipping malformed EC line: {s!r}")
                raise
            ec_to_txs.append(tx_idx)

    all_read_names = []
    all_rnames = []
    compat_list = []
    p_nm = []
    weights = []

    for global_idx, ec_id in enumerate(ec_ids):
        tx_idx_list = ec_to_txs[ec_id]
        if not tx_idx_list:
            continue

        # number of transcripts in this EC
        k = len(tx_idx_list)
        # EC multiplicity
        c = counts[global_idx]
        ec_label = f"EC_{ec_id}"

        for tx_idx in tx_idx_list:
            all_read_names.append(ec_label)
            all_rnames.append(tx_names[tx_idx])
            compat_list.append(1.0 / k)
            p_nm.append(1.0)
            weights.append(c)

    all_read_names = np.array(all_read_names, dtype=object)
    all_rnames = np.array(all_rnames, dtype=object)
    compat_list = np.array(compat_list, dtype=float)
    p_nm = np.array(p_nm, dtype=float)
    weights = np.array(weights, dtype=float)

    return all_read_names, all_rnames, compat_list, p_nm, weights