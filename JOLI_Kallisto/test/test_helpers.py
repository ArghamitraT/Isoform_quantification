"""
test_helpers.py
===============
Shared utilities for all JOLI-Kallisto test files.

Provides:
  - make_sample_dir() : create a temporary directory with synthetic bustools output files
  - check()           : print PASS / FAIL for a single assertion and return bool

Inputs:  (called by test files, not run directly)
Outputs: (helpers only, no outputs)
"""

import json
import os
import tempfile

import numpy as np


def make_sample_dir(
    ec_transcripts: list,
    ec_counts: list,
    transcript_names: list,
    n_pseudoaligned: int = None,
    n_targets: int = None,
    include_run_info: bool = True,
) -> str:
    """
    Create a temporary directory with synthetic bustools output files.

    Creates the following files matching the real kallisto/bustools format:
      - matrix.ec       : "ec_id<TAB>tx_idx1,tx_idx2,..."
      - transcripts.txt : one transcript name per line
      - count.mtx       : MatrixMarket sparse (1 row = 1 sample, cols = ECs, 1-indexed)
      - run_info.json   : JSON with n_processed, n_pseudoaligned, n_targets

    Args:
        ec_transcripts   : list[list[int]] -- transcript indices per EC (0-indexed).
        ec_counts        : list[int]       -- read count per EC.
        transcript_names : list[str]       -- transcript name per index.
        n_pseudoaligned  : int             -- written to run_info.json (default: sum of counts).
        n_targets        : int             -- written to run_info.json (default: n transcripts).
        include_run_info : bool            -- whether to write run_info.json.

    Returns:
        str -- path to temporary directory. Caller is responsible for cleanup via shutil.rmtree().
    """
    tmpdir = tempfile.mkdtemp(prefix="joli_test_")
    n_ecs = len(ec_transcripts)
    n_tx  = len(transcript_names)
    if n_pseudoaligned is None:
        n_pseudoaligned = sum(ec_counts)
    if n_targets is None:
        n_targets = n_tx

    # matrix.ec
    with open(os.path.join(tmpdir, "matrix.ec"), "w") as fh:
        for ec_id, txs in enumerate(ec_transcripts):
            fh.write(f"{ec_id}\t{','.join(str(t) for t in txs)}\n")

    # transcripts.txt
    with open(os.path.join(tmpdir, "transcripts.txt"), "w") as fh:
        for name in transcript_names:
            fh.write(f"{name}\n")

    # count.mtx (MatrixMarket: 1 row = 1 sample, cols = ECs, 1-indexed)
    nonzero_entries = [(1, ec_id + 1, cnt)
                       for ec_id, cnt in enumerate(ec_counts) if cnt > 0]
    with open(os.path.join(tmpdir, "count.mtx"), "w") as fh:
        fh.write("%%MatrixMarket matrix coordinate real general\n")
        fh.write("%\n")
        fh.write(f"1 {n_ecs} {len(nonzero_entries)}\n")
        for row, col, val in nonzero_entries:
            fh.write(f"{row} {col} {val}\n")

    # run_info.json
    if include_run_info:
        info = {
            "n_targets":      n_targets,
            "n_bootstraps":   0,
            "n_processed":    n_pseudoaligned + 1000,
            "n_pseudoaligned": n_pseudoaligned,
            "n_unique":       n_pseudoaligned // 2,
            "p_pseudoaligned": 68.8,
            "kallisto_version": "0.51.1",
        }
        with open(os.path.join(tmpdir, "run_info.json"), "w") as fh:
            json.dump(info, fh)

    return tmpdir


def check(condition: bool, test_name: str, detail: str = "") -> bool:
    """
    Print PASS or FAIL for a single assertion.

    Args:
        condition : bool -- True if the assertion passed.
        test_name : str  -- Short description of what is being checked.
        detail    : str  -- Optional extra info printed on failure.

    Returns:
        bool -- same as condition (for chaining with &=).
    """
    if condition:
        print(f"  PASS  {test_name}")
    else:
        print(f"  FAIL  {test_name}" + (f": {detail}" if detail else ""))
    return condition
