"""
run_batch_comparison.py
=======================
Run compare_abundance_files.py for every common sample across two experiment
folders. Samples are auto-discovered — any subfolder present in BOTH experiment
directories is compared (non-sample dirs like code_snapshot are ignored because
they contain no abundance.tsv).

Usage:
    # Pass experiment folder names as CLI args (recommended):
    python analysis/run_batch_comparison.py exprmnt_2026_03_28__00_29_23 exprmnt_2026_03_28__00_38_39

    # Or set EXP1 / EXP2 in the CONFIG section and run with no args:
    python analysis/run_batch_comparison.py

Inputs:
    - EXP1, EXP2: two timestamped experiment folder names (CLI args or CONFIG).
    - Samples are auto-discovered as subfolders present in both experiments that
      contain an abundance.tsv file. No hardcoded sample list needed.
    - RESULTS_BASE: base directory containing experiment folders.

Outputs:
    For each sample, compare_abundance_files.py is called, which produces:
      {RESULTS_BASE}/{EXP1}/comparison_vs_{EXP2}__{timestamp}.txt
      {RESULTS_BASE}/{EXP2}/comparison_vs_{EXP1}__{timestamp}.txt
      {RESULTS_BASE}/{EXP1}/comparison_detail_vs_{EXP2}__{timestamp}/
          merged_abundance.tsv
          top_differences_raw.tsv
          top_differences_log1p.tsv
          summary_metrics.tsv
          summary_metrics.json
    A combined summary table is printed and saved to:
      {RESULTS_BASE}/{EXP1}/batch_comparison_summary__{timestamp}.tsv
"""

import os
import subprocess
import sys
import time
from pathlib import Path


# ============================================================
# CONFIG — edit these variables before running; do not edit below
# ============================================================

# Timestamped experiment folder names — overridden by CLI args when provided
EXP1 = ""   # e.g. "exprmnt_2026_03_28__00_29_23"
EXP2 = ""   # e.g. "exprmnt_2026_03_28__00_38_39"

# Base directory containing the experiment folders
RESULTS_BASE = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

# Abundance filename within each sample subfolder
ABUNDANCE_FILE = "abundance.tsv"

# Number of most-different transcripts to write per sample comparison
TOP_N = 50

# ============================================================
# END CONFIG
# ============================================================


THIS_DIR   = Path(os.path.dirname(os.path.realpath(__file__)))
COMPARE_PY = THIS_DIR / "compare_abundance_files.py"
TIMESTAMP  = time.strftime("%Y_%m_%d__%H_%M_%S")



def discover_common_samples(exp1_dir: Path, exp2_dir: Path, abundance_file: str) -> list:
    """
    Find sample names present in both experiment folders that contain abundance.tsv.

    Args:
        exp1_dir       : Path -- first experiment directory.
        exp2_dir       : Path -- second experiment directory.
        abundance_file : str  -- abundance filename to look for inside each subfolder.

    Returns:
        list[str] -- sorted list of common sample names.
    """
    def sample_subdirs(exp_dir: Path) -> set:
        return {
            d.name for d in exp_dir.iterdir()
            if d.is_dir() and (d / abundance_file).exists()
        }

    samples1 = sample_subdirs(exp1_dir)
    samples2 = sample_subdirs(exp2_dir)
    common   = sorted(samples1 & samples2)

    only1 = sorted(samples1 - samples2)
    only2 = sorted(samples2 - samples1)
    if only1:
        print(f"  [INFO] Samples only in EXP1 (skipped): {only1}")
    if only2:
        print(f"  [INFO] Samples only in EXP2 (skipped): {only2}")

    return common


def main() -> None:
    """
    Parse CLI args, auto-discover common samples, and run compare_abundance_files.py
    for each. Tracks successes/failures and exits non-zero if any comparison failed.

    Returns:
        None
    """
    import argparse
    parser = argparse.ArgumentParser(
        description="Batch comparison of two experiment folders."
    )
    parser.add_argument("exp1", nargs="?", default=EXP1,
                        help="First experiment folder name (overrides CONFIG EXP1).")
    parser.add_argument("exp2", nargs="?", default=EXP2,
                        help="Second experiment folder name (overrides CONFIG EXP2).")
    parser.add_argument("--results-base", default=RESULTS_BASE,
                        help=f"Base results directory. Default: {RESULTS_BASE}")
    parser.add_argument("--abundance-file", default=ABUNDANCE_FILE,
                        help=f"Abundance filename. Default: {ABUNDANCE_FILE}")
    parser.add_argument("--top-n", type=int, default=TOP_N,
                        help=f"Top-N differences per comparison. Default: {TOP_N}")
    args = parser.parse_args()

    exp1         = args.exp1
    exp2         = args.exp2
    results_base = args.results_base
    abund_file   = args.abundance_file
    top_n        = args.top_n

    if not exp1 or not exp2:
        print("ERROR: provide EXP1 and EXP2 as CLI args or set them in CONFIG.")
        sys.exit(1)

    base = Path(results_base)

    # Sanity check: both experiment folders must exist
    for exp in [exp1, exp2]:
        exp_dir = base / exp
        if not exp_dir.is_dir():
            print(f"ERROR: experiment folder not found: {exp_dir}")
            sys.exit(1)

    # Auto-discover common samples
    samples = discover_common_samples(base / exp1, base / exp2, abund_file)

    if not samples:
        print(f"ERROR: no common samples with {abund_file} found in both experiments.")
        sys.exit(1)

    print("=" * 70)
    print(f"Batch comparison — {TIMESTAMP}")
    print(f"  EXP1   : {exp1}")
    print(f"  EXP2   : {exp2}")
    print(f"  Samples: {samples}")
    print(f"  Base   : {results_base}")
    print("=" * 70)

    succeeded = []
    failed    = []

    for sample in samples:
        cmd = [
            sys.executable, str(COMPARE_PY),
            exp1, exp2,
            "--results-base",   results_base,
            "--sample",         sample,
            "--abundance-file", abund_file,
            "--top-n",          str(top_n),
        ]
        print(f"\n{'=' * 70}")
        print(f"Sample: {sample}")
        print(f"  Command: {' '.join(cmd)}")
        print(f"{'=' * 70}")

        rc = subprocess.run(cmd, text=True).returncode
        if rc == 0:
            succeeded.append(sample)
        else:
            print(f"[ERROR] compare_abundance_files.py returned {rc} for {sample}")
            failed.append(sample)

    # Print final summary
    print("\n" + "=" * 70)
    print(f"Batch comparison complete — {TIMESTAMP}")
    print(f"  Succeeded : {len(succeeded)} / {len(samples)}")
    if succeeded:
        for s in succeeded:
            print(f"    [OK]   {s}")
    if failed:
        print(f"  Failed / skipped : {len(failed)}")
        for s in failed:
            print(f"    [FAIL] {s}")
    print("=" * 70)

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
