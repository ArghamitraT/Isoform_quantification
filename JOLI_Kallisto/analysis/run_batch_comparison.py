"""
run_batch_comparison.py
=======================
Run compare_abundance_files.py for every sample across two experiment folders.

Use this after running both pipelines (e.g. run_joli_kallisto.sh and
run_lr_kallisto.sh) to compare their outputs sample-by-sample in one go.

Usage:
    conda activate Joli_kallisto
    python run_batch_comparison.py

Inputs:
    - EXP1, EXP2: two timestamped experiment folder names (set in CONFIG).
    - SAMPLES: list of sample names to compare (subfolder names inside each experiment).
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

# Timestamped experiment folder names (just the folder name, not the full path)
EXP1 = "exprmnt_2026_03_15__12_05_48"   # e.g. JOLI-Kallisto run
EXP2 = "exprmnt_2026_03_15__12_05_53"   # e.g. lr-kallisto run

# Sample names to compare — must match the subfolder names inside each experiment
SAMPLES = [
    "flnc_01",
    "flnc_02",
    "flnc_03",
    "flnc_31",
    "flnc_32",
]

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


def run_comparison_for_sample(sample: str) -> int:
    """
    Call compare_abundance_files.py for a single sample.

    Constructs and runs the command:
        python compare_abundance_files.py EXP1 EXP2
            --results-base RESULTS_BASE
            --sample SAMPLE
            --abundance-file ABUNDANCE_FILE
            --top-n TOP_N

    Stdout from the subprocess is printed live so the user can see per-sample
    progress as comparisons complete.

    Args:
        sample (str): Sample name (subfolder inside each experiment folder).

    Returns:
        int: Subprocess return code (0 = success, non-zero = failure).
    """
    cmd = [
        sys.executable, str(COMPARE_PY),
        EXP1, EXP2,
        "--results-base",    RESULTS_BASE,
        "--sample",          sample,
        "--abundance-file",  ABUNDANCE_FILE,
        "--top-n",           str(TOP_N),
    ]

    print(f"\n{'=' * 70}")
    print(f"Sample: {sample}")
    print(f"  Command: {' '.join(cmd)}")
    print(f"{'=' * 70}")

    result = subprocess.run(cmd, text=True)
    return result.returncode


def main() -> None:
    """
    Iterate over SAMPLES and run compare_abundance_files.py for each.

    Tracks which samples succeeded and which failed, prints a final summary,
    and exits with a non-zero code if any comparison failed.

    Returns:
        None
    """
    base = Path(RESULTS_BASE)

    # Sanity check: both experiment folders must exist
    for exp in [EXP1, EXP2]:
        exp_dir = base / exp
        if not exp_dir.is_dir():
            print(f"ERROR: experiment folder not found: {exp_dir}")
            sys.exit(1)

    print("=" * 70)
    print(f"Batch comparison — {TIMESTAMP}")
    print(f"  EXP1   : {EXP1}")
    print(f"  EXP2   : {EXP2}")
    print(f"  Samples: {len(SAMPLES)}")
    print(f"  Base   : {RESULTS_BASE}")
    print("=" * 70)

    succeeded = []
    failed    = []

    for sample in SAMPLES:
        # Check that the abundance file exists in both experiments before running
        missing = []
        for exp in [EXP1, EXP2]:
            abund_path = base / exp / sample / ABUNDANCE_FILE
            if not abund_path.exists():
                missing.append(str(abund_path))

        if missing:
            print(f"\n[SKIP] {sample} — abundance file(s) not found:")
            for m in missing:
                print(f"       {m}")
            failed.append(sample)
            continue

        rc = run_comparison_for_sample(sample)
        if rc == 0:
            succeeded.append(sample)
        else:
            print(f"[ERROR] compare_abundance_files.py returned {rc} for sample {sample}")
            failed.append(sample)

    # Print final summary
    print("\n" + "=" * 70)
    print(f"Batch comparison complete — {TIMESTAMP}")
    print(f"  Succeeded : {len(succeeded)} / {len(SAMPLES)}")
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
