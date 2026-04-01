"""
utility.py — Shared helpers for the RNA-seq simulation pipeline.

What it does:
    Provides three functions used by every pipeline module:
    - create_run_dir: creates a timestamped output folder under files/results/
    - save_runtime:   writes elapsed wall-clock time to runtime.txt
    - save_code_snapshot: copies all code files into run_dir/code_snapshot/

Inputs:  none (imported as a library)
Outputs: none (side effects only — creates dirs/files on disk)
"""

import shutil
import time
from datetime import datetime
from pathlib import Path


def create_run_dir(base: str = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results") -> Path:
    """
    Create a timestamped results folder and return its path.

    Args:
        base (str): Parent directory under which the run folder is created.

    Returns:
        Path: Path to the newly created folder, e.g. .../exprmnt_2026_03_30__14_05_00
    """
    timestamp = datetime.now().strftime("%Y_%m_%d__%H_%M_%S")
    run_dir = Path(base) / f"exprmnt_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=True)
    print(f"[utility] Run dir created: {run_dir}")
    return run_dir


def save_runtime(run_dir: Path, elapsed_seconds: float) -> None:
    """
    Write elapsed wall-clock time to runtime.txt inside run_dir.

    Args:
        run_dir (Path):         Path to the timestamped results folder.
        elapsed_seconds (float): Total elapsed time in seconds.

    Returns:
        None
    """
    out = run_dir / "runtime.txt"
    out.write_text(f"{elapsed_seconds:.2f} seconds\n")
    print(f"[utility] Runtime saved: {elapsed_seconds:.2f}s → {out}")


def save_code_snapshot(run_dir: Path, src_root: Path) -> None:
    """
    Copy all code files from src_root into run_dir/code_snapshot/, preserving
    relative directory structure.

    Copies files with extensions: .py .sh .txt .yml .yaml .md

    Args:
        run_dir  (Path): Path to the timestamped results folder.
        src_root (Path): Root of the code directory to snapshot (e.g. code/Simulations/).

    Returns:
        None
    """
    extensions = {".py", ".sh", ".txt", ".yml", ".yaml", ".md"}
    snapshot_dir = run_dir / "code_snapshot"
    snapshot_dir.mkdir(parents=True, exist_ok=True)

    copied = 0
    for src_file in src_root.rglob("*"):
        if src_file.is_file() and src_file.suffix in extensions:
            rel = src_file.relative_to(src_root)
            dst = snapshot_dir / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src_file, dst)
            copied += 1

    print(f"[utility] Code snapshot saved: {copied} files → {snapshot_dir}")


# ── Quick self-test ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    import tempfile

    print("=== utility.py self-test ===")

    with tempfile.TemporaryDirectory() as tmp:
        # Test create_run_dir
        run_dir = create_run_dir(base=tmp)
        assert run_dir.exists(), "run_dir was not created"
        assert run_dir.name.startswith("exprmnt_"), "run_dir name format wrong"

        # Test save_runtime
        save_runtime(run_dir, 123.456)
        rt = (run_dir / "runtime.txt").read_text()
        assert "123.46" in rt, f"Unexpected runtime content: {rt}"

        # Test save_code_snapshot (snapshot this file itself)
        src_root = Path(__file__).parent
        save_runtime_start = time.time()
        save_code_snapshot(run_dir, src_root)
        snapshot_dir = run_dir / "code_snapshot"
        assert snapshot_dir.exists(), "code_snapshot dir not created"
        snapped = list(snapshot_dir.rglob("*.py"))
        assert len(snapped) > 0, "No .py files copied into snapshot"

    print("=== All checks passed ===")
