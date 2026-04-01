"""
test_utility.py — Tests for src/utility.py

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/Simulations
    conda activate lrgsp_simulation
    python tests/test_utility.py
"""

import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from utility import create_run_dir, save_runtime, save_code_snapshot


def test_create_run_dir():
    with tempfile.TemporaryDirectory() as tmp:
        run_dir = create_run_dir(base=tmp)
        assert run_dir.exists(), "run_dir not created"
        assert run_dir.name.startswith("exprmnt_"), f"Bad name format: {run_dir.name}"
        print("[PASS] create_run_dir")


def test_save_runtime():
    with tempfile.TemporaryDirectory() as tmp:
        run_dir = create_run_dir(base=tmp)
        save_runtime(run_dir, 99.5)
        content = (run_dir / "runtime.txt").read_text()
        assert "99.50" in content, f"Unexpected content: {content}"
        print("[PASS] save_runtime")


def test_save_code_snapshot():
    with tempfile.TemporaryDirectory() as tmp:
        run_dir = create_run_dir(base=tmp)
        src_root = Path(__file__).parent.parent / "src"
        save_code_snapshot(run_dir, src_root)
        snapshot_dir = run_dir / "code_snapshot"
        assert snapshot_dir.exists(), "code_snapshot dir missing"
        files = list(snapshot_dir.rglob("*"))
        assert len(files) > 0, "No files copied into snapshot"
        print(f"[PASS] save_code_snapshot ({len(files)} files)")


if __name__ == "__main__":
    print("=== test_utility.py ===")
    test_create_run_dir()
    test_save_runtime()
    test_save_code_snapshot()
    print("=== All tests passed ===")
