"""
test_main_multisample_joli.py
=============================
Tests for JOLI_Kallisto/main_multisample_joli.py.

Covers:
  - resolve_sample_dirs() with explicit SAMPLE_DIRS
  - resolve_sample_dirs() with SAMPLES_FOLDER auto-discovery
  - create_run_dir() produces a timestamped directory
  - save_experiment_description() writes correct file with all CONFIG keys
  - save_runtime() writes runtime.txt
  - save_code_snapshot() copies .py/.sh files into code_snapshot/

Note: full end-to-end main() is NOT tested here (that requires real TCC data
and would be an integration test). These tests cover the helper functions only.

Run with:
  conda activate NanoCount_5
  cd JOLI_Kallisto/
  python test/test_main_multisample_joli.py
"""

import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import main_multisample_joli as mmj
from test_helpers import check


# ============================================================
# Test 1 — resolve_sample_dirs with explicit SAMPLE_DIRS
# ============================================================

def test_resolve_sample_dirs_explicit():
    """SAMPLE_DIRS list is returned as absolute paths when SAMPLES_FOLDER is empty."""
    print("\n[test_resolve_sample_dirs_explicit]")
    passed = True

    # Temporarily patch CONFIG values
    original_dirs   = mmj.SAMPLE_DIRS
    original_folder = mmj.SAMPLES_FOLDER

    dir_a = tempfile.mkdtemp(prefix="joli_test_sa_")
    dir_b = tempfile.mkdtemp(prefix="joli_test_sb_")
    try:
        mmj.SAMPLE_DIRS    = [dir_a, dir_b]
        mmj.SAMPLES_FOLDER = ""

        result = mmj.resolve_sample_dirs()
        passed &= check(len(result) == 2,
                        "resolve_sample_dirs() returns 2 entries",
                        f"got {len(result)}")
        passed &= check(all(os.path.isabs(p) for p in result),
                        "all paths are absolute")
        passed &= check(os.path.abspath(dir_a) in result,
                        "dir_a in result")
        passed &= check(os.path.abspath(dir_b) in result,
                        "dir_b in result")
    finally:
        mmj.SAMPLE_DIRS    = original_dirs
        mmj.SAMPLES_FOLDER = original_folder
        shutil.rmtree(dir_a)
        shutil.rmtree(dir_b)
    return passed


# ============================================================
# Test 2 — resolve_sample_dirs with SAMPLES_FOLDER auto-discovery
# ============================================================

def test_resolve_sample_dirs_folder():
    """SAMPLES_FOLDER auto-discovers all subdirectories."""
    print("\n[test_resolve_sample_dirs_folder]")
    passed = True

    original_dirs   = mmj.SAMPLE_DIRS
    original_folder = mmj.SAMPLES_FOLDER

    parent = tempfile.mkdtemp(prefix="joli_test_parent_")
    try:
        # Create 3 subdirectories inside parent
        for name in ["sample_1", "sample_2", "sample_3"]:
            os.makedirs(os.path.join(parent, name))
        # Also create a file (should NOT be included)
        open(os.path.join(parent, "readme.txt"), "w").close()

        mmj.SAMPLE_DIRS    = []
        mmj.SAMPLES_FOLDER = parent

        result = mmj.resolve_sample_dirs()
        passed &= check(len(result) == 3,
                        "auto-discovers 3 subdirs (file excluded)",
                        f"got {len(result)}")
        passed &= check(all(os.path.isdir(p) for p in result),
                        "all discovered paths are directories")
        passed &= check(result == sorted(result),
                        "results are sorted")
    finally:
        mmj.SAMPLE_DIRS    = original_dirs
        mmj.SAMPLES_FOLDER = original_folder
        shutil.rmtree(parent)
    return passed


# ============================================================
# Test 3 — create_run_dir creates timestamped directory
# ============================================================

def test_create_run_dir():
    """create_run_dir() creates a new directory named exprmnt_{timestamp}."""
    print("\n[test_create_run_dir]")
    passed = True

    original_base  = mmj.RESULTS_BASE
    tmp_base       = tempfile.mkdtemp(prefix="joli_test_base_")
    try:
        mmj.RESULTS_BASE = tmp_base
        run_dir = mmj.create_run_dir()

        passed &= check(os.path.isdir(run_dir),
                        "run_dir exists and is a directory")
        passed &= check(os.path.basename(run_dir).startswith("exprmnt_"),
                        "run_dir name starts with 'exprmnt_'",
                        f"got {os.path.basename(run_dir)}")
        passed &= check(run_dir.startswith(tmp_base),
                        "run_dir is inside RESULTS_BASE")
    finally:
        mmj.RESULTS_BASE = original_base
        shutil.rmtree(tmp_base)
    return passed


# ============================================================
# Test 4 — save_experiment_description writes log file
# ============================================================

def test_save_experiment_description():
    """experiment_description.log is written with all CONFIG keys and sample list."""
    print("\n[test_save_experiment_description]")
    passed = True

    run_dir = tempfile.mkdtemp(prefix="joli_test_run_")
    sample_dirs = ["/fake/sample_1", "/fake/sample_2"]
    try:
        mmj.save_experiment_description(run_dir, sample_dirs)

        log_path = os.path.join(run_dir, "experiment_description.log")
        passed &= check(os.path.exists(log_path),
                        "experiment_description.log exists")

        with open(log_path) as fh:
            content = fh.read()

        for key in ["READ_TYPE", "EFF_LEN_MODE", "CONVERGENCE_MODE",
                    "MAX_EM_ROUNDS", "MAX_GD_ROUNDS", "GD_LR",
                    "ALPHA_INITIAL", "GD_CONVERGENCE_TOL"]:
            passed &= check(key in content,
                            f"'{key}' present in experiment_description.log")

        passed &= check("/fake/sample_1" in content,
                        "sample_1 path listed in log")
        passed &= check("/fake/sample_2" in content,
                        "sample_2 path listed in log")
    finally:
        shutil.rmtree(run_dir)
    return passed


# ============================================================
# Test 5 — save_runtime writes runtime.txt
# ============================================================

def test_save_runtime():
    """runtime.txt is written with elapsed time in seconds."""
    print("\n[test_save_runtime]")
    passed = True

    run_dir = tempfile.mkdtemp(prefix="joli_test_run_")
    try:
        mmj.save_runtime(run_dir, 123.45)

        rt_path = os.path.join(run_dir, "runtime.txt")
        passed &= check(os.path.exists(rt_path), "runtime.txt exists")

        with open(rt_path) as fh:
            content = fh.read().strip()
        passed &= check("123.45" in content,
                        "runtime.txt contains elapsed time '123.45'",
                        f"got: '{content}'")
    finally:
        shutil.rmtree(run_dir)
    return passed


# ============================================================
# Test 6 — save_code_snapshot copies source files
# ============================================================

def test_save_code_snapshot():
    """code_snapshot/ is created and contains .py files from JOLI_Kallisto/."""
    print("\n[test_save_code_snapshot]")
    passed = True

    run_dir = tempfile.mkdtemp(prefix="joli_test_run_")
    try:
        mmj.save_code_snapshot(run_dir)

        snapshot_dir = os.path.join(run_dir, "code_snapshot")
        passed &= check(os.path.isdir(snapshot_dir),
                        "code_snapshot/ directory exists")

        files = os.listdir(snapshot_dir)
        py_files = [f for f in files if f.endswith(".py")]
        passed &= check(len(py_files) > 0,
                        f"code_snapshot/ contains .py files",
                        f"got {len(py_files)}")
        passed &= check("main_multisample_joli.py" in py_files,
                        "main_multisample_joli.py is in snapshot")
        passed &= check("main_joli.py" in py_files,
                        "main_joli.py is in snapshot")
    finally:
        shutil.rmtree(run_dir)
    return passed


# ============================================================
# Runner
# ============================================================

def run_all():
    tests = [
        test_resolve_sample_dirs_explicit,
        test_resolve_sample_dirs_folder,
        test_create_run_dir,
        test_save_experiment_description,
        test_save_runtime,
        test_save_code_snapshot,
    ]

    results = []
    for t in tests:
        try:
            results.append(t())
        except Exception as e:
            print(f"  ERROR in {t.__name__}: {e}")
            import traceback; traceback.print_exc()
            results.append(False)

    print("\n" + "=" * 50)
    passed = sum(results)
    total  = len(results)
    print(f"Results: {passed}/{total} test functions passed")
    if passed == total:
        print("ALL PASSED")
    else:
        print("SOME FAILED — see details above")
    return passed == total


if __name__ == "__main__":
    success = run_all()
    sys.exit(0 if success else 1)
