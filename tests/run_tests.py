"""Minimal test runner (pytest is not installed in ast192-hj; files are pytest-compatible).

    python tests/run_tests.py
"""
import importlib
import pathlib
import sys
import traceback

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

failed = passed = 0
for f in sorted((ROOT / "tests").glob("test_*.py")):
    mod = importlib.import_module(f"tests.{f.stem}")
    for name in sorted(n for n in dir(mod) if n.startswith("test_")):
        try:
            getattr(mod, name)()
            passed += 1
            print(f"PASS {f.stem}.{name}")
        except Exception:
            failed += 1
            print(f"FAIL {f.stem}.{name}")
            traceback.print_exc()
print(f"\n{passed} passed, {failed} failed")
sys.exit(1 if failed else 0)
