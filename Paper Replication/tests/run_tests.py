"""Minimal stand-in for pytest: runs every test_* function in tests/test_*.py.
Usage:  python tests/run_tests.py"""

import importlib.util
import sys
import time
import traceback
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
failed = 0
for path in sorted(Path(__file__).parent.glob("test_*.py")):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for name in sorted(n for n in dir(module) if n.startswith("test_")):
        t0 = time.time()
        try:
            getattr(module, name)()
            print(f"PASS  {path.name}::{name}  ({time.time() - t0:.1f}s)")
        except Exception:
            failed += 1
            print(f"FAIL  {path.name}::{name}")
            traceback.print_exc()
print(f"\n{failed} failed" if failed else "\nall tests passed")
sys.exit(1 if failed else 0)
