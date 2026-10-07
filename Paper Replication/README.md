# Morton & Winn (2014) replication with corrected inclination inference

Full plan: `PROJECT_SPEC.md`. Papers live in `../../Papers/`.

```
src/        physics (importable, tested)
  inclination.py   per-star cos i posteriors: correct (Masuda & Winn 2020 Eq. 9) and MW14 (Eq. 13)
  obliquity.py     Fisher distribution, f_cosI(z|κ), κ prior
notebooks/  one notebook per stage; they orchestrate and plot only
tests/      test_*.py (pytest-compatible); run with  python tests/run_tests.py
data/raw/        untouched downloads (cached, never re-queried)
data/processed/  derived tables
data/morton_winn_table1.csv   Table 1 of MW14, transcribed from the PDF (scripts/)
figures/    saved figures
```

Environment: `conda activate ast192` (pytest is not in it, hence `tests/run_tests.py`).

## Status

| Stage | Status |
|---|---|
| 0 — inference engine validation | done, all checks pass (`notebooks/stage0_validation.ipynb`) |
| 1 — reproduce MW14 | not started |
| 2 — corrected inference | not started |
| 3 — extended sample | not started |
| 4 — significance and controls | not started |
