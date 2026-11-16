# ARES+MOOG [Fe/H] on the cluster

The equivalent-width method used by SWEET-Cat (Sousa 2014; Sousa et al. 2021): ARES measures Fe I and Fe II
equivalent widths, then MOOG and Kurucz ATLAS9 models are iterated to excitation and ionisation balance. It runs
through FASMA (Andreasen et al. 2017). The download and RV-alignment code is shared with `metallicity_pipeline.py`.

| File | What it does |
|---|---|
| `ares_moog/setup_ares_moog.sh` | One-time install (conda env in `$ARESMOOG_HOME/env`, FASMA + MOOG, ARES, Kurucz models, line lists), then a solar test |
| `ares_moog/patch_fasma.py` | Three one-line Python 3 fixes to FASMA (run by the setup script) |
| `ares_moog/validation_targets.csv` | The 4 validation stars, their TIC priors, and the published values to reproduce |
| `ew_pipeline.py` | Stages: `manifest` → (download via `metallicity_pipeline.py`) → `coadd` → `fit` → `collect`; plus `sun-test` |
| `ew.sbatch` | `coadd fit collect` as one Slurm job |

## 0. Copy the code to the cluster

```bash
# on your laptop
rsync -av ~/Desktop/Programming/Astro192/AST_192_Fall26/hpc/ <user>@login.pax.tufts.edu:~/Astro192/hpc/
```
The `coadd` stage also needs `eblm_feh.py`, `eso_coverage.py`, `sophie_coverage.py` and
`eblm_metallicity_data/phoenix_grid.npz` in `~/Astro192`. These are already there from the main pipeline.

## 1. Install (once, on a login node, ~20–40 min, mostly conda)

```bash
cd ~/Astro192
export ARESMOOG_HOME=/cluster/tufts/martinlab/ndahle01/ares_moog
bash hpc/ares_moog/setup_ares_moog.sh
```
It finishes with a solar test: ARES+MOOG on the HARPS Ganymede spectrum that ships with ARES. You want
`Sun: Teff ~5777, logg ~4.44, [Fe/H] ~+0.00 -> PASS`. If it says `CHECK` or `FAILED`, stop there and look at
`$ARESMOOG_HOME/test_root/ew_work/_sun_test/{fasma.log,captain.log,result.json}`.

If conda comes from a module, set `MODULES="anaconda/..."` first. MOOG can only store a 60-character path to its data
files, so the script links `~/.moog_fasma` to them. Don't delete that link.

## 2. Validation run (the 4 stars)

Use a separate `METAL_ROOT` so nothing touches the main metallicity output. The manifest stage refuses to run on the
main directory.

```bash
cd ~/Astro192
conda activate ast192
export METAL_ROOT=/cluster/tufts/martinlab/ndahle01/ew_validation
export ARESMOOG_HOME=/cluster/tufts/martinlab/ndahle01/ares_moog

# login node (needs internet): list and fetch the HARPS spectra (264 spectra, ~1.3 GB, ~10-20 min)
python hpc/ew_pipeline.py manifest --targets hpc/ares_moog/validation_targets.csv
python hpc/metallicity_pipeline.py download --task 0 --ntasks 1

# compute node: co-add, ARES+MOOG, and the comparison table (~30-60 min)
COMPARE=hpc/ares_moog/validation_targets.csv sbatch hpc/ew.sbatch
```

The comparison table is printed at the end of `ew_<jobid>.log` and saved to `$METAL_ROOT/results/ew_validation.csv`.

What to expect:
* **TOI-426** (Castro-González+2026, +0.16 ± 0.04) and **TOI-2211** (Egger+2026, −0.13 ± 0.04): the same archive
  spectra and the same recipe, so the differences should be ≲ 0.03–0.05 dex. Our manifest finds 73 of TOI-426's
  78 HARPS spectra and all 19 of TOI-2211's.
* **HD 48265** (SWEET-Cat, +0.39 ± 0.03) and **HD 176986** (Tsantaki+2013, +0.03 ± 0.05; SWEET-Cat's HARPS-N value is
  +0.01): older analyses on possibly different spectra, so ~0.05 dex differences are acceptable.
* HD 48265 also gets a second co-add from its single HARPS spectrum at R = 80,000 (EGGS mode). The comparison uses
  the higher-S/N co-add.

## 3. Debugging one star

Every co-add has its own working directory, `$METAL_ROOT/ew_work/<group>/`, containing:

| File | Contents |
|---|---|
| `result.json` | Status, parameters, line counts, run time |
| `fasma.log` | ARES and FASMA output, including the iteration log |
| `captain.log` | FASMA's own log |
| `linelist/star*.moog` | The measured EWs (`*_outlier.moog` after 3σ clipping) |
| `summary.out` | MOOG's last abfind output |

To rerun only the fit for one star:

```bash
STAGES="fit collect" ONLY=TOI-426 COMPARE=hpc/ares_moog/validation_targets.csv sbatch hpc/ew.sbatch
```

Statuses:

| Status | Meaning |
|---|---|
| `ok` | Converged |
| `not_converged` | FASMA hit 150 iterations; the parameters are still written |
| `poor_coverage` | The co-add covers less than 60% of the line list |
| `ares_failed` | ARES produced no output |
| `failed: ...` | A Python exception |

## 4. Main sample (after validation passes)

Point `METAL_ROOT` at the main metallicity directory. The main manifest and the `raw/` spectra are reused, and new
output goes only to `coadd_ew/`, `ew_work/`, `ew_fits/` and `results/ew_fits.csv`.

```bash
export METAL_ROOT=/cluster/tufts/martinlab/ndahle01/metallicity
sbatch -c 32 -t 2-00:00:00 hpc/ew.sbatch
```

* Re-co-adding all 830 groups over the wider range takes a few hours.
* Expect ~1–5 min per star for the fit on one core.
* CORALIE (the EBLM-only spectra from the earlier local analysis) isn't in the cluster manifest. Adding it is a
  separate step.

## After the main run: rescue fits and adopt one [Fe/H] per star

```bash
STAGES="rescue collect adopt" sbatch hpc/ew.sbatch      # ~375 fits, about an hour on 16 cores
```

* `rescue` re-fits, with logg fixed at the TIC value, every co-add whose free fit failed or ran to logg = 5, plus all
  converged HARPS co-adds (`role=calib`) to measure what fixing logg does to [Fe/H] -> `results/ew_fits_fixlogg.csv`.
* `adopt` measures each instrument's [Fe/H] offset from HARPS two ways, from free-logg fits: (a) stars with a HARPS
  fit too, and (b) SWEET-Cat's homogeneous values (`eblm_metallicity_data/sweetcat.csv`, `SWFlag = 1`, 5″ match) for
  slow rotators (vbroad < 10 km/s), relative to HARPS's own offset from SWEET-Cat. Modes with < 5 stars are pooled
  per instrument. It subtracts the inverse-variance mean of the two when significant (> 2σ), adds the intrinsic
  scatter as a systematic, adds a 0.03 dex method floor, and picks the candidate with the smallest total error for
  each star -> `results/feh_adopted.csv` (`feh_candidates.csv` = every usable fit, `feh_offsets.csv` = the offset
  table, with both estimates in `offset_pairs` and `offset_sc`). A note is printed when the two disagree by > 2σ.
  Flags: `fixed_logg`, `sophie_he`, `snr_lt_50`, `teff_outside_4500_6500`, `vt_ge_2.5`, `n_fe2_lt_5`,
  `vbroad_ge_10` (fast rotator: ARES+MOOG reads ~0.2–0.4 dex low against SWEET-Cat above 10 km/s),
  `instrument_tension` (another instrument disagrees by > 3 sigma). vbroad (also a column in `feh_adopted.csv`)
  comes from the PHOENIX fits in `results/metallicity_final.csv`; without that file there is no flag and no rotation
  cut. Excluded stars (e.g. SB2s) are listed in `EXCLUDE` in `ew_pipeline.py`.
* Why the SWEET-Cat estimate: with HARPS pairs alone, SOPHIE HE got −0.095 dex from 10 stars; SWEET-Cat (51 fits)
  gives −0.120, and the combined −0.116 brings SOPHIE HE into line with SWEET-Cat (`Metallicity_ARES_MOOG.ipynb`).
  `adopt` runs locally as well as on the cluster (only pandas and numpy needed).
* Why SOPHIE HE is treated this way: `hpc/ares_moog/sophie_test.py` (EWs ~12% low from an additive pedestal plus
  a low continuum; not resolution).
