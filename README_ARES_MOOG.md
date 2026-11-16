# Running the ARES+MOOG [Fe/H] analysis on the Tufts cluster

This guide reruns the ARES+MOOG metallicity analysis on the spectra that are **already downloaded on the cluster**,
then brings the results back to a laptop for `Metallicity_ARES_MOOG.ipynb`.

ARES measures Fe I and Fe II equivalent widths. MOOG, with Kurucz ATLAS9 models and run through FASMA, then iterates
Teff, logg, microturbulence and [Fe/H] until the Fe lines are in excitation and ionisation balance. This is the
SWEET-Cat recipe (Sousa 2014; Sousa et al. 2021). Every stage reuses the downloaded spectra, so none of the steps
below needs internet access.

More detail: [`hpc/ares_moog/README.md`](hpc/ares_moog/README.md) (installation, validation, debugging one star)
and [`hpc/README.md`](hpc/README.md) (the download pipeline that produced the spectra).

## What is where

The paths below are the ones used so far. Replace `ndahle01` and `martinlab` with your own user and lab directory.

| Location on the cluster | Contents |
|---|---|
| `~/Astro192/` | The project code. `hpc/` plus `eblm_feh.py`, `eso_coverage.py`, `sophie_coverage.py` and `eblm_metallicity_data/` (`phoenix_grid.npz`, `sweetcat.csv`) |
| `/cluster/tufts/martinlab/ndahle01/metallicity/` (`METAL_ROOT`) | The main sample: `catalogs/` (targets and manifest), `raw/` (downloaded spectra, ~50 GB), and the outputs below |
| `/cluster/tufts/martinlab/ndahle01/ares_moog/` (`ARESMOOG_HOME`) | The ARES, MOOG and FASMA installation and its conda env |
| `/cluster/tufts/martinlab/ndahle01/ew_validation/` | The 4-star validation run (separate `METAL_ROOT`) |

Keep all of these out of `$HOME`, which has a 30 GB quota.

## 0. Update the code on the cluster

On your laptop:

```bash
rsync -av ~/Desktop/Programming/Astro192/AST_192_Fall26/hpc/ ndahle01@login.pax.tufts.edu:~/Astro192/hpc/
```

Also copy `eblm_feh.py`, `eso_coverage.py`, `sophie_coverage.py` or `eblm_metallicity_data/sweetcat.csv` if you
changed them.

## 1. One-time install (skip if `ARESMOOG_HOME` already exists)

On a login node, this takes ~20–40 min:

```bash
cd ~/Astro192
export ARESMOOG_HOME=/cluster/tufts/martinlab/ndahle01/ares_moog
bash hpc/ares_moog/setup_ares_moog.sh
```

It ends with a solar test. You want `Sun: Teff ~5777, logg ~4.44, [Fe/H] ~+0.00 -> PASS`. The script links
`~/.moog_fasma` to MOOG's data files, because MOOG can only store a 60-character path. Don't delete that link.

## 2. Set up the shell

Do this every time you log in and before every `sbatch`, from the project root:

```bash
cd ~/Astro192
conda activate ast192
export METAL_ROOT=/cluster/tufts/martinlab/ndahle01/metallicity
export ARESMOOG_HOME=/cluster/tufts/martinlab/ndahle01/ares_moog
```

`hpc/ew.sbatch` stops immediately if either variable is missing. It switches between the `ast192` env (co-adding)
and `$ARESMOOG_HOME/env` (ARES, MOOG, FASMA) by itself. If conda comes from a module, also
`export MODULES="anaconda/..."`.

## 3. Run the analysis

There are two jobs. Each prints its log to `ew_<jobid>.log` in the directory you submitted from.

**Job 1: co-add, fit, collect.** This takes a few hours on 32 cores.

```bash
sbatch -c 32 -t 2-00:00:00 hpc/ew.sbatch
```

| Stage | What it does | Output (in `$METAL_ROOT`) |
|---|---|---|
| `coadd` | For each (star, instrument, resolution) group: RV-align and co-add the spectra in `raw/` over 4450–6900 Å. All 830 groups take a few hours. | `coadd_ew/` |
| `fit` | ARES equivalent widths, then FASMA/MOOG iterates to the best Teff, logg, ξ and [Fe/H]. ~1–5 min per star per core. | `ew_work/<group>/`, `ew_fits/` |
| `collect` | Gathers every fit into one table. | `results/ew_fits.csv` |

**Job 2: rescue fits, then one [Fe/H] per star.** About an hour on 16 cores.

```bash
STAGES="rescue collect adopt" sbatch hpc/ew.sbatch
```

- **`rescue`:** refits every co-add whose free fit failed or ran to logg = 5, this time with logg fixed at the TIC
  value. It also refits all converged HARPS co-adds the same way, as a calibration of what fixing logg does to
  [Fe/H]. Output: `results/ew_fits_fixlogg.csv`.
- **`adopt`:** picks one value per star and writes `results/feh_adopted.csv`, `feh_candidates.csv` and
  `feh_offsets.csv`. For each instrument it:
  - measures the [Fe/H] offset from HARPS, using both stars that also have HARPS fits and SWEET-Cat's slow rotators;
  - subtracts the offset if it is significant, adds systematic errors, and chooses the candidate with the smallest
    total error;
  - sets flags such as `sophie_he`, `fixed_logg` and `vbroad_ge_10`.

  `adopt` reads vbroad from the PHOENIX results in `results/metallicity_final.csv`. Without that file there is no
  rotation flag. `adopt` needs only pandas and numpy, so it also runs on a laptop
  (`python hpc/ew_pipeline.py adopt`).

You can submit job 2 as soon as job 1 is queued:
`STAGES="rescue collect adopt" sbatch -d afterok:<job 1 id> hpc/ew.sbatch`.

**Useful variations:**

```bash
ONLY=WASP-43 STAGES="coadd fit collect" sbatch hpc/ew.sbatch   # one star, as a quick test
STAGES="fit collect" sbatch -c 32 hpc/ew.sbatch                 # refit without re-co-adding
```

Every stage can be rerun. To redo a stage from scratch, delete its output folder first (`coadd_ew/`, `ew_fits/` or
`ew_fits_fixlogg/`).

## 4. Check the run

- **The end of each job log:** `tail -50 ew_<jobid>.log`.
- **Fit statuses:** `ok` means converged. `not_converged` means FASMA hit 150 iterations (the parameters are still
  written). `poor_coverage` means the co-add covers less than 60% of the line list. `ares_failed` means ARES produced
  no output.
- **To debug one star:** look in `$METAL_ROOT/ew_work/<group>/`. It holds `result.json`, `fasma.log`,
  `captain.log`, the measured EWs in `linelist/`, and MOOG's `summary.out`.

## 5. Bring the results back and look at them

On your laptop, from the project root, check first with `-n` (a dry run that copies nothing):

```bash
rsync -avn ndahle01@login.pax.tufts.edu:/cluster/tufts/martinlab/ndahle01/metallicity/results/ metallicity_hpc_data/results/
rsync -avn ndahle01@login.pax.tufts.edu:/cluster/tufts/martinlab/ndahle01/ew_validation/results/ metallicity_hpc_data/ew_validation_results/
```

Then run the same commands without `-n`. Copying over `results/` replaces the local tables. If you have rerun
`adopt` or `combine` locally since the last cluster run, the local files are newer, so back them up first or add
`-u` (skip files that are newer locally).

Then open `Metallicity_ARES_MOOG.ipynb` and run all cells. It reads:

- `metallicity_hpc_data/results/`: `ew_fits.csv`, `ew_fits_fixlogg.csv`, `feh_candidates.csv`, `feh_adopted.csv`,
  `feh_offsets.csv` and `metallicity_final.csv`;
- `metallicity_hpc_data/catalogs/targets.csv`;
- `metallicity_hpc_data/ew_validation_results/ew_validation.csv`;
- `eblm_metallicity_data/` (`sweetcat.csv`, `eblm_literature_feh.csv`).

## Known limitations

- **CORALIE is missing.** CORALIE spectra aren't in the ESO or SOPHIE archives, so ~100 EBLMs that the PHOENIX
  analysis has are missing from this one.
- **Fast rotators (vbroad ≥ 10 km/s) read 0.2–0.4 dex low** compared with SWEET-Cat. They are flagged
  `vbroad_ge_10`.
- **About 29% of free fits fail** and are rescued with logg fixed. Those stars are flagged `fixed_logg`.
- **SOPHIE HE needs a correction.** Its equivalent widths come out ~12% low, so `adopt` subtracts an offset
  (−0.116 dex), and the stars are flagged `sophie_he`. A small S/N-dependent residual (about −0.04 dex for the
  faintest stars) remains.
