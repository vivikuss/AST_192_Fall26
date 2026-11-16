# Metallicity pipeline on the cluster

Homogeneous spectroscopic [Fe/H] for every hot-Jupiter host and EBLM primary with public ESO (HARPS, FEROS, UVES,
ESPRESSO) or SOPHIE spectra, plus an independent PAWS/iSpec re-analysis as a check.

| Script | What it does |
|---|---|
| `metallicity_pipeline.py` | manifest → download → co-add → fit (the `eblm_feh` PHOENIX fit) → combine |
| `paws_pipeline.py` | prepare → run PAWS headless (iSpec EW + synthesis) → collect |
| `submit_all.sh`, `submit_paws.sh` | submit each pipeline as a chain of dependent Slurm jobs |
| `metallicity.sbatch`, `paws.sbatch` | each pipeline as a single `sbatch` job (all stages in one job) |
| `run_stage.sbatch` | the generic job script both chains use |
| `setup_paws_env.sh` | one-time install of iSpec v2023.08.04 + PAWS (pinned) in a `paws` conda env |
| `ew_pipeline.py`, `ew.sbatch`, `ares_moog/` | ARES+MOOG (SWEET-Cat recipe, via FASMA) [Fe/H]; see `ares_moog/README.md` |

Results go to `$METAL_ROOT/results/`; `Metallicity_Survey.ipynb` (project root) plots them.

## 1. One-time set-up

```bash
# on your laptop: copy the project, including the untracked data the pipeline needs
rsync -av --exclude 'metallicity_hpc_data' --exclude '*.fits' --exclude '.DS_Store' \
    ~/Desktop/Programming/Astro192/AST_192_Fall26/ <user>@login.pax.tufts.edu:~/Astro192/
#   needed from eblm_metallicity_data/: phoenix_grid.npz, solar_line_mask.npz, sweetcat.csv,
#   nasa_pscomppars_P10.csv, eblm_literature_feh.csv, eblm_targets_tic.csv, sophie_footprint.csv;
#   and eblm_master_table.csv at the top level

# on the cluster
cd ~/Astro192
conda env create -f environment_ast192.yaml        # once; add `requests pyvo` if they are missing
export METAL_ROOT=/cluster/tufts/<lab>/<user>/metallicity   # ~60 GB free; not $HOME
export PAWS_HOME=/cluster/tufts/<lab>/<user>/paws_software  # ~4 GB
bash hpc/setup_paws_env.sh                          # only needed for the PAWS check
```

If conda comes from a module, set e.g. `export MODULES="anaconda/2023.07"` (check `module avail`) before both
scripts. The partition defaults to `batch`; override with `PARTITION=...`.

## 2. Run

Simplest: one job per pipeline, all stages inside it (run from the project root, with `METAL_ROOT` and `PAWS_HOME`
exported):

```bash
sbatch hpc/metallicity.sbatch                            # 16 CPUs, 2 days max; prints the job id
sbatch -d afterok:<that job id> hpc/paws.sbatch          # 32 CPUs, 2 days max; starts when the first one ends
ONLY=WASP-43 STAGES="manifest download coadd fit" sbatch hpc/metallicity.sbatch   # quick test on one star
```

Options are environment variables listed at the top of each file (`STAGES`, `ONLY`, `N_DOWNLOAD`, `PAWS_TARGETS`,
`PAWS_ARGS`, …); resources can be overridden on the command line (`sbatch -p <partition> -c 32 -t 3-00:00:00 ...`).
The job log goes to `metallicity_<jobid>.log` / `paws_<jobid>.log` in the directory you submitted from.

Or, to spread the work over many smaller jobs (job arrays, which usually start sooner on a busy cluster):

```bash
bash hpc/submit_all.sh                           # manifest -> download -> coadd -> fit -> combine
AFTER=<download array job id> bash hpc/submit_paws.sh   # PAWS can start once the downloads finish
```

Expected size (Nov 2026 archive contents): **13,600 spectra, ~50 GB, 830 co-add groups, 606 stars** (365 hot-Jupiter
hosts, 241 EBLMs). The download dominates the wall time (hours, set by the archives, not by CPUs). Co-adding and
fitting take ~2 s per group on one core. PAWS takes ~0.5–3 h per group on one core, for 87 groups (published-PAWS
EBLMs) or ~505 groups (`PAWS_TARGETS=all`).

**Internet on compute nodes.** The `manifest` and `download` stages need outbound HTTPS to archive.eso.org,
dataportal.eso.org, atlas.obs-hp.fr, mast.stsci.edu and exoplanetarchive.ipac.caltech.edu. If the compute nodes
cannot reach them, point `DL_PARTITION` at a partition that can, or run those two stages on a login/transfer node:

```bash
python hpc/metallicity_pipeline.py manifest
for i in 0 1 2 3; do python hpc/metallicity_pipeline.py download --task $i --ntasks 4 & done; wait
```

Every stage can be re-run: downloads, co-adds and cached catalogs are kept. Delete `$METAL_ROOT/fits/` to refit,
`$METAL_ROOT/coadd/` to re-co-add, or `$METAL_ROOT/catalogs/manifest.csv` to rebuild the download list. `--only
<target_id,...>` restricts any stage to a few stars, for tests.

## 3. What each stage does

* **manifest**: targets are the 696 hot-Jupiter hosts (`eblm_feh.is_hot_jupiter`) and the 794 EBLMs in
  `eblm_master_table.csv`, with TIC v8 Teff/logg priors and NASA-archive transit ephemerides. Spectra must cover
  5150–6100 Å, be public, and have S/N ≥ 10. The match radius is 10″ (hot Jupiters, SOPHIE) or 20″ (ESO EBLM pointings
  are often 10–20″ off). Exposures taken in transit are dropped (they carry the RM distortion). Each (star,
  instrument, resolution) group keeps its 50 best exposures (15 for ESPRESSO) and is kept only if the co-added
  S/N ≥ 50.
* **co-add**: each exposure is read in its own format: ESPRESSO vacuum → air, UVES topocentric (handled by the RV
  shift), and SOPHIE images. Each is normalised, its RV is measured by cross-correlation against a PHOENIX template at
  the TIC Teff, and it is shifted to the rest frame. Exposures with in-transit header times or CCF peak < 0.3 are
  rejected. The rest are combined (inverse-variance, 5σ clip) on a 0.01 Å grid in the format `eblm_feh.read_spectrum`
  reads.
* **fit**: `eblm_feh.fit_spectrum` with the group's resolution and TIC Teff/logg priors. The PHOENIX grid is
  interpolated linearly in Teff and logg and with a natural cubic spline in [Fe/H] (Nov 2026 change: linear
  interpolation across the 0.5 dex [Fe/H] steps trapped ~13% of fits exactly on the 0.0 grid point). The fit stage
  also runs off the cluster: copy `coadd/` (~400 MB) into a local `METAL_ROOT` and run
  `python hpc/metallicity_pipeline.py fit`; paths in the co-add index are re-rooted to the local `coadd/` folder.
* **combine**: first adds the CORALIE fits from the earlier local analysis (`eblm_metallicity_data/
  eblm_spectroscopic_fits.csv`, CORALIE rows only; CORALIE is not in the ESO or SOPHIE archives), matched to EBLM
  targets by the coordinates in the co-add file names (≤ 30″; `results/coralie_matches.csv`). Nov 2026: 375 of 381
  co-adds matched, adding ~100 EBLMs after the cuts; CORALIE agrees with HARPS to +0.00/+0.03 dex (pre/post-2014).
  `combine --no-coralie` leaves them out. Then `eblm_feh.quality_cuts` (plus: fits with a NaN [Fe/H] error are
  rejected) and `combine_per_star`. The
  zero point is fitted against ARES+MOOG/SME values (SWEET-Cat homogeneous hot-Jupiter hosts plus the EBLM
  literature). It depends on Teff, so the adopted one is a robust fit of a line rising up to a fitted break
  temperature and constant above it (Nov 2026 run, spline fits: +0.235 dex above 5660 K, +0.039 dex per 100 K below), held at the
  edge value outside the calibrators' Teff range (`teff_outside_calib`). Outputs: `zero_point_teff.csv` (the curve
  and its bootstrap error), `zero_point_teff_bins.csv` (binned medians, a check), `zero_points.json` (constant zero
  points, pooled and per source, plus the Teff fit). `feh_constzp` keeps the constant-zero-point value. Instrument
  offsets are computed relative to HARPS.

Tested locally end to end on three stars (WASP-43, EBLM J1934-42, EBLM J0040+01). On the same HARPS data the new
co-adds reproduce the existing `eblm_spectroscopic_fits.csv` to 0.01–0.02 dex in [Fe/H] and 10 K in Teff. The PAWS
wrapper's plumbing (stub GUI, file discovery, iSpec patches) was tested with a stand-in iSpec. iSpec itself has not
been run yet, so check the first few `paws/<group>/paws.log` files before trusting the whole batch.

## PAWS notes

`paws_pipeline.py` runs PAWS's own `coadding.py`, `pipeline.py` and `synth.py` from commit `8ec055c`. A stand-in
replaces the ipywidgets GUI, using G2 starting values, 420–680 nm and WIDTH for the EW stage (the BEBOP V set-up). Where
the code and the paper disagree:

* **Per-exposure S/N cut**: the code uses 10, the paper 20. `prepare --snr-min 20` (the default) applies the paper's
  value.
* **[M/H] in the synthesis stage**: the code leaves it free, the paper fixes it at the EW value. The default follows
  the code; `PAWS_ARGS=--fix-mh` follows the paper.
* **Missing line mask for R ≤ 50,000**: the synthesis stage asks for a file iSpec doesn't have
  (`47000_GES/spectrum_ew_…`). It is redirected to `47000_GES/spectrum_synth_good_for_params_all_extended.txt`.
