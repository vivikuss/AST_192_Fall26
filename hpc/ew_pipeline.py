"""[Fe/H] by the ARES+MOOG equivalent-width method (the SWEET-Cat recipe: Sousa 2014; Sousa et al. 2021), run through
FASMA (Andreasen et al. 2017), on archival spectra. Uses the download / RV-alignment / co-add machinery of
metallicity_pipeline.py and writes everything under $METAL_ROOT, next to that pipeline's output.

Stages:

  manifest  ESO spectra of a custom target list (e.g. hpc/ares_moog/validation_targets.csv)
            -> catalogs/targets.csv + catalogs/manifest.csv                        (ast192 env; needs internet)
            then fetch them with:  python hpc/metallicity_pipeline.py download     (same METAL_ROOT)
            (for the main HJ/EBLM sample, skip this: the main pipeline's manifest and raw/ spectra are reused)
  coadd     per (star, instrument, resolution): RV-align and co-add 4450-6900 A, 50 A continuum -> coadd_ew/
                                                                                   (ast192 env)
  fit       ARES equivalent widths -> MOOG excitation/ionisation balance (FASMA) -> ew_fits/
                                                                                   (aresmoog env, $ARESMOOG_HOME)
  collect   ew_fits/*.csv -> results/ew_fits.csv; with --compare <csv>, results/ew_validation.csv
  sun-test  the fit stage on ARES's solar (Ganymede) HARPS spectrum: checks the installation   (aresmoog env)
  rescue    re-fit with logg fixed at the TIC value: every co-add whose free fit failed or ran to logg = 5, plus
            the converged HARPS co-adds as a calibration -> ew_fits_fixlogg/ (then collect again)  (aresmoog env)
  adopt     one [Fe/H] per star: instrument offsets vs HARPS, fixed-logg offset, error floor, flags
            -> results/feh_adopted.csv, feh_candidates.csv, feh_offsets.csv

    python hpc/ew_pipeline.py manifest --targets hpc/ares_moog/validation_targets.csv
    python hpc/ew_pipeline.py coadd --only TOI-426
    python hpc/ew_pipeline.py fit --task 3 --ntasks 8          # or let Slurm array variables set these
    python hpc/ew_pipeline.py collect --compare hpc/ares_moog/validation_targets.csv

Method choices (all FASMA / SWEET-Cat defaults unless noted):
  * line list: Sousa et al. (2008) Fe I + Fe II, with log gf for Kurucz models; below 5200 K FASMA drops the lines
    not in the Tsantaki et al. (2013) cool-star list ('teffrange'). A target can ask for the Tsantaki list from
    the start (column `linelist` in targets.csv), as the TOI-2211 and HD 176986 analyses did.
  * ARES: local continuum from its own S/N estimate (rejt from the noise in 6 line-free windows), lines with
    EW > 200 mA dropped; lines falling on missing data in a co-add are removed from the list first.
  * MOOG 2014 (abfind) on interpolated Kurucz ATLAS9 models; FASMA's minimiser iterates Teff (excitation balance),
    logg (Fe I = Fe II), vt (no trend with reduced EW) to convergence, then: 3-sigma outliers removed once and
    re-converged ('outlier:allOnce'), vt fixed if it runs to the edge ('autofixvt'), stricter re-convergence
    ('refine'). Errors are FASMA's (slope errors propagated, as in Gonzalez & Vanture 1998 / Sousa 2014).
"""
import argparse
import contextlib
import glob
import json
import multiprocessing as mproc
import os
import shutil
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.abspath(os.environ.get('METAL_ROOT', os.path.join(HERE, 'metallicity_hpc_data')))
CAT, RAW, RES = (os.path.join(ROOT, d) for d in ('catalogs', 'raw', 'results'))
EWC = os.path.join(ROOT, 'coadd_ew')       # wide, gently normalised co-adds
WORK = os.path.join(ROOT, 'ew_work')       # one FASMA working directory per co-add
FITS_EW = os.path.join(ROOT, 'ew_fits')    # per-task result tables
FITS_FIX = os.path.join(ROOT, 'ew_fits_fixlogg')   # rescue stage: per-task result tables (logg fixed)

EW_LO, EW_HI = 4450.0, 6900.0   # co-add / ARES range (air, A): covers the whole Sousa+2008 list (4508-6864 A)
NORM_WIDTH = 50.0               # A; running 95th-percentile pre-normalisation (ARES then fits the local continuum)
COVER_EM = (4.5e-7, 6.85e-7)    # ObsCore em_min / em_max (m) a spectrum must span to be useful
MATCH_ARCSEC = 15.0             # cone radius for custom targets (nearby stars move several arcsec over 20 yr)
MAX_PER_GROUP = 100             # highest-S/N exposures co-added per group (custom targets)
LINE_PAD = 3.0                  # A of valid data required on each side of a line
MIN_LINE_FRACTION = 0.6         # skip co-adds covering less than this fraction of the line list

LINELISTS = {'sousa': 'Sousa2007_opt_kurucz.lst', 'tsantaki': 'Tsantaki2013_cool_opt.lst'}
FASMA_SUPPORT = ('Sousa2007_opt_kurucz.lst', 'Tsantaki2013_cool_opt.lst', 'coolNormalDiff.lines')
FASMA_OPTIONS = {'sousa': 'teffrange,autofixvt,refine,outlier:allOnce',
                 'tsantaki': 'autofixvt,refine,outlier:allOnce'}
ARES_OPTIONS = f'lambdai:{EW_LO},lambdaf:{EW_HI},rvmask:0,force,output:star.ares'
SUN = {'teff': 5777, 'logg': 4.44, 'feh': 0.00}

# rescue / adopt
LOGG_RUNAWAY = 4.99             # free-logg fits at the grid edge (logg -> 5) are not usable
REF_INST = 'HARPS_R115000'      # reference for instrument offsets (validated against published ARES+MOOG values)
INST_PRIORITY = ['HARPS_R115000', 'ESPRESSO', 'UVES_R107000', 'HARPS_R80000', 'UVES', 'SOPHIE_HR', 'FEROS',
                 'SOPHIE_HE']   # tie-break only (prefix match); the adopted value is the one with the smallest error
METHOD_FLOOR = 0.03             # dex, added in quadrature to every adopted [Fe/H] (method / line-list systematics)
DEFAULT_INST_SYS = 0.05         # dex, systematic for an instrument with too little overlap with REF_INST to measure
MIN_OVERLAP = 5                 # stars in common with REF_INST needed to measure (and correct) an instrument offset
EXCLUDE = {'HD210027': 'iota Peg, double-lined binary: companion lines in the spectrum'}
# External check of the instrument offsets: SWEET-Cat's homogeneous ARES+MOOG values (SWFlag = 1). Only slow rotators
# are used, because ARES+MOOG [Fe/H] reads low for broad lines (Metallicity_ARES_MOOG.ipynb: -0.18 dex at
# vbroad 10-15 km/s, -0.4 dex above 15 km/s, vs SWEET-Cat), and line broadening is taken from the PHOENIX fits.
SWEETCAT = os.path.join(HERE, 'eblm_metallicity_data', 'sweetcat.csv')
SC_MATCH_ARCSEC = 5.0
VBROAD_FAST = 10.0              # km/s; at or above this a star is flagged (vbroad_ge_10) and kept out of the SWEET-Cat offsets


def _dirs():
    for d in (CAT, RAW, RES, EWC, os.path.join(EWC, 'index'), WORK, FITS_EW, FITS_FIX):
        os.makedirs(d, exist_ok=True)


def _task_args(a):
    task = a.task if a.task is not None else int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
    ntasks = a.ntasks if a.ntasks is not None else int(os.environ.get('SLURM_ARRAY_TASK_COUNT', 1))
    return task, ntasks


def _workers():
    return int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count() or 1))


def slug(name):
    """Same as metallicity_pipeline.slug (copied so the light stages don't import that module)."""
    return ''.join(c if c.isalnum() or c in '+-.' else '_' for c in str(name).replace('EBLM ', 'EBLM_')).strip('_')


def _targets():
    return pd.read_csv(os.path.join(CAT, 'targets.csv')).set_index('target_id', drop=False)


# ----------------------------------------------------------------------------- manifest (custom targets)

def run_manifest(targets_csv, instruments=('HARPS',), max_per_group=MAX_PER_GROUP, force=False):
    """ESO Phase 3 spectra of the stars in `targets_csv` -> catalogs/targets.csv and catalogs/manifest.csv."""
    import pyvo
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import metallicity_pipeline as mp

    tpath, mpath = os.path.join(CAT, 'targets.csv'), os.path.join(CAT, 'manifest.csv')
    if os.path.exists(tpath) and not force:
        old = pd.read_csv(tpath)
        if 'sample' in old and old['sample'].isin(['HJ', 'EBLM']).any():
            sys.exit(f'{tpath} holds the main HJ/EBLM sample. Point METAL_ROOT at a separate directory for custom '
                     f'targets (the main sample needs no manifest stage here), or pass --force.')
    t = pd.read_csv(targets_csv, comment='#')
    if 'target_id' not in t:
        t['target_id'] = t['name'].map(slug)
    t['sample'] = t['sample'] if 'sample' in t else 'VAL'
    for c in ('pl_name', 'pl_tranmid', 'pl_orbper', 'pl_trandur'):   # no ephemeris: nothing rejected as in-transit
        if c not in t:
            t[c] = np.nan

    svc = pyvo.dal.TAPService(mp.ec.ESO_TAP)
    inst = ', '.join(f"'{i}'" for i in instruments)
    parts = []
    for r in t.itertuples():
        dra, ddec = MATCH_ARCSEC / 3600 / np.cos(np.radians(r.dec)), MATCH_ARCSEC / 3600   # (no RA=0 wrap handling)
        q = (f"select dp_id, s_ra, s_dec, obs_collection, snr, em_min, em_max, em_res_power, t_min, access_estsize, "
             f"obs_release_date from ivoa.ObsCore where dataproduct_type = 'spectrum' and obs_collection in ({inst}) "
             f"and s_ra between {r.ra - dra} and {r.ra + dra} and s_dec between {r.dec - ddec} and {r.dec + ddec}")
        e = svc.search(q, maxrec=200_000).to_table().to_pandas()
        sep = SkyCoord(e.s_ra.values, e.s_dec.values, unit='deg').separation(SkyCoord(r.ra, r.dec, unit='deg'))
        released = pd.to_datetime(e.obs_release_date, utc=True, format='ISO8601', errors='coerce') \
            < pd.Timestamp.now(tz='UTC')
        keep = ((sep.arcsec <= MATCH_ARCSEC) & (e.em_min < COVER_EM[0]) & (e.em_max > COVER_EM[1])
                & (e.snr >= mp.SNR_SPEC_MIN) & released).values
        print(f'{r.name}: {len(e)} ESO spectra in the box, {keep.sum()} usable', flush=True)
        parts.append(e[keep].assign(target_id=r.target_id, sep_arcsec=sep.arcsec[keep]))
    e = pd.concat(parts, ignore_index=True)
    man = pd.DataFrame({
        'target_id': e.target_id, 'archive': 'ESO', 'instrument': e.obs_collection,
        'resolution': (e.em_res_power.fillna(0) / 1000).round().astype(int) * 1000,
        'file_id': e.dp_id, 'url': mp.ESO_FILE + e.dp_id, 'snr': e.snr, 'size_mb': e.access_estsize / 1024,
        'sep_arcsec': e.sep_arcsec})
    man = man[man.resolution > 0]
    man['group'] = man.target_id + '__' + man.instrument + '_R' + man.resolution.astype(str)
    man = man.sort_values('snr', ascending=False)
    man = man[man.groupby('group').cumcount() < max_per_group]
    man['path'] = [os.path.join(RAW, 'ESO', r.target_id, r.file_id.replace(':', '_') + '.fits')
                   for r in man.itertuples()]
    man = man.sort_values(['target_id', 'group', 'snr'], ascending=[True, True, False]).reset_index(drop=True)
    _dirs()
    t.to_csv(tpath, index=False)
    man.to_csv(mpath, index=False)
    summary = man.groupby('group').agg(n=('file_id', 'size'), snr_best=('snr', 'max'),
                                       snr_coadd=('snr', lambda x: np.sqrt(np.sum(x ** 2))))
    print(summary.round(0).to_string())
    print(f'{len(man)} spectra (~{man.size_mb.sum() / 1024:.1f} GB) in {man.group.nunique()} groups -> {mpath}')
    missing = set(t.target_id) - set(man.target_id)
    if missing:
        print(f'no usable spectra for: {", ".join(sorted(missing))}')


# ----------------------------------------------------------------------------- co-add (wide range)

def _coadd_worker(args):
    import metallicity_pipeline as mp
    rows, target = args
    try:
        return mp.coadd_group(rows, target, mp._grid(), lo=EW_LO, hi=EW_HI, norm_width=NORM_WIDTH, out_root=EWC,
                              note=f'Rest-frame co-add for ARES+MOOG, {NORM_WIDTH:.0f} A pre-normalisation; '
                                   'ew_pipeline.py')
    except Exception as e:
        return {'group': rows.iloc[0].group, 'target_id': rows.iloc[0].target_id, 'file': None, 'error': repr(e)}


def run_coadd(task, ntasks, only=None, instruments=None):
    import metallicity_pipeline as mp   # noqa: F401  (imported here: it needs the ast192 env)
    _dirs()
    man = pd.read_csv(os.path.join(CAT, 'manifest.csv'))
    t = _targets()
    if only:
        man = man[man.target_id.isin(only)]
    if instruments:
        man = man[man.instrument.isin(instruments)]
    groups = sorted(man.group.unique())[task::ntasks]
    jobs = [(man[man.group == g], t.loc[g.split('__')[0]]) for g in groups]
    with mproc.Pool(_workers()) as pool:
        out = pool.map(_coadd_worker, jobs, chunksize=1)
    pd.DataFrame(out).to_csv(os.path.join(EWC, 'index', f'coadd_{task:03d}.csv'), index=False)
    ok = [o for o in out if o.get('file') is not None]
    print(f'task {task}/{ntasks}: {len(ok)}/{len(out)} groups co-added')
    for o in out:
        print(f"  {o['group']}: {o.get('n_used', 0)}/{o.get('n_files', '?')} exposures used, "
              f"S/N {o.get('snr_coadd', float('nan')):.0f}" + (f"  ERROR {o['error']}" if o.get('error') else ''))


def _coadd_index():
    files = sorted(glob.glob(os.path.join(EWC, 'index', 'coadd_*.csv')))
    if not files:
        sys.exit(f'no co-add index in {EWC}/index; run the coadd stage first')
    idx = pd.concat([pd.read_csv(f) for f in files], ignore_index=True).drop_duplicates('group', keep='last')
    idx = idx[idx.file.notna()].copy()
    # re-root paths written on another machine (e.g. co-adds copied from the cluster to a laptop)
    idx['file'] = [p if os.path.exists(p) else os.path.join(EWC, p[p.rfind('/coadd_ew/') + 10:])
                   if '/coadd_ew/' in p else p for p in idx.file]
    return idx


# ----------------------------------------------------------------------------- fit (ARES + MOOG through FASMA)

def _aresmoog_paths():
    home = os.environ.get('ARESMOOG_HOME')
    if not home:
        sys.exit('set ARESMOOG_HOME (the directory hpc/ares_moog/setup_ares_moog.sh installed into)')
    fasma, bindir = os.path.join(home, 'FASMA'), os.path.join(home, 'bin')
    for p in (os.path.join(bindir, 'ARES'), os.path.join(bindir, 'MOOGSILENT'), os.path.join(fasma, 'ewDriver.py'),
              os.path.join(fasma, 'models', 'kurucz95')):
        if not os.path.exists(p):
            sys.exit(f'{p} missing; run hpc/ares_moog/setup_ares_moog.sh')
    return home, fasma, bindir


def _covered_linelist(spec, src, dst, pad=LINE_PAD):
    """Copy line list `src` to `dst`, keeping only lines with valid data within +-pad A in spectrum `spec`.

    ARES crashes or measures garbage on lines that fall on gaps (zeros in a co-add, the ends of a spectrum)."""
    from astropy.io import fits
    with fits.open(spec) as h:
        flux = np.asarray(h[0].data, float).ravel()
        w0, dw = h[0].header['CRVAL1'], h[0].header['CDELT1']
    good = np.isfinite(flux) & (flux > 0)
    cum = np.concatenate([[0], np.cumsum(~good)])   # number of bad pixels before each index
    with open(src) as f:
        lines = f.readlines()
    head, body = lines[:2], [ln for ln in lines[2:] if ln.strip()]
    keep = []
    for ln in body:
        wl = float(ln.split()[0])
        i0, i1 = int(np.floor((wl - pad - w0) / dw)), int(np.ceil((wl + pad - w0) / dw)) + 1
        if i0 >= 0 and i1 <= flux.size and cum[i1] - cum[i0] == 0:
            keep.append(ln)
    with open(dst, 'w') as f:
        f.writelines(head + keep)
    return len(body), len(keep)


def _initial_vt(teff, logg):
    """FASMA's starting microturbulence: Tsantaki+2013 (dwarfs) / Adibekyan+2015 (giants)."""
    vt = 6.932e-4 * teff - 0.348 * logg - 1.437 if logg >= 3.95 else 2.72 - 0.457 * logg
    return float(np.clip(vt, 0.3, 3.0))


def _count_fe(path):
    d = np.loadtxt(path, skiprows=1, usecols=(1,), ndmin=1)
    return int(np.sum(np.isclose(d, 26.0))), int(np.sum(np.isclose(d, 26.1)))


def fit_one(job):
    """ARES then FASMA's EW minimisation for one spectrum, in its own working directory (FASMA works in the cwd).

    Runs in a fresh worker process (maxtasksperchild=1): FASMA keeps module-level state and logging handlers."""
    t0 = time.time()
    home, fasma, bindir = job['home'], job['fasma'], job['bindir']
    wd, spec, key = job['workdir'], job['file'], job['linelist']
    ares_opts = job.get('ares_options') or ARES_OPTIONS          # overrides: hpc/ares_moog/sophie_test.py
    fasma_opts = job.get('fasma_options') or FASMA_OPTIONS[key]
    res = {k: job.get(k) for k in ('group', 'target_id', 'file', 'snr_coadd', 'n_used', 'variant', 'role')}
    res.update(linelist=key, status='ok', ares_options=ares_opts, fasma_options=fasma_opts)
    os.environ['PATH'] = bindir + os.pathsep + os.environ.get('PATH', '')
    os.environ['OMP_NUM_THREADS'] = '1'
    try:
        shutil.rmtree(wd, ignore_errors=True)
        for d in ('linelist', 'results', 'rawLinelist'):
            os.makedirs(os.path.join(wd, d))
        os.symlink(os.path.join(fasma, 'models'), os.path.join(wd, 'models'))
        for f in FASMA_SUPPORT:   # FASMA reads these from ./rawLinelist (log gf update, cool-star switch)
            shutil.copy(os.path.join(fasma, 'rawLinelist', f), os.path.join(wd, 'rawLinelist', f))
        n_list, n_cov = _covered_linelist(spec, os.path.join(fasma, 'rawLinelist', LINELISTS[key]),
                                          os.path.join(wd, 'rawLinelist', 'star_lines.lst'))
        res.update(n_lines_list=n_list, n_lines_covered=n_cov)
        if n_cov < MIN_LINE_FRACTION * n_list:
            res['status'] = 'poor_coverage'
            return res

        os.chdir(wd)
        sys.path.insert(0, fasma)
        import aresDriver
        import ewDriver
        with open('fasma.log', 'w') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            with open('StarMe_ares.cfg', 'w') as f:
                f.write(f'star_lines.lst {spec} {ares_opts}\n')
            res['snr_ares'] = aresDriver.aresdriver('StarMe_ares.cfg')
            if not os.path.exists('linelist/star.moog'):
                res['status'] = 'ares_failed'
                return res
            res['n_fe1_ares'], res['n_fe2_ares'] = _count_fe('linelist/star.moog')

            teff0 = job['teff0'] if np.isfinite(job['teff0']) else SUN['teff']
            logg0 = job['logg0'] if np.isfinite(job['logg0']) else SUN['logg']
            if 'logg' not in fasma_opts.split(','):   # a fixed logg ('logg' option) is used as given
                teff0, logg0 = float(np.clip(teff0, 4000, 7000)), float(np.clip(logg0, 3.0, 4.9))
            teff0 = float(np.clip(teff0, 4000, 7000))
            res.update(teff0=teff0, logg0=logg0)
            with open('StarMe_ew.cfg', 'w') as f:
                f.write(f'star.moog {int(round(teff0))} {logg0:.2f} 0.00 {_initial_vt(teff0, logg0):.2f} '
                        f'{fasma_opts}\n')
            ewDriver.EWmethod(cfgfile='StarMe_ew.cfg', overwrite=True).ewdriver()

        out = pd.read_csv('EWresults.dat', sep='\t').iloc[-1]
        res.update(teff=out.teff, teff_err=out.tefferr, logg=out.logg, logg_err=out.loggerr, feh=out.feh,
                   feh_err=out.feherr, vt=out.vt, vt_err=out.vterr, converged=str(out.convergence) == 'True',
                   final_linelist=out.linelist)
        final = os.path.join('linelist', out.linelist)
        if os.path.exists(final):
            res['n_fe1'], res['n_fe2'] = _count_fe(final)
        if not res['converged']:
            res['status'] = 'not_converged'
    except Exception as e:
        res['status'] = f'failed: {e!r}'[:300]
    finally:
        res['seconds'] = round(time.time() - t0, 1)
        with open(os.path.join(wd, 'result.json'), 'w') as f:
            json.dump({k: (v.item() if hasattr(v, 'item') else v) for k, v in res.items()}, f, indent=1)
    return res


def _run_fits(jobs):
    ctx = mproc.get_context('fork' if sys.platform.startswith('linux') else 'spawn')
    with ctx.Pool(max(1, min(_workers(), len(jobs))), maxtasksperchild=1) as pool:
        return pool.map(fit_one, jobs, chunksize=1)


def run_fit(task, ntasks, only=None, linelist=None):
    home, fasma, bindir = _aresmoog_paths()
    _dirs()
    idx = _coadd_index()
    t = _targets()
    if only:
        idx = idx[idx.target_id.isin(only)]
    idx = idx.sort_values('group').iloc[task::ntasks]
    jobs = []
    for r in idx.itertuples():
        tg = t.loc[r.target_id]
        key = linelist or (tg.get('linelist') if isinstance(tg.get('linelist'), str) else 'sousa')
        if key not in LINELISTS:
            sys.exit(f'{r.target_id}: unknown linelist {key!r} (choose from {", ".join(LINELISTS)})')
        jobs.append(dict(group=r.group, target_id=r.target_id, file=r.file, snr_coadd=r.snr_coadd, n_used=r.n_used,
                         teff0=float(tg.get('tic_teff', np.nan)), logg0=float(tg.get('tic_logg', np.nan)),
                         linelist=key, workdir=os.path.join(WORK, r.group), home=home, fasma=fasma, bindir=bindir))
    print(f'task {task}/{ntasks}: fitting {len(jobs)} co-adds with {min(_workers(), max(len(jobs), 1))} workers',
          flush=True)
    out = pd.DataFrame(_run_fits(jobs)) if jobs else pd.DataFrame()
    out.to_csv(os.path.join(FITS_EW, f'ew_fits_{task:03d}.csv'), index=False)
    if len(out):
        print(out.status.value_counts().to_string())
        cols = [c for c in ('group', 'status', 'teff', 'logg', 'feh', 'feh_err', 'vt', 'n_fe1', 'n_fe2', 'seconds')
                if c in out]
        print(out[cols].to_string(index=False))


def _inst(group):
    return group.str.split('__').str[1]


def run_rescue(task, ntasks, only=None, calib=True):
    """Re-fit with logg FIXED at the TIC value (FASMA option 'logg'; Teff, vt, [Fe/H] free):
      role=rescue  every co-add whose free fit did not converge or ran to logg = 5 (results/ew_fits.csv)
      role=calib   every converged REF_INST (HARPS) fit, to measure what fixing logg does to [Fe/H] (adopt stage)
    Same co-adds, line lists and options as the fit stage otherwise. -> ew_fits_fixlogg/ (collect merges them)."""
    home, fasma, bindir = _aresmoog_paths()
    _dirs()
    path = os.path.join(RES, 'ew_fits.csv')
    if not os.path.exists(path):
        sys.exit(f'{path} not found: run fit + collect first')
    f = pd.read_csv(path)
    st = f.status.astype(str)
    bad = ((st == 'not_converged') | st.str.startswith('failed') | ((st == 'ok') & (f.logg >= LOGG_RUNAWAY)))
    good_ref = (st == 'ok') & (f.logg < LOGG_RUNAWAY) & (_inst(f.group) == REF_INST)
    f['role'] = np.where(bad, 'rescue', np.where(good_ref & calib, 'calib', ''))
    f = f[f.role != '']
    if only:
        f = f[f.target_id.isin(only)]
    t = _targets()
    f = f.assign(tic_logg=f.target_id.map(t['tic_logg']), tic_teff=f.target_id.map(t['tic_teff']))
    no_logg = ~np.isfinite(f.tic_logg)
    if no_logg.any():
        print(f'{int(no_logg.sum())} co-adds skipped: no TIC logg ({", ".join(f[no_logg].group.head(5))}...)')
    f = f[~no_logg].sort_values('group').iloc[task::ntasks]
    jobs = []
    for r in f.itertuples():
        key = r.linelist if isinstance(r.linelist, str) and r.linelist in LINELISTS else 'sousa'
        jobs.append(dict(group=r.group, target_id=r.target_id, file=r.file, snr_coadd=r.snr_coadd, n_used=r.n_used,
                         variant='fixlogg', role=r.role, teff0=float(r.tic_teff),
                         logg0=float(np.clip(r.tic_logg, 3.0, 5.0)), linelist=key,
                         fasma_options=FASMA_OPTIONS[key] + ',logg',
                         workdir=os.path.join(WORK, r.group + '__fixlogg'), home=home, fasma=fasma, bindir=bindir))
    print(f'task {task}/{ntasks}: {len(jobs)} fixed-logg fits '
          f'({sum(j["role"] == "rescue" for j in jobs)} rescue, {sum(j["role"] == "calib" for j in jobs)} calib), '
          f'{min(_workers(), max(len(jobs), 1))} workers', flush=True)
    out = pd.DataFrame(_run_fits(jobs)) if jobs else pd.DataFrame()
    os.makedirs(FITS_FIX, exist_ok=True)
    out.to_csv(os.path.join(FITS_FIX, f'ew_fits_{task:03d}.csv'), index=False)
    if len(out):
        print(out.groupby('role').status.apply(lambda x: x.astype(str).str.split(':').str[0].value_counts().to_dict())
              .to_string())


def run_sun_test():
    """ARES+MOOG on the solar spectrum shipped with ARES (HARPS, via Ganymede). Pass = solar parameters back."""
    home, fasma, bindir = _aresmoog_paths()
    _dirs()
    job = dict(group='SUN_harps_ganymede', target_id='Sun', file=os.path.join(home, 'ARES', 'sun_harps_ganymede.fits'),
               snr_coadd=np.nan, n_used=1, teff0=5600.0, logg0=4.30, linelist='sousa',
               workdir=os.path.join(WORK, '_sun_test'), home=home, fasma=fasma, bindir=bindir)
    r = _run_fits([job])[0]
    print(json.dumps({k: (v.item() if hasattr(v, 'item') else v) for k, v in r.items()}, indent=1, default=str))
    if r['status'] != 'ok':
        sys.exit(f"sun test FAILED ({r['status']}); see {job['workdir']}/fasma.log and captain.log")
    d = {k: r[k] - SUN[k] for k in SUN}
    ok = abs(d['teff']) <= 60 and abs(d['logg']) <= 0.10 and abs(d['feh']) <= 0.03
    print(f"Sun: Teff {r['teff']:.0f} ({d['teff']:+.0f} K), logg {r['logg']:.2f} ({d['logg']:+.2f}), "
          f"[Fe/H] {r['feh']:+.3f} ({d['feh']:+.3f}), vt {r['vt']:.2f} -> {'PASS' if ok else 'CHECK'}")
    if not ok:
        sys.exit('sun test: parameters off the solar values by more than 60 K / 0.1 dex / 0.03 dex')


# ----------------------------------------------------------------------------- collect / compare

def run_collect(compare=None):
    files = sorted(glob.glob(os.path.join(FITS_EW, 'ew_fits_*.csv')))
    if not files:
        sys.exit(f'no fit results in {FITS_EW}; run the fit stage first')
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True).drop_duplicates('group', keep='last')
    t = _targets()
    df = df.merge(t[['target_id', 'name', 'sample']].reset_index(drop=True), on='target_id', how='left')
    os.makedirs(RES, exist_ok=True)
    df.to_csv(os.path.join(RES, 'ew_fits.csv'), index=False)
    print(f'{len(df)} co-adds -> {os.path.join(RES, "ew_fits.csv")}')
    print(df.status.str.split(':').str[0].value_counts().to_string())
    fix = sorted(glob.glob(os.path.join(FITS_FIX, 'ew_fits_*.csv')))
    if fix:
        dx = pd.concat([pd.read_csv(f) for f in fix], ignore_index=True).drop_duplicates('group', keep='last')
        dx = dx.merge(t[['target_id', 'name', 'sample']].reset_index(drop=True), on='target_id', how='left')
        dx.to_csv(os.path.join(RES, 'ew_fits_fixlogg.csv'), index=False)
        print(f'{len(dx)} fixed-logg fits -> {os.path.join(RES, "ew_fits_fixlogg.csv")}')
        print(dx.groupby('role').status.apply(lambda x: x.astype(str).str.split(':').str[0].value_counts().to_dict())
              .to_string())
    if not compare:
        return

    c = pd.read_csv(compare, comment='#')
    if 'target_id' not in c:
        c['target_id'] = c['name'].map(slug)
    ok = df[df.status == 'ok'].sort_values('snr_coadd', ascending=False)
    best = ok.drop_duplicates('target_id')   # highest-S/N converged co-add per star
    v = c.merge(best, on='target_id', how='left', suffixes=('', '_fit'))
    v['d_feh'] = v.feh - v.pub_feh
    v['d_feh_sigma'] = v.d_feh / np.hypot(v.feh_err, v.pub_feh_err)
    for k in ('teff', 'logg'):
        if f'pub_{k}' in v:
            v[f'd_{k}'] = v[k] - v[f'pub_{k}']
    v.to_csv(os.path.join(RES, 'ew_validation.csv'), index=False)
    cols = ['name', 'group', 'n_used', 'snr_coadd', 'linelist', 'teff', 'logg', 'vt', 'feh', 'feh_err', 'pub_feh',
            'pub_feh_err', 'd_feh', 'd_feh_sigma', 'd_teff', 'd_logg', 'n_fe1', 'n_fe2']
    print(v[[k for k in cols if k in v]].round(3).to_string(index=False))
    d = v.d_feh.dropna()
    if len(d):
        print(f'[Fe/H] ours - published: mean {d.mean():+.3f}, rms {np.sqrt(np.mean(d ** 2)):.3f} dex (N={len(d)}); '
              f'max |diff| / combined error = {np.nanmax(np.abs(v.d_feh_sigma)):.1f}')
    print(f'-> {os.path.join(RES, "ew_validation.csv")}')


# ----------------------------------------------------------------------------- adopt (one [Fe/H] per star)

def _priority(inst):
    for i, p in enumerate(INST_PRIORITY):
        if str(inst).startswith(p):
            return i
    return len(INST_PRIORITY)


def _offset(d, err2=None):
    """Median offset, robust scatter, its error, and the intrinsic scatter (robust scatter minus the quoted errors)."""
    d = np.asarray(d, float)
    n = len(d)
    if not n:
        return dict(n=0, offset=np.nan, scatter=np.nan, offset_err=np.nan, sys=np.nan)
    med = float(np.median(d))
    sc = float(1.4826 * np.median(np.abs(d - med))) if n > 1 else np.nan
    sysv = sc if err2 is None else float(np.sqrt(max(sc ** 2 - np.median(err2), 0.0))) if n > 1 else np.nan
    return dict(n=n, offset=med, scatter=sc, offset_err=1.2533 * sc / np.sqrt(n) if n > 1 else np.nan, sys=sysv)


def _radec_deg(ra, dec):
    """SWEET-Cat's sexagesimal 'hh mm ss.s' / '+dd mm ss.s' -> degrees (numpy only: the aresmoog env has no astropy)."""
    def one(s, scale):
        p = str(s).split()
        sign = -1.0 if p[0].startswith('-') else 1.0
        v = abs(float(p[0])) + float(p[1]) / 60 + float(p[2]) / 3600
        return sign * v * scale
    return (np.array([one(s, 15.0) for s in ra]), np.array([one(s, 1.0) for s in dec]))


def _unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def _sweetcat_match(t):
    """Homogeneous SWEET-Cat [Fe/H] for each target within SC_MATCH_ARCSEC (NaN where none)."""
    out = pd.DataFrame(index=t.index, columns=['sc_feh', 'sc_feh_err'], dtype=float)
    if not os.path.exists(SWEETCAT):
        print(f'note: {SWEETCAT} not found: instrument offsets from HARPS pairs only')
        return out
    sc = pd.read_csv(SWEETCAT).dropna(subset=['RA', 'DEC', '[Fe/H]'])
    sc = sc[sc.SWFlag == 1]
    u_sc = _unit(*_radec_deg(sc.RA, sc.DEC))
    u_t = _unit(t.ra.values, t.dec.values)
    best = np.empty(len(t), int)
    cosd = np.empty(len(t))
    for i in range(0, len(t), 500):                   # chunks keep the dot-product matrix small
        dots = u_t[i:i + 500] @ u_sc.T
        best[i:i + 500] = dots.argmax(axis=1)
        cosd[i:i + 500] = dots.max(axis=1)
    hit = np.degrees(np.arccos(np.clip(cosd, -1, 1))) * 3600 < SC_MATCH_ARCSEC
    out.loc[hit, 'sc_feh'] = sc['[Fe/H]'].values[best[hit]]
    out.loc[hit, 'sc_feh_err'] = sc['e[Fe/H]'].values[best[hit]]
    return out


def _vbroad():
    """Line broadening (rotation + macroturbulence, km/s) per star from the PHOENIX fits (metallicity_pipeline.py)."""
    p = os.path.join(RES, 'metallicity_final.csv')
    if not os.path.exists(p):
        print(f'note: {p} not found: no vbroad, so no fast-rotator flag and no rotation cut in the SWEET-Cat offsets')
        return pd.Series(dtype=float)
    return pd.read_csv(p).drop_duplicates('target_id').set_index('target_id').vbroad


def _combine(a, b):
    """Inverse-variance mean of two offset estimates (dicts from _offset); either may be missing (None)."""
    est = [e for e in (a, b) if e is not None and np.isfinite(e['offset_err']) and e['offset_err'] > 0]
    if not est:
        return np.nan, np.nan
    w = np.array([1 / e['offset_err'] ** 2 for e in est])
    return float(np.sum(w * [e['offset'] for e in est]) / w.sum()), float(1 / np.sqrt(w.sum()))


def run_adopt():
    """One [Fe/H] per star from all converged fits (free logg, and fixed logg from the rescue stage).

    1. Usable fits: status ok, finite [Fe/H] and error; free-logg fits must not have run to logg = 5. Stars in
       EXCLUDE (e.g. double-lined binaries) are dropped.
    2. Instrument offsets vs REF_INST (HARPS), measured two independent ways from free-logg fits:
       (a) pairs: stars with fits from both this instrument and HARPS, median difference;
       (b) SWEET-Cat: median (this work - SWEET-Cat homogeneous) for this instrument minus the same for HARPS,
           slow rotators only (vbroad < VBROAD_FAST). This is needed because the pairs can be few (10 for SOPHIE HE,
           which gave -0.095 dex, against -0.15 from 43 SWEET-Cat stars) and SOPHIE's overlap stars are hotter.
       An estimate counts if it has >= MIN_OVERLAP stars; the offset is the inverse-variance mean of those that
       count, applied (subtracted) when |offset| > 2 x its error. The intrinsic scatter of the pairs (robust scatter
       minus quoted errors in quadrature) is the instrument's systematic error, or the SWEET-Cat one if there are
       too few pairs (DEFAULT_INST_SYS if neither has enough stars; 0 for REF_INST).
    3. Fixed-logg offset: fixed minus free [Fe/H] on the same HARPS co-adds (rescue stage, role=calib); applied the
       same way, its scatter added as a systematic to every fixed-logg value.
    4. Error = sqrt(FASMA error^2 + instrument sys^2 + fixed-logg sys^2 + METHOD_FLOOR^2); the adopted value is the
       candidate with the smallest error (ties: INST_PRIORITY, then S/N).
    -> results/feh_adopted.csv (one row per star), feh_candidates.csv (every usable fit), feh_offsets.csv"""
    t = _targets()
    main_p, fix_p = os.path.join(RES, 'ew_fits.csv'), os.path.join(RES, 'ew_fits_fixlogg.csv')
    if not os.path.exists(main_p):
        sys.exit(f'{main_p} not found: run collect first')
    free = pd.read_csv(main_p).assign(method='free')
    parts = [free]
    if os.path.exists(fix_p):
        parts.append(pd.read_csv(fix_p).assign(method='fixlogg'))
    else:
        print(f'note: {fix_p} not found (rescue stage not run): free-logg fits only')
    a = pd.concat(parts, ignore_index=True)
    a['inst'] = _inst(a.group)
    a['role'] = a['role'] if 'role' in a else np.nan
    usable = ((a.status.astype(str) == 'ok') & np.isfinite(a.feh) & np.isfinite(a.feh_err) & (a.feh_err > 0)
              & ((a.method == 'fixlogg') | (a.logg < LOGG_RUNAWAY)))
    a = a[usable].copy()
    excluded = a.target_id.isin(EXCLUDE)
    a = a[~excluded]

    # ---- instrument offsets (free-logg fits only): (a) HARPS pairs, (b) SWEET-Cat relative to HARPS
    vb = _vbroad()
    scm = _sweetcat_match(t)
    fr = a[a.method == 'free'].copy()
    fr['vbroad'] = fr.target_id.map(vb)
    fr = fr.join(scm, on='target_id')
    ref = fr[fr.inst == REF_INST].set_index('target_id')

    def _vs_ref(g):
        m = g.set_index('target_id').join(ref[['feh', 'feh_err']], rsuffix='_ref', how='inner')
        return _offset(m.feh - m.feh_ref, m.feh_err ** 2 + m.feh_err_ref ** 2)

    def _vs_sc(g):            # this work - SWEET-Cat, slow rotators with a homogeneous SWEET-Cat value
        s = g[np.isfinite(g.sc_feh) & (g.vbroad < VBROAD_FAST)]
        return _offset(s.feh - s.sc_feh, s.feh_err ** 2 + s.sc_feh_err.fillna(0) ** 2)

    sc_ref = _vs_sc(fr[fr.inst == REF_INST])
    rows = [dict(kind='instrument', inst=REF_INST, basis='reference', n=len(ref), offset=0.0, scatter=0.0,
                 offset_err=0.0, sys=0.0, applied=False, offset_pairs=0.0, offset_pairs_err=0.0,
                 n_sc=sc_ref['n'], offset_sc=0.0, offset_sc_err=0.0)]
    print(f"{REF_INST} - SWEET-Cat (slow rotators): {sc_ref['offset']:+.3f} +- {sc_ref['offset_err']:.3f} dex "
          f"from {sc_ref['n']} stars (the zero point of the SWEET-Cat offsets below)")

    others = fr[fr.inst != REF_INST].assign(family=lambda d: d.inst.str.split('_R').str[0])
    for inst, g in others.groupby('inst'):
        fam = others[others.family == g.family.iloc[0]]
        pooled = fam.inst.nunique() > 1
        o, basis = _vs_ref(g), 'pairs: this mode'
        if o['n'] < MIN_OVERLAP and pooled:    # too few stars in common: pool the instrument's other modes (e.g. all UVES)
            o, basis = _vs_ref(fam), f'pairs: all {g.family.iloc[0]} modes'
        s, sbasis = _vs_sc(g), 'SWEET-Cat: this mode'
        if s['n'] < MIN_OVERLAP and pooled:
            s, sbasis = _vs_sc(fam), f'SWEET-Cat: all {g.family.iloc[0]} modes'
        s_rel = None
        if s['n'] >= MIN_OVERLAP and sc_ref['n'] >= MIN_OVERLAP:
            s_rel = dict(s, offset=s['offset'] - sc_ref['offset'], offset_err=float(np.hypot(s['offset_err'], sc_ref['offset_err'])))
        p_ok = o if o['n'] >= MIN_OVERLAP else None
        off_v, off_e = _combine(p_ok, s_rel)
        used = [b for b, e in ((basis, p_ok), (sbasis, s_rel)) if e is not None]
        if p_ok is not None and s_rel is not None:
            tension = abs(p_ok['offset'] - s_rel['offset']) / np.hypot(p_ok['offset_err'], s_rel['offset_err'])
            if tension > 2:
                print(f'note: {inst}: pairs ({p_ok["offset"]:+.3f}) and SWEET-Cat ({s_rel["offset"]:+.3f}) offsets '
                      f'differ by {tension:.1f} sigma')
        sys_v = o['sys'] if p_ok is not None else (s['sys'] if s_rel is not None else DEFAULT_INST_SYS)
        rows.append(dict(kind='instrument', inst=inst, basis=' + '.join(used) if used else 'too few (default sys)',
                         n=o['n'], offset=off_v, offset_err=off_e, scatter=o['scatter'], sys=sys_v,
                         applied=bool(used and abs(off_v) > 2 * off_e),
                         offset_pairs=o['offset'], offset_pairs_err=o['offset_err'], n_sc=s['n'],
                         offset_sc=s_rel['offset'] if s_rel else np.nan, offset_sc_err=s_rel['offset_err'] if s_rel else np.nan))
    # ---- fixed-logg offset (same HARPS co-adds, fixed minus free)
    cal = a[(a.method == 'fixlogg') & (a.role == 'calib')].set_index('group')
    m = cal.join(fr.set_index('group')[['feh']], rsuffix='_free', how='inner')
    o = _offset(m.feh - m.feh_free)
    o.update(kind='fixed_logg', inst='(any)', basis='HARPS fixed - free', applied=bool(o['n'] >= MIN_OVERLAP and abs(o['offset']) > 2 * o['offset_err']),
             sys=o['scatter'] if o['n'] >= MIN_OVERLAP else DEFAULT_INST_SYS)
    rows.append(o)
    off = pd.DataFrame(rows)[['kind', 'inst', 'basis', 'n', 'offset', 'offset_err', 'scatter', 'sys', 'applied',
                              'offset_pairs', 'offset_pairs_err', 'n_sc', 'offset_sc', 'offset_sc_err']]
    off.to_csv(os.path.join(RES, 'feh_offsets.csv'), index=False)

    # ---- corrected values and errors for every candidate (rescue fits only; calib fits duplicate a free fit)
    c = a[(a.method == 'free') | (a.role == 'rescue')].copy()
    io = off[off.kind == 'instrument'].set_index('inst')
    c['inst_offset'] = c.inst.map(io.offset.where(io.applied, 0.0)).fillna(0.0)
    c['inst_sys'] = c.inst.map(io.sys).fillna(DEFAULT_INST_SYS)
    fx = off[off.kind == 'fixed_logg'].iloc[0]
    is_fix = c.method == 'fixlogg'
    c['logg_offset'] = np.where(is_fix & fx.applied, fx.offset, 0.0)
    c['logg_sys'] = np.where(is_fix, fx.sys, 0.0)
    c['feh_raw'] = c.feh
    c['feh'] = c.feh_raw - c.inst_offset - c.logg_offset
    c['feh_err_total'] = np.sqrt(c.feh_err ** 2 + c.inst_sys ** 2 + c.logg_sys ** 2 + METHOD_FLOOR ** 2)
    c['priority'] = c.inst.map(_priority)
    c = c.sort_values(['target_id', 'feh_err_total', 'priority', 'snr_coadd'], ascending=[True, True, True, False])
    c.to_csv(os.path.join(RES, 'feh_candidates.csv'), index=False)

    best = c.drop_duplicates('target_id').set_index('target_id')
    # tension: another candidate disagrees with the adopted value by > 3 sigma
    j = c.join(best[['feh', 'feh_err_total', 'group']], on='target_id', rsuffix='_adopt')
    j = j[j.group != j.group_adopt]
    tens = j[np.abs(j.feh - j.feh_adopt) > 3 * np.hypot(j.feh_err_total, j.feh_err_total_adopt)].target_id.unique()
    n_cand = c.groupby('target_id').size()

    def flags(r):
        f = []
        if r.method == 'fixlogg':
            f.append('fixed_logg')
        if r.inst.startswith('SOPHIE_HE'):
            f.append('sophie_he')
        if not (r.snr_coadd >= 50):
            f.append('snr_lt_50')
        if not (4500 <= r.teff <= 6500):
            f.append('teff_outside_4500_6500')
        if r.vt >= 2.5:
            f.append('vt_ge_2.5')
        if getattr(r, 'n_fe2', np.nan) < 5:
            f.append('n_fe2_lt_5')
        if r.vbroad >= VBROAD_FAST:      # broad lines: ARES+MOOG reads low (NaN vbroad -> not flagged)
            f.append('vbroad_ge_10')
        if r.Index in tens:
            f.append('instrument_tension')
        return ';'.join(f)

    best['vbroad'] = vb.reindex(best.index).values if len(vb) else np.nan
    best['flags'] = [flags(r) for r in best.itertuples()] if len(best) else []
    best['n_candidates'] = n_cand.reindex(best.index)
    best = best.join(t[['tic_teff', 'tic_logg']], how='left')
    cols = ['name', 'sample', 'feh', 'feh_err_total', 'feh_err', 'feh_raw', 'inst', 'method', 'group', 'teff',
            'teff_err', 'logg', 'logg_err', 'vt', 'vt_err', 'snr_coadd', 'n_used', 'n_fe1', 'n_fe2', 'linelist',
            'inst_offset', 'logg_offset', 'vbroad', 'n_candidates', 'flags', 'tic_teff', 'tic_logg']
    best = best[[k for k in cols if k in best]].rename(columns={'feh_err_total': 'feh_err', 'feh_err': 'feh_err_fasma'})
    best.reset_index().to_csv(os.path.join(RES, 'feh_adopted.csv'), index=False)

    # ---- summary
    print('Offsets vs', REF_INST, '(applied = subtracted; sys = intrinsic scatter used as systematic error):')
    print(off.round(3).to_string(index=False))
    allt = t[['sample']].copy()
    fitted = set(pd.read_csv(main_p).target_id)
    allt['status'] = np.where(allt.index.isin(best.index), 'adopted',
                     np.where(allt.index.isin(EXCLUDE), 'excluded',
                     np.where(allt.index.isin(fitted), 'no usable fit', 'no spectra')))
    print('\nStars by sample:')
    print(allt.groupby(['sample', 'status']).size().unstack(fill_value=0).to_string())
    print('\nAdopted values by instrument and method:')
    print(best.groupby(['sample', 'inst', 'method']).size().unstack(fill_value=0).to_string())
    fl = best['flags'].str.split(';').explode()
    print('\nFlags (stars):')
    print(fl[fl != ''].groupby(best['sample'].reindex(fl[fl != ''].index)).value_counts().unstack(0, fill_value=0)
          .to_string())
    print(f'\nmedian adopted [Fe/H] error: {best.feh_err.median():.3f} dex '
          f'(free logg {best[best.method == "free"].feh_err.median():.3f}, '
          f'fixed logg {best[best.method == "fixlogg"].feh_err.median():.3f})')
    for k, why in EXCLUDE.items():
        print(f'excluded {k}: {why}')
    print(f'-> {RES}/feh_adopted.csv, feh_candidates.csv, feh_offsets.csv')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('stage', choices=['manifest', 'coadd', 'fit', 'collect', 'sun-test', 'rescue', 'adopt'])
    p.add_argument('--no-calib', action='store_true',
                   help='rescue: skip the fixed-logg calibration fits on converged HARPS co-adds')
    p.add_argument('--targets', help='manifest: CSV of custom targets (name, ra, dec, tic_teff, tic_logg, ...)')
    p.add_argument('--instruments', default=None,
                   help='comma-separated instruments (manifest default: HARPS; coadd default: all in the manifest)')
    p.add_argument('--max-per-group', type=int, default=MAX_PER_GROUP)
    p.add_argument('--force', action='store_true', help='manifest: allow overwriting an existing targets.csv')
    p.add_argument('--task', type=int, default=None, help='array task index (default: $SLURM_ARRAY_TASK_ID or 0)')
    p.add_argument('--ntasks', type=int, default=None, help='array size (default: $SLURM_ARRAY_TASK_COUNT or 1)')
    p.add_argument('--only', default=None, help='comma-separated target_ids to restrict a stage to')
    p.add_argument('--linelist', choices=sorted(LINELISTS), default=None,
                   help='fit: force one line list for every star (default: targets.csv column `linelist`, else sousa)')
    p.add_argument('--compare', default=None, help='collect: CSV with pub_feh (etc.) per target to compare against')
    a = p.parse_args()
    only = a.only.split(',') if a.only else None
    inst = a.instruments.split(',') if a.instruments else None
    if a.stage == 'manifest':
        if not a.targets:
            p.error('manifest needs --targets')
        run_manifest(a.targets, inst or ('HARPS',), a.max_per_group, a.force)
    elif a.stage == 'coadd':
        run_coadd(*_task_args(a), only=only, instruments=inst)
    elif a.stage == 'fit':
        run_fit(*_task_args(a), only=only, linelist=a.linelist)
    elif a.stage == 'collect':
        run_collect(a.compare)
    elif a.stage == 'rescue':
        run_rescue(*_task_args(a), only=only, calib=not a.no_calib)
    elif a.stage == 'adopt':
        run_adopt()
    else:
        run_sun_test()


if __name__ == '__main__':
    main()
