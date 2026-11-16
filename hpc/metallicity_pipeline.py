"""Homogeneous spectroscopic [Fe/H] for every hot-Jupiter host and EBLM primary with public archival spectra.

Runs the eblm_feh PHOENIX fit on ESO (HARPS, FEROS, UVES, ESPRESSO) and SOPHIE spectra, on a Slurm cluster or locally.
Stages (each one reads the previous stage's output under $METAL_ROOT):

  manifest  targets + TIC priors + ephemerides, and the list of spectra to fetch       (one job, needs internet)
  download  fetch the spectra                                                         (array job, needs internet)
  coadd     per (star, instrument, resolution): RV-align, drop in-transit exposures, co-add   (array job)
  fit       eblm_feh.fit_spectrum on every co-add                                     (array job)
  combine   quality cuts, per-star mean, zero point, instrument offsets -> results/   (one job)

    python hpc/metallicity_pipeline.py manifest
    python hpc/metallicity_pipeline.py download --task 3 --ntasks 8     # or let Slurm array variables set these
    python hpc/metallicity_pipeline.py coadd --only HD_209458,WASP-43    # restrict any stage to some targets

See hpc/README.md for the cluster set-up and hpc/submit_all.sh for the chained Slurm submission.
"""
import argparse
import glob
import json
import os
import sys
import time
import urllib.parse
import warnings
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import numpy as np
import pandas as pd
from astropy.io import fits

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
os.chdir(HERE)   # eblm_feh and the coverage scripts use paths relative to the project root
warnings.filterwarnings('ignore', category=RuntimeWarning)   # all-NaN pixel columns at the co-add edges

import eblm_feh as ef          # noqa: E402
import eso_coverage as ec      # noqa: E402
import sophie_coverage as sc   # noqa: E402

ROOT = os.path.abspath(os.environ.get('METAL_ROOT', os.path.join(HERE, 'metallicity_hpc_data')))
CAT, RAW, COADD, FITS, RES = (os.path.join(ROOT, d) for d in ('catalogs', 'raw', 'coadd', 'fits', 'results'))

ESO_FILE = 'https://dataportal.eso.org/dataportal_new/file/'
SOPHIE_FILE = 'http://atlas.obs-hp.fr/sophie/sophie.cgi?c=i&a=mime:application/fits&o=sophie:[s1d,{seq}]'
SOPHIE_R = {'HR': 75000, 'HE': 39000}
VACUUM = {'ESPRESSO'}           # ESPRESSO Phase 3 WAVE is in vacuum (checked: +83 km/s offset against an air model)

MATCH_ARCSEC = {'HJ': 10.0, 'EBLM': 20.0}   # ESO pointings of EBLMs are often 10-20" off (see eblm_rm.eso_spectra)
SNR_SPEC_MIN = 10        # drop single exposures below this S/N before download
SNR_GROUP_MIN = 50       # co-added S/N a (star, instrument) group needs to be worth fetching
MAX_PER_GROUP = 50       # highest-S/N exposures kept per group; co-added S/N saturates long before 400 RM exposures
MAX_PER_GROUP_INST = {'ESPRESSO': 15}   # ESPRESSO 1-D products are ~40 MB each
TRANSIT_PAD_H = 0.5      # hours added on each side of the transit when rejecting in-transit (RM-distorted) exposures

WAVE_LO, WAVE_HI, DWAVE = 5000.0, 6300.0, 0.01   # co-add grid (air, Angstrom); covers every eblm_feh window
CCF_WINDOWS = (1, 2, 4)  # eblm_feh windows used to measure RVs
CCF_MIN = 0.3            # exposures whose CCF peak is below this are rejected (wrong star, no signal)


def _dirs():
    for d in (CAT, RAW, COADD, FITS, RES, os.path.join(COADD, 'index')):
        os.makedirs(d, exist_ok=True)


def _task_args(a):
    task = a.task if a.task is not None else int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
    ntasks = a.ntasks if a.ntasks is not None else int(os.environ.get('SLURM_ARRAY_TASK_COUNT', 1))
    return task, ntasks


def _workers():
    return int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count() or 1))


def slug(name):
    return ''.join(c if c.isalnum() or c in '+-.' else '_' for c in str(name).replace('EBLM ', 'EBLM_')).strip('_')


# ----------------------------------------------------------------------------- manifest

def _tic_by_id(ids, chunk=200):
    from astroquery.mast import Catalogs
    rows = []
    ids = sorted({int(i) for i in ids if np.isfinite(i)})
    for k in range(0, len(ids), chunk):
        t = Catalogs.query_criteria(catalog='Tic', ID=ids[k:k + chunk]).to_pandas()
        rows.append(t[['ID', 'Teff', 'e_Teff', 'logg', 'e_logg', 'Vmag', 'Tmag', 'GAIA']])
    t = pd.concat(rows, ignore_index=True)
    t['ID'] = t.ID.astype(int)
    return t.rename(columns={'ID': 'tic', 'Teff': 'tic_teff', 'e_Teff': 'tic_teff_err', 'logg': 'tic_logg',
                             'e_logg': 'tic_logg_err', 'Vmag': 'tic_vmag', 'Tmag': 'tmag', 'GAIA': 'gaia_dr2'})


def _hj_ephemerides():
    q = ('select pl_name, hostname, tic_id, pl_tranmid, pl_orbper, pl_trandur from pscomppars where pl_orbper < 10')
    return pd.read_csv('https://exoplanetarchive.ipac.caltech.edu/TAP/sync?format=csv&query=' + urllib.parse.quote(q))


def build_targets(path=None):
    """Hot-Jupiter hosts and EBLMs, with TIC Teff/logg priors and (for hot Jupiters) the transit ephemeris."""
    path = path or os.path.join(CAT, 'targets.csv')
    if os.path.exists(path):
        return pd.read_csv(path)
    hj = ec.hot_jupiter_targets()
    eph = _hj_ephemerides()
    hj = hj.merge(eph[['pl_name', 'tic_id', 'pl_tranmid', 'pl_orbper', 'pl_trandur']], on='pl_name', how='left')
    hj['tic'] = pd.to_numeric(hj.tic_id.astype(str).str.replace('TIC', '').str.strip(), errors='coerce')
    hj = hj.assign(sample='HJ', vmag=hj.sy_vmag)
    eb = ec.eblm_targets().rename(columns={'TIC': 'tic', 'Vmag': 'vmag'})
    eb = eb.assign(sample='EBLM', pl_name=np.nan, st_teff=np.nan, st_met=np.nan, disc_facility='WASP (EBLM)',
                   pl_tranmid=np.nan, pl_orbper=np.nan, pl_trandur=np.nan)
    cols = ['name', 'sample', 'pl_name', 'ra', 'dec', 'vmag', 'tic', 'st_teff', 'st_met', 'disc_facility',
            'pl_tranmid', 'pl_orbper', 'pl_trandur']
    t = pd.concat([hj[cols], eb[cols + ['have_spectra']]], ignore_index=True)
    t['target_id'] = t.name.map(slug)
    dup = t.target_id.duplicated(keep=False)   # a few EBLM names appear twice in the master table with different TICs
    t.loc[dup, 'target_id'] = t.target_id[dup] + '_TIC' + t.tic[dup].astype('Int64').astype(str)
    t = t.merge(_tic_by_id(t.tic.dropna()), on='tic', how='left')
    t.to_csv(path, index=False)
    return t


def eso_footprint_full(path=None):
    """Like eso_coverage.eso_footprint, plus what is needed to download and group each spectrum."""
    path = path or os.path.join(CAT, 'eso_footprint_full.csv')
    if os.path.exists(path):
        return pd.read_csv(path)
    import pyvo
    svc = pyvo.dal.TAPService(ec.ESO_TAP)
    parts = []
    for inst in ec.INSTRUMENTS:
        q = (f"select dp_id, s_ra, s_dec, obs_collection, snr, em_min, em_max, em_res_power, t_min, access_estsize, "
             f"obs_release_date from ivoa.ObsCore where obs_collection = '{inst}' and dataproduct_type = 'spectrum'")
        parts.append(svc.run_async(q, maxrec=1_000_000).to_table().to_pandas())
        print(f'ESO {inst}: {len(parts[-1])} spectra', flush=True)
    fp = pd.concat(parts, ignore_index=True)
    fp.to_csv(path, index=False)
    return fp


def _match(targets, ra, dec, radius_by_sample):
    """(target index, footprint index) pairs; each spectrum goes to its nearest target only."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c_fp = SkyCoord(ra, dec, unit='deg')
    c_t = SkyCoord(targets.ra.values, targets.dec.values, unit='deg')
    it, ifp, sep, _ = c_fp.search_around_sky(c_t, max(radius_by_sample.values()) * u.arcsec)
    m = pd.DataFrame({'it': it, 'ifp': ifp, 'sep': sep.arcsec})
    radius = targets['sample'].map(radius_by_sample).values[m.it.values]
    m = m[m.sep.values <= radius]
    return m.sort_values('sep').drop_duplicates('ifp')


def in_transit(jd, t):
    """True where `jd` falls within the transit (plus TRANSIT_PAD_H) of target row `t`; False without an ephemeris."""
    jd = np.asarray(jd, float)
    if not (np.isfinite(t.pl_tranmid) and np.isfinite(t.pl_orbper) and np.isfinite(t.pl_trandur)):
        return np.zeros(jd.shape, bool)
    ph = (jd - t.pl_tranmid) / t.pl_orbper
    dt_h = np.abs(ph - np.round(ph)) * t.pl_orbper * 24
    return dt_h < 0.5 * t.pl_trandur + TRANSIT_PAD_H


def build_manifest():
    """One row per spectrum to download, grouped into (target, instrument, resolution) co-add groups."""
    t = build_targets()
    rows = []

    fp = eso_footprint_full()
    fp = fp[(fp.em_min < ec.WAVE_MIN) & (fp.em_max > ec.WAVE_MAX) & (fp.snr >= SNR_SPEC_MIN)
            & (pd.to_datetime(fp.obs_release_date, utc=True, format='ISO8601', errors='coerce')
               < pd.Timestamp.now(tz='UTC'))]
    fp = fp.dropna(subset=['s_ra', 's_dec']).reset_index(drop=True)
    m = _match(t, fp.s_ra.values, fp.s_dec.values, MATCH_ARCSEC)
    e = fp.iloc[m.ifp.values].reset_index(drop=True)
    e['target_id'] = t.target_id.values[m.it.values]
    e['sep_arcsec'] = m.sep.values
    tt = t.set_index('target_id')
    e['in_transit'] = [bool(in_transit(r.t_min + 2400000.5, tt.loc[r.target_id])) for r in e.itertuples()]
    e = e[~e.in_transit]
    rows.append(pd.DataFrame({
        'target_id': e.target_id, 'archive': 'ESO', 'instrument': e.obs_collection,
        'resolution': (e.em_res_power.fillna(0) / 1000).round().astype(int) * 1000,
        'file_id': e.dp_id, 'url': ESO_FILE + e.dp_id, 'snr': e.snr, 'size_mb': e.access_estsize / 1024,
        'sep_arcsec': e.sep_arcsec}))

    sfp = sc.sophie_footprint()
    sfp = sfp[sfp.public & (sfp.snr >= SNR_SPEC_MIN) & sfp['mode'].isin(SOPHIE_R)].reset_index(drop=True)
    m = _match(t, sfp.s_ra.values, sfp.s_dec.values, {'HJ': ec.MATCH_ARCSEC, 'EBLM': ec.MATCH_ARCSEC})
    s = sfp.iloc[m.ifp.values].reset_index(drop=True)
    rows.append(pd.DataFrame({
        'target_id': t.target_id.values[m.it.values], 'archive': 'SOPHIE', 'instrument': 'SOPHIE_' + s['mode'],
        'resolution': s['mode'].map(SOPHIE_R), 'file_id': 'SOPHIE_' + s.seq.astype(str),
        'url': [SOPHIE_FILE.format(seq=q) for q in s.seq], 'snr': s.snr, 'size_mb': 1.2,
        'sep_arcsec': m.sep.values}))

    man = pd.concat(rows, ignore_index=True)
    man = man[man.resolution > 0]
    man['group'] = man.target_id + '__' + man.instrument + '_R' + man.resolution.astype(str)
    man = man.sort_values('snr', ascending=False)
    rank = man.groupby('group').cumcount()
    man = man[rank < man.instrument.map(MAX_PER_GROUP_INST).fillna(MAX_PER_GROUP)]
    gsnr = man.groupby('group').snr.apply(lambda x: np.sqrt(np.sum(x ** 2)))
    man = man[man.group.map(gsnr) >= SNR_GROUP_MIN]
    man['path'] = [os.path.join(RAW, r.archive, r.target_id, r.file_id.replace(':', '_') + '.fits')
                   for r in man.itertuples()]
    man = man.sort_values(['target_id', 'group', 'snr'], ascending=[True, True, False]).reset_index(drop=True)
    man.to_csv(os.path.join(CAT, 'manifest.csv'), index=False)

    g = man.groupby('target_id')
    print(f'{len(man)} spectra in {man.group.nunique()} co-add groups for {g.ngroups} targets '
          f'({t[t.target_id.isin(g.groups)]["sample"].value_counts().to_dict()}); '
          f'~{man.size_mb.sum() / 1024:.1f} GB to download')
    print(man.groupby('instrument').agg(n_spectra=('file_id', 'size'), n_groups=('group', 'nunique'),
                                        gb=('size_mb', lambda x: x.sum() / 1024)).round(1).to_string())
    return man


# ----------------------------------------------------------------------------- download

def _download(url, path, retries=5):
    import requests
    if os.path.exists(path) and os.path.getsize(path) > 50_000:
        return 'cached'
    os.makedirs(os.path.dirname(path), exist_ok=True)
    for k in range(retries):
        try:
            r = requests.get(url, timeout=600)
            if r.status_code == 200 and r.content[:6] == b'SIMPLE':
                with open(path + '.part', 'wb') as f:
                    f.write(r.content)
                os.replace(path + '.part', path)
                return 'ok'
        except requests.RequestException:
            pass
        time.sleep(15 * (k + 1))
    return 'failed'


def run_download(task, ntasks, only=None, threads=4):
    man = pd.read_csv(os.path.join(CAT, 'manifest.csv'))
    if only:
        man = man[man.target_id.isin(only)]
    mine = man.iloc[task::ntasks]
    with ThreadPoolExecutor(threads) as ex:
        status = list(ex.map(lambda r: _download(r[0], r[1]), zip(mine.url, mine.path)))
    out = mine.assign(status=status)[['file_id', 'path', 'status']]
    out.to_csv(os.path.join(CAT, f'download_status_{task:03d}.csv'), index=False)
    print(f'task {task}/{ntasks}: {pd.Series(status).value_counts().to_dict()}')


# ----------------------------------------------------------------------------- reading and co-adding

def read_raw(path, archive, instrument):
    """Air wavelength [A], flux, flux error, JD of mid-exposure and header S/N of one archival 1-D spectrum."""
    with fits.open(path) as h:
        hd = h[0].header
        if archive == 'SOPHIE':
            flux = h[0].data.astype(float)
            wave = hd['CRVAL1'] + hd['CDELT1'] * np.arange(flux.size)
            jd = hd.get('HIERARCH OHP DRS BJD', np.nan)
            snr = hd.get('HIERARCH OHP DRS CAL EXT SN26', np.nan)
            err = None
        else:
            t = h[1].data
            names = [c.upper() for c in t.columns.names]
            col = lambda k: np.asarray(t[t.columns.names[names.index(k)]][0], float)
            fk, ek = ('FLUX_REDUCED', 'ERR_REDUCED') if 'FLUX_REDUCED' in names else ('FLUX', 'ERR')
            wave, flux = col('WAVE'), col(fk)
            err = col(ek) if ek in names else None
            if wave.mean() < 1000:
                wave = wave * 10
            if instrument in VACUUM:
                wave = ef.vac_to_air(wave)
            jd = hd.get('MJD-OBS', np.nan) + 2400000.5 + hd.get('EXPTIME', 0) / 2 / 86400
            snr = hd.get('SNR', np.nan)
    if err is None or not np.isfinite(err).any() or np.nanmax(err) <= 0:
        # photon-noise shape, scaled to the header S/N near 550 nm
        err = np.sqrt(np.clip(flux, 0, None) + 1e-3 * np.nanmedian(np.abs(flux)))
        near = (wave > 5450) & (wave < 5550)
        if np.isfinite(snr) and snr > 0 and near.any():
            err *= np.nanmedian(flux[near] / err[near]) / snr
    return wave, flux, err, jd, snr


def _normalise(wave, flux, err, width_a=20.0):
    """Divide by a running 95th-percentile pseudo-continuum (removes blaze and flux-calibration differences)."""
    n = max(5, int(width_a / np.median(np.diff(wave))))
    cont = pd.Series(flux).rolling(n, center=True, min_periods=n // 3).quantile(0.95).to_numpy()
    cont = pd.Series(cont).rolling(n, center=True, min_periods=1).mean().to_numpy()
    return flux / cont, err / cont


class RVTemplate:
    """PHOENIX template (from the eblm_feh grid) on a log-lambda grid, for cross-correlation RVs."""

    def __init__(self, grid, teff=5800.0, logg=4.4, feh=0.0):
        teff = float(np.clip(teff if np.isfinite(teff) else 5800.0, grid.teff[0], grid.teff[-1]))
        logg = float(np.clip(logg if np.isfinite(logg) else 4.4, grid.logg[0], grid.logg[-1]))
        full = grid.interpolate(teff, logg, feh)
        self.w = [grid.waves[i] for i in CCF_WINDOWS]
        self.f = [full[grid.slices[i]] - 1 for i in CCF_WINDOWS]

    def ccf(self, wave, flux, vgrid):
        out = np.zeros(vgrid.size)
        for mw, mf in zip(self.w, self.f):
            s = (wave > mw[0] - 8) & (wave < mw[-1] + 8) & np.isfinite(flux)
            if s.sum() < 100:
                continue
            for k, v in enumerate(vgrid):
                o = np.interp(mw * (1 + v / ef.C_KMS), wave[s], flux[s]) - 1
                out[k] += np.dot(o - o.mean(), mf - mf.mean()) / (np.std(o) * np.std(mf) * o.size + 1e-30)
        return out / len(self.w)

    def rv(self, wave, flux):
        coarse = np.arange(-300, 300.1, 2.0)
        c = self.ccf(wave, flux, coarse)
        v0 = coarse[np.argmax(c)]
        fine = np.arange(v0 - 4, v0 + 4.01, 0.1)
        cf = self.ccf(wave, flux, fine)
        return fine[np.argmax(cf)], cf.max()


def coadd_group(rows, target, grid, lo=WAVE_LO, hi=WAVE_HI, norm_width=20.0, out_root=None,
                note='Rest-frame co-add (air wavelengths), normalised; metallicity_pipeline.py'):
    """Co-add the exposures of one group in the stellar rest frame; returns the index row (and writes the FITS).

    The defaults are this pipeline's co-adds; ew_pipeline.py calls it with a wider range (lo, hi), a gentler
    normalisation (norm_width, in A) and its own output folder (out_root)."""
    r0 = rows.iloc[0]
    out_path = os.path.join(out_root or COADD, r0.target_id, f'{r0.instrument}_R{r0.resolution}.fits')
    info = {'group': r0.group, 'target_id': r0.target_id, 'instrument': f'{r0.instrument}_R{r0.resolution}',
            'archive': r0.archive, 'resolution': r0.resolution, 'file': out_path, 'n_files': len(rows)}
    tmpl = RVTemplate(grid, target.tic_teff, target.tic_logg)
    x = np.arange(lo, hi, DWAVE)
    stack, ivar, used, rvs, n_transit, n_badccf, n_missing = [], [], [], [], 0, 0, 0
    for r in rows.itertuples():
        if not os.path.exists(r.path):
            n_missing += 1
            continue
        try:
            w, f, e, jd, snr = read_raw(r.path, r.archive, r.instrument)
        except Exception:
            n_missing += 1
            continue
        if in_transit(jd, target):
            n_transit += 1
            continue
        s = (w > lo - 30) & (w < hi + 30) & np.isfinite(f) & np.isfinite(e) & (e > 0)
        if s.sum() < 1000:
            n_missing += 1
            continue
        w, f, e = w[s], f[s], e[s]
        fn, en = _normalise(w, f, e, width_a=norm_width)
        rv, peak = tmpl.rv(w, fn)
        if peak < CCF_MIN:
            n_badccf += 1
            continue
        wr = w / (1 + rv / ef.C_KMS)
        stack.append(np.interp(x, wr, fn, left=np.nan, right=np.nan))
        ivar.append(1 / np.interp(x, wr, en, left=np.nan, right=np.nan) ** 2)
        used.append(os.path.basename(r.path)); rvs.append(rv)
    info.update(n_used=len(used), n_in_transit=n_transit, n_bad_ccf=n_badccf, n_missing=n_missing,
                rv_median=np.median(rvs) if rvs else np.nan, rv_std=np.std(rvs) if len(rvs) > 1 else np.nan)
    if not used:
        info['file'] = None
        return info
    F, W = np.array(stack), np.nan_to_num(np.array(ivar))
    F = np.where(np.isfinite(F), F, 0.0); W = np.where(F > 0, W, 0.0)
    if len(used) >= 3:   # clip cosmics / bad pixels against the per-pixel median
        med = np.nanmedian(np.where(W > 0, F, np.nan), axis=0)
        mad = 1.4826 * np.nanmedian(np.abs(np.where(W > 0, F, np.nan) - med), axis=0)
        W = np.where(np.abs(F - med) > 5 * np.maximum(mad, 1 / np.sqrt(W.max(axis=0) + 1e-30)), 0.0, W)
    wsum = W.sum(axis=0)
    flux = np.where(wsum > 0, (W * F).sum(axis=0) / np.where(wsum > 0, wsum, 1), np.nan)
    info['snr_coadd'] = float(np.nanmedian(np.sqrt(wsum[wsum > 0]) * flux[wsum > 0]))
    hdr = fits.Header()
    hdr['CRVAL1'], hdr['CDELT1'], hdr['CTYPE1'], hdr['CUNIT1'] = lo, DWAVE, 'AWAV', 'Angstrom'
    # FITS headers must be ASCII; some catalogue names carry a Unicode minus (e.g. 'CD−31 3271')
    obj = str(target['name']).replace('−', '-').encode('ascii', 'replace').decode()
    hdr['OBJECT'], hdr['INSTRUME'], hdr['SPEC_RES'] = obj, r0.instrument, int(r0.resolution)
    hdr['NCOMBINE'], hdr['SNR'], hdr['RVMED'], hdr['RVSTD'] = len(used), info['snr_coadd'], info['rv_median'], \
        np.nan_to_num(info['rv_std'])
    hdr['COMMENT'] = note
    for k, u in enumerate(used[:999]):
        hdr[f'FILE{k:03d}'] = u
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fits.PrimaryHDU(np.nan_to_num(flux, nan=0.0).astype(np.float32), hdr).writeto(out_path, overwrite=True)
    return info


def _coadd_worker(args):
    rows, target = args
    grid = _grid()
    try:
        return coadd_group(rows, target, grid)
    except Exception as e:
        return {'group': rows.iloc[0].group, 'target_id': rows.iloc[0].target_id, 'file': None, 'error': repr(e)}


_G = {}


def _grid():
    if 'grid' not in _G:
        _G['grid'] = ef.ModelGrid()
    return _G['grid']


def run_coadd(task, ntasks, only=None):
    man = pd.read_csv(os.path.join(CAT, 'manifest.csv'))
    t = pd.read_csv(os.path.join(CAT, 'targets.csv')).set_index('target_id', drop=False)
    if only:
        man = man[man.target_id.isin(only)]
    groups = sorted(man.group.unique())[task::ntasks]
    jobs = [(man[man.group == g], t.loc[g.split('__')[0]]) for g in groups]
    with ProcessPoolExecutor(_workers()) as ex:
        out = list(ex.map(_coadd_worker, jobs))
    pd.DataFrame(out).to_csv(os.path.join(COADD, 'index', f'coadd_{task:03d}.csv'), index=False)
    print(f'task {task}/{ntasks}: {sum(o.get("file") is not None for o in out)}/{len(out)} groups co-added')


# ----------------------------------------------------------------------------- fitting

def _fit_worker(args):
    path, label, resolution, teff0, logg = args
    ef.RESOLUTION[label] = resolution
    if 'mask' not in _G:
        z = np.load(os.path.join(ef.D, 'solar_line_mask.npz'))
        _G['mask'] = [z[f'w{i}'] for i in range(len(ef.WINDOWS))]
    try:
        res, _ = ef.fit_spectrum(_grid(), path, label, _G['mask'], teff0, logg)
    except Exception as e:
        res = {'error': repr(e)}
    return {'file': path, **res}


def _coadd_index():
    files = sorted(glob.glob(os.path.join(COADD, 'index', 'coadd_*.csv')))
    idx = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    return idx.drop_duplicates('group', keep='last')


def _local_coadd_path(p):
    """The co-add index stores absolute paths from the machine that ran the co-add stage. If the co-adds were copied
    elsewhere (e.g. from the cluster to a laptop), find the same file under this METAL_ROOT's coadd/ folder."""
    if os.path.exists(p):
        return p
    i = str(p).replace('\\', '/').rfind('/coadd/')
    return os.path.join(COADD, str(p)[i + len('/coadd/'):]) if i >= 0 else p


def run_fit(task, ntasks, only=None):
    idx = _coadd_index().dropna(subset=['file'])
    idx['file'] = idx.file.map(_local_coadd_path)
    t = pd.read_csv(os.path.join(CAT, 'targets.csv')).set_index('target_id')
    if only:
        idx = idx[idx.target_id.isin(only)]
    mine = idx.sort_values('group').iloc[task::ntasks]
    jobs = [(r.file, r.instrument, int(r.resolution), t.loc[r.target_id, 'tic_teff'], t.loc[r.target_id, 'tic_logg'])
            for r in mine.itertuples()]
    with ProcessPoolExecutor(_workers()) as ex:
        out = list(ex.map(_fit_worker, jobs))
    res = mine.merge(pd.DataFrame(out), on='file').rename(columns={'target_id': 'name'})
    res.to_csv(os.path.join(FITS, f'fits_{task:03d}.csv'), index=False)
    print(f'task {task}/{ntasks}: {res.feh.notna().sum() if "feh" in res else 0}/{len(res)} co-adds fitted')


# ----------------------------------------------------------------------------- combining

def calibrators(t):
    """[Fe/H] from ARES+MOOG/SME-type analyses, keyed by target_id: SWEET-Cat homogeneous hot-Jupiter hosts and the
    EBLM literature compilation."""
    hj = ef.hot_jupiter_hosts(ef.load_hot_jupiters(), ef.load_sweetcat())
    hj = hj[hj.sc_homogeneous & hj.sc_feh.notna()]
    ids = t[t['sample'] == 'HJ'].set_index('name').target_id
    sw = pd.DataFrame({'name': hj.hostname.map(ids), 'feh': hj.sc_feh, 'method': 'ARES+MOOG', 'source': 'SWEET-Cat'})
    lit = pd.read_csv(os.path.join(ef.D, 'eblm_literature_feh.csv'), dtype={'arxiv_id': str})
    ids = t[t['sample'] == 'EBLM'].assign(short=lambda d: d.name.str.replace('EBLM ', '')).drop_duplicates('short')
    lit = lit.assign(name=lit.name.map(ids.set_index('short').target_id), source='EBLM literature')
    return pd.concat([sw.dropna(subset=['name']), lit.dropna(subset=['name'])[['name', 'feh', 'method', 'source']]],
                     ignore_index=True)


def instrument_offsets(fits_ok, ref='HARPS_R115000'):
    """Median [Fe/H] difference of each instrument against `ref`, over stars measured by both."""
    p = fits_ok.pivot_table(index='name', columns='instrument', values='feh')
    rows = []
    for col in p.columns:
        if col == ref:
            continue
        d = (p[col] - p[ref]).dropna() if ref in p else pd.Series(dtype=float)
        rows.append({'instrument': col, 'vs': ref, 'n': len(d), 'median_diff': d.median() if len(d) else np.nan,
                     'robust_rms': 1.4826 * np.median(np.abs(d - d.median())) if len(d) > 2 else np.nan})
    return pd.DataFrame(rows)


ZP_PER_BIN = 20                            # calibrators per Teff bin (bins are for display/diagnostics only)
ZP_BREAKS = np.arange(5000.0, 6505.0, 10.0)  # break temperatures searched by the zero-point fit


def _huber_line(x, y, scale=None, n_iter=8, k=1.345):
    """Robust (Huber IRLS) fit of y = a + b x. Returns (a, b, Huber loss, robust scatter of the residuals). With
    `scale` given, weights and loss use that fixed scale, so losses of different models are comparable."""
    w = np.ones_like(y)
    A = np.column_stack([np.ones_like(x), x])
    for _ in range(n_iter):
        sw = np.sqrt(w)
        a, b = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)[0]
        r = y - a - b * x
        s = scale or max(1.4826 * np.median(np.abs(r - np.median(r))), 1e-3)
        u = np.abs(r) / (k * s)
        w = np.where(u <= 1, 1.0, 1.0 / u)
    u = np.abs(r) / s
    loss = np.where(u <= k, 0.5 * u ** 2, k * u - 0.5 * k ** 2).sum()
    return a, b, loss, 1.4826 * np.median(np.abs(r - np.median(r)))


def _fit_broken_line(teff, d, breaks=ZP_BREAKS):
    """zp(T) = zp_hot + slope * (min(T, T_break) - T_break): rises linearly up to T_break, constant above.
    T_break is chosen on a grid by the lowest Huber loss, all candidates scored with one fixed scale (the robust
    scatter about a first-pass fit), so the losses are comparable."""
    mid = breaks[len(breaks) // 2]
    scale = max(_huber_line(np.minimum(teff, mid) - mid, d)[3], 1e-3)
    best = None
    for tb in breaks:
        a, b, loss, _ = _huber_line(np.minimum(teff, tb) - tb, d, scale=scale)
        if best is None or loss < best[3]:
            best = (tb, a, b, loss)
    return {'teff_break': best[0], 'zp_hot': best[1], 'slope': best[2]}


def _broken_line(teff, p):
    return p['zp_hot'] + p['slope'] * (np.minimum(teff, p['teff_break']) - p['teff_break'])


def teff_zero_point(stars, cal, n_boot=300, seed=1, per_bin=ZP_PER_BIN):
    """Zero point as a smooth function of Teff, fitted to (literature - this work) of the calibrators: a straight line
    rising up to a break temperature, constant above it (robust Huber fit; the break is fitted too). Outside the
    calibrators' Teff range it is held at its value at the range edge. The error at each Teff is the bootstrap scatter
    of the whole fit. Returns a dict: params (dict), grid (teff, zp, zp_err on a 10 K grid over the calibrated range),
    bins (binned medians, for plots), resid (the calibrators' residuals and Teff), and zp(teff) / zp_err(teff)."""
    _, _, d = ef.zero_point(stars, cal)
    df = pd.DataFrame({'d': d, 'teff': stars.set_index('name').teff.reindex(d.index)}).dropna().sort_values('teff')
    T, y = df.teff.values, df.d.values
    lo, hi = T.min(), T.max()
    p = _fit_broken_line(T, y)
    grid = np.arange(np.floor(lo / 10) * 10, hi + 10, 10.0)
    rng = np.random.default_rng(seed)
    boot = []
    for _ in range(n_boot):
        i = rng.integers(0, len(T), len(T))
        boot.append(_broken_line(grid, _fit_broken_line(T[i], y[i])))
    g = pd.DataFrame({'teff': grid, 'zp': _broken_line(grid, p), 'zp_err': np.std(boot, axis=0)})

    def zp(teff):
        return np.interp(np.clip(teff, lo, hi), g.teff, g.zp)

    def zp_err(teff):
        return np.interp(np.clip(teff, lo, hi), g.teff, g.zp_err)

    rows = []
    for idx in np.array_split(np.arange(len(df)), max(1, len(df) // per_bin)):
        c = df.iloc[idx]
        bb = np.median(rng.choice(c.d.values, (1000, len(c))), axis=1)
        rows.append({'teff': c.teff.median(), 'teff_lo': c.teff.min(), 'teff_hi': c.teff.max(), 'n': len(c),
                     'zp': c.d.median(), 'zp_err': bb.std()})
    p = {k: float(v) for k, v in p.items()} | {'teff_min': float(lo), 'teff_max': float(hi)}
    return {'params': p, 'grid': g, 'bins': pd.DataFrame(rows), 'resid': df, 'zp': zp, 'zp_err': zp_err}


# CORALIE (Euler 1.2 m; not in the ESO or SOPHIE archives): fits from the earlier local analysis
# (eblm_feh.fit_all on the EBLM_DATA co-adds), added to the archive fits at the combine stage. Only CORALIE rows are
# used; that file's HARPS rows come from the same ESO data the pipeline co-adds itself.
CORALIE_FITS = os.path.join(HERE, 'eblm_metallicity_data', 'eblm_spectroscopic_fits.csv')
CORALIE_MATCH_ARCSEC = 30.0   # file-name coordinates are truncated to 1 s of RA (up to 15") and 1" of Dec


def coralie_fits(t, path=CORALIE_FITS):
    """CORALIE fits mapped onto EBLM target_ids, in the pipeline's fit-table format.

    Matched by the coordinates encoded in the co-add file names to the nearest EBLM target within
    CORALIE_MATCH_ARCSEC, falling back to the name (J0008+02 -> EBLM J0008+02) when a row has no coordinates or no
    target is close enough. Returns (fits, match table); the match table lists every CORALIE row with its target_id
    (NaN if unmatched), separation and how it was matched."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    if not os.path.exists(path):
        return pd.DataFrame(), pd.DataFrame()
    lf = pd.read_csv(path)
    lf = lf[lf.instrument.isin(['CORALIE07', 'CORALIE14'])].reset_index(drop=True)
    eb = t[(t['sample'] == 'EBLM') & t.ra.notna() & t.dec.notna()].reset_index(drop=True)
    by_name = (eb.assign(short=eb.name.str.replace('EBLM ', '', regex=False))
               .drop_duplicates('short', keep=False).set_index('short').target_id)   # ambiguous names left out
    tid = pd.Series([None] * len(lf), dtype=object)
    sep = np.full(len(lf), np.nan)
    how = np.array(['unmatched'] * len(lf), dtype=object)
    has = (lf.ra_str.notna() & lf.dec_str.notna()).values
    if has.any():
        c = SkyCoord(lf.ra_str[has].values, lf.dec_str[has].values)
        idx, s, _ = c.match_to_catalog_sky(SkyCoord(eb.ra.values * u.deg, eb.dec.values * u.deg))
        s = s.arcsec
        close = s < CORALIE_MATCH_ARCSEC
        rows = np.flatnonzero(has)
        tid.iloc[rows[close]] = eb.target_id.values[idx[close]]
        sep[rows] = s
        how[rows[close]] = 'coordinates'
    byn = lf.name.map(by_name)
    fill = tid.isna() & byn.notna()
    tid[fill] = byn[fill]
    how[fill.values] = 'name'
    match = pd.DataFrame({'coralie_name': lf.name, 'instrument': lf.instrument, 'target_id': tid,
                          'sep_arcsec': sep, 'matched_by': how,
                          'name_agrees': (byn == tid) | byn.isna()})
    res = lf.instrument.map(ef.RESOLUTION).astype(int)
    out = lf.assign(name=tid, instrument=lf.instrument + '_R' + res.astype(str), resolution=res,
                    archive='CORALIE (local)', n_files=lf.n_exposures, n_used=lf.n_exposures)
    out['group'] = out.name + '__' + out.instrument
    out = out[out.name.notna()].drop(columns=['ra_str', 'dec_str', 'n_exposures'])
    return out.reset_index(drop=True), match


def run_combine(coralie=True):
    t = pd.read_csv(os.path.join(CAT, 'targets.csv'))
    fits_tab = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(FITS, 'fits_*.csv')))],
                         ignore_index=True)
    fits_tab = fits_tab.drop_duplicates('file', keep='last')
    if 'error' in fits_tab:
        fits_tab = fits_tab[fits_tab.error.isna()]
    if coralie:
        cf, match = coralie_fits(t)
        if not match.empty:
            match.to_csv(os.path.join(RES, 'coralie_matches.csv'), index=False)
            new = set(cf.name) - set(fits_tab.name)
            print(f'CORALIE: {len(cf)}/{len(match)} co-adds matched to {cf.name.nunique()} EBLM targets '
                  f'({len(new)} not in the archive sample); {(match.matched_by == "unmatched").sum()} unmatched, '
                  f'{(~match.name_agrees).sum()} where name and coordinates disagree (see coralie_matches.csv)')
            fits_tab = pd.concat([fits_tab, cf], ignore_index=True)
    tq = t.rename(columns={'name': 'display_name', 'target_id': 'name'})
    fits_tab = ef.quality_cuts(fits_tab, tq)
    # a fit whose covariance failed (NaN [Fe/H] error) would give its star NaN [Fe/H] and infinite error in the
    # inverse-variance mean (WASP-98 and WASP-183, UVES R=57000, in the Nov 2026 run): treat it as not ok
    fits_tab['ok'] &= np.isfinite(fits_tab.feh_err)
    fits_tab.to_csv(os.path.join(RES, 'all_fits.csv'), index=False)
    print(f'{fits_tab.ok.sum()}/{len(fits_tab)} co-adds pass the cuts')
    if not fits_tab.ok.any():
        print('no co-add passed the quality cuts (see all_fits.csv): nothing to combine')
        return
    stars = ef.combine_per_star(fits_tab)

    cal = calibrators(t)
    zps = {}
    for label, c in [('pooled', cal), ('SWEET-Cat', cal[cal.source == 'SWEET-Cat']),
                     ('EBLM literature', cal[cal.source == 'EBLM literature'])]:
        if c.empty:
            continue
        zp, zp_err, d = ef.zero_point(stars, c)
        zps[label] = {'zp': float(zp), 'zp_err': float(zp_err), 'n': int(d.size), 'scatter': float(d.std(ddof=1))}
    z = zps['pooled']
    if z['n'] < 3:
        print(f'WARNING: only {z["n"]} calibration stars; the zero point is not meaningful yet')
    # The offset depends on Teff (~0 dex at 5000 K, ~+0.23 dex above ~5600 K), so the adopted zero point is a smooth
    # fit in Teff (rising line + plateau); the constant one is kept as feh_constzp for comparison.
    stars['feh_constzp'] = stars.feh_raw + z['zp']
    zt = teff_zero_point(stars, cal)
    zt['grid'].to_csv(os.path.join(RES, 'zero_point_teff.csv'), index=False)
    zt['bins'].to_csv(os.path.join(RES, 'zero_point_teff_bins.csv'), index=False)
    resid = zt['resid']
    post = resid.d - zt['zp'](resid.teff.values)
    zp_p = zt['params']
    zps['pooled, Teff-dependent'] = {'zp': zp_p['zp_hot'], 'zp_err': float(zt['grid'].zp_err.mean()),
                                     'n': int(resid.d.size), 'scatter': float(post.std(ddof=1)),
                                     'teff_break': zp_p['teff_break'], 'slope_per_100K': 100 * zp_p['slope'],
                                     'teff_min': zp_p['teff_min'], 'teff_max': zp_p['teff_max']}
    has_t = stars.teff.notna()
    tf = stars.teff.fillna(0).values
    stars['zp'] = np.where(has_t, zt['zp'](tf), z['zp'])
    stars['zp_err'] = np.where(has_t, zt['zp_err'](tf), z['zp_err'])
    stars['teff_outside_calib'] = has_t & ((stars.teff < zp_p['teff_min']) | (stars.teff > zp_p['teff_max']))
    stars['feh'] = stars.feh_raw + stars.zp
    stars['feh_err'] = np.sqrt(stars.feh_stat ** 2 + stars.zp_err ** 2 + post.std(ddof=1) ** 2)
    print(f"Teff-dependent zero point: {zp_p['zp_hot']:+.3f} dex above {zp_p['teff_break']:.0f} K, "
          f"{100 * zp_p['slope']:+.3f} dex per 100 K below; calibrated {zp_p['teff_min']:.0f}-{zp_p['teff_max']:.0f} K; "
          f"scatter {z['scatter']:.3f} -> {post.std(ddof=1):.3f} dex")
    hq = cal[cal.method.str.startswith(ef.CALIB_METHODS)]
    cal_mean = hq.groupby('name').agg(cal_feh=('feh', 'mean'), cal_source=('source', 'first'))
    paws = cal[cal.method.str.startswith('PAWS')].groupby('name').feh.mean().rename('paws_lit_feh')
    stars = (stars.merge(tq, on='name', how='left').merge(cal_mean, left_on='name', right_index=True, how='left')
             .merge(paws, left_on='name', right_index=True, how='left'))
    stars = stars.rename(columns={'name': 'target_id', 'display_name': 'name'})
    stars.to_csv(os.path.join(RES, 'metallicity_final.csv'), index=False)
    instrument_offsets(fits_tab[fits_tab.ok]).to_csv(os.path.join(RES, 'instrument_offsets.csv'), index=False)
    with open(os.path.join(RES, 'zero_points.json'), 'w') as f:
        json.dump(zps, f, indent=2)
    print(f'{len(stars)} stars ({stars["sample"].value_counts().to_dict()})')
    print('zero points:', json.dumps(zps, indent=1))


# ----------------------------------------------------------------------------- CLI

def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('stage', choices=['manifest', 'download', 'coadd', 'fit', 'combine'])
    p.add_argument('--task', type=int, default=None, help='array task index (default: $SLURM_ARRAY_TASK_ID or 0)')
    p.add_argument('--ntasks', type=int, default=None, help='array size (default: $SLURM_ARRAY_TASK_COUNT or 1)')
    p.add_argument('--only', default=None, help='comma-separated target_ids to restrict a stage to')
    p.add_argument('--no-coralie', action='store_true',
                   help='combine: leave out the local CORALIE fits (eblm_metallicity_data/eblm_spectroscopic_fits.csv)')
    a = p.parse_args()
    _dirs()
    task, ntasks = _task_args(a)
    only = a.only.split(',') if a.only else None
    {'manifest': lambda: build_manifest(),
     'download': lambda: run_download(task, ntasks, only),
     'coadd': lambda: run_coadd(task, ntasks, only),
     'fit': lambda: run_fit(task, ntasks, only),
     'combine': lambda: run_combine(coralie=not a.no_coralie)}[a.stage]()


if __name__ == '__main__':
    main()
