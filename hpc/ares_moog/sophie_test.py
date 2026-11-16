"""Why do SOPHIE HE (R = 39,000) ARES+MOOG fits come out ~0.11 dex metal-poor relative to HARPS, and fail more often?

Uses the stars that have a converged HARPS (R = 115,000) fit AND a SOPHIE HE co-add from the main EW run, and re-fits
them with one change at a time:

  SOPHIE HE co-add
    he_base        as in the main run (ARES estimates the S/N itself; continuum window space = 2 A)
    he_snr_x1.5    ARES continuum level (rejt = 1 - 1/(S/N)) from 1.5 x the S/N ARES estimated in the main run
    he_snr_x2.5    ... 2.5 x                     -> is the continuum placed too low?
    he_space3      local-continuum window +-3 A (ARES's own default)
    he_space4      +-4 A                         -> do neighbouring line wings at low R drag the continuum down?
    he_fixlogg     logg fixed at the TIC (Gaia-based) value -> does the logg runaway cause the offset?
  HARPS co-add (controls)
    harps_deg39k        convolved to R = 39,000, noiseless      -> is it resolution itself?
    harps_deg39k_noise  ... plus white noise matching the SOPHIE co-add's S/N (as measured by ARES)
    harps_snr_x0.67     continuum level from 0.67 x ARES's S/N  -> can a low continuum alone make -0.1 dex?
    harps_fixlogg       logg fixed at the TIC value             -> what fixing logg does on good spectra
  The reference is each star's HARPS fit from the main run (results/ew_fits.csv).

Stages (aresmoog env; METAL_ROOT = the main metallicity directory, ARESMOOG_HOME set):
    python hpc/ares_moog/sophie_test.py select     # pick the stars           -> $METAL_ROOT/sophie_test/stars.csv
    python hpc/ares_moog/sophie_test.py degrade    # make the degraded HARPS spectra
    python hpc/ares_moog/sophie_test.py run        # all variants (Slurm: hpc/ares_moog/sophie_test.sbatch)
    python hpc/ares_moog/sophie_test.py report     # summary table           -> $METAL_ROOT/sophie_test/report.csv
    python hpc/ares_moog/sophie_test.py pedestal   # additive light in SOPHIE vs HARPS_39k; writes he_pedcorr spectra
    python hpc/ares_moog/sophie_test.py run --variants he_pedcorr   # re-fit the corrected SOPHIE spectra
    python hpc/ares_moog/sophie_test.py lines      # line-by-line EWs vs HARPS (seconds; login node is fine)
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import ew_pipeline as ew   # noqa: E402

OUT = os.path.join(ew.ROOT, 'sophie_test')
HE, HARPS = 'SOPHIE_HE_R39000', 'HARPS_R115000'
R_HE, R_HARPS = 39000.0, 115000.0
SEED = 192


def _ares(**kw):
    """ARES options for FASMA: the main-run options plus overrides (snr -> rejt = 1 - 1/snr; space in A)."""
    return ew.ARES_OPTIONS + ''.join(f',{k}:{v}' for k, v in kw.items())


def _main_fits():
    f = pd.read_csv(os.path.join(ew.RES, 'ew_fits.csv'))
    f['inst'] = f.group.str.split('__').str[1]
    return f


# ----------------------------------------------------------------------------- select

def run_select():
    f = _main_fits()
    idx = ew._coadd_index()
    good_h = f[(f.inst == HARPS) & (f.status == 'ok') & (f.logg < 4.99)].set_index('target_id')
    he = idx[idx.group.str.endswith('__' + HE)].set_index('target_id')
    ids = sorted(set(good_h.index) & set(he.index))
    fh = f[f.inst == HE].set_index('target_id')
    t = ew._targets()
    rows = []
    for i in ids:
        rows.append(dict(
            target_id=i, sample=t.loc[i, 'sample'] if i in t.index else None,
            tic_teff=t.loc[i, 'tic_teff'], tic_logg=t.loc[i, 'tic_logg'],
            harps_file=good_h.loc[i, 'file'], harps_snr=good_h.loc[i, 'snr_coadd'],
            harps_snr_ares=good_h.loc[i, 'snr_ares'], harps_teff=good_h.loc[i, 'teff'],
            harps_logg=good_h.loc[i, 'logg'], harps_feh=good_h.loc[i, 'feh'], harps_feh_err=good_h.loc[i, 'feh_err'],
            he_file=he.loc[i, 'file'], he_snr=he.loc[i, 'snr_coadd'],
            he_snr_ares=fh.loc[i, 'snr_ares'] if i in fh.index else np.nan,
            he_status_main=fh.loc[i, 'status'] if i in fh.index else 'not fitted',
            he_feh_main=fh.loc[i, 'feh'] if i in fh.index else np.nan))
    s = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    s.to_csv(os.path.join(OUT, 'stars.csv'), index=False)
    print(f'{len(s)} stars with a converged HARPS fit and a SOPHIE HE co-add '
          f'(SOPHIE HE in the main run: {s.he_status_main.value_counts().to_dict()})')
    print(s[['target_id', 'sample', 'harps_snr', 'he_snr', 'harps_feh', 'he_feh_main', 'he_status_main']]
          .round(2).to_string(index=False))


# ----------------------------------------------------------------------------- degrade HARPS to SOPHIE HE

def degrade(wave, flux, r_from=R_HARPS, r_to=R_HE, chunk=200.0):
    """Gaussian convolution from resolution r_from to r_to, in ~200 A chunks (kernel width at each chunk's centre).

    Missing pixels (flux <= 0) are excluded by normalised convolution and stay missing, with the mask grown by the
    kernel so no degraded pixel mixes data and gaps."""
    from scipy.ndimage import gaussian_filter1d, binary_dilation
    good = np.isfinite(flux) & (flux > 0)
    f0 = np.where(good, flux, 0.0)
    m0 = good.astype(float)
    dw = np.median(np.diff(wave))
    out = np.zeros_like(flux)
    edges = np.arange(wave[0], wave[-1] + chunk, chunk)
    pad = 60.0   # A of overlap so chunk edges are convolved with real neighbours
    for lo, hi in zip(edges[:-1], edges[1:]):
        lam = 0.5 * (lo + hi)
        sig_pix = lam * np.sqrt(1 / r_to ** 2 - 1 / r_from ** 2) / 2.3548 / dw
        s = (wave >= lo - pad) & (wave < hi + pad)
        num = gaussian_filter1d(f0[s], sig_pix, mode='nearest')
        den = gaussian_filter1d(m0[s], sig_pix, mode='nearest')
        core = (wave[s] >= lo) & (wave[s] < hi)
        out[np.where(s)[0][core]] = (num / np.where(den > 0, den, 1))[core]
    sig_max = wave[-1] * np.sqrt(1 / r_to ** 2 - 1 / r_from ** 2) / 2.3548 / dw
    bad = binary_dilation(~good, iterations=int(np.ceil(3 * sig_max)))
    return np.where(bad, 0.0, out)


def run_degrade():
    from astropy.io import fits
    s = pd.read_csv(os.path.join(OUT, 'stars.csv'))
    rng = np.random.default_rng(SEED)
    for r in s.itertuples():
        with fits.open(r.harps_file) as h:
            flux = np.asarray(h[0].data, float)
            hdr = h[0].header.copy()
        wave = hdr['CRVAL1'] + hdr['CDELT1'] * np.arange(flux.size)
        deg = degrade(wave, flux)
        d = os.path.join(OUT, 'spectra', r.target_id)
        os.makedirs(d, exist_ok=True)
        hdr['SPEC_RES'] = int(R_HE)
        hdr['COMMENT'] = f'HARPS co-add convolved to R={R_HE:.0f} (sophie_test.py)'
        fits.PrimaryHDU(deg.astype(np.float32), hdr).writeto(os.path.join(d, 'harps_deg39k.fits'), overwrite=True)
        # white noise so ARES measures the same S/N as on this star's SOPHIE HE co-add (fallback: its co-add S/N)
        target = r.he_snr_ares if np.isfinite(r.he_snr_ares) and r.he_snr_ares > 0 else r.he_snr
        noisy = np.where(deg > 0, deg + rng.normal(0, 1, deg.size) * deg / target, 0.0)
        hdr['COMMENT'] = f'... plus white noise, S/N {target:.0f}'
        fits.PrimaryHDU(noisy.astype(np.float32), hdr).writeto(os.path.join(d, 'harps_deg39k_noise.fits'),
                                                                overwrite=True)
        print(f'{r.target_id}: degraded (noise S/N {target:.0f})')


# ----------------------------------------------------------------------------- run variants

def _jobs(s, variants=None):
    home, fasma, bindir = ew._aresmoog_paths()
    base = ew.FASMA_OPTIONS['sousa']
    fix = base + ',logg'
    jobs = []
    for r in s.itertuples():
        spec = os.path.join(OUT, 'spectra', r.target_id)
        he_sn = r.he_snr_ares if np.isfinite(r.he_snr_ares) and r.he_snr_ares > 0 else r.he_snr
        h_sn = r.harps_snr_ares if np.isfinite(r.harps_snr_ares) and r.harps_snr_ares > 0 else r.harps_snr
        logg_fix = r.tic_logg if np.isfinite(r.tic_logg) else 4.4
        v = {
            'he_base': (r.he_file, ew.ARES_OPTIONS, base, r.tic_logg),
            'he_snr_x1.5': (r.he_file, _ares(snr=f'{1.5 * he_sn:.1f}'), base, r.tic_logg),
            'he_snr_x2.5': (r.he_file, _ares(snr=f'{2.5 * he_sn:.1f}'), base, r.tic_logg),
            'he_space3': (r.he_file, _ares(space=3.0), base, r.tic_logg),
            'he_space4': (r.he_file, _ares(space=4.0), base, r.tic_logg),
            'he_fixlogg': (r.he_file, ew.ARES_OPTIONS, fix, logg_fix),
            'harps_deg39k': (os.path.join(spec, 'harps_deg39k.fits'), ew.ARES_OPTIONS, base, r.tic_logg),
            'harps_deg39k_noise': (os.path.join(spec, 'harps_deg39k_noise.fits'), ew.ARES_OPTIONS, base, r.tic_logg),
            'harps_snr_x0.67': (r.harps_file, _ares(snr=f'{0.67 * h_sn:.1f}'), base, r.tic_logg),
            'harps_fixlogg': (r.harps_file, ew.ARES_OPTIONS, fix, logg_fix),
            'he_pedcorr': (os.path.join(spec, 'he_pedcorr.fits'), ew.ARES_OPTIONS, base, r.tic_logg),
        }
        if variants is None:
            v.pop('he_pedcorr')   # only on request (made by the pedestal stage)
        else:
            v = {k: x for k, x in v.items() if k in variants}
        for name, (f, ares, fopts, logg0) in v.items():
            jobs.append(dict(group=f'{r.target_id}__{name}', target_id=r.target_id, variant=name, file=f,
                             snr_coadd=np.nan, n_used=np.nan, teff0=float(r.tic_teff), logg0=float(logg0),
                             linelist='sousa', ares_options=ares, fasma_options=fopts,
                             workdir=os.path.join(OUT, 'work', r.target_id, name),
                             home=home, fasma=fasma, bindir=bindir))
    return jobs


def run_variants(only=None, variants=None):
    s = pd.read_csv(os.path.join(OUT, 'stars.csv'))
    if only:
        s = s[s.target_id.isin(only)]
    jobs = _jobs(s, variants)
    missing = [j['file'] for j in jobs if not os.path.exists(j['file'])]
    if missing:
        sys.exit(f'{len(missing)} input spectra missing (run the degrade stage first?), e.g. {missing[0]}')
    print(f'{len(jobs)} fits ({len(s)} stars), {ew._workers()} workers', flush=True)
    out = pd.DataFrame(ew._run_fits(jobs))
    path = os.path.join(OUT, 'results.csv')
    if os.path.exists(path):   # a partial rerun replaces only the fits it redid
        old = pd.read_csv(path)
        out = pd.concat([old[~old.group.isin(set(out.group))], out], ignore_index=True)
    out.to_csv(path, index=False)
    print(out.groupby('variant').status.apply(lambda x: (x == 'ok').sum()).to_string())


# ----------------------------------------------------------------------------- report

ORDER = ['he_base', 'he_snr_x1.5', 'he_snr_x2.5', 'he_space3', 'he_space4', 'he_fixlogg',
         'harps_deg39k', 'harps_deg39k_noise', 'harps_snr_x0.67', 'harps_fixlogg', 'he_pedcorr']


def run_report():
    s = pd.read_csv(os.path.join(OUT, 'stars.csv')).set_index('target_id')
    r = pd.read_csv(os.path.join(OUT, 'results.csv'))
    r = r.join(s[['harps_feh', 'harps_teff', 'harps_logg', 'tic_logg']], on='target_id')
    r['ok'] = (r.status == 'ok') & (r.logg < 4.99)
    r['d_feh'] = r.feh - r.harps_feh
    r['d_teff'] = r.teff - r.harps_teff
    r['d_logg'] = r.logg - r.harps_logg
    rows = []
    for v in ORDER:
        x = r[r.variant == v]
        if not len(x):
            continue
        y = x[x.ok]
        mad = lambda z: 1.4826 * np.median(np.abs(z - np.median(z))) if len(z) else np.nan
        rows.append(dict(variant=v, n=len(x), n_ok=len(y),
                         d_feh_median=np.median(y.d_feh) if len(y) else np.nan, d_feh_scatter=mad(y.d_feh),
                         d_teff_median=np.median(y.d_teff) if len(y) else np.nan,
                         d_logg_median=np.median(y.d_logg) if len(y) else np.nan,
                         vt_median=np.median(y.vt) if len(y) else np.nan))
    rep = pd.DataFrame(rows)
    rep.to_csv(os.path.join(OUT, 'report.csv'), index=False)
    wide = r.pivot_table(index='target_id', columns='variant', values='d_feh')[[v for v in ORDER if v in set(r.variant)]]
    wide.to_csv(os.path.join(OUT, 'report_per_star.csv'))
    print('Differences from each star\'s HARPS fit in the main run (converged fits only; scatter = 1.48 x MAD):')
    print(rep.round({'d_feh_median': 3, 'd_feh_scatter': 3, 'd_teff_median': 0, 'd_logg_median': 2,
                     'vt_median': 2}).to_string(index=False))
    print('\nper star, [Fe/H] minus HARPS (NaN = did not converge or ran to logg = 5):')
    print(wide.round(2).to_string())
    print(f'\n-> {OUT}/report.csv, report_per_star.csv, results.csv')


# ----------------------------------------------------------------------------- pedestal (additive light in SOPHIE HE)

def _read_spec(path):
    from astropy.io import fits
    with fits.open(path) as h:
        flux = np.asarray(h[0].data, float)
        hdr = h[0].header.copy()
    return hdr['CRVAL1'] + hdr['CDELT1'] * np.arange(flux.size), flux, hdr


def _robust_line(x, y, n_iter=4, clip=3.0):
    """y = a x + b by least squares with sigma clipping; returns a, b, rms, n."""
    ok = np.ones(x.size, bool)
    for _ in range(n_iter):
        a, b = np.polyfit(x[ok], y[ok], 1)
        res = y - (a * x + b)
        rms = 1.4826 * np.median(np.abs(res[ok]))
        new = np.abs(res) < clip * rms
        if new.sum() == ok.sum():
            break
        ok = new
    return a, b, rms, int(ok.sum())


def run_pedestal(chunk=100.0, lo=ew.EW_LO, hi=ew.EW_HI):
    """Measure an additive light level in each SOPHIE HE co-add against the same star's HARPS co-add at R = 39,000.

    If SOPHIE = scale x (HARPS + p) / (1 + p) (a featureless pedestal p, in units of the continuum: scattered light,
    background, sky), then per chunk SOPHIE = a x HARPS_39k + b with p = b / a, whatever the local continuum scale.
    Lines are made shallower by 1 / (1 + p), so EW ratio = 1 / (1 + p). Writes spectra/<star>/he_pedcorr.fits =
    SOPHIE - b(wavelength) (a smooth fit to the chunk values), for the he_pedcorr variant."""
    from astropy.io import fits
    s = pd.read_csv(os.path.join(OUT, 'stars.csv'))
    rows, summary = [], []
    for r in s.itertuples():
        ws, fs, hdr = _read_spec(r.he_file)
        degp = os.path.join(OUT, 'spectra', r.target_id, 'harps_deg39k.fits')
        if not os.path.exists(degp):
            print(f'{r.target_id}: no degraded HARPS spectrum (run the degrade stage)')
            continue
        wh, fh, _ = _read_spec(degp)
        okh = fh > 0

        def harps_at(v):   # HARPS_39k at SOPHIE's wavelengths, shifted by v km/s; NaN where missing
            w = ws / (1 + v / 299792.458)
            f = np.interp(w, wh, np.where(okh, fh, np.nan), left=np.nan, right=np.nan)
            return f

        # residual RV mismatch between the two rest-frame co-adds (both corrected by their own CCF RVs)
        sel = (ws > 5000) & (ws < 6500) & (fs > 0)
        best = (-np.inf, 0.0)
        for v in np.arange(-4.0, 4.001, 0.1):
            h = harps_at(v)
            m = sel & np.isfinite(h)
            c = np.corrcoef(fs[m], h[m])[0, 1]
            best = max(best, (c, v))
        v = best[1]
        h = harps_at(v)
        good = (fs > 0) & np.isfinite(h)

        ch = []
        for c0 in np.arange(lo, hi, chunk):
            m = good & (ws >= c0) & (ws < c0 + chunk)
            if m.sum() < 0.5 * chunk / np.median(np.diff(ws)) or np.std(h[m]) < 0.02:
                continue
            a, b, rms, n = _robust_line(h[m], fs[m])
            ch.append(dict(target_id=r.target_id, wave=c0 + chunk / 2, a=a, b=b, p=b / a, rms=rms, n=n))
        if len(ch) < 5:
            print(f'{r.target_id}: too few usable chunks')
            continue
        c = pd.DataFrame(ch)
        rows.append(c)
        # smooth b(wavelength): quadratic, clipped
        x = (c.wave.values - 5500.0) / 1000.0
        keep = np.ones(len(c), bool)
        for _ in range(3):
            coef = np.polyfit(x[keep], c.b.values[keep], 2)
            res = c.b.values - np.polyval(coef, x)
            sig = 1.4826 * np.median(np.abs(res[keep])) + 1e-4
            keep = np.abs(res) < 3 * sig
        bfit = np.polyval(coef, (ws - 5500.0) / 1000.0)
        corr = np.where(fs > 0, fs - bfit, 0.0)
        corr = np.where(corr > 0, corr, 0.0)
        hdr['COMMENT'] = f'SOPHIE HE co-add minus fitted additive pedestal (sophie_test.py pedestal, v={v:+.1f} km/s)'
        out = os.path.join(OUT, 'spectra', r.target_id, 'he_pedcorr.fits')
        fits.PrimaryHDU(corr.astype(np.float32), hdr).writeto(out, overwrite=True)
        blue, red = c[c.wave < 5000].p.median(), c[c.wave > 6400].p.median()
        summary.append(dict(target_id=r.target_id, v_shift=v, corr=best[0], p_median=c.p.median(),
                            p_blue=blue, p_red=red, predicted_ew_ratio=1 / (1 + c.p.median()),
                            fit_rms=c.rms.median(), n_chunks=len(c)))
    pd.concat(rows, ignore_index=True).to_csv(os.path.join(OUT, 'pedestal_chunks.csv'), index=False)
    t = pd.DataFrame(summary)
    lp = os.path.join(OUT, 'lines_per_star.csv')
    if os.path.exists(lp):
        l = pd.read_csv(lp)
        t = t.merge(l[l.variant == 'he_base'][['target_id', 'ew_ratio']].rename(columns={'ew_ratio': 'measured_ew_ratio'}),
                    on='target_id', how='left')
    t = t.merge(s[['target_id', 'he_snr']], on='target_id', how='left')
    t.to_csv(os.path.join(OUT, 'pedestal.csv'), index=False)
    print('Additive pedestal p (fraction of the continuum) in each SOPHIE HE co-add, from SOPHIE = a HARPS_39k + b:')
    print(t.round({'v_shift': 1, 'corr': 4, 'p_median': 3, 'p_blue': 3, 'p_red': 3, 'predicted_ew_ratio': 3,
                   'fit_rms': 4, 'measured_ew_ratio': 3, 'he_snr': 0}).to_string(index=False))
    print(f'\n-> {OUT}/pedestal.csv, pedestal_chunks.csv, spectra/<star>/he_pedcorr.fits')


# ----------------------------------------------------------------------------- lines (line-by-line EWs)

def _read_ares(wd):
    """ARES output of one fit (wavelength, depth, FWHM, EW in mA); falls back to FASMA's linelist/star.moog (EW only)."""
    import glob
    for p in sorted(glob.glob(os.path.join(wd, '*.ares')) + glob.glob(os.path.join(wd, '*', '*.ares'))):
        a = np.loadtxt(p, ndmin=2)
        if len(a):
            return pd.DataFrame(dict(wave=a[:, 0], depth=a[:, 2], fwhm=a[:, 3], ew=a[:, 4]))
    p = os.path.join(wd, 'linelist', 'star.moog')
    if os.path.exists(p):
        a = np.loadtxt(p, skiprows=1, ndmin=2)
        return pd.DataFrame(dict(wave=a[:, 0], depth=np.nan, fwhm=np.nan, ew=a[:, -1]))
    return None


def _match(a, b, tol=0.05):
    a, b = a.sort_values('wave'), b.sort_values('wave')
    m = pd.merge_asof(a, b, on='wave', direction='nearest', tolerance=tol, suffixes=('', '_ref'))
    return m[(m.ew > 0) & (m.ew_ref > 0) & np.isfinite(m.ew_ref)]


def run_lines(variants=('he_base', 'he_pedcorr', 'harps_deg39k', 'harps_deg39k_noise')):
    """Compare each variant's ARES EWs line by line with the same star's main-run HARPS EWs.

    The pattern says what kind of problem it is:
      EW ratio < 1 for all lines alike, depth ratio ~ EW ratio -> an additive pedestal (scattered light / background
                                                                   / sky) or a continuum normalised too low
      ratio falls with line strength                          -> saturated/blended lines (resolution, line wings)
      ratio changes with wavelength                           -> order merging / blaze / normalisation"""
    import warnings
    warnings.filterwarnings('ignore', 'All-NaN slice', RuntimeWarning)   # depth/FWHM absent when only star.moog exists
    s = pd.read_csv(os.path.join(OUT, 'stars.csv'))
    rows, pooled = [], []
    for r in s.itertuples():
        ref = _read_ares(os.path.join(ew.WORK, f'{r.target_id}__{HARPS}'))
        if ref is None:
            print(f'{r.target_id}: no main-run HARPS ARES output in {ew.WORK}')
            continue
        for v in variants:
            x = _read_ares(os.path.join(OUT, 'work', r.target_id, v))
            if x is None:
                continue
            m = _match(x, ref)
            m = m[(m.ew_ref >= 5) & (m.ew_ref <= 150)]
            if len(m) < 20:
                continue
            lr = np.log10(m.ew / m.ew_ref)
            slope_ew = np.polyfit(np.log10(m.ew_ref), lr, 1)[0]
            slope_wl = np.polyfit(m.wave / 1000.0, lr, 1)[0]
            rows.append(dict(target_id=r.target_id, variant=v, n=len(m),
                             ew_ratio=10 ** np.median(lr),
                             depth_ratio=np.nanmedian(m.depth / m.depth_ref),
                             fwhm_ratio=np.nanmedian(m.fwhm / m.fwhm_ref),
                             weak_ratio=10 ** np.median(lr[m.ew_ref < 30]) if (m.ew_ref < 30).sum() > 5 else np.nan,
                             strong_ratio=10 ** np.median(lr[m.ew_ref > 70]) if (m.ew_ref > 70).sum() > 5 else np.nan,
                             blue_ratio=10 ** np.median(lr[m.wave < 5500]) if (m.wave < 5500).sum() > 5 else np.nan,
                             red_ratio=10 ** np.median(lr[m.wave > 6100]) if (m.wave > 6100).sum() > 5 else np.nan,
                             dlogratio_dlogew=slope_ew, dlogratio_per_kA=slope_wl))
            pooled.append(m.assign(target_id=r.target_id, variant=v, lr=lr))
    if not rows:
        print('no matched line lists found')
        return
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'lines_per_star.csv'), index=False)
    p = pd.concat(pooled, ignore_index=True)
    p.to_csv(os.path.join(OUT, 'lines_matched.csv'), index=False)
    print('Line-by-line EWs vs the main-run HARPS EWs (lines with 5-150 mA in HARPS; ratios are medians):')
    print('\nmedian over stars:')
    print(t.groupby('variant')[[c for c in t.columns if c not in ('target_id', 'variant')]].median()
          .reindex([v for v in variants if v in set(t.variant)]).round(3).to_string())
    print('\nper star (EW ratio / depth ratio / weak-line ratio / strong-line ratio):')
    w = t.pivot(index='target_id', columns='variant', values=['ew_ratio', 'depth_ratio', 'weak_ratio', 'strong_ratio'])
    print(w.round(3).to_string())
    print('\npooled EW ratio by HARPS EW (mA) and by wavelength (A):')
    p['ew_bin'] = pd.cut(p.ew_ref, [5, 15, 30, 50, 70, 100, 150])
    p['wl_bin'] = pd.cut(p.wave, [4400, 5000, 5500, 6000, 6500, 7000])
    for b in ('ew_bin', 'wl_bin'):
        g = p.groupby(['variant', b], observed=True).lr.median().unstack(0)
        print((10 ** g)[[v for v in variants if v in g.columns]].round(3).to_string(), '\n')
    print(f'-> {OUT}/lines_per_star.csv, lines_matched.csv')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('stage', choices=['select', 'degrade', 'run', 'report', 'pedestal', 'lines'])
    p.add_argument('--only', default=None, help='run: comma-separated target_ids')
    p.add_argument('--variants', default=None, help=f'run: comma-separated subset of {",".join(ORDER)}')
    a = p.parse_args()
    if a.stage == 'select':
        run_select()
    elif a.stage == 'degrade':
        run_degrade()
    elif a.stage == 'run':
        run_variants(a.only.split(',') if a.only else None, a.variants.split(',') if a.variants else None)
    elif a.stage == 'pedestal':
        run_pedestal()
    elif a.stage == 'lines':
        run_lines()
    else:
        run_report()


if __name__ == '__main__':
    main()
