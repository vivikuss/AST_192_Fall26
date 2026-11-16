"""Spectroscopic [Fe/H] for EBLM primaries from co-added CORALIE/HARPS s1d spectra, plus hot-Jupiter host samples.

The fit is a direct comparison of each co-added spectrum against interpolated PHOENIX-ACES synthetic spectra
(Husser et al. 2013) in a few Fe-line-rich windows, with Teff/logg priors from the TESS Input Catalog.
"""
import glob
import io
import os
import re
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.ndimage import maximum_filter1d
from scipy.optimize import least_squares

D = 'eblm_metallicity_data'
SPEC_DIRS = {'CORALIE07': 'EBLM_DATA/combined_s1d_eblm07', 'CORALIE14': 'EBLM_DATA/combined_s1d_eblm14',
             'HARPS': 'EBLM_DATA/combined_s1d_harps'}
RESOLUTION = {'CORALIE07': 55000, 'CORALIE14': 60000, 'HARPS': 115000}
C_KMS = 299792.458

PHX_URL = 'https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/'
TEFFS = np.r_[np.arange(4200, 7001, 100), 7200, 7400].astype(float)   # PHOENIX steps by 200 K above 7000 K
LOGGS = np.array([3.5, 4.0, 4.5, 5.0])
FEHS = np.array([-1.5, -1.0, -0.5, 0.0, 0.5])
# Fe-line-rich windows (air, Angstrom) that avoid Mg b, Na D, H-alpha and the strongest telluric bands
WINDOWS = [(5040, 5140), (5215, 5320), (5370, 5480), (5560, 5660), (6000, 6100), (6130, 6260)]
PAD = 4.0          # Angstrom of model padding on each side of a window, for broadening and RV shifts
DV_MODEL = 1.0     # km/s per pixel of the log-lambda model grid
POLY_DEG = 3       # per-window continuum polynomial


# ----------------------------------------------------------------------------- targets and spectra

def jname(path):
    """Short EBLM designation (e.g. J0008+02) and full coordinates (if the filename has them) of a spectrum file."""
    b = os.path.basename(path).replace('_combined.fits', '')
    m = re.match(r'(\d{2})(\d{2})(\d{2})([NS])(\d{2})(\d{2})(\d{2})$', b)
    if m:
        sign = '+' if m.group(4) == 'N' else '-'
        ra = f'{m.group(1)}h{m.group(2)}m{m.group(3)}s'
        dec = f'{sign}{m.group(5)}d{m.group(6)}m{m.group(7)}s'
        return f'J{m.group(1)}{m.group(2)}{sign}{m.group(5)}', ra, dec
    return b, None, None


def spectrum_files():
    """One row per (target, instrument) spectrum."""
    rows = []
    for inst, d in SPEC_DIRS.items():
        for f in sorted(glob.glob(f'{d}/*_combined.fits')):
            name, ra, dec = jname(f)
            hdr = fits.getheader(f)
            rows.append({'name': name, 'instrument': inst, 'file': f, 'ra_str': ra, 'dec_str': dec,
                         'n_exposures': sum(k.startswith('FILE') for k in hdr)})
    return pd.DataFrame(rows)


def read_spectrum(path):
    h = fits.open(path)[0]
    flux = h.data.astype(float)
    wave = h.header['CRVAL1'] + h.header['CDELT1'] * np.arange(flux.size)
    return wave, flux


# ----------------------------------------------------------------------------- PHOENIX grid

def vac_to_air(w):
    """Vacuum to air wavelengths (Angstrom), Morton (2000) / VALD."""
    s2 = (1e4 / w) ** 2
    return w / (1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2))


def model_wavelengths():
    """Log-lambda (DV_MODEL km/s) air-wavelength grid of each padded window."""
    return [np.exp(np.arange(np.log(lo - PAD), np.log(hi + PAD), DV_MODEL / C_KMS)) for lo, hi in WINDOWS]


def pseudo_continuum(w, f):
    """Upper envelope (2 A running max, smoothed by a cubic) used to normalise model spectra."""
    env = maximum_filter1d(f, size=max(3, int(2.0 / np.median(np.diff(w)))))
    x = (w - w.mean()) / np.ptp(w)
    keep = np.ones_like(w, bool)
    for _ in range(5):
        c = np.polyfit(x[keep], env[keep], 3)
        keep = env >= np.polyval(c, x) * 0.995
    return np.polyval(c, x)


def _phoenix_name(teff, logg, feh):
    z = '-0.0' if feh == 0 else f'{feh:+.1f}'
    return f'PHOENIX-ACES-AGSS-COND-2011/Z{z}/lte{int(teff):05d}-{logg:.2f}{z}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'


def _fetch(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read()


def build_phoenix_grid(path=f'{D}/phoenix_grid.npz', workers=8):
    """Download PHOENIX HiRes spectra once and keep only the normalised fitting windows (~100 MB on disk)."""
    if os.path.exists(path):
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    wvac = fits.getdata(io.BytesIO(_fetch(PHX_URL + 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')))
    sel = (wvac > 5000) & (wvac < 6300)
    wair = vac_to_air(wvac[sel])
    mw = model_wavelengths()
    grid = np.full((TEFFS.size, LOGGS.size, FEHS.size, sum(w.size for w in mw)), np.nan, np.float32)

    def one(idx):
        i, j, k = idx
        f = fits.getdata(io.BytesIO(_fetch(PHX_URL + _phoenix_name(TEFFS[i], LOGGS[j], FEHS[k]))))[sel]
        out = []
        for w in mw:
            fw = np.interp(w, wair, f)
            out.append(fw / pseudo_continuum(w, fw))
        grid[i, j, k] = np.concatenate(out)

    idx = [(i, j, k) for i in range(TEFFS.size) for j in range(LOGGS.size) for k in range(FEHS.size)]
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(one, idx))
    if np.isnan(grid).any():
        raise RuntimeError(f'{np.isnan(grid[..., 0]).sum()} PHOENIX models failed to download; re-run')
    np.savez_compressed(path, grid=grid, teff=TEFFS, logg=LOGGS, feh=FEHS,
                        sizes=np.array([w.size for w in mw]))


def _natural_spline_m2(y, x, axis):
    """Second derivatives of the natural cubic spline through y(x) along `axis` (zero at both ends), for every
    other index at once. x may be non-uniform."""
    y = np.moveaxis(np.asarray(y, float), axis, -1)
    n = x.size
    m2 = np.zeros_like(y)
    if n > 2:
        h = np.diff(x)
        A = np.zeros((n - 2, n - 2))
        for r in range(n - 2):
            A[r, r] = 2 * (h[r] + h[r + 1])
            if r > 0:
                A[r, r - 1] = h[r]
            if r < n - 3:
                A[r, r + 1] = h[r + 1]
        d = np.diff(y, axis=-1) / h
        rhs = 6 * (d[..., 1:] - d[..., :-1])
        m2[..., 1:-1] = rhs @ np.linalg.inv(A).T
    return np.moveaxis(m2, -1, axis)


class ModelGrid:
    """PHOENIX grid lookup. Teff and logg are interpolated linearly. [Fe/H] is interpolated with a natural cubic
    spline through all grid points by default: the grid steps 0.5 dex in [Fe/H], and linear interpolation there
    puts a kink in the model at every grid point, which traps fits exactly on a grid point ([Fe/H] = 0.000 for ~13% of
    the Nov 2026 fits). `feh_interp='linear'` restores the old behaviour."""

    def __init__(self, path=f'{D}/phoenix_grid.npz', feh_interp='spline'):
        z = np.load(path)
        self.grid, self.teff, self.logg, self.feh = z['grid'], z['teff'], z['logg'], z['feh']
        self.waves = model_wavelengths()
        self.slices = np.split(np.arange(self.grid.shape[-1]), np.cumsum(z['sizes'])[:-1])
        self.bounds = ([self.teff[0], self.logg[0], self.feh[0]], [self.teff[-1], self.logg[-1], self.feh[-1]])
        if feh_interp not in ('spline', 'linear'):
            raise ValueError(feh_interp)
        self.feh_interp = feh_interp
        if feh_interp == 'spline':
            self.m2 = _natural_spline_m2(self.grid, self.feh, axis=2).astype(np.float32)

    def interpolate(self, teff, logg, feh):
        """Normalised model flux: linear in Teff and logg; spline (default) or linear in [Fe/H]."""
        idx, wts = [], []
        for v, ax in ((teff, self.teff), (logg, self.logg), (feh, self.feh)):
            i = int(np.clip(np.searchsorted(ax, v) - 1, 0, ax.size - 2))
            t = (v - ax[i]) / (ax[i + 1] - ax[i])
            idx.append(i); wts.append(t)
        (i, j, k), (a, b, c) = idx, wts
        cube = self.grid[i:i + 2, j:j + 2, k:k + 2].astype(float)
        # [Fe/H] first (the interpolation is linear in the grid values, so the order of the axes does not matter)
        f = cube[:, :, 0] * (1 - c) + cube[:, :, 1] * c
        if self.feh_interp == 'spline':
            m = self.m2[i:i + 2, j:j + 2, k:k + 2].astype(float)
            h2 = (self.feh[k + 1] - self.feh[k]) ** 2 / 6
            f += h2 * (((1 - c) ** 3 - (1 - c)) * m[:, :, 0] + (c ** 3 - c) * m[:, :, 1])
        f = f[0] * (1 - a) + f[1] * a
        return f[0] * (1 - b) + f[1] * b


def broadening_kernel(vbroad, resolution, eps=0.6):
    """Rotational profile (Gray 2005; vbroad stands in for vsini + macroturbulence) convolved with the instrumental Gaussian."""
    sig = C_KMS / resolution / 2.3548
    half = int(np.ceil(vbroad + 5 * sig)) + 1
    v = np.arange(-half, half + 1) * DV_MODEL
    g = np.exp(-0.5 * (v / sig) ** 2)
    if vbroad > 0.5 * DV_MODEL:
        x = np.clip(1 - (v / vbroad) ** 2, 0, None)
        rot = (2 * (1 - eps) * np.sqrt(x) + 0.5 * np.pi * eps * x)
        g = np.convolve(g, rot / rot.sum(), mode='same')
    return g / g.sum()


# ----------------------------------------------------------------------------- fitting

PARAMS = ['teff', 'logg', 'feh', 'vbroad', 'rv']
SUN = [5772.0, 4.438, 0.0]   # IAU 2015 nominal Teff, logg


class SpectrumFit:
    """Least-squares fit of (a subset of) Teff, logg, [Fe/H], vbroad and RV to one co-added spectrum.

    `line_mask` (one boolean array per window, on the model wavelength grid) restricts the fit to pixels where
    PHOENIX reproduces the solar spectrum, which makes the analysis differential with respect to the Sun.
    """

    def __init__(self, grid, wave, flux, resolution, line_mask=None):
        self.g, self.R = grid, resolution
        self.w, self.f, self.x, self.good = [], [], [], []
        for i, (lo, hi) in enumerate(WINDOWS):
            s = (wave > lo) & (wave < hi) & np.isfinite(flux) & (flux > 0)
            if s.sum() < 100:
                raise ValueError(f'window {lo}-{hi} not covered')
            ww, ff = wave[s], flux[s] / np.median(flux[s])
            self.w.append(ww); self.f.append(ff); self.x.append((ww - ww.mean()) / np.ptp(ww))
            good = np.ones(ww.size, bool)
            if line_mask is not None:
                good = np.interp(ww, grid.waves[i], line_mask[i].astype(float)) > 0.5
            self.good.append(good)
        self.mask = [g.copy() for g in self.good]
        self.sigma = 0.02

    def model(self, p):
        teff, logg, feh, vb, rv = p
        full = self.g.interpolate(teff, logg, feh)
        kern = broadening_kernel(vb, self.R)
        out = []
        for mw, sl, ww, ff, xx, m in zip(self.g.waves, self.g.slices, self.w, self.f, self.x, self.mask):
            fm = np.convolve(full[sl], kern, mode='same')
            fm = np.interp(ww, mw * (1 + rv / C_KMS), fm)
            A = np.vander(xx, POLY_DEG + 1) * fm[:, None]   # continuum polynomial times model, linear in coefficients
            coef = np.linalg.lstsq(A[m], ff[m], rcond=None)[0]
            out.append(A @ coef)
        return out

    def _full(self, q):
        p = self.p.copy()
        p[self.free] = q
        return p

    def residuals(self, q):
        return np.concatenate([((ff - mm) / self.sigma)[m]
                               for ff, mm, m in zip(self.f, self.model(self._full(q)), self.mask)])

    def fit(self, p0, free=(0, 1, 2, 3, 4), n_iter=3, clip=4.0):
        free = list(free)
        lo = np.array(self.g.bounds[0] + [0.0, -15.0])
        hi = np.array(self.g.bounds[1] + [80.0, 15.0])
        self.p, self.free = np.clip(np.asarray(p0, float), lo + 1e-3, hi - 1e-3), free
        q = self.p[free]
        for _ in range(n_iter):
            sol = least_squares(self.residuals, q, bounds=(lo[free], hi[free]),
                                x_scale=np.array([100, 0.1, 0.05, 1.0, 0.5])[free],
                                diff_step=np.array([1e-3, 1e-3, 1e-2, 1e-2, 1e-2])[free])
            q = sol.x
            # re-estimate the per-pixel noise and clip outliers (cosmics, tellurics, badly modelled lines)
            mod = self.model(self._full(q))
            res = np.concatenate([(ff - mm)[m] for ff, mm, m in zip(self.f, mod, self.mask)])
            self.sigma = 1.4826 * np.median(np.abs(res - np.median(res)))
            self.mask = [g & (np.abs(ff - mm) < clip * self.sigma) for g, ff, mm in zip(self.good, self.f, mod)]
        p = self._full(q)
        err = np.full(p.size, np.nan)
        try:
            err[free] = np.sqrt(np.diag(np.linalg.inv(sol.jac.T @ sol.jac)))
        except np.linalg.LinAlgError:
            pass
        n = sum(m.sum() for m in self.mask)
        return dict(zip(PARAMS, p), **dict(zip([k + '_err' for k in PARAMS], err)),
                    snr=1 / self.sigma, chi2_red=2 * sol.cost / max(n - len(free), 1), n_pix=n,
                    logg_fixed=1 not in free, success=sol.success)


def load_solar_atlas(path=f'{D}/kurucz2005_solar_windows.npz'):
    """Kitt Peak solar flux atlas (Kurucz 2005): normalised flux, air wavelengths, cut to the fitting windows."""
    if not os.path.exists(path):
        txt = _fetch('http://kurucz.harvard.edu/sun/fluxatlas2005/solarfluxintwl.asc').decode().splitlines()
        wf = np.array([l.split() for l in txt if re.match(r'\s*\d+\.\d+\s+-?\d+\.\d+\s*$', l)], float)
        w, f = wf[:, 0] * 10, wf[:, 1]
        keep = np.zeros(w.size, bool)
        for lo, hi in WINDOWS:
            keep |= (w > lo - PAD) & (w < hi + PAD)
        np.savez_compressed(path, wave=w[keep], flux=f[keep])
    z = np.load(path)
    return z['wave'], z['flux']


def build_line_mask(grid, tol=0.03, grow_kms=4.0):
    """Pixels (on the model grid) where the best-fitting solar PHOENIX model matches the solar atlas to within `tol`.

    Disagreements are smoothed to ~instrumental resolution and grown by `grow_kms` on each side, so that a
    badly modelled line is dropped as a whole rather than only at its core.
    """
    from scipy.ndimage import gaussian_filter1d, binary_dilation
    w, f = load_solar_atlas()
    sf = SpectrumFit(grid, w, f, resolution=500000)
    sun = sf.fit(SUN + [3.0, 0.0], free=[3, 4], clip=1e9)
    mod = sf.model([sun[k] for k in PARAMS])
    mask = []
    for i, (ww, ff, mm) in enumerate(zip(sf.w, sf.f, mod)):
        step = C_KMS * np.median(np.diff(ww) / ww[1:])            # km/s per atlas pixel
        diff = gaussian_filter1d(np.abs(ff - mm), 2.5 / step)     # ~ CORALIE resolution element
        bad = binary_dilation(diff > tol, iterations=max(1, int(grow_kms / step)))
        mask.append(np.interp(grid.waves[i], ww, (~bad).astype(float), left=0, right=0) > 0.5)
    return mask, sun


def fit_spectrum(grid, path, instrument, line_mask, teff0=np.nan, logg=np.nan):
    """Fit Teff, [Fe/H], vbroad, RV with logg fixed (Gaia/TIC) when available; best of two starting [Fe/H]."""
    wave, flux = read_spectrum(path)
    t0 = teff0 if np.isfinite(teff0) else 6000.0
    free = [0, 2, 3, 4] if np.isfinite(logg) else [0, 1, 2, 3, 4]
    g0 = float(np.clip(logg, LOGGS[0], LOGGS[-1])) if np.isfinite(logg) else 4.3
    best = None
    for feh0 in (-0.3, 0.2):
        sf = SpectrumFit(grid, wave, flux, RESOLUTION[instrument], line_mask)
        r = sf.fit([float(np.clip(t0, TEFFS[0] + 50, TEFFS[-1] - 50)), g0, feh0, 8.0, 0.0], free=free)
        if best is None or r['chi2_red'] < best['chi2_red']:
            best, best_fit = r, sf
    return best, best_fit


_WORKER = {}


def _fit_one(args):
    path, instrument, teff0, logg, mask_path = args
    if 'grid' not in _WORKER:
        _WORKER['grid'] = ModelGrid()
        z = np.load(mask_path)
        _WORKER['mask'] = [z[f'w{i}'] for i in range(len(WINDOWS))]
    try:
        res, _ = fit_spectrum(_WORKER['grid'], path, instrument, _WORKER['mask'], teff0, logg)
    except Exception as e:   # keep going; the failure is recorded in the table
        res = {'error': repr(e)}
    return {'file': path, **res}


def fit_all(spec, targets, line_mask, path=f'{D}/eblm_spectroscopic_fits.csv', workers=8):
    """Fit every spectrum in `spec` in parallel; cached to `path`. Delete the file to refit."""
    if os.path.exists(path):
        return pd.read_csv(path)
    mask_path = f'{D}/solar_line_mask.npz'
    np.savez_compressed(mask_path, **{f'w{i}': m for i, m in enumerate(line_mask)})
    t = targets.set_index('name')
    jobs = [(s.file, s.instrument, t.loc[s['name'], 'tic_teff'], t.loc[s['name'], 'tic_logg'], mask_path)
            for _, s in spec.iterrows()]
    from concurrent.futures import ProcessPoolExecutor
    from tqdm.auto import tqdm
    with ProcessPoolExecutor(workers) as ex:
        out = list(tqdm(ex.map(_fit_one, jobs, chunksize=2), total=len(jobs)))
    res = spec.merge(pd.DataFrame(out), on='file')
    res.to_csv(path, index=False)
    return res


# ----------------------------------------------------------------------------- catalogues

def eso_coordinates(path):
    """Pointing of the first exposure of a combined HARPS file, from the ESO archive (for files named only Jhhmm+dd)."""
    import urllib.parse
    dp_id = fits.getheader(path)['FILE000'].replace('_s1d_A.fits', '')
    q = f"select object, ra, dec from dbo.raw where dp_id='{dp_id}'"
    url = 'https://archive.eso.org/tap_obs/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=' + urllib.parse.quote(q)
    r = pd.read_csv(url)
    return r.iloc[0] if len(r) else None


def query_tic(targets):
    """TIC v8 Teff/logg (Gaia-based): brightest TIC star (T < 14) within 20 arcsec.

    The coordinates in the CORALIE filenames are truncated to whole seconds, so they can be ~15 arcsec off."""
    from astroquery.mast import Catalogs
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    rows = []
    for _, t in targets.iterrows():
        c = SkyCoord(t.ra_str, t.dec_str) if isinstance(t.ra_str, str) else SkyCoord(t.ra, t.dec, unit='deg')
        tab = Catalogs.query_region(c, radius=20 * u.arcsec, catalog='TIC').to_pandas()
        tab = tab[(tab.objType == 'STAR') & (tab.Tmag < 14)].sort_values('Tmag')
        if len(tab):
            s = tab.iloc[0]
            rows.append({'name': t['name'], 'tic': s.ID, 'gaia_dr2': s.GAIA, 'ra': s.ra, 'dec': s.dec, 'Tmag': s.Tmag,
                         'Vmag': s.Vmag, 'sep_arcsec': s.dstArcSec, 'tic_teff': s.Teff, 'tic_teff_err': s.e_Teff,
                         'tic_logg': s.logg, 'tic_logg_err': s.e_logg, 'tic_rad': s.rad, 'tic_mass': s.mass})
    return pd.DataFrame(rows)


def load_targets(path=f'{D}/eblm_targets_tic.csv'):
    """One row per EBLM target: coordinates (filename or ESO archive) and TIC priors. Cached after the first call."""
    if os.path.exists(path):
        return pd.read_csv(path)
    spec = spectrum_files()
    t = spec.sort_values('ra_str').drop_duplicates('name')[['name', 'ra_str', 'dec_str']].reset_index(drop=True)
    t['ra'], t['dec'], t['eso_object'] = np.nan, np.nan, None
    for i, r in t[t.ra_str.isna()].iterrows():
        e = eso_coordinates(spec[(spec.name == r['name']) & (spec.instrument == 'HARPS')].file.iloc[0])
        if e is not None:
            t.loc[i, ['ra', 'dec', 'eso_object']] = [e.ra, e.dec, e.object]
    tic = query_tic(t)
    t = t.merge(tic, on='name', how='left', suffixes=('_file', ''))
    t.to_csv(path, index=False)
    return t


def query_hot_jupiters():
    """NASA Exoplanet Archive composite parameters for confirmed planets with P < 10 d."""
    import urllib.parse
    q = ('select pl_name,hostname,pl_orbper,pl_bmassj,pl_bmassprov,pl_radj,pl_rvamp,st_mass,discoverymethod,'
         'disc_facility,st_teff,st_logg,st_met,st_meterr1,st_meterr2,st_metratio,st_met_reflink,sy_vmag,ra,dec '
         'from pscomppars where pl_orbper < 10')
    url = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?format=csv&query=' + urllib.parse.quote(q)
    return pd.read_csv(url)


def reliable_mass(df, tol=3.0):
    """True where the archive mass agrees (within a factor `tol`) with the mass implied by the measured RV
    semi-amplitude. Some composite-table entries carry Earth masses in the Jupiter-mass column (e.g. HD 40307 c,
    listed as 5.9 MJ with K = 2.5 m/s); those come out ~318x too heavy and fail this check."""
    k_pred = 28.4329 * df.pl_bmassj * (df.pl_orbper / 365.25) ** (-1 / 3) * df.st_mass.fillna(1.0) ** (-2 / 3)
    ratio = k_pred / df.pl_rvamp
    return ratio.between(1 / tol, tol)


def is_hot_jupiter(df, mmin=0.3, mmax=13.0, rmin=0.7, pmax=10.0):
    """P < pmax days and a giant: a mass of mmin-mmax MJ that is consistent with the RV amplitude, or, failing
    that, a radius above rmin RJ (and no implausible >13 MJ mass)."""
    ok = reliable_mass(df)
    by_mass = ok & df.pl_bmassj.between(mmin, mmax, inclusive='left')
    by_radius = ~ok & (df.pl_radj > rmin) & ~(df.pl_bmassj >= mmax)
    return (df.pl_orbper < pmax) & (by_mass | by_radius)


def load_sweetcat(path=f'{D}/sweetcat.csv'):
    """SWEET-Cat (Santos et al. 2013; Sousa et al. 2021): homogeneous ARES+MOOG parameters of planet hosts."""
    if not os.path.exists(path):
        urllib.request.urlretrieve('https://sweetcat.iastro.pt/catalog/SWEETCAT_Dataframe.csv', path)
    return pd.read_csv(path)


# ----------------------------------------------------------------------------- combining, calibration, samples

INSTRUMENT_FLOOR = 0.04   # dex; robust scatter between CORALIE07, CORALIE14 and HARPS fits of the same star
CALIB_METHODS = ('ARES+MOOG', 'SME')   # literature analyses used to set the zero point (SWEET-Cat uses ARES+MOOG)


def quality_cuts(fits_tab, targets, snr_min=20, vbroad_max=40.0, dteff_max=500.0, tic_teff_max=7300.0):
    """Flag spectra whose fit is trustworthy.

    Enough S/N; lines not rotationally blended; Teff and [Fe/H] off the grid edges; the star not hotter than the
    grid (TIC Teff); and spectroscopic Teff within `dteff_max` (3x the robust scatter) of the photometric TIC
    Teff, which rejects the occasional fit that slides down the Teff-[Fe/H] degeneracy into a false minimum.
    """
    f = fits_tab.copy()
    tic = f.name.map(targets.set_index('name').tic_teff)
    f['ok'] = ((f.snr >= snr_min) & (f.vbroad <= vbroad_max) &
               (f.teff > TEFFS[0] + 1) & (f.teff < TEFFS[-1] - 1) & (f.feh > FEHS[0] + 0.01) &
               ~(tic > tic_teff_max) & ~((f.teff - tic).abs() > dteff_max))
    return f


def combine_per_star(fits_tab):
    """Inverse-variance mean over the good spectra of each star; per-spectrum error = formal (+) instrument floor."""
    g = fits_tab[fits_tab.ok].copy()
    g['e'] = np.hypot(g.feh_err, INSTRUMENT_FLOOR)
    rows = []
    for name, s in g.groupby('name'):
        w = 1 / s.e ** 2
        rows.append({'name': name, 'feh_raw': np.sum(w * s.feh) / w.sum(), 'feh_stat': 1 / np.sqrt(w.sum()),
                     'teff': np.sum(w * s.teff) / w.sum(), 'logg': s.logg.iloc[0], 'logg_fixed': s.logg_fixed.all(),
                     'vbroad': np.median(s.vbroad), 'snr_best': s.snr.max(), 'n_spectra': len(s),
                     'instruments': '+'.join(sorted(s.instrument))})
    return pd.DataFrame(rows)


def zero_point(stars, lit, methods=CALIB_METHODS, n_boot=10000, seed=1):
    """Median (literature - this work) over stars with high-quality literature analyses, with a bootstrap error."""
    hq = lit[lit.method.str.startswith(methods)].groupby('name').feh.mean()
    d = (hq - stars.set_index('name').feh_raw).dropna()
    rng = np.random.default_rng(seed)
    boot = np.median(rng.choice(d.values, (n_boot, d.size)), axis=1)
    return np.median(d), np.std(boot), d


def hot_jupiter_hosts(nasa, sweetcat, match_arcsec=5.0):
    """One row per hot-Jupiter host: archive [Fe/H] and, where available, the homogeneous SWEET-Cat value."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    hj = nasa[is_hot_jupiter(nasa)].sort_values('pl_bmassj', ascending=False).drop_duplicates('hostname')
    hj = hj.assign(wasp=hj.disc_facility.str.contains('WASP', na=False)).reset_index(drop=True)
    sc = sweetcat.dropna(subset=['RA_NASA', 'DEC_NASA', '[Fe/H]'])
    c1 = SkyCoord(hj.ra.values, hj.dec.values, unit='deg')
    c2 = SkyCoord(sc.RA_NASA.values, sc.DEC_NASA.values, unit='deg')
    idx, sep, _ = c1.match_to_catalog_sky(c2)
    ok = sep.arcsec < match_arcsec
    for col, src in [('sc_feh', '[Fe/H]'), ('sc_feh_err', 'e[Fe/H]'), ('sc_teff', 'Teff'), ('sc_homogeneous', 'SWFlag')]:
        hj[col] = np.where(ok, sc[src].values[idx], np.nan)
    hj['sc_homogeneous'] = hj.sc_homogeneous == 1
    return hj


def load_hot_jupiters(path=f'{D}/nasa_pscomppars_P10.csv'):
    """Cached NASA Exoplanet Archive query (all P < 10 d planets; `is_hot_jupiter` selects the giants)."""
    if not os.path.exists(path):
        query_hot_jupiters().to_csv(path, index=False)
    return pd.read_csv(path)


# ----------------------------------------------------------------------------- obliquities

def _sys_key(s):
    k = re.sub(r'[\s_\-\*]', '', str(s)).lower()
    return re.sub(r'(?<![0-9])0+(?=[0-9])', '', k)   # TEPCat zero-pads: HAT-P-01, HD_003167b


def hot_jupiter_obliquities(hj_hosts):
    """TEPCat preferred obliquities of hot Jupiters (0.3-13 MJ, P < 10 d), with host [Fe/H].

    [Fe/H] is the homogeneous SWEET-Cat value where one exists (the scale the EBLM values are tied to),
    otherwise TEPCat's compiled value; `feh_source` says which.
    """
    import obliquity_load as ol
    import obliquity_stats as ost
    t = ost.load_tepcat_hosts()
    t = t[ol.is_hot_jupiter(t.Mp_MJ, t.P_d)].copy()
    sc = hj_hosts[hj_hosts.sc_homogeneous].copy()
    by_key = {}
    for _, r in sc.iterrows():
        k = _sys_key(r.hostname)
        by_key.setdefault(k, r)
        by_key.setdefault(k[:-1] if k.endswith('a') and k[-2].isdigit() else k, r)   # 'KELT-2 A' -> kelt2
    hits = []
    for s in t.System:
        k = _sys_key(s)
        cands = [k, k[:-1] if k[-1] in 'bcdef' else None]
        hits.append(next((by_key[c] for c in cands if c and c in by_key), None))
    t['sc_feh'] = [h.sc_feh if h is not None else np.nan for h in hits]
    t['sc_feh_err'] = [h.sc_feh_err if h is not None else np.nan for h in hits]
    t['feh'] = t.sc_feh.fillna(t.FeH)
    t['feh_err'] = t.sc_feh_err.where(t.sc_feh.notna(), (t.FeH_up + t.FeH_lo) / 2)
    t['feh_source'] = np.where(t.sc_feh.notna(), 'SWEET-Cat', np.where(t.FeH.notna(), 'TEPCat', ''))
    return t[t.feh.notna()].reset_index(drop=True)


def _num(v):
    """Leading number of a table entry such as '0.07+-0.08 [spectroscopic]' (NaN if none)."""
    m = re.match(r'\s*([-+]?\d*\.?\d+)', str(v))
    return float(m.group(1)) if m else np.nan


def eblm_obliquities(stars, path=f'{D}/eblm_obliquity_literature.csv'):
    """One row per low-mass eclipsing binary with a published obliquity, with primary-star [Fe/H].

    |lambda| comes from the preferred RM row, psi from any row that has it (e.g. J0608-59: Kunovac Hodzic+2020 RM,
    Marcussen & Albrecht 2022 psi). [Fe/H] is this work's calibrated value when the star is in the sample,
    otherwise the value quoted in the obliquity paper (`feh_source`).
    Wells+2026 values could not be checked against the paper (`verified = False`).
    """
    import obliquity_stats as ost
    o = pd.read_csv(path, dtype={'arxiv_id': str})
    rows = []
    for name, g in o.groupby('name', sort=False):
        lam = g[g.lambda_deg.notna() & g.notes.fillna('').str.contains('preferred')]
        lam = lam.iloc[0] if len(lam) else (g[g.lambda_deg.notna()].iloc[0] if g.lambda_deg.notna().any() else None)
        psi = g[g.psi_deg.notna()]
        psi = psi.iloc[0] if len(psi) else None
        first = g.iloc[0]
        rows.append({'name': name, 'alt_name': first.alt_name, 'is_eblm': bool(first.is_eblm),
                     'period_d': first.period_d, 'm2_msun': first.m2_msun, 'teff_paper': _num(first.teff),
                     'abs_lambda': ost.fold(lam.lambda_deg) if lam is not None else np.nan,
                     'lambda_up': lam.lambda_err_up if lam is not None else np.nan,
                     'lambda_lo': lam.lambda_err_lo if lam is not None else np.nan,
                     'psi': psi.psi_deg if psi is not None else np.nan,
                     'psi_up': psi.psi_err_up if psi is not None else np.nan,
                     'psi_lo': psi.psi_err_lo if psi is not None else np.nan,
                     'feh_paper': _num(first.feh), 'reference': '; '.join(g.reference.unique()),
                     'verified': not g.reference.str.contains('Wells\\+2026').any()})
    out = pd.DataFrame(rows)
    s = stars.set_index('name')
    out['feh'] = out.name.map(s.feh).fillna(out.feh_paper)
    out['feh_err'] = out.name.map(s.feh_err).fillna(0.1)
    out['teff'] = out.name.map(s.teff).fillna(out.teff_paper)
    out['feh_source'] = np.where(out.name.isin(s.index), 'this work', np.where(out.feh_paper.notna(), 'paper', ''))
    return out


# ----------------------------------------------------------------------------- additional planets

def load_all_planets(path=f'{D}/nasa_pscomppars_all.csv'):
    """Every confirmed planet in the NASA Exoplanet Archive (composite parameters), cached."""
    if not os.path.exists(path):
        import urllib.parse
        q = ('select pl_name,hostname,ra,dec,sy_pnum,sy_snum,cb_flag,pl_orbper,pl_bmassj,pl_bmassprov,pl_radj,'
             'pl_rvamp,st_mass,discoverymethod,disc_facility,disc_year,pl_orbsmax from pscomppars')
        pd.read_csv('https://exoplanetarchive.ipac.caltech.edu/TAP/sync?format=csv&query=' + urllib.parse.quote(q)) \
          .to_csv(path, index=False)
    return pd.read_csv(path)


def _describe(planets):
    ok = reliable_mass(planets) if len(planets) else []
    return '; '.join(f'{r.pl_name} (P={r.pl_orbper:.4g} d, ' + (f'M={r.pl_bmassj:.3g} MJ' if good else
                     f'R={r.pl_radj:.2g} RJ' if np.isfinite(r.pl_radj) else 'M uncertain') + ')'
                     for r, good in zip(planets.itertuples(), ok))


def eblm_planets(targets, planets, match_arcsec=10.0):
    """Known (circumbinary) planets around each EBLM target, by position."""
    from astropy.coordinates import SkyCoord
    c1 = SkyCoord(targets.ra.values, targets.dec.values, unit='deg')
    c2 = SkyCoord(planets.ra.values, planets.dec.values, unit='deg')
    rows = []
    for name, c in zip(targets.name, c1):
        p = planets[c2.separation(c).arcsec < match_arcsec]
        rows.append({'name': name, 'n_planets': len(p), 'planets': _describe(p)})
    return pd.DataFrame(rows)


def hot_jupiter_companions(hj_hosts, planets):
    """Other known planets in each hot-Jupiter system (same archive host name, excluding the hot Jupiter itself)."""
    rows = []
    planets = planets.assign(giant=(reliable_mass(planets) & (planets.pl_bmassj > 0.3)) |
                                   (~reliable_mass(planets) & (planets.pl_radj > 0.7)))
    for r in hj_hosts.itertuples():
        others = planets[(planets.hostname == r.hostname) & (planets.pl_name != r.pl_name)]
        inner = others[others.pl_orbper < r.pl_orbper]
        rows.append({'hostname': r.hostname, 'n_other': len(others), 'has_companion': len(others) > 0,
                     'n_inner': len(inner), 'n_outer': len(others) - len(inner),
                     'n_outer_giant': ((others.pl_orbper > r.pl_orbper) & others.giant).sum(),
                     'nearby': ((others.pl_orbper / r.pl_orbper).between(1 / 4, 4)).any(),   # within a factor 4 in P
                     'companions': _describe(others)})
    return hj_hosts.merge(pd.DataFrame(rows), on='hostname')


def binomial_ci(k, n, level=0.68):
    """Clopper-Pearson interval on a fraction k/n."""
    from scipy.stats import beta
    a = (1 - level) / 2
    lo = beta.ppf(a, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - a, k + 1, n - k) if k < n else 1.0
    return lo, hi


# ----------------------------------------------------------------------------- unconfirmed candidates

CANDIDATE_TABLES = {
    'TOI': ('toi', 'toi,tid,tfopwg_disp,pl_orbper,pl_rade,ra,dec,toi_created', "tfopwg_disp in ('PC','APC')"),
    'KOI': ('cumulative', 'kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_prad,koi_score,ra,dec',
            "koi_disposition = 'CANDIDATE'"),
    'K2': ('k2pandc', 'pl_name,epic_hostname,disposition,pl_orbper,pl_rade,ra,dec', "disposition = 'CANDIDATE' and default_flag = 1"),
}


def load_candidates(path=f'{D}/archive_candidates.csv'):
    """Unconfirmed transit candidates from the NASA Exoplanet Archive: TESS TOIs (TFOPWG disposition PC or APC),
    Kepler KOIs (CANDIDATE) and K2 candidates (CANDIDATE). One row per candidate signal, cached."""
    if os.path.exists(path):
        return pd.read_csv(path)
    import urllib.parse
    out = []
    for survey, (table, cols, where) in CANDIDATE_TABLES.items():
        q = f'select {cols} from {table} where {where}'
        t = pd.read_csv('https://exoplanetarchive.ipac.caltech.edu/TAP/sync?format=csv&query=' + urllib.parse.quote(q))
        name = {'TOI': 'TOI-' + t.get('toi', pd.Series(dtype=str)).astype(str),
                'KOI': t.get('kepoi_name'), 'K2': t.get('pl_name')}[survey]
        out.append(pd.DataFrame({'survey': survey, 'candidate': name,
                                 'tic': t.tid if survey == 'TOI' else np.nan,
                                 'disposition': t.iloc[:, 2],
                                 'period': t.pl_orbper if 'pl_orbper' in t else t.koi_period,
                                 'radius_re': t.pl_rade if 'pl_rade' in t else t.koi_prad,
                                 'ra': t.ra, 'dec': t.dec}))
    c = pd.concat(out, ignore_index=True)
    c.to_csv(path, index=False)
    return c


def _same_period(p, known, tol=0.01):
    """p matches any known period (or 1/2x, 2x of it) to within a fractional `tol`."""
    known = np.asarray([k for k in known if np.isfinite(k)])
    if not np.isfinite(p) or known.size == 0:
        return False
    return bool(np.any(np.abs(p / np.concatenate([known, known / 2, known * 2]) - 1) < tol))


def match_candidates(systems, cands, known_periods, match_arcsec=5.0):
    """Candidates within `match_arcsec` of each system, minus those whose period matches a known body.

    systems : DataFrame with 'key', 'ra', 'dec'; known_periods : dict key -> list of periods already accounted for
    (confirmed planets, or the eclipsing binary itself).
    """
    from astropy.coordinates import SkyCoord
    cs = SkyCoord(cands.ra.values, cands.dec.values, unit='deg')
    rows = []
    for r in systems.itertuples():
        if not (np.isfinite(r.ra) and np.isfinite(r.dec)):
            continue
        near = cands[cs.separation(SkyCoord(r.ra, r.dec, unit='deg')).arcsec < match_arcsec]
        for c in near.itertuples():
            rows.append({'key': r.key, 'survey': c.survey, 'candidate': c.candidate, 'disposition': c.disposition,
                         'period': c.period, 'radius_re': c.radius_re,
                         'is_known_body': _same_period(c.period, known_periods.get(r.key, []))})
    return pd.DataFrame(rows, columns=['key', 'survey', 'candidate', 'disposition', 'period', 'radius_re',
                                       'is_known_body'])
