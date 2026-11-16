"""Which hot-Jupiter hosts and EBLM primaries have public SOPHIE (OHP 1.93 m) spectra, and the combined ESO + SOPHIE
coverage. Run eso_coverage.py first.

The SOPHIE archive (http://atlas.obs-hp.fr/sophie/) has no TAP service, but its CGI returns the whole spectra table
for an RA range, so the full archive is pulled in 24 one-hour slices and cross-matched locally, as for ESO.
"""
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

import eso_coverage as ec

D = 'eblm_metallicity_data'
SOPHIE_CGI = 'http://atlas.obs-hp.fr/sophie/sophie.cgi'
MODES = ('HR', 'HE')   # high resolution (R ~ 75000) and high efficiency (R ~ 40000); both cover 387-694 nm
COLS = ['objname', 'j2000', 'seq', 'slen', 'date', 'mode', 'fiber_b', 'exptime', 'sn26', 'get_spec']
J2000 = re.compile(r'^J\d{6}(\.\d)?[+-]\d{6}$')


def _parse(line):
    """Fields of one data row, located relative to the J2000 designation (some object names contain tabs)."""
    f = [x.strip() for x in line.split('\t')]
    j = next((i for i, x in enumerate(f) if J2000.match(x)), None)
    if j is None or len(f) < j + 12:
        return None
    # after J2000: Simbad link, E link, seq, slen, date, mode, fiber_b, exptime, sn26, view_spec, view_head, get_spec
    return [' '.join(f[:j]), f[j], *f[j + 3:j + 10], f[j + 12] if len(f) > j + 12 else '']


def _slice(h, cache=f'{D}/sophie_slices'):
    """All spectra with B1950 RA in [h, h+1) hours (the CGI's range search is in B1950; the slices tile the sky)."""
    os.makedirs(cache, exist_ok=True)
    path = f'{cache}/ra{h:02d}.csv'
    if os.path.exists(path):
        return pd.read_csv(path, dtype=str)
    r = f'[{h:02d} 00 00 -90],[{h + 1:02d} 00 00 +90]'
    url = f'{SOPHIE_CGI}?n=sophie&a=csv&ob=ra,seq&r={urllib.parse.quote(r)}'
    for attempt in range(3):
        try:
            txt = urllib.request.urlopen(url, timeout=300).read().decode('latin-1')
            break
        except OSError:
            time.sleep(10)
    else:
        raise RuntimeError(f'SOPHIE query failed for RA slice {h}')
    rows = [_parse(l) for l in txt.splitlines() if l and not l.startswith(('#', '$link'))]
    n_expected = int(re.search(r'# (\d+) lines were successfully retrieved', txt).group(1))
    df = pd.DataFrame([r for r in rows if r is not None], columns=COLS)
    if len(df) != n_expected:
        print(f'RA {h:02d}h: parsed {len(df)} of {n_expected} rows')
    df.to_csv(path, index=False)
    return df


def sophie_footprint(path=f'{D}/sophie_footprint.csv'):
    """Position, mode, S/N and public status of every SOPHIE spectrum. Cached after the first call."""
    if os.path.exists(path):
        return pd.read_csv(path)
    with ThreadPoolExecutor(4) as ex:
        parts = list(ex.map(_slice, range(24)))
    fp = pd.concat(parts, ignore_index=True)
    j = fp.j2000.str.extract(r'^J(\d\d)(\d\d)(\d\d(?:\.\d)?)([+-])(\d\d)(\d\d)(\d\d)$').astype({3: str})
    n = j.drop(columns=3).astype(float)
    fp['s_ra'] = 15 * (n[0] + n[1] / 60 + n[2] / 3600)
    fp['s_dec'] = np.where(j[3] == '-', -1, 1) * (n[4] + n[5] / 60 + n[6] / 3600)
    fp['public'] = fp.get_spec == 'get_spec'
    fp['snr'] = pd.to_numeric(fp.sn26, errors='coerce')
    fp = fp[['objname', 's_ra', 's_dec', 'seq', 'date', 'mode', 'exptime', 'snr', 'public']]
    fp.to_csv(path, index=False)
    return fp


def coverage(targets, fp, match_arcsec=ec.MATCH_ARCSEC, snr_min=ec.SNR_MIN):
    """Per target and SOPHIE mode: public spectra, those above `snr_min`, co-added S/N; plus proprietary counts."""
    c_fp = SkyCoord(fp.s_ra.values, fp.s_dec.values, unit='deg')
    c_t = SkyCoord(targets.ra.values, targets.dec.values, unit='deg')
    it, ifp, _, _ = c_fp.search_around_sky(c_t, match_arcsec * u.arcsec)
    m = fp.iloc[ifp].assign(name=targets.name.values[it]).reset_index(drop=True)
    m['snr'] = m.snr.fillna(0)
    out = targets.set_index('name')
    out['n_SOPHIE_private'] = m[~m.public].groupby('name').size()
    pub = m[m.public]
    g = pub.groupby(['name', 'mode'])
    tab = pd.DataFrame({'n': g.size(), 'n_good': g.snr.apply(lambda s: (s >= snr_min).sum()),
                        'snr_coadd': g.snr.apply(lambda s: np.sqrt(np.sum(s ** 2)))}).unstack('mode')
    tab.columns = [f'{a}_SOPHIE_{b}' for a, b in tab.columns]
    out = out.join(tab).reset_index()
    out['n_SOPHIE_private'] = out.n_SOPHIE_private.fillna(0).astype(int)
    for mode in MODES:
        for k in ('n', 'n_good'):
            col = f'{k}_SOPHIE_{mode}'
            out[col] = out[col].fillna(0).astype(int) if col in out else 0
        if f'snr_coadd_SOPHIE_{mode}' not in out:
            out[f'snr_coadd_SOPHIE_{mode}'] = np.nan
    out['sophie_snr_coadd'] = out[[f'snr_coadd_SOPHIE_{m}' for m in MODES]].max(axis=1)
    out['sophie_usable'] = out.sophie_snr_coadd >= 50
    return out


def summarize(cov, eso, label):
    """SOPHIE coverage, and what it adds to ESO."""
    # name + position, because eblm_master_table.csv repeats a few names with different TIC entries
    keep = ['name', 'ra', 'dec', 'usable', 'best_instrument']
    cov = cov.round({'ra': 6, 'dec': 6})
    c = cov.merge(eso[keep].rename(columns={'usable': 'eso_usable'}).round({'ra': 6, 'dec': 6}),
                  on=['name', 'ra', 'dec'], how='left', validate='1:1')
    c['eso_usable'] = c.eso_usable.fillna(False).astype(bool)
    north = c.dec > 25
    print(f'\n=== {label}: {len(c)} stars ({north.sum()} at dec > +25) ===')
    for mode in MODES:
        print(f'  SOPHIE {mode}  public spectra: {(c[f"n_SOPHIE_{mode}"] > 0).sum():4d}   '
              f'with >= 1 spectrum S/N>={ec.SNR_MIN}: {(c[f"n_good_SOPHIE_{mode}"] > 0).sum():4d}')
    print(f'  only proprietary SOPHIE data:    {((c.n_SOPHIE_private > 0) & c.sophie_snr_coadd.isna()).sum():4d}')
    print(f'  SOPHIE usable (co-add S/N >= 50): {c.sophie_usable.sum():4d}')
    both = c.eso_usable | c.sophie_usable
    print(f'  usable in ESO only / SOPHIE only / both: {(c.eso_usable & ~c.sophie_usable).sum()} / '
          f'{(~c.eso_usable & c.sophie_usable).sum()} / {(c.eso_usable & c.sophie_usable).sum()}')
    print(f'  ESO + SOPHIE usable:            {both.sum():4d} / {len(c)} ({both.mean():.0%})')
    print(f'    north (dec > +25):            {(both & north).sum():4d} / {north.sum()}')
    print(f'    south (dec < +25):            {(both & ~north).sum():4d} / {(~north).sum()}')
    return c


if __name__ == '__main__':
    fp = sophie_footprint()
    print(f'{len(fp)} SOPHIE spectra, {fp.public.mean():.0%} public')
    hj = summarize(coverage(ec.hot_jupiter_targets(), fp), pd.read_csv(f'{D}/eso_coverage_hot_jupiters.csv'),
                   'Hot-Jupiter hosts')
    hj.to_csv(f'{D}/sophie_coverage_hot_jupiters.csv', index=False)
    eb = summarize(coverage(ec.eblm_targets(), fp), pd.read_csv(f'{D}/eso_coverage_eblm.csv'), 'EBLMs (master table)')
    eb.to_csv(f'{D}/sophie_coverage_eblm.csv', index=False)
    for flag, lab in [(True, 'already have CORALIE/HARPS co-adds'), (False, 'no spectra yet')]:
        s = eb[eb.have_spectra == flag]
        print(f'  {lab}: {len(s)} stars, {s.sophie_usable.sum()} SOPHIE-usable, '
              f'{(s.sophie_usable | s.eso_usable).sum()} usable in ESO or SOPHIE')
