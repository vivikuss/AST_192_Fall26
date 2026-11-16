"""Which hot-Jupiter hosts and EBLM primaries have public reduced high-resolution spectra in the ESO archive.

The full Phase 3 footprint of HARPS, FEROS, UVES and ESPRESSO (~570k spectra) is pulled once from the ESO TAP
service and cross-matched locally, which is far faster than one cone search per star.
"""
import os

import numpy as np
import pandas as pd
import pyvo
from astropy.coordinates import SkyCoord
import astropy.units as u

import eblm_feh as ef

D = 'eblm_metallicity_data'
ESO_TAP = 'https://archive.eso.org/tap_obs'
INSTRUMENTS = ('HARPS', 'FEROS', 'UVES', 'ESPRESSO')
MATCH_ARCSEC = 10.0   # generous: ESO coordinates are pointing positions, and proper motion moves stars over 20 yr
SNR_MIN = 20          # per-spectrum S/N for a spectrum to count as usable on its own
# most of the eblm_feh fitting windows (5040-6260 A) must be covered; drops e.g. UVES blue-arm-only and NIR setups
WAVE_MIN, WAVE_MAX = 5.15e-7, 6.10e-7


def eso_footprint(path=f'{D}/eso_footprint.csv'):
    """Position, S/N and release date of every Phase 3 spectrum from INSTRUMENTS. Cached after the first call."""
    if os.path.exists(path):
        return pd.read_csv(path)
    svc = pyvo.dal.TAPService(ESO_TAP)
    parts = []
    for inst in INSTRUMENTS:
        q = (f"select s_ra, s_dec, obs_collection, snr, em_min, em_max, obs_release_date from ivoa.ObsCore "
             f"where obs_collection = '{inst}' and dataproduct_type = 'spectrum'")
        parts.append(svc.run_async(q, maxrec=1_000_000).to_table().to_pandas())
        print(f'{inst}: {len(parts[-1])} spectra')
    fp = pd.concat(parts, ignore_index=True)
    fp.to_csv(path, index=False)
    return fp


def coverage(targets, fp, match_arcsec=MATCH_ARCSEC, snr_min=SNR_MIN):
    """Per target and instrument: number of spectra, number above `snr_min`, and the S/N of a co-add of all of them."""
    fp = fp.dropna(subset=['s_ra', 's_dec'])
    fp = fp[(fp.em_min < WAVE_MIN) & (fp.em_max > WAVE_MAX)]
    c_fp = SkyCoord(fp.s_ra.values, fp.s_dec.values, unit='deg')
    c_t = SkyCoord(targets.ra.values, targets.dec.values, unit='deg')
    it, ifp, _, _ = c_fp.search_around_sky(c_t, match_arcsec * u.arcsec)
    m = fp.iloc[ifp].assign(name=targets.name.values[it]).reset_index(drop=True)
    m['snr'] = m.snr.fillna(0)
    g = m.groupby(['name', 'obs_collection'])
    tab = pd.DataFrame({'n': g.size(), 'n_good': g.snr.apply(lambda s: (s >= snr_min).sum()),
                        'snr_coadd': g.snr.apply(lambda s: np.sqrt(np.sum(s ** 2)))}).unstack('obs_collection')
    tab.columns = [f'{a}_{b}' for a, b in tab.columns]
    out = targets.set_index('name').join(tab).reset_index()
    for inst in INSTRUMENTS:
        for k in ('n', 'n_good'):
            col = f'{k}_{inst}'
            out[col] = out[col].fillna(0).astype(int) if col in out else 0
        if f'snr_coadd_{inst}' not in out:
            out[f'snr_coadd_{inst}'] = np.nan
    snr = out[[f'snr_coadd_{i}' for i in INSTRUMENTS]]
    has = snr.notna().any(axis=1)
    out['best_instrument'] = None
    out.loc[has, 'best_instrument'] = snr[has].idxmax(axis=1).str.replace('snr_coadd_', '')
    out['best_snr_coadd'] = snr.max(axis=1)
    out['usable'] = out.best_snr_coadd >= 50   # co-added S/N good enough for a ~0.05 dex [Fe/H]
    return out


def hot_jupiter_targets():
    nasa = ef.load_hot_jupiters()
    hj = nasa[ef.is_hot_jupiter(nasa)].sort_values('pl_bmassj', ascending=False).drop_duplicates('hostname')
    return hj.rename(columns={'hostname': 'name'})[['name', 'pl_name', 'ra', 'dec', 'sy_vmag', 'st_teff',
                                                    'st_met', 'disc_facility']].reset_index(drop=True)


def eblm_targets():
    """All EBLMs in the master table, flagged by whether we already have CORALIE/HARPS co-adds for them."""
    t = pd.read_csv('eblm_master_table.csv').rename(columns={'RA': 'ra', 'Dec': 'dec'})
    have = ef.load_targets().dropna(subset=['tic']).tic.astype(int)
    t['have_spectra'] = t.TIC.isin(have)
    return t[['name', 'TIC', 'ra', 'dec', 'Vmag', 'have_spectra']]


def summarize(cov, label):
    n = len(cov)
    south = cov.dec < 25
    print(f'\n=== {label}: {n} stars ({south.sum()} at dec < +25) ===')
    for inst in INSTRUMENTS:
        k = (cov[f'n_{inst}'] > 0).sum()
        print(f'  {inst:9s} any spectra: {k:4d}   with >= 1 spectrum S/N>={SNR_MIN}: {(cov[f"n_good_{inst}"] > 0).sum():4d}')
    anyc = cov.best_snr_coadd.notna()
    print(f'  any ESO spectrum:             {anyc.sum():4d} ({anyc.mean():.0%})')
    print(f'  co-added S/N >= 50 (usable):  {cov.usable.sum():4d} ({cov.usable.mean():.0%})')
    print(f'    of which south (dec < +25): {(cov.usable & south).sum():4d} / {south.sum()}')
    print(f'    of which north (dec > +25): {(cov.usable & ~south).sum():4d} / {(~south).sum()}')


if __name__ == '__main__':
    fp = eso_footprint()
    hj = coverage(hot_jupiter_targets(), fp)
    hj.to_csv(f'{D}/eso_coverage_hot_jupiters.csv', index=False)
    summarize(hj, 'Hot-Jupiter hosts')
    eb = coverage(eblm_targets(), fp)
    eb.to_csv(f'{D}/eso_coverage_eblm.csv', index=False)
    summarize(eb, 'EBLMs (master table)')
    for flag, lab in [(True, 'already have CORALIE/HARPS co-adds'), (False, 'no spectra yet')]:
        s = eb[eb.have_spectra == flag]
        print(f'  {lab}: {len(s)} stars, {s.usable.sum()} usable in ESO archive')
