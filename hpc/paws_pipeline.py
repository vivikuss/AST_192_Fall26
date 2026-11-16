"""Headless re-run of PAWS (Freckelton et al. 2024, BEBOP V; github.com/alixviolet/PAWS) on the archival spectra,
as an independent check on the PHOENIX-fit metallicities from metallicity_pipeline.py.

PAWS is an ipywidgets GUI around iSpec (EW stage: ARES-style Gaussian EWs + WIDTH, Fe I/II, ATLAS9; synthesis stage:
SPECTRUM, ATLAS9, Teff/logg/[M/H]/vsini free). Its analysis functions read their settings from the GUI widgets, so
`run` executes PAWS's own coadding.py / pipeline.py / synth.py unchanged, with stand-in widget objects that carry the
settings, one co-add group per subprocess.

Stages (after metallicity_pipeline.py manifest + download):
  prepare   pick groups; write each exposure in iSpec text format (nm) into $METAL_ROOT/paws/<group>/   (ast192 env)
  run       PAWS co-add -> EW -> synthesis for each group                                              (paws env)
  collect   gather results -> $METAL_ROOT/results/paws_results.csv                                     (ast192 env)

Differences between the PAWS code (commit 8ec055c) and the paper, and how they are handled:
  * the code co-adds exposures with S/N > 10, the paper says S/N >= 20: `prepare --snr-min` (default 20) applies the
    paper's cut before PAWS sees the data;
  * the paper fixes [M/H] at the EW value in the synthesis stage, the code leaves it free: default = the code as
    published, `run --fix-mh` = the paper's description;
  * for R <= 50000 (SOPHIE HE) the synthesis stage asks iSpec for 47000_GES/spectrum_ew_ispec_good_for_params_all_
    extended.txt, which does not exist in iSpec v2023.08.04; it is redirected to spectrum_synth_good_for_params_all_
    extended.txt (the GES synthesis mask for SPECTRUM).
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import time
import types

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.abspath(os.environ.get('METAL_ROOT', os.path.join(HERE, 'metallicity_hpc_data')))
PAWS_ROOT = os.path.join(ROOT, 'paws')
PAWS_HOME = os.path.abspath(os.environ.get('PAWS_HOME', os.path.join(ROOT, 'paws_software')))
ISPEC_DIR = os.environ.get('ISPEC_DIR', os.path.join(PAWS_HOME, 'iSpec_v20230804'))
PAWS_DIR = os.environ.get('PAWS_DIR', os.path.join(PAWS_HOME, 'PAWS'))

WAVE_NM = (420.0, 680.0)            # PAWS / BEBOP V default range
SPECTRAL_TYPE = 'G2'                # BEBOP V used the G2 starting point for all stars
DEFAULT_INSTRUMENTS = ('SOPHIE_HE', 'SOPHIE_HR', 'HARPS')   # the instruments PAWS was built and validated for
GROUP_TIMEOUT_H = 10.0


# ----------------------------------------------------------------------------- prepare (ast192 env)

def published_paws_targets(t):
    """target_ids of EBLMs with a published PAWS [Fe/H] in eblm_literature_feh.csv."""
    lit = pd.read_csv(os.path.join(HERE, 'eblm_metallicity_data', 'eblm_literature_feh.csv'))
    names = set(lit[lit.method.str.startswith('PAWS')].name)
    eb = t[t['sample'] == 'EBLM']
    return set(eb[eb.name.str.replace('EBLM ', '').isin(names)].target_id)


def run_prepare(which, instruments, snr_min, only=None):
    sys.path.insert(0, os.path.join(HERE, 'hpc'))
    import metallicity_pipeline as mp
    man = pd.read_csv(os.path.join(mp.CAT, 'manifest.csv'))
    t = pd.read_csv(os.path.join(mp.CAT, 'targets.csv')).set_index('target_id', drop=False)
    man = man[man.instrument.isin(instruments) & (man.snr >= snr_min)]
    if which == 'published':
        man = man[man.target_id.isin(published_paws_targets(t.reset_index(drop=True)))]
    if only:
        man = man[man.target_id.isin(only)]
    rows = []
    for g, m in man.groupby('group'):
        d = os.path.join(PAWS_ROOT, g)
        os.makedirs(d, exist_ok=True)
        target = t.loc[m.target_id.iloc[0]]
        n = 0
        for r in m.itertuples():
            out = os.path.join(d, os.path.basename(r.path).replace('.fits', '_s1d.txt.gz'))
            if os.path.exists(out):
                n += 1
                continue
            if not os.path.exists(r.path):
                continue
            w, f, e, jd, _ = mp.read_raw(r.path, r.archive, r.instrument)
            if mp.in_transit(jd, target):
                continue
            wn = w / 10
            s = (wn > WAVE_NM[0] - 1) & (wn < WAVE_NM[1] + 1) & np.isfinite(f) & np.isfinite(e) & (e > 0)
            pd.DataFrame({'waveobs': wn[s], 'flux': f[s], 'err': e[s]}).to_csv(
                out, sep='\t', index=False, float_format='%.7g', compression='gzip')
            n += 1
        rows.append({'group': g, 'target_id': m.target_id.iloc[0], 'instrument': m.instrument.iloc[0],
                     'resolution': int(m.resolution.iloc[0]), 'n_spectra': n, 'dir': d})
    groups = pd.DataFrame(rows)
    groups = groups[groups.n_spectra > 0]
    groups.to_csv(os.path.join(PAWS_ROOT, 'groups.csv'), index=False)
    print(f'{len(groups)} PAWS groups ({groups.instrument.value_counts().to_dict()}) for '
          f'{groups.target_id.nunique()} targets -> {PAWS_ROOT}')


# ----------------------------------------------------------------------------- run (paws env)

class _Widget:
    """Stand-in for any ipywidgets object: holds .value, ignores callbacks and layout, works as a context manager."""

    def __init__(self, *args, value=None, **kw):
        self.value = value
        self.style = types.SimpleNamespace()

    def on_click(self, *a, **k):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _install_stubs(group_dir, resolution):
    widgets = types.ModuleType('ipywidgets')

    def widget_class(name):            # ipywidgets.Button, .Text, ... (but not __file__ etc.)
        if name.startswith('_'):
            raise AttributeError(name)
        return _Widget

    widgets.__getattr__ = widget_class
    widgets.HBox = widgets.VBox = _Widget
    ipy = types.ModuleType('IPython')
    disp = types.ModuleType('IPython.display')
    disp.display = lambda *a, **k: None
    ipy.display = disp
    setup = types.ModuleType('setup')   # replaces PAWS's GUI setup.py
    for k, v in dict(datapathWidget=group_dir, number='Single', coadd='Yes', res_widget=int(resolution),
                     ispecpathWidget=ISPEC_DIR, min_widget=WAVE_NM[0], max_widget=WAVE_NM[1], s1d_widget='s1d',
                     unit='Nanometres', input_widget=SPECTRAL_TYPE, osx='Ubuntu/Linux').items():
        setattr(setup, k, _Widget(value=v))
    sys.modules.update({'ipywidgets': widgets, 'IPython': ipy, 'IPython.display': disp, 'setup': setup})


def _patch_ispec(ispec, fix_mh):
    orig_regions = ispec.read_line_regions

    def read_line_regions(path, *a, **k):
        if not os.path.exists(path) and 'spectrum_ew_ispec_good_for_params_all_extended' in path:
            path = path.replace('spectrum_ew_ispec_good_for_params_all_extended',
                                'spectrum_synth_good_for_params_all_extended')
        return orig_regions(path, *a, **k)

    ispec.read_line_regions = read_line_regions
    if fix_mh:
        orig_model = ispec.model_spectrum

        def model_spectrum(*a, **k):
            a = [[p for p in x if p != 'mh'] if isinstance(x, list) and 'teff' in x else x for x in a]
            return orig_model(*a, **k)

        ispec.model_spectrum = model_spectrum


def _params(path):
    if not os.path.exists(path):
        return {}
    p = pd.read_csv(path)
    return {f'{r.Parameter}': float(r.Value) for r in p.itertuples()} | \
           {f'{r.Parameter}_err': float(r.Errors) for r in p.itertuples()}


def run_one(group_dir, resolution, fix_mh):
    """PAWS on one group folder (runs inside its own subprocess)."""
    import matplotlib
    matplotlib.use('Agg')
    _install_stubs(group_dir, resolution)
    sys.path[:0] = [PAWS_DIR, ISPEC_DIR]
    os.chdir(PAWS_DIR)                  # PAWS opens its solar/ line masks by relative path
    import ispec
    _patch_ispec(ispec, fix_mh)
    t0 = time.time()
    status = {}
    import coadding
    coadding.coadd_spectra(None)
    status['coadd_s'] = time.time() - t0
    if not os.path.exists(os.path.join(group_dir, 'prepared_1_.fits')):
        status['error'] = 'no spectrum passed the PAWS co-add (all S/N <= 10?)'
    else:
        import pipeline
        pipeline.parameters_from_ew(None)
        status['ew_s'] = time.time() - t0 - status['coadd_s']
        import synth
        synth.merged_synthesis(None)
        status['synth_s'] = time.time() - t0 - status['coadd_s'] - status['ew_s']
    res = {'ew_' + k: v for k, v in _params(os.path.join(group_dir, 'params.csv')).items()}
    res |= {'synth_' + k: v for k, v in _params(os.path.join(group_dir, 'params_synth_pipeline.csv')).items()}
    with open(os.path.join(group_dir, 'paws_result.json'), 'w') as f:
        json.dump(res | status | {'fix_mh': fix_mh}, f, indent=1)


def run_tasks(task, ntasks, fix_mh, only=None):
    groups = pd.read_csv(os.path.join(PAWS_ROOT, 'groups.csv'))
    if only:
        groups = groups[groups.target_id.isin(only)]
    for g in groups.sort_values('group').iloc[task::ntasks].itertuples():
        if os.path.exists(os.path.join(g.dir, 'paws_result.json')):
            print(f'{g.group}: done already')
            continue
        t0 = time.time()
        cmd = [sys.executable, os.path.abspath(__file__), '_one', g.dir, str(g.resolution)] + \
              (['--fix-mh'] if fix_mh else [])
        with open(os.path.join(g.dir, 'paws.log'), 'w') as log:
            try:
                rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=GROUP_TIMEOUT_H * 3600).returncode
            except subprocess.TimeoutExpired:
                rc = 'timeout'
        if rc != 0:
            with open(os.path.join(g.dir, 'paws_result.json'), 'w') as f:
                json.dump({'error': f'PAWS exited with {rc}; see paws.log'}, f)
        print(f'{g.group}: rc={rc} in {(time.time() - t0) / 60:.1f} min', flush=True)


# ----------------------------------------------------------------------------- collect (ast192 env)

def run_collect():
    groups = pd.read_csv(os.path.join(PAWS_ROOT, 'groups.csv'))
    rows = []
    for g in groups.itertuples():
        p = os.path.join(g.dir, 'paws_result.json')
        r = json.load(open(p)) if os.path.exists(p) else {'error': 'not run'}
        rows.append({'group': g.group, 'target_id': g.target_id, 'instrument': g.instrument,
                     'resolution': g.resolution, 'n_spectra': g.n_spectra, **r})
    out = pd.DataFrame(rows)
    os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
    out.to_csv(os.path.join(ROOT, 'results', 'paws_results.csv'), index=False)
    ok = out.get('synth_MH', pd.Series(dtype=float)).notna()
    print(f'{ok.sum()}/{len(out)} groups with a PAWS synthesis result; '
          f'errors: {out.get("error", pd.Series(dtype=str)).value_counts().head(5).to_dict()}')


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('stage', choices=['prepare', 'run', 'collect', '_one'])
    p.add_argument('args', nargs='*', help=argparse.SUPPRESS)
    p.add_argument('--targets', choices=['published', 'all'], default='published',
                   help='published: EBLMs with a published PAWS [Fe/H]; all: every group of --instruments')
    p.add_argument('--instruments', default=','.join(DEFAULT_INSTRUMENTS))
    p.add_argument('--snr-min', type=float, default=20.0, help='per-exposure S/N cut (BEBOP V Sect. 4.4)')
    p.add_argument('--fix-mh', action='store_true', help='fix [M/H] at the EW value in synthesis (as in the paper)')
    p.add_argument('--task', type=int, default=None)
    p.add_argument('--ntasks', type=int, default=None)
    p.add_argument('--only', default=None, help='comma-separated target_ids')
    a = p.parse_args()
    only = a.only.split(',') if a.only else None
    task = a.task if a.task is not None else int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
    ntasks = a.ntasks if a.ntasks is not None else int(os.environ.get('SLURM_ARRAY_TASK_COUNT', 1))
    if a.stage == 'prepare':
        run_prepare(a.targets, a.instruments.split(','), a.snr_min, only)
    elif a.stage == 'run':
        run_tasks(task, ntasks, a.fix_mh, only)
    elif a.stage == 'collect':
        run_collect()
    else:
        run_one(a.args[0], int(a.args[1]), a.fix_mh)


if __name__ == '__main__':
    main()
