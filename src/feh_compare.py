"""[Fe/H] of the 186 shared stars (107 WASP hosts, 79 EBLM primaries): teammate values vs Gaia DR3 GSP-Spec.
EXPLORATORY, not pre-registered (DECISIONS E65).

    python -m src.feh_compare            # GSP-Spec query (ARI, cached) + offsets + group comparison + plot

Inputs
  - Teammate file (PRIVATE, CORALIE-derived values): ~/data/feh_teammate/hj_eblm_feh.csv (override $FEH_TEAMMATE_DIR).
  - The shared star list (public ids): ~/data/feh_teammate/hj_eblm_star_list.csv (E64 sha256).
Outputs
  - data/private/feh/  (gitignored): per-star merged table, aggregates json, figure. Nothing per-star is printed
    except counts and aggregates.
  - data/raw/gaia/feh/gspspec_186.*  (public Gaia values, gitignored cache).

No velocity, RV or dispersion is read: every table read is column-restricted and passes a velocity-name guard.
`vbroad_phoenix` (spectral line broadening) is read on purpose, for the broadening test.

GSP-Spec (Recio-Blanco et al. 2023, A&A 674, A29):
  - Calibration: Sect. 9.1.1 eq. 2 with Table 3 literature row (third order) on the UNCALIBRATED logg_gspspec:
    [M/H]_cal = mh_gspspec + 0.274 - 0.1373 g - 0.0050 g^2 + 0.0048 g^3.  The paper gives no logg validity range.
  - Quality: primary = "medium quality" (Sect. 9.1.1): flags 1-13 <= 1, except fluxNoise (7) <= 3 and
    KMgiantPar (13) = 0.  Variant = "best quality": flags 1-13 all 0.
"""
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table
from scipy import stats

from . import gaia_tap
from .config import ROOT

TEAM = Path(os.environ.get("FEH_TEAMMATE_DIR", Path.home() / "data" / "feh_teammate"))
FEH_FILE, LIST_FILE = TEAM / "hj_eblm_feh.csv", TEAM / "hj_eblm_star_list.csv"
LIST_SHA256 = "0d9451ac3ef71dd51c872b83e423cf392156e40caac3a26f6f084e6563ed7230"   # E64
OUT = ROOT / "data" / "private" / "feh"
FORBIDDEN = re.compile(r"radial_velocity|^rv|_rv|radv|^[UVW]$|v_tot|velocit|gamma|disp", re.I)
TEAM_COLS = ["sample", "gaia_dr3_source_id", "feh", "feh_err", "flags", "feh_phoenix", "feh_phoenix_err",
             "vbroad_phoenix", "phoenix_instruments", "status"]
P_MH = (0.274, -0.1373, -0.0050, 0.0048)          # Table 3, [M/H] literature calibration, p0..p3
SEED = 20261117
N_BOOT = 10000


def _guard(cols):
    bad = [c for c in cols if FORBIDDEN.search(c)]
    assert not bad, f"velocity-like columns requested: {bad}"
    return cols


def load_team():
    import hashlib
    assert hashlib.sha256(LIST_FILE.read_bytes()).hexdigest() == LIST_SHA256, "star list changed since E64"
    ids = pd.read_csv(LIST_FILE, usecols=_guard(["sample", "gaia_dr3_source_id"]), dtype=str)
    t = pd.read_csv(FEH_FILE, usecols=_guard(TEAM_COLS), dtype={"gaia_dr3_source_id": str, "flags": str})
    assert len(t) == 186 and t.gaia_dr3_source_id.is_unique
    assert (t[["sample", "gaia_dr3_source_id"]].values == ids.values).all(), "rows differ from the shared list"
    return t


def gspspec(dr3_ids, refresh=False):
    """GSP-Spec columns for the given DR3 ids. ESA Gaia archive (an ARI async job stayed QUEUED > 600 s on
    2026-11-17); a plain IN list, no upload. The ids are the public ones of the shared list (E64)."""
    cols = _guard(["source_id", "teff_gspspec", "logg_gspspec", "mh_gspspec", "mh_gspspec_lower",
                   "mh_gspspec_upper", "flags_gspspec"])
    ids = sorted(int(i) for i in dr3_ids)
    q = (f"SELECT {', '.join(cols)} FROM gaiadr3.astrophysical_parameters "
         f"WHERE source_id IN ({', '.join(map(str, ids))})")
    g = gaia_tap.query(q, "feh/gspspec_186", refresh=refresh).to_pandas()
    g.columns = [c.lower() for c in g.columns]
    g["gaia_dr3_source_id"] = g.pop("source_id").astype(np.int64).astype(str)
    for c in ("mh_gspspec", "logg_gspspec", "teff_gspspec", "mh_gspspec_lower", "mh_gspspec_upper"):
        g[c] = pd.to_numeric(g[c], errors="coerce")       # masked -> NaN
    g["flags_gspspec"] = g["flags_gspspec"].astype(str)
    return g


def quality(flags):
    """(medium, best) booleans from the 41-character flags_gspspec chain (Recio-Blanco+2023 Sect. 9.1.1)."""
    if not isinstance(flags, str) or len(flags) < 13 or not flags[:13].isdigit():
        return False, False
    f = [int(c) for c in flags[:13]]
    best = all(v == 0 for v in f)
    medium = all(v <= 1 for i, v in enumerate(f) if i not in (6, 12)) and f[6] <= 3 and f[12] == 0
    return medium, best


def calibrate(mh, logg):
    return mh + sum(p * logg ** i for i, p in enumerate(P_MH))


def robust(x, rng=None):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return {"n": 0}
    out = {"n": int(len(x)), "median": float(np.median(x)),
           "robust_sd": float(1.4826 * np.median(np.abs(x - np.median(x)))),
           "iqr": float(np.subtract(*np.percentile(x, [75, 25])))}
    if rng is not None and len(x) >= 5:
        b = np.median(rng.choice(x, (N_BOOT, len(x))), axis=1)
        out["median_ci68"] = [float(v) for v in np.percentile(b, [16, 84])]
    return out


def group_diff(a, b, rng):
    a, b = (np.asarray(v, float) for v in (a, b))
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    d = np.median(rng.choice(a, (N_BOOT, len(a))), 1) - np.median(rng.choice(b, (N_BOOT, len(b))), 1)
    return {"wasp_minus_eblm_median": float(np.median(a) - np.median(b)),
            "ci95": [float(v) for v in np.percentile(d, [2.5, 97.5])],
            "mannwhitney_p": float(stats.mannwhitneyu(a, b).pvalue), "ks_p": float(stats.ks_2samp(a, b).pvalue)}


def main(refresh=False):
    rng = np.random.default_rng(SEED)
    t = load_team()
    g = gspspec(t.gaia_dr3_source_id, refresh)
    d = t.merge(g, on="gaia_dr3_source_id", how="left", validate="1:1")
    d["in_ap"] = d.gaia_dr3_source_id.isin(g.gaia_dr3_source_id)
    q = d.flags_gspspec.map(quality)
    d["gsp_medium"], d["gsp_best"] = [x[0] for x in q], [x[1] for x in q]
    d["mh_cal"] = calibrate(d.mh_gspspec, d.logg_gspspec)
    has_mh = d.mh_gspspec.notna()
    d["gsp"] = d.mh_cal.where(has_mh & d.gsp_medium)               # primary GSP-Spec value
    d["gsp_best"] = d.mh_cal.where(has_mh & d.gsp_best)
    d["feh_noflag"] = d.feh.where(d["flags"].isna())
    # README_feh says feh_phoenix is fitted to ESO/SOPHIE spectra, but phoenix_instruments lists CORALIE for most
    # EBLMs: split by whether any non-CORALIE spectrum went into the PHOENIX fit.
    inst = d["phoenix_instruments"].fillna("")
    d["phx_has_coralie"] = inst.str.contains("CORALIE")
    d["phx_has_eso_sophie"] = inst.str.contains("HARPS|ESPRESSO|UVES|FEROS|SOPHIE")
    d["feh_phoenix_eso"] = d.feh_phoenix.where(d.phx_has_eso_sophie)
    d["feh_phoenix_coralie_only"] = d.feh_phoenix.where(d.phx_has_coralie & ~d.phx_has_eso_sophie)
    R = {"status": "EXPLORATORY", "n_rows": len(d), "seed": SEED, "n_boot": N_BOOT,
         "calibration": "Recio-Blanco+2023 eq.2, Table 3 literature row, uncalibrated logg_gspspec",
         "quality_primary": "medium (flags 1-13 <=1; fluxNoise<=3; KMgiantPar=0)",
         "quality_variant": "best (flags 1-13 all 0)"}

    # coverage
    cov = {}
    for s, x in [("all", d)] + list(d.groupby("sample")):
        cov[s] = {"n": len(x), "in_ap_table": int(x.in_ap.sum()),
                  "mh_gspspec": int(x.mh_gspspec.notna().sum()), "medium": int(x.gsp.notna().sum()),
                  "best": int(x.gsp_best.notna().sum()), "feh": int(x.feh.notna().sum()),
                  "feh_noflag": int(x.feh_noflag.notna().sum()), "feh_phoenix": int(x.feh_phoenix.notna().sum()),
                  "feh_phoenix_with_coralie": int((x.feh_phoenix.notna() & x.phx_has_coralie).sum()),
                  "feh_phoenix_coralie_only": int(x.feh_phoenix_coralie_only.notna().sum()),
                  "logg_range_gsp_medium": [float(x.loc[x.gsp.notna(), "logg_gspspec"].min()),
                                            float(x.loc[x.gsp.notna(), "logg_gspspec"].max())]
                  if x.gsp.notna().any() else None}
    R["coverage"] = cov

    # offsets (teammate - GSP-Spec)
    off = {}
    for ref in ("gsp", "gsp_best", "mh_gspspec_medium_uncal"):
        refv = d.mh_gspspec.where(has_mh & d.gsp_medium) if ref == "mh_gspspec_medium_uncal" else d[ref]
        for col in ("feh_noflag", "feh", "feh_phoenix", "feh_phoenix_eso", "feh_phoenix_coralie_only"):
            diff = d[col] - refv
            off[f"{col}-{ref}"] = {s: robust(diff[d["sample"].eq(s)] if s != "all" else diff, rng)
                                   for s in ("all", "WASP", "EBLM")}
    R["offsets"] = off

    # dependence of the offset on line broadening
    br = {}
    for col in ("feh_noflag", "feh", "feh_phoenix"):
        diff = d[col] - d["gsp"]
        m = diff.notna() & d.vbroad_phoenix.notna()
        res = {"n": int(m.sum())}
        if m.sum() >= 5:
            r = stats.spearmanr(d.vbroad_phoenix[m], diff[m])
            res.update(spearman_rho=float(r.statistic), spearman_p=float(r.pvalue),
                       theil_sen_slope_per_kms=float(stats.theilslopes(diff[m], d.vbroad_phoenix[m]).slope))
            for lo, hi in ((0, 5), (5, 10), (10, np.inf)):
                k = m & d.vbroad_phoenix.ge(lo) & d.vbroad_phoenix.lt(hi)
                res[f"vbroad_{lo}_{hi}"] = robust(diff[k])
        br[col] = res
    R["broadening"] = br

    # group comparison on one method
    grp = {}
    for col in ("feh_phoenix", "feh_phoenix_eso", "gsp", "gsp_best"):
        w, e = d.loc[d["sample"].eq("WASP"), col], d.loc[d["sample"].eq("EBLM"), col]
        grp[col] = {"WASP": robust(w, rng), "EBLM": robust(e, rng)}
        if w.notna().sum() >= 5 and e.notna().sum() >= 5:
            grp[col]["difference"] = group_diff(w, e, rng)
    R["groups"] = grp

    OUT.mkdir(parents=True, exist_ok=True)
    os.chmod(OUT, 0o700)
    keep = ["sample", "gaia_dr3_source_id", "feh", "feh_noflag", "feh_phoenix", "phoenix_instruments", "vbroad_phoenix", "mh_gspspec",
            "logg_gspspec", "teff_gspspec", "flags_gspspec", "mh_cal", "gsp", "gsp_best"]
    d[keep].to_csv(OUT / "feh_gspspec_merged.csv", index=False)
    (OUT / "feh_summary.json").write_text(json.dumps(R, indent=1))
    figure(d, OUT / "feh_wasp_vs_eblm.png")
    for f in OUT.iterdir():
        os.chmod(f, 0o600)
    return R


def figure(d, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, (col, lab) in zip(axes[:2], [("feh_phoenix", "teammate PHOENIX [Fe/H]"),
                                         ("gsp", "GSP-Spec [M/H] (calibrated, medium quality)")]):
        for s, c in (("WASP", "C0"), ("EBLM", "C3")):
            x = np.sort(d.loc[d["sample"].eq(s), col].dropna().values)
            if len(x):
                ax.step(x, np.arange(1, len(x) + 1) / len(x), where="post", color=c,
                        label=f"{s} (N={len(x)}, median {np.median(x):+.2f})")
        ax.set_xlabel(lab); ax.set_ylabel("cumulative fraction"); ax.legend(fontsize=8)
    ax = axes[2]
    for s, c in (("WASP", "C0"), ("EBLM", "C3")):
        k = d["sample"].eq(s)
        ax.scatter(d.vbroad_phoenix[k], (d.feh_phoenix - d.gsp)[k], s=12, color=c, label=f"{s}: PHOENIX − GSP")
        ax.scatter(d.vbroad_phoenix[k], (d.feh_noflag - d.gsp)[k], s=14, marker="x", color=c,
                   label=f"{s}: feh (no flags) − GSP")
    ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("vbroad_phoenix (km/s)"); ax.set_ylabel("Δ[Fe/H] (dex)")
    ax.legend(fontsize=7)
    fig.text(0.99, 0.01, "EXPLORATORY; PRIVATE (CORALIE-derived values); do not distribute", ha="right", fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


if __name__ == "__main__":
    import sys
    R = main(refresh="--refresh" in sys.argv)
    print(json.dumps({k: R[k] for k in ("coverage",)}, indent=1))
