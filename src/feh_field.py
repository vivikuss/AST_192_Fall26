"""EXPLORATORY (not pre-registered; DECISIONS E67). Field kinematics vs GSP-Spec [M/H], and the dispersion change
that the WASP hosts' (and, caveated, the EBLMs') metallicity alone would predict.

    python -m src.feh_field field     # host ranges, pool union, field selection (no GSP-Spec yet)
    python -m src.feh_field gsp       # GSP-Spec for the selected field (ESA archive, upload, cached)
    python -m src.feh_field run       # [M/H]-binned dispersions, reweighting predictions, figures

Hard rule: NO velocity, RV or dispersion of any host or EBLM is read or computed. Hosts contribute only their
Method A membership (names from data/processed/dr3_exploratory/host_steps_wasp107.json, no RV used) and their
d, G, (BP-RP)_0 (from the cached DR3 rows, reading only non-velocity columns), plus [M/H]. Velocities used here
are those of FIELD stars (Kayla's field file, NEA hosts removed at load as in Method A).

Field sets:
  - "pool": union of the Method A matching pools of the 97 Method A WASP hosts (<= 10 deg, -100..+50 pc,
    -1..+0.1 mag in G, |d colour| <= 0.2).  Primary.
  - "box": field stars inside the hosts' overall range of (BP-RP)_0, G and d (min..max).  Secondary.
[M/H]: GSP-Spec, calibration and medium-quality cut exactly as E66 (src.feh_compare.calibrate/quality).
Private: Nick's per-star values are only read from ~/data/feh_teammate; outputs that hold them stay in
data/private/feh_field/ (gitignored).
"""
import json
import re
import sys

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

from . import dr3x_method_a as ma
from . import gaia_tap
from .config import ROOT, RAW
from .feh_compare import calibrate, quality, FEH_FILE, load_team

PROC = ROOT / "data" / "processed" / "feh_field"
PRIV = ROOT / "data" / "private" / "feh_field"
HOST_STEPS = ROOT / "data" / "processed" / "dr3_exploratory" / "host_steps_wasp107.json"
HOST_FITS = RAW / "ari" / "dr3x_hosts_wasp107.fits"
HOST_SAFE_COLS = ["source_id", "ra", "dec", "parallax", "phot_g_mean_mag", "bp_rp", "ag_gspphot",
                  "ebpminrp_gspphot", "parallax_over_error", "visibility_periods_used", "astrometric_chi2_al",
                  "astrometric_n_good_obs_al"]
FORBIDDEN = re.compile(r"radial_velocity|^rv|_rv|radv|^[UVW]$|v_tot|velocit|gamma|disp", re.I)
BIN = 0.1
SEED = 20261118
N_BOOT = 1000
WASP_FEH_PHOENIX_OFFSET = None      # filled from data/private/feh/feh_summary.json (E66, WASP feh_phoenix - GSP)


def _guard(cols):
    bad = [c for c in cols if FORBIDDEN.search(c)]
    assert not bad, bad
    return cols


# ----------------------------------------------------------------------------- hosts (no kinematics)
def method_a_hosts():
    """The 97 Method A WASP hosts: name, dr3, ra, dec, d, G, color0. No velocity column is touched."""
    steps = json.loads(HOST_STEPS.read_text())["method_A_main_with_companion_exclusions"]
    removed = {n for s in steps for n in s["removed"]}
    n_final = steps[-1]["n_after"]
    w = ma.wasp_hosts()
    w = w[~w["name"].isin(removed)].reset_index(drop=True)
    assert len(w) == n_final == 97, (len(w), n_final)
    with fits.open(HOST_FITS, memmap=True) as hd:
        data = hd[1].data
        cols = {c.lower(): c for c in data.columns.names}
        h = pd.DataFrame({c: np.asarray(data[cols[c]]).astype(float if c != "source_id" else np.int64)
                          for c in _guard(HOST_SAFE_COLS)})
    h = h.rename(columns={"source_id": "dr3"})
    w = w.merge(h, on="dr3", how="left", validate="1:1")
    w = ma.photometry(w)
    # consistency with the logged steps (astrometry, MS + colour, d < 600)
    ok = ma.astrometry_ok(w) & (w["zams_offset"] > -ma.MS_MAX_ABOVE) & w["color0"].between(*ma.COLOUR_WINDOW) \
        & (w["d"] < ma.D_MAX)
    assert ok.all(), f"{(~ok).sum()} host(s) disagree with the logged Method A steps"
    return w[["name", "dr3", "ra", "dec", "d", "phot_g_mean_mag", "color0"]]


# ----------------------------------------------------------------------------- field selection
def select_field():
    hosts = method_a_hosts()
    field, n_removed = ma.load_field()
    pools = ma.pools(hosts, field)
    in_pool = np.zeros(len(field), bool)
    for p in pools:
        in_pool[p] = True
    rng = {c: [float(hosts[c].min()), float(hosts[c].max())] for c in ("color0", "phot_g_mean_mag", "d")}
    in_box = np.ones(len(field), bool)
    for c, (lo, hi) in rng.items():
        in_box &= field[c].between(lo, hi).to_numpy()
    sel = field[in_pool | in_box].copy()
    sel["in_pool"], sel["in_box"] = in_pool[in_pool | in_box], in_box[in_pool | in_box]
    PROC.mkdir(parents=True, exist_ok=True)
    sel.to_csv(PROC / "field_selected.csv.gz", index=False)
    meta = {"n_field_loaded": len(field), "n_nea_hosts_removed": n_removed, "host_ranges": rng,
            "n_pool_union": int(in_pool.sum()), "n_box": int(in_box.sum()), "n_selected": len(sel),
            "pool_sizes": {"min": int(min(map(len, pools))), "median": float(np.median(list(map(len, pools)))),
                           "max": int(max(map(len, pools)))}}
    (PROC / "field_selection.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))


# ----------------------------------------------------------------------------- GSP-Spec for the field
GSP_ADQL = ("SELECT g.source_id, a.teff_gspspec, a.logg_gspspec, a.mh_gspspec, a.flags_gspspec "
            "FROM gaiadr3.gaia_source AS g JOIN gaiadr3.astrophysical_parameters AS a ON a.source_id = g.source_id "
            "WHERE {range} AND a.mh_gspspec IS NOT NULL AND g.parallax > 1.5 "
            "AND g.phot_g_mean_mag BETWEEN 8.0 AND 13.2")


def gsp_field():
    """Every GSP-Spec star with parallax > 1.5 mas and 8.0 <= G <= 13.2 (covers the pool union and the box:
    d <= 667 pc, G 8.17-13.10), one ARI async job per HEALPix level-1 range. ESA ID uploads stalled and ESA chunked jobs got HTTP 429 (E67)."""
    _guard(re.findall(r"\b\w+_gspspec\b", GSP_ADQL))
    from . import ari_tap
    tbl, failed = ari_tap.query_ranges(GSP_ADQL, "feh/field_gspspec", level=1)
    g = tbl.to_pandas()
    g.columns = [c.lower() for c in g.columns]
    PROC.mkdir(parents=True, exist_ok=True)
    g.to_csv(PROC / "field_gspspec_all.csv.gz", index=False)
    (PROC / "field_gspspec_query.json").write_text(json.dumps({"adql": GSP_ADQL, "rows": len(g),
                                                                "failed_pixels": failed}, indent=1))
    print(len(g), "rows; failed pixels:", len(failed))


# ----------------------------------------------------------------------------- analysis
MH_EDGES = np.round(np.arange(-1.0, 0.8001, BIN), 2)          # 0.1 dex bins, -1.0 .. +0.8
MIN_BIN = 30
COMPONENTS = ("sigma_U", "sigma_V", "sigma_W", "sigma_tot")


def disp(uvw, w=None):
    """Weighted (or plain) std per component (ddof 0) and sigma_tot = sqrt(sum var). FIELD stars only."""
    w = np.ones(len(uvw)) if w is None else w
    m = np.average(uvw, axis=0, weights=w)
    var = np.average((uvw - m) ** 2, axis=0, weights=w)
    return np.r_[np.sqrt(var), np.sqrt(var.sum())]


def load_field_mh():
    sel = pd.read_csv(PROC / "field_selected.csv.gz")
    g = pd.read_csv(PROC / "field_gspspec_all.csv.gz", dtype={"flags_gspspec": str})
    f = sel.merge(g, on="source_id", how="left", validate="1:1")
    q = f["flags_gspspec"].map(quality)
    f["medium"] = [x[0] for x in q]
    f["mh"] = calibrate(f["mh_gspspec"], f["logg_gspspec"]).where(f["mh_gspspec"].notna() & f["medium"])
    return f


def host_mh(hosts):
    """[M/H] samples for the reweighting targets. Private inputs (Nick); only aggregates leave this function."""
    t = load_team()
    m = pd.read_csv(ROOT / "data" / "private" / "feh" / "feh_gspspec_merged.csv", dtype={"gaia_dr3_source_id": str},
                    usecols=["gaia_dr3_source_id", "gsp"])
    t = t.merge(m, on="gaia_dr3_source_id", how="left", validate="1:1")
    off = json.loads((ROOT / "data" / "private" / "feh" / "feh_summary.json").read_text())["offsets"]
    o_feh, o_phx = off["feh-gsp"]["WASP"]["median"], off["feh_phoenix-gsp"]["WASP"]["median"]
    ids = set(hosts["dr3"].astype(str))
    wa = t[t["gaia_dr3_source_id"].isin(ids)]
    eb = t[t["sample"].eq("EBLM")]
    targets = {
        "WASP_gsp": wa["gsp"].dropna().to_numpy(),
        "WASP_nick_feh_on_gsp_scale": (wa["feh"] - o_feh).dropna().to_numpy(),
        "WASP_nick_feh_noflag_on_gsp_scale": (wa["feh"].where(wa["flags"].isna()) - o_feh).dropna().to_numpy(),
        "EBLM_feh_phoenix_shifted_ASSUMED": (eb["feh_phoenix"] - o_phx).dropna().to_numpy(),
    }
    meta = {"n_method_a_hosts": len(hosts), "offsets_used": {"feh-gsp WASP": o_feh, "feh_phoenix-gsp WASP": o_phx},
            "n": {k: int(len(v)) for k, v in targets.items()}}
    return targets, meta


def weights(field_mh, target, edges):
    fb = np.digitize(field_mh, edges) - 1
    tb = np.digitize(target, edges) - 1
    nb = len(edges) - 1
    fh = np.bincount(fb[(fb >= 0) & (fb < nb)], minlength=nb).astype(float)
    th = np.bincount(tb[(tb >= 0) & (tb < nb)], minlength=nb).astype(float)
    lost = th[fh == 0].sum() + ((tb < 0) | (tb >= nb)).sum()           # target mass with no field stars
    ratio = np.divide(th / max(th.sum(), 1), fh / max(fh.sum(), 1), out=np.zeros(nb), where=fh > 0)
    w = np.where((fb >= 0) & (fb < nb), ratio[np.clip(fb, 0, nb - 1)], 0.0)
    return w, float(lost / max(len(target), 1))


def predict(f, target, rng, edges=MH_EDGES):
    """Delta sigma (reweighted - unweighted) for each component; bootstrap over field stars AND target stars."""
    uvw, mh = f[["U", "V", "W"]].to_numpy(), f["mh"].to_numpy()
    base = disp(uvw)
    w, lost = weights(mh, target, edges)
    point = disp(uvw, w) - base
    bs = np.empty((N_BOOT, 4))
    for b in range(N_BOOT):
        i = rng.integers(0, len(f), len(f))
        tt = target[rng.integers(0, len(target), len(target))]
        wb, _ = weights(mh[i], tt, edges)
        bs[b] = disp(uvw[i], wb) - disp(uvw[i])
    return {"unweighted": dict(zip(COMPONENTS, base.round(3).tolist())),
            "delta": dict(zip(COMPONENTS, point.round(3).tolist())),
            "delta_ci68": {c: np.percentile(bs[:, k], [16, 84]).round(3).tolist() for k, c in enumerate(COMPONENTS)},
            "delta_ci95": {c: np.percentile(bs[:, k], [2.5, 97.5]).round(3).tolist() for k, c in enumerate(COMPONENTS)},
            "relative_delta_pct": dict(zip(COMPONENTS, (100 * point / base).round(2).tolist())),
            "target_mass_without_field_stars": round(lost, 3), "n_field": len(f), "n_target": int(len(target)),
            "field_median_mh": round(float(np.median(mh)), 3), "target_median_mh": round(float(np.median(target)), 3),
            "eff_n_field": round(float(w.sum() ** 2 / (w ** 2).sum()), 1)}


def binned(f, rng):
    uvw, mh = f[["U", "V", "W"]].to_numpy(), f["mh"].to_numpy()
    rows = []
    for lo, hi in zip(MH_EDGES[:-1], MH_EDGES[1:]):
        k = np.flatnonzero((mh >= lo) & (mh < hi))
        r = {"mh_lo": float(lo), "mh_hi": float(hi), "n": int(len(k))}
        if len(k) >= MIN_BIN:
            v = uvw[k]
            r.update(dict(zip(COMPONENTS, disp(v).round(3).tolist())))
            bs = np.array([disp(v[rng.integers(0, len(k), len(k))]) for _ in range(N_BOOT)])
            r.update({c + "_err": round(float(bs[:, j].std()), 3) for j, c in enumerate(COMPONENTS)})
        rows.append(r)
    return pd.DataFrame(rows)


def run():
    rng = np.random.default_rng(SEED)
    f = load_field_mh()
    hosts = method_a_hosts()
    out = {"status": "EXPLORATORY", "seed": SEED, "n_boot": N_BOOT, "bin_dex": BIN, "min_bin_n": MIN_BIN,
           "statistic": "np.std per component (ddof 0); sigma_tot = sqrt(sum var); FIELD stars only"}
    cov = {}
    for name, m in (("pool", f["in_pool"]), ("box", f["in_box"])):
        x = f[m]
        gb = pd.cut(x["phot_g_mean_mag"], [8, 10, 11, 12, 12.5, 13.2])
        cov[name] = {"n": len(x), "any_mh_gspspec": int(x["mh_gspspec"].notna().sum()),
                     "medium": int(x["mh"].notna().sum()),
                     "medium_frac_by_G": {str(k): round(float(v), 3) for k, v in
                                          x["mh"].notna().groupby(gb, observed=False).mean().items()},
                     "medium_frac_by_color0": {str(k): round(float(v), 3) for k, v in x["mh"].notna().groupby(
                         pd.cut(x["color0"], [0.5, 0.7, 0.9, 1.1, 1.4]), observed=False).mean().items()}}
    out["coverage"] = cov
    targets, tmeta = host_mh(hosts)
    out["targets"] = tmeta
    out["target_medians"] = {k: round(float(np.median(v)), 3) for k, v in targets.items()}
    PROC.mkdir(parents=True, exist_ok=True)
    for name in ("pool", "box"):
        x = f[f["in_" + name] & f["mh"].notna()].reset_index(drop=True)
        b = binned(x, rng)
        b.to_csv(PROC / f"binned_dispersions_{name}.csv", index=False)
        out[f"predictions_{name}"] = {k: predict(x, v, rng) for k, v in targets.items()}
        out[f"predictions_{name}_bin0.2"] = {"WASP_gsp": predict(x, targets["WASP_gsp"], rng,
                                                                 np.round(np.arange(-1.0, 0.8001, 0.2), 2))}
        if name == "pool":
            fig_binned(b, out, PROC.parent.parent.parent / "figures" / "feh_field")
            fig_uv(x, PROC.parent.parent.parent / "figures" / "feh_field")
    (PROC / "feh_field_summary.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("coverage", "targets", "target_medians")}, indent=1))


def fig_binned(b, out, d):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d.mkdir(parents=True, exist_ok=True)
    ok = b.dropna(subset=["sigma_U"])
    x = (ok["mh_lo"] + ok["mh_hi"]) / 2
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
    for ax, c in zip(axes.flat, COMPONENTS):
        ax.errorbar(x, ok[c], yerr=ok[c + "_err"], fmt="o", ms=4, color="k")
        for k, col in (("WASP_gsp", "C0"), ("EBLM_feh_phoenix_shifted_ASSUMED", "C3")):
            ax.axvline(out["target_medians"][k], color=col, ls="--", lw=1,
                       label={"WASP_gsp": "WASP hosts median (GSP-Spec)",
                              "EBLM_feh_phoenix_shifted_ASSUMED": "EBLM median (PHOENIX, shifted; assumed)"}[k])
        ax.set_ylabel(c.replace("sigma_", "σ_") + " (km/s)")
        ax.set_title(f"Field (Method A pool union), {c}", fontsize=9)
    for ax in axes[1]:
        ax.set_xlabel("GSP-Spec [M/H] (calibrated, medium quality)")
    axes[0, 0].legend(fontsize=7)
    fig.text(0.99, 0.005, f"EXPLORATORY. Field stars only; bins of {BIN} dex with N >= {MIN_BIN}; "
             f"bootstrap errors ({N_BOOT})", ha="right", fontsize=7)
    fig.tight_layout()
    fig.savefig(d / "field_sigma_vs_mh.png", dpi=150)
    plt.close(fig)


def fig_uv(x, d):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    bins = [(-np.inf, -0.3), (-0.3, -0.1), (-0.1, 0.1), (0.1, 0.3), (0.3, np.inf)]
    fig, axes = plt.subplots(1, len(bins), figsize=(18, 4), sharex=True, sharey=True)
    for ax, (lo, hi) in zip(axes, bins):
        k = (x["mh"] >= lo) & (x["mh"] < hi)
        h = ax.hist2d(x.loc[k, "U"], x.loc[k, "V"], bins=[np.linspace(-150, 150, 61), np.linspace(-120, 60, 37)],
                      density=True, norm=LogNorm(vmin=1e-6, vmax=3e-3), cmap="viridis")
        lab = f"[M/H] < {hi}" if lo == -np.inf else (f"[M/H] >= {lo}" if hi == np.inf else f"{lo} <= [M/H] < {hi}")
        ax.set_title(f"{lab} (N={int(k.sum())})", fontsize=9)
        ax.set_xlabel("U (km/s)")
    axes[0].set_ylabel("V (km/s)")
    fig.colorbar(h[3], ax=axes, label="density (per km²/s²)", shrink=0.9)
    fig.text(0.99, 0.005, "EXPLORATORY. Field stars only (Method A pool union, GSP-Spec medium quality); "
             "UVW as in Kayla's field file", ha="right", fontsize=7)
    fig.savefig(d / "field_uv_by_mh.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    {"field": select_field, "gsp": gsp_field,
     "run": lambda: run()}[sys.argv[1]]()
