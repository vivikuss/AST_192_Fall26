"""EXPLORATORY (not pre-registered; DECISIONS E69). Metallicity-matched Method A.

    python -m src.dr3x_feh_matched

Hosts: the 97 Stage 2 Method A hosts with their Stage 2 kinematic values, rebuilt through the exact Stage 2 code
path (no new host quantity); guarded against the logged E61 host sigma_W. Field: Method A field with GSP-Spec
[M/H] (E68). Each variant restricts each host's Method A pool, then runs Method A exactly as in Stage 2
(A.pools, A.draws with 10,000 draws, empty pools dropped, A.one_sided_sigma).
Outputs: data/processed/dr3_exploratory/feh_matched_results.json (aggregates); host names dropped for empty pools
go to data/private/feh_matched/ only.
"""
import json

import numpy as np
import pandas as pd

from . import config, dr3x_hosts
from . import dr3x_method_a as A
from .feh_compare import calibrate, quality, load_team
from .feh_field import PROC as FEH_PROC

E61_SIGMA_W_HOST = 18.009676119939627          # stage2_results.json, method_A main (E61)
PRIV = config.ROOT / "data" / "private" / "feh_matched"
OUT = A.OUT / "feh_matched_results.json"
NICK_SHIFT = None                               # E66 WASP median of feh(no flags) - GSP-Spec, read at run time

STATS = {"sigma_U": lambda v: np.std(v[..., 0], axis=-1), "sigma_V": lambda v: np.std(v[..., 1], axis=-1),
         "sigma_W": A.sigma_w, "sigma_tot": A.sigma_tot, "hs19_eq1": A.eq1}


def stage2_hosts():
    cfg = config.load()
    base, _, nea_rv = dr3x_hosts.bases(cfg)
    _, _, h = dr3x_hosts.method_a(base, nea_rv, exclude_companions=True)
    h = h.copy()
    h["rv"] = h["dr3"].map(nea_rv).fillna(h["radial_velocity"])
    h = A.add_uvw(h).reset_index(drop=True)
    sw = float(A.sigma_w(h[["U", "V", "W"]].to_numpy()))
    assert len(h) == 97 and abs(sw - E61_SIGMA_W_HOST) < 1e-9, f"host rows differ from Stage 2: N {len(h)}, {sw}"
    return h


def field_with_mh():
    fld, _ = A.load_field()
    g = pd.read_csv(FEH_PROC / "field_gspspec_pool.csv.gz", dtype={"flags_gspspec": str})
    med = g["flags_gspspec"].map(lambda f: quality(f)[0])
    g["mh"] = calibrate(g["mh_gspspec"], g["logg_gspspec"]).where(g["mh_gspspec"].notna() & med)
    return fld.merge(g[["source_id", "mh"]], on="source_id", how="left", validate="1:1")


def host_mh(h):
    global NICK_SHIFT
    t = load_team()
    m = pd.read_csv(config.ROOT / "data" / "private" / "feh" / "feh_gspspec_merged.csv",
                    usecols=["gaia_dr3_source_id", "gsp"], dtype={"gaia_dr3_source_id": str})
    t = t.merge(m, on="gaia_dr3_source_id", how="left", validate="1:1").set_index("gaia_dr3_source_id")
    off = json.loads((config.ROOT / "data" / "private" / "feh" / "feh_summary.json").read_text())["offsets"]
    NICK_SHIFT = off["feh_noflag-gsp"]["WASP"]["median"]
    key = h["dr3"].astype(np.int64).astype(str)
    gsp = key.map(t["gsp"])
    nick = key.map(t["feh"].where(t["flags"].isna())) - NICK_SHIFT
    return gsp.to_numpy(float), nick.to_numpy(float)


def run_variant(h, fld, label, host_mh_vals=None, window=None, require_mh=False, subset=None):
    """Method A on hosts `subset` (bool mask); pools optionally restricted by GSP-Spec [M/H] window."""
    hs = h[subset].reset_index(drop=True)
    hm = None if host_mh_vals is None else host_mh_vals[subset]
    pl = A.pools(hs, fld)
    fmh = fld["mh"].to_numpy()
    if require_mh or window is not None:
        pl = [p[np.isfinite(fmh[p])] for p in pl]
    if window is not None:
        pl = [p[np.abs(fmh[p] - hm[i]) <= window + 1e-9] for i, p in enumerate(pl)]
    empty = [hs["name"].iat[i] for i, p in enumerate(pl) if len(p) == 0]
    keep = np.array([len(p) > 0 for p in pl])
    hs, pl = hs[keep].reset_index(drop=True), [p for p in pl if len(p)]
    saved = A.STATS
    A.STATS = STATS                                  # same draws; extra per-component statistics
    try:
        mc = A.draws(pl, fld[["U", "V", "W"]].to_numpy(), label)
    finally:
        A.STATS = saved
    sizes = [len(p) for p in pl]
    res = {"label": label, "n_hosts_in": int(subset.sum()), "n_hosts": len(hs), "n_dropped_empty_pool": len(empty),
           "pool_size_median": float(np.median(sizes)), "pool_size_min": int(min(sizes)),
           "window_dex": window, "pool_requires_gsp_mh": bool(require_mh or window is not None)}
    hv = hs[["U", "V", "W"]].to_numpy()
    for s, fn in STATS.items():
        host = float(fn(hv))
        p, sig = A.one_sided_sigma(host, mc[s])
        q = np.percentile(mc[s], [2.5, 16, 50, 84, 97.5])
        res[s] = {"host": host, "field_median": float(q[2]), "one_sigma": [float(q[1]), float(q[3])],
                  "two_sigma": [float(q[0]), float(q[4])], "tail_fraction": p, "one_sided_sigma": sig}
    return res, empty


def main():
    h = stage2_hosts()
    fld = field_with_mh()
    gsp, nick = host_mh(h)
    has_gsp, has_nick = np.isfinite(gsp), np.isfinite(nick)
    R = {"header": A.HEADER, "what": "Metallicity-matched Method A (DECISIONS E69)", "n_draws": A.N_DRAWS,
         "hosts_stage2": len(h), "hosts_with_gsp_medium": int(has_gsp.sum()),
         "hosts_with_nick_feh_noflag": int(has_nick.sum()), "nick_shift_dex": NICK_SHIFT,
         "field_rows": len(fld), "field_with_gsp_medium_mh": int(fld["mh"].notna().sum())}
    dropped = {}
    runs = [("dr3x_A_feh_matched", dict(host_mh_vals=gsp, window=0.1, subset=has_gsp)),
            ("dr3x_A_feh_matched_w02", dict(host_mh_vals=gsp, window=0.2, subset=has_gsp)),
            ("dr3x_A_feh_matched_nick", dict(host_mh_vals=nick, window=0.1, subset=has_nick)),
            ("dr3x_A_feh_subset", dict(subset=has_gsp)),
            ("dr3x_A_feh_gsponly", dict(require_mh=True, subset=has_gsp)),
            ("dr3x_A_feh_subset_nick", dict(subset=has_nick))]
    for label, kw in runs:
        R[label], dropped[label] = run_variant(h, fld, label, **kw)
        r = R[label]
        print(f"{label}: N {r['n_hosts']}/{r['n_hosts_in']}, pool median {r['pool_size_median']:.0f}, "
              f"sigma_W {r['sigma_W']['host']:.2f} vs {r['sigma_W']['field_median']:.2f} "
              f"({r['sigma_W']['one_sided_sigma']:+.2f} sigma)", flush=True)
    A.write_json(OUT, R)
    PRIV.mkdir(parents=True, exist_ok=True)
    (PRIV / "dropped_empty_pool.json").write_text(json.dumps(dropped, indent=1))
    return R


if __name__ == "__main__":
    main()
