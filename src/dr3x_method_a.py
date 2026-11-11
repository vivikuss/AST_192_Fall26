"""EXPLORATORY (not pre-registered; DECISIONS E52, E53). Method A: Hu et al. (2026) as implemented by Kayla.

Our own implementation, written from reading Kayla's script as text. Her script was not imported, run or
copied in:
    ~/data/hu2026_dr3_Kayla/code/hu2026_replication.py
    sha256 6f7bd10a539e13de70b3b8a73f21720aaafbb2ee12c3550f7b58d81057057887 (DECISIONS E50)
Settings follow that script as received. The field is her `field_main_sequence_uvw.csv`, whose build settings
were verified from her source text (cache key, E51).

    python -m src.dr3x_method_a hosts         # host step counts (no kinematics)
    python -m src.dr3x_method_a validate      # Hu et al. HS19-all and TESS d<600 (needs the published tables)
    python -m src.dr3x_method_a wasp --stage2 # WASP sigma_W etc. Refuses without --stage2 (Stage 2 only)

Steps (her logic):
  - Hosts: Gaia DR3 rows (her host columns + ruwe); App. A cuts 1-3 (parallax_over_error > 10,
    visibility_periods_used > 8, UWE < 1.2 max(1, exp(-0.2 (G - 19.5)))); RV = NEA st_radv, else Gaia DR3 RV;
    d = 1000/parallax; extinction = GSP-Phot (A_G, E(BP-RP)), else Bayestar19 x 0.884 with CV18 coefficients,
    else zero ONLY where the map has no value by design (dec < -30 deg, or closer than its first distance bin);
    main sequence = < 1 mag above the HS19 Pleiades ZAMS polynomial and 0.5 <= (BP-RP)_0 <= 1.5; d < 600 pc.
  - A missing Bayestar19 file is a hard stop (RuntimeError). It is never replaced by zero.
  - Field: her file; every NEA host with a DR3 id is removed at load.
  - Matching: per-host pool (separation <= 10 deg, -100 <= d_FS - d_HJ <= +50 pc, -1 <= G_FS - G_HJ <= +0.1,
    |Delta (BP-RP)_0| <= 0.2, all inclusive). 10,000 draws of one pool star per host, with replacement. A host
    with an empty pool is dropped (named).
  - Statistics on the same draws: sigma_W = std(W) (ddof 0; primary), sigma_tot, HS19 eq. 1.
  - Significance: tail fraction on the hosts' side of the control median -> one-sided Gaussian sigma (floor
    1/N); 1 sigma = [16, 84], 2 sigma = [2.5, 97.5] percentiles of the controls.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import norm

from . import seeds
from .config import ROOT

HEADER = "EXPLORATORY — not pre-registered"
KAYLA = Path.home() / "data" / "hu2026_dr3_Kayla"
FIELD_CSV = KAYLA / "field_main_sequence_uvw.csv"
NEA_CSV = KAYLA / "nea_pscomppars.csv"
BAYESTAR_H5 = Path.home() / "data" / "bayestar19" / "bayestar2019.h5"   # location to be agreed (E53)
CROSSMATCH = ROOT / "data" / "processed" / "dr3_crossmatch" / "wasp107_dr2_dr3_summary.csv"
OUT = ROOT / "data" / "processed" / "dr3_exploratory"

# values as in her script (E50); names are ours
R_G, R_BP, R_RP = 2.740, 3.374, 2.035
BAYESTAR_TO_EBV = 0.884
ZAMS = [1.07857, 6.23258, -10.85944, 21.65561, -20.27879, 9.78665, -2.52543, 0.33004, -0.0170556]
MS_MAX_ABOVE = 1.0
COLOUR_WINDOW = (0.5, 1.5)
D_MAX = 600.0
POOL_HJ = {"sep_deg": 10.0, "dd": (-100.0, 50.0), "dg": (-1.0, 0.1), "dcol": 0.2}
N_DRAWS = 10_000
HOST_COLUMNS = ["source_id", "ra", "dec", "parallax", "pmra", "pmdec", "radial_velocity", "phot_g_mean_mag",
                "bp_rp", "ag_gspphot", "ebpminrp_gspphot", "parallax_over_error", "visibility_periods_used",
                "astrometric_chi2_al", "astrometric_n_good_obs_al", "ruwe"]
FIELD_COLUMNS = ["source_id", "ra", "dec", "d", "phot_g_mean_mag", "color0", "ext_source", "U", "V", "W"]


# ----------------------------------------------------------------------------- output
def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"header": HEADER, **obj}, indent=1, default=float))


def write_csv(path, df):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(f"# {HEADER}\n")
        df.to_csv(fh, index=False)


# ----------------------------------------------------------------------------- HEALPix (NESTED)
def _interleave(v):
    """Bits b0 b1 b2 ... -> b0 0 b1 0 b2 ... (v < 2**16)."""
    v = np.asarray(v, dtype=np.int64)
    out = np.zeros_like(v)
    for bit in range(16):
        out |= ((v >> bit) & 1) << (2 * bit)
    return out


def ang2pix_nest(nside, lon_deg, lat_deg):
    """NESTED HEALPix pixel of (lon, lat); the ang2pix_nest algorithm of Gorski et al. (2005) / chealpix."""
    lon_deg, lat_deg = np.atleast_1d(lon_deg).astype(float), np.atleast_1d(lat_deg).astype(float)
    z = np.sin(np.radians(lat_deg))
    za = np.abs(z)
    tt = np.mod(np.radians(lon_deg), 2 * np.pi) / (np.pi / 2)
    face, ix, iy = (np.zeros(len(z), dtype=np.int64) for _ in range(3))
    eq = za <= 2.0 / 3.0
    if eq.any():
        a = nside * (0.5 + tt[eq])
        bz = nside * z[eq] * 0.75
        jp, jm = (a - bz).astype(np.int64), (a + bz).astype(np.int64)
        fp, fm = jp // nside, jm // nside
        face[eq] = np.where(fp == fm, fp | 4, np.where(fp < fm, fp, fm + 8))
        ix[eq] = jm & (nside - 1)
        iy[eq] = nside - (jp & (nside - 1)) - 1
    cap = ~eq
    if cap.any():
        ntt = np.minimum(tt[cap].astype(np.int64), 3)
        tp = tt[cap] - ntt
        r = nside * np.sqrt(3 * (1 - za[cap]))
        jp = np.minimum((tp * r).astype(np.int64), nside - 1)
        jm = np.minimum(((1 - tp) * r).astype(np.int64), nside - 1)
        north = z[cap] >= 0
        face[cap] = np.where(north, ntt, ntt + 8)
        ix[cap] = np.where(north, nside - jm - 1, jp)
        iy[cap] = np.where(north, nside - jp - 1, jm)
    return face * nside * nside + _interleave(ix) + (_interleave(iy) << 1)


# ----------------------------------------------------------------------------- extinction
def galactic(ra, dec):
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    g = SkyCoord(ra=np.asarray(ra) * u.deg, dec=np.asarray(dec) * u.deg).galactic
    return g.l.deg, g.b.deg


def bayestar_ebv(ra, dec, d_pc, path=BAYESTAR_H5):
    """E(B-V) from Bayestar19 (median of the samples in each distance bin, x 0.884), interpolated linearly in
    distance modulus and held at the last bin beyond it. Returns (ebv, why) with why in
    {'bayestar19', 'zero: outside map', 'zero: closer than first bin'}. Missing file -> RuntimeError."""
    import h5py
    path = Path(path)
    if not path.is_file():
        raise RuntimeError(f"Bayestar19 file not found at {path}: Method A stops here (no zero fallback).")
    n = len(ra)
    ebv, why = np.zeros(n), np.full(n, "zero: outside map", dtype=object)
    if n == 0:
        return ebv, why
    l, b = galactic(ra, dec)
    with h5py.File(path, "r") as f:
        info = f["pixel_info"][:]
        dm_grid = np.asarray(f["pixel_info"].attrs["DM_bin_edges"], dtype=float)
        row = np.full(n, -1, dtype=np.int64)
        for nside in np.unique(info["nside"]):
            sel = np.flatnonzero(info["nside"] == nside)
            hp = info["healpix_index"][sel].astype(np.int64)
            o = np.argsort(hp)
            hp, sel = hp[o], sel[o]
            pix = ang2pix_nest(int(nside), l, b)
            k = np.clip(np.searchsorted(hp, pix), 0, len(hp) - 1)
            hit = hp[k] == pix
            row[hit] = sel[k[hit]]
        found = np.flatnonzero(row >= 0)
        if len(found) == 0:
            return ebv, why
        rows = np.unique(row[found])
        prof = np.median(f["samples"][rows.tolist()], axis=1)       # (n_rows, n_dm)
    p = prof[np.searchsorted(rows, row[found])]
    dm = 5 * np.log10(np.asarray(d_pc, dtype=float)[found]) - 5
    j = np.clip(np.searchsorted(dm_grid, dm), 1, len(dm_grid) - 1)
    w = np.clip((dm - dm_grid[j - 1]) / (dm_grid[j] - dm_grid[j - 1]), 0, 1)
    k = np.arange(len(found))
    val = ((1 - w) * p[k, j - 1] + w * p[k, j]) * BAYESTAR_TO_EBV
    close = dm < dm_grid[0]
    ebv[found] = np.where(close, 0.0, val)
    why[found] = np.where(close, "zero: closer than first bin", "bayestar19")
    return ebv, why


def photometry(h, bayestar_path=BAYESTAR_H5):
    """d, A_G, E(BP-RP), (BP-RP)_0, M_G0, ZAMS offset, ext_source. Requires parallax > 0."""
    h = h.copy()
    h["d"] = 1000.0 / h["parallax"]
    gsp = (h["ag_gspphot"].notna() & h["ebpminrp_gspphot"].notna()).to_numpy()
    ebv, why = np.zeros(len(h)), np.full(len(h), "gspphot", dtype=object)
    if (~gsp).any():
        e, w = bayestar_ebv(h["ra"].to_numpy()[~gsp], h["dec"].to_numpy()[~gsp], h["d"].to_numpy()[~gsp],
                            bayestar_path)
        ebv[~gsp], why[~gsp] = e, w
    h["A_G"] = np.where(gsp, h["ag_gspphot"], R_G * ebv)
    h["E_BPRP"] = np.where(gsp, h["ebpminrp_gspphot"], (R_BP - R_RP) * ebv)
    h["ext_source"] = why
    h["color0"] = h["bp_rp"] - h["E_BPRP"]
    h["mg0"] = h["phot_g_mean_mag"] - 5 * np.log10(h["d"] / 10.0) - h["A_G"]
    h["zams_offset"] = h["mg0"] - np.polynomial.polynomial.polyval(h["color0"], ZAMS)
    return h


# ----------------------------------------------------------------------------- hosts
def astrometry_ok(h):
    uwe = np.sqrt(h["astrometric_chi2_al"] / (h["astrometric_n_good_obs_al"] - 5))
    lim = 1.2 * np.maximum(1.0, np.exp(-0.2 * (h["phot_g_mean_mag"] - 19.5)))
    return (h["parallax_over_error"] > 10) & (h["visibility_periods_used"] > 8) & (uwe < lim)


def nea_radv(dr3_ids, nea_csv=NEA_CSV):
    """NEA st_radv per DR3 id: first non-null over the host's HJ rows sorted by (disc_year, pl_name), as her
    host table does (pandas 'first' skips nulls); hosts without an HJ row use all their rows."""
    n = pd.read_csv(nea_csv, usecols=["pl_name", "hostname", "pl_orbper", "pl_bmassj", "disc_year", "st_radv",
                                      "gaia_dr3_id"], dtype={"gaia_dr3_id": str})
    n["dr3"] = n["gaia_dr3_id"].str.extract(r"(\d{6,})$")[0].astype("Int64")
    n = n[n["dr3"].isin(set(dr3_ids))]
    hj = (n["pl_orbper"] < 10.0) & (n["pl_bmassj"] > 0.1)
    n = pd.concat([n[hj], n[~hj & ~n["dr3"].isin(n.loc[hj, "dr3"])]])
    n = n.sort_values(["disc_year", "pl_name"])
    return n.groupby("dr3")["st_radv"].first()


def host_steps(h, nea_rv, bayestar_path=BAYESTAR_H5):
    """Apply her host steps in order; returns (survivors, steps) with the names removed at each step.
    Only the presence of an RV is used here, never its value."""
    steps = []

    def cut(name, keep):
        nonlocal h
        keep = np.asarray(keep, dtype=bool)
        steps.append({"step": name, "removed": h.loc[~keep, "name"].tolist(), "n_after": int(keep.sum())})
        h = h[keep].reset_index(drop=True)

    steps.append({"step": "start", "removed": [], "n_after": len(h)})
    cut("no Gaia DR3 row", h["parallax"].notna())
    cut("astrometric cuts 1-3", astrometry_ok(h))
    has_rv = h["dr3"].map(nea_rv).notna() | h["radial_velocity"].notna()
    h["rv_source"] = np.where(h["dr3"].map(nea_rv).notna(), "NEA st_radv",
                              np.where(h["radial_velocity"].notna(), "Gaia DR3", "none"))
    cut("no RV (NEA or Gaia DR3)", has_rv)
    h = photometry(h, bayestar_path)
    cut("main sequence + 0.5<=(BP-RP)_0<=1.5",
        (h["zams_offset"] > -MS_MAX_ABOVE) & h["color0"].between(*COLOUR_WINDOW))
    cut("d < 600 pc", h["d"] < D_MAX)
    return h, steps


# ----------------------------------------------------------------------------- field, pools, draws
def load_field(field_csv=FIELD_CSV, nea_csv=NEA_CSV):
    f = pd.read_csv(field_csv, usecols=FIELD_COLUMNS, dtype={"source_id": np.int64})
    ids = pd.read_csv(nea_csv, usecols=["gaia_dr3_id"], dtype=str)["gaia_dr3_id"] \
        .str.extract(r"(\d{6,})$")[0].dropna().astype(np.int64)
    keep = ~f["source_id"].isin(set(ids))
    return f[keep].reset_index(drop=True), int((~keep).sum())


def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def pools(hosts, field, crit=POOL_HJ):
    tree = cKDTree(unit(field["ra"].to_numpy(), field["dec"].to_numpy()))
    near = tree.query_ball_point(unit(hosts["ra"].to_numpy(), hosts["dec"].to_numpy()),
                                 2 * np.sin(np.radians(crit["sep_deg"]) / 2))
    fd, fg, fc = (field[c].to_numpy() for c in ("d", "phot_g_mean_mag", "color0"))
    out = []
    for i, cand in enumerate(near):
        c = np.asarray(cand, dtype=np.int64)
        dd, dg = fd[c] - hosts["d"].iat[i], fg[c] - hosts["phot_g_mean_mag"].iat[i]
        ok = ((dd >= crit["dd"][0]) & (dd <= crit["dd"][1]) & (dg >= crit["dg"][0]) & (dg <= crit["dg"][1])
              & (np.abs(fc[c] - hosts["color0"].iat[i]) <= crit["dcol"]))
        out.append(c[ok])
    return out


def sigma_w(v):
    return np.std(v[..., 2], axis=-1)


def sigma_tot(v):
    return np.sqrt(np.var(v, axis=-2).sum(axis=-1))


def eq1(v):
    return np.mean(np.linalg.norm(v - v.mean(axis=-2, keepdims=True), axis=-1), axis=-1)


STATS = {"sigma_W": sigma_w, "sigma_tot": sigma_tot, "hs19_eq1": eq1}


def draws(pool_list, field_uvw, label, n=N_DRAWS, block=500):
    """n control samples (one pool star per host, with replacement); all statistics on the same draws."""
    rng = np.random.default_rng(np.random.SeedSequence(list(seeds.stream_entropy(label))))
    seeds.label_code(label)                                   # registered-label check
    assert label in seeds.REGISTRY, label
    cnt = np.array([len(p) for p in pool_list])
    flat, off = np.concatenate(pool_list), np.concatenate([[0], np.cumsum(cnt)[:-1]])
    out = {k: np.empty(n) for k in STATS}
    for i0 in range(0, n, block):
        k = min(block, n - i0)
        v = field_uvw[flat[off + (rng.random((k, len(pool_list))) * cnt).astype(np.int64)]]
        for s, fn in STATS.items():
            out[s][i0:i0 + k] = fn(v)
    return out


def one_sided_sigma(host, mc):
    below = host < np.median(mc)
    p = np.mean(mc <= host) if below else np.mean(mc >= host)
    s = norm.isf(max(p, 1.0 / len(mc)))
    return float(p), float(-s if below else s)


def compare(hosts, field, label, crit=POOL_HJ):
    pl = pools(hosts, field, crit)
    empty = [hosts["name"].iat[i] for i, p in enumerate(pl) if len(p) == 0]
    keep = np.array([len(p) > 0 for p in pl])
    hosts, pl = hosts[keep].reset_index(drop=True), [p for p in pl if len(p)]
    mc = draws(pl, field[["U", "V", "W"]].to_numpy(), label)
    res = {"n_hosts": len(hosts), "dropped_empty_pool": empty,
           "pool_size_median": float(np.median([len(p) for p in pl])), "pool_size_min": int(min(map(len, pl)))}
    hv = hosts[["U", "V", "W"]].to_numpy()
    for s, fn in STATS.items():
        p, sig = one_sided_sigma(float(fn(hv)), mc[s])
        q = np.percentile(mc[s], [2.5, 16, 50, 84, 97.5])
        res[s] = {"host": float(fn(hv)), "field_median": float(q[2]), "one_sigma": [q[1], q[3]],
                  "two_sigma": [q[0], q[4]], "tail_fraction": p, "one_sided_sigma": sig}
    return res


# ----------------------------------------------------------------------------- Gaia DR3 host rows
def query_hosts(dr3_ids, tag):
    from . import gaia_tap
    ids = ", ".join(str(int(i)) for i in sorted(set(dr3_ids)))
    q = f"SELECT {', '.join('g.' + c for c in HOST_COLUMNS)} FROM gaiadr3.gaia_source AS g WHERE g.source_id IN ({ids})"
    t = gaia_tap.query(q, f"dr3x_hosts_{tag}").to_pandas()
    t.columns = [c.lower() for c in t.columns]
    return t.rename(columns={"source_id": "dr3"})


def wasp_hosts():
    s = pd.read_csv(CROSSMATCH, usecols=["archive_name_2019", "closest_dr3_source_id", "routes"],
                    dtype={"closest_dr3_source_id": "Int64"})
    assert len(s) == 107 and s["routes"].str.startswith("agree").all()
    return s.rename(columns={"archive_name_2019": "name", "closest_dr3_source_id": "dr3"})[["name", "dr3"]]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["hosts", "validate", "wasp"])
    ap.add_argument("--stage2", action="store_true")
    a = ap.parse_args(argv)
    if a.what == "wasp" and not a.stage2:
        raise SystemExit("WASP kinematics are Stage 2 (DECISIONS E53); refusing without --stage2")
    if not BAYESTAR_H5.is_file():
        raise SystemExit(f"STOP: Bayestar19 file missing ({BAYESTAR_H5}). Method A needs it (no zero fallback).")
    if a.what == "hosts":
        w = wasp_hosts()
        g = query_hosts(w["dr3"], "wasp107")
        h = w.merge(g, on="dr3", how="left", validate="one_to_one")
        surv, steps = host_steps(h, nea_radv(w["dr3"]))
        write_json(OUT / "method_a_wasp_host_steps.json",
                   {"steps": steps, "ext_source": surv["ext_source"].value_counts().to_dict(),
                    "rv_source": surv["rv_source"].value_counts().to_dict()})
        for s in steps:
            print(f"{s['step']}: {s['n_after']}" + (f"  removed {s['removed']}" if s["removed"] else ""))
    else:
        raise SystemExit(f"'{a.what}' is not wired up in Stage 1 (needs the published tables / Stage 2)")


if __name__ == "__main__":
    main(sys.argv[1:])
