"""Gaia DR2 -> DR3 crossmatch of the 107 WASP hosts (frozen amendment 1 membership). Phase 2 groundwork.

    python -m src.dr3_crossmatch           # -> data/processed/dr3_crossmatch/  (gitignored)

IDs, positions and magnitudes only (DECISIONS E48). Every table is read with an explicit column list, so no
velocity, radial-velocity or kinematic column of any host is loaded:
  - membership: hosts_hj_table1_final.csv (archive_name_2019, hostname, source_id, ra, dec, phot_g_mean_mag) plus
    the archive WASP facility label joined on the DR2 source_id (hu2026_checks.facility_map, E32/E36), rebuilt
    here exactly as wasp_amendment1.samples() does, but without loading the other columns. Must give 107 hosts.
  - route 1 (primary): host DR2 source_id -> NEA gaia_dr2_id -> NEA gaia_dr3_id, from Kayla's nea_pscomppars.csv
    (~/data/hu2026_dr3_Kayla/), columns pl_name, hostname, gaia_dr2_id, gaia_dr3_id only.
  - route 2 (check): gaiadr3.dr2_neighbourhood for every host DR2 source_id, all candidates kept, joined to
    gaiadr3.gaia_source for ra, dec, phot_g_mean_mag of each candidate only. The raw answer is cached.
"""
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from . import config, hu2026_checks
from .config import ROOT
from .hosts import RAW, _dr2_ids, select_planets

OUT = ROOT / "data" / "processed" / "dr3_crossmatch"
NEA_KAYLA = Path.home() / "data" / "hu2026_dr3_Kayla" / "nea_pscomppars.csv"
ARCHIVE = RAW / "exoarchive" / "pscomppars_2026-10-28.csv"     # hosts.archive_table snapshot
TAP = "Gaia archive via src.gaia_tap (astroquery, async)"
HOST_COLS = ["archive_name_2019", "hostname", "source_id", "ra", "dec", "phot_g_mean_mag"]
FORBIDDEN = re.compile(r"radial_velocity|^rv|_rv|radv|^[UVW]$|v_tot|velocit", re.I)


def _cols_ok(df):
    bad = [c for c in df.columns if FORBIDDEN.search(c)]
    assert not bad, f"velocity-like columns loaded: {bad}"
    return df


def membership(cfg):
    """The 107 WASP hosts of amendment 1 (source_id join), IDs/positions/magnitudes only."""
    fin = _cols_ok(pd.read_csv(config.out_dir(cfg) / "hosts_hj_table1_final.csv", usecols=HOST_COLS,
                               dtype={"source_id": "Int64"}))
    pl = _cols_ok(pd.read_csv(ARCHIVE, usecols=lambda c: not FORBIDDEN.search(c)))
    hj = select_planets(pl, cfg, "hj")
    disc = hu2026_checks.disc_table(cfg).set_index("pl_name")
    fac = hu2026_checks.facility_map(hj, disc, "source_id")
    fin["facility"] = fin["source_id"].map(fac).fillna("").map(hu2026_checks.facility_label)
    w = fin[fin["facility"].notna()].reset_index(drop=True)
    assert len(w) == 107 and w["source_id"].is_unique and w["source_id"].notna().all(), len(w)
    return w.rename(columns={"source_id": "dr2_source_id", "ra": "dr2_ra", "dec": "dr2_dec",
                             "phot_g_mean_mag": "dr2_g"})


def _ids(series, tag):
    return series.astype(str).str.extract(rf"Gaia {tag} (\d+)")[0].astype("Int64")


def route_nea(w):
    """Route 1: DR2 id -> every NEA row with that gaia_dr2_id -> its gaia_dr3_id(s)."""
    n = pd.read_csv(NEA_KAYLA, usecols=["pl_name", "hostname", "gaia_dr2_id", "gaia_dr3_id"], dtype=str)
    n = n.assign(dr2=_ids(n["gaia_dr2_id"], "DR2"), dr3=_ids(n["gaia_dr3_id"], "DR3")).dropna(subset=["dr2"])
    rows = []
    for _, h in w.iterrows():
        m = n[n["dr2"] == h["dr2_source_id"]]
        d3 = sorted(set(m["dr3"].dropna().astype(np.int64)))
        status = ("no NEA row with this gaia_dr2_id" if m.empty else "NEA row has no gaia_dr3_id" if not d3
                  else "one DR3 id" if len(d3) == 1 else "several DR3 ids")
        rows.append({"dr2_source_id": h["dr2_source_id"], "nea_hostnames": "; ".join(sorted(set(m["hostname"]))),
                     "nea_n_planet_rows": len(m), "nea_dr3_ids": ";".join(map(str, d3)),
                     "nea_dr3_source_id": d3[0] if len(d3) == 1 else pd.NA, "nea_status": status})
    out = pd.DataFrame(rows)
    out["nea_dr3_source_id"] = out["nea_dr3_source_id"].astype("Int64")
    return out


def neighbourhood_query(ids):
    return ("SELECT n.dr2_source_id, n.dr3_source_id, n.angular_distance, n.magnitude_difference, "
            "n.proper_motion_propagation, g.ra AS dr3_ra, g.dec AS dr3_dec, g.phot_g_mean_mag AS dr3_g "
            "FROM gaiadr3.dr2_neighbourhood AS n JOIN gaiadr3.gaia_source AS g ON g.source_id = n.dr3_source_id "
            f"WHERE n.dr2_source_id IN ({', '.join(str(int(i)) for i in ids)})")


def route_neighbourhood(w, refresh=False):
    """Route 2: all dr2_neighbourhood candidates per host (cached raw answer + query + time)."""
    raw, meta = OUT / "neighbourhood_raw.csv", OUT / "neighbourhood_query.json"
    q = neighbourhood_query(sorted(w["dr2_source_id"].astype(np.int64)))
    if refresh or not raw.exists() or json.loads(meta.read_text())["query"] != q:
        from . import gaia_tap                     # async job, retries, cache in data/raw/gaia (A5)
        t = gaia_tap.query(q, "dr3_crossmatch_wasp107_neighbourhood", refresh=refresh).to_pandas()
        t.columns = [c.lower() for c in t.columns]
        _cols_ok(t)
        OUT.mkdir(parents=True, exist_ok=True)
        t.to_csv(raw, index=False)
        meta.write_text(json.dumps({"tap": TAP, "query": q, "utc": datetime.now(timezone.utc).isoformat(
            timespec="seconds"), "rows": len(t), "query_sha256": hashlib.sha256(q.encode()).hexdigest()}, indent=1))
    assert json.loads(meta.read_text())["query"] == q, "cached neighbourhood answer is for a different host list"
    t = _cols_ok(pd.read_csv(raw, dtype={"dr2_source_id": "Int64", "dr3_source_id": "Int64"}))
    assert set(t["dr2_source_id"]) <= set(w["dr2_source_id"])
    return t


def build(cfg, refresh=False):
    w = membership(cfg)
    a = route_nea(w)
    t = route_neighbourhood(w, refresh)
    t["rank_by_distance"] = t.groupby("dr2_source_id")["angular_distance"].rank(method="first").astype(int)
    t["n_candidates"] = t.groupby("dr2_source_id")["dr3_source_id"].transform("size")
    host = w[["archive_name_2019", "hostname", "facility", "dr2_source_id", "dr2_ra", "dr2_dec", "dr2_g"]] \
        .merge(a, on="dr2_source_id", how="left", validate="one_to_one")
    long = host.merge(t, on="dr2_source_id", how="left", validate="one_to_many")
    long["n_candidates"] = long["n_candidates"].fillna(0).astype(int)
    long["is_nea_dr3"] = long["dr3_source_id"].eq(long["nea_dr3_source_id"]).fillna(False).astype(bool)
    long["same_id_dr2_dr3"] = long["dr3_source_id"].eq(long["dr2_source_id"]).fillna(False).astype(bool)

    s = []
    for k, g in long.groupby("dr2_source_id", sort=False):
        h = g.iloc[0]
        c = g.dropna(subset=["dr3_source_id"]).sort_values("angular_distance")
        closest = c.iloc[0]["dr3_source_id"] if len(c) else pd.NA
        nea = h["nea_dr3_source_id"]
        if pd.isna(nea) and not len(c):
            agree = "no DR3 match on either route"
        elif pd.isna(nea):
            agree = "NEA route gives no single DR3 id"
        elif not len(c):
            agree = "no neighbourhood candidate"
        elif nea == closest:
            agree = "agree (NEA id = closest candidate)"
        elif (c["dr3_source_id"] == nea).any():
            agree = "DISAGREE: NEA id is a candidate but not the closest"
        else:
            agree = "DISAGREE: NEA id not among candidates"
        s.append({"dr2_source_id": k, "n_candidates": len(c), "closest_dr3_source_id": closest,
                  "closest_angular_distance_mas": c.iloc[0]["angular_distance"] if len(c) else np.nan,
                  "closest_magnitude_difference": c.iloc[0]["magnitude_difference"] if len(c) else np.nan,
                  "routes": agree})
    s = host.merge(pd.DataFrame(s), on="dr2_source_id", validate="one_to_one")
    s["closest_dr3_source_id"] = s["closest_dr3_source_id"].astype("Int64")
    return w, long, s


def main(argv=()):
    cfg = config.load()
    w, long, s = build(cfg, refresh="--refresh" in argv)
    OUT.mkdir(parents=True, exist_ok=True)
    long.to_csv(OUT / "wasp107_dr2_dr3_candidates.csv", index=False)
    s.to_csv(OUT / "wasp107_dr2_dr3_summary.csv", index=False)
    print(f"hosts {len(s)}; neighbourhood rows {long['dr3_source_id'].notna().sum()}")
    print("NEA route:", s["nea_status"].value_counts().to_dict())
    print("candidates per host:", s["n_candidates"].value_counts().sort_index().to_dict())
    print("routes:", s["routes"].value_counts().to_dict())
    print("DR3 id == DR2 id for the chosen match:",
          int((s["closest_dr3_source_id"] == s["dr2_source_id"]).sum()))
    d = long.dropna(subset=["dr3_source_id"])
    print("angular distance (mas) all candidates: min %.3f median %.3f max %.3f" %
          tuple(d["angular_distance"].quantile([0, .5, 1])))
    print("|magnitude difference| all candidates: min %.3f median %.3f max %.3f" %
          tuple(d["magnitude_difference"].abs().quantile([0, .5, 1])))
    for k in ("wasp107_dr2_dr3_candidates.csv", "wasp107_dr2_dr3_summary.csv", "neighbourhood_raw.csv"):
        print(k, hashlib.sha256((OUT / k).read_bytes()).hexdigest())


if __name__ == "__main__":
    main(sys.argv[1:])
