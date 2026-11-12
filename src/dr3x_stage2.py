"""EXPLORATORY (not pre-registered; DECISIONS E52-E61). Stage 2 of notes/dr3_exploratory_plan.md.

    python -m src.dr3x_stage2 field   # Method B DR3 field from ~/data/dr3x_field_b -> field_b_ms.csv.gz
    python -m src.dr3x_stage2 run     # all rows -> data/processed/dr3_exploratory/stage2_results.json

Rows:
  - Method A (Hu et al. as implemented by Kayla; src.dr3x_method_a): main (close-companion exclusions) and
    without exclusions; sigma_W (primary), sigma_tot, HS19 eq. 1; one stream, label dr3x_A_wasp (E61 item 3).
  - Method B (frozen DR2 primary ported to DR3, plan section 3): eq. 1, z matching, 40,000 controls,
    Delta = sigma - median(control) with the 40,000-resample percentile bootstrap (frozen test B); label
    dr3x_B_wasp. The frozen six-case rule is reported "for comparison only" (E54 item 9).
  - Method B on DR2 for the same hosts (plan section 3.6): frozen Phase 1 rows and field; label dr3x_B_dr2match.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import (ari_tap, config, cuts, dr3x_hosts, eblm_kinematics as ek, extinction, kinematics, main_sequence,
               pipeline, wasp_kinematics as wk)
from . import dr3x_method_a as A

FIELD_DIR = Path.home() / "data" / "dr3x_field_b"
FIELD_MS = A.OUT / "field_b_ms.csv.gz"
CKPT = A.OUT / "mc_checkpoints"
FIELD_CUTS = [1, 3, 4, 5, 6, 7, 8, 9, 10]                  # cut 2 dropped (E54 item 5); cut 8 at 0.015 (item 6)
HOST_LISTS = ("hj_table1", "hj", "longp")                 # D17 lists (E61 item 2)
FIELD_KEEP = ["source_id", "l", "b", "dist_pc", "ebv", "bp_rp0", "M_G0", "U", "V", "W", "z"]
N_LEAVES = 48
FROZEN_A1 = config.ROOT / "data" / "processed" / "wasp_amendment1" / "results.json"
HU_WASP_D600 = {"n": 102, "sigma_W_hosts": 17.7, "sigma_W_field": 18.4}   # Hu et al., as given by the user


def log(m):
    print(time.strftime("%H:%M:%S"), m, flush=True)


# ----------------------------------------------------------------------------- planet-host removal (D17)
def planet_host_dr3_ids(cfg):
    base = config.out_dir(cfg)
    rows = pd.concat([pd.read_csv(base / f"hosts_{k}_after_cuts.csv", usecols=["hostname", "source_id", "gaia_dr3_id"],
                                  dtype={"source_id": "Int64", "gaia_dr3_id": str}).assign(list=k)
                      for k in HOST_LISTS], ignore_index=True)
    rows["dr3"] = rows["gaia_dr3_id"].str.extract(r"Gaia DR3 (\d+)")[0].astype("Int64")
    miss = rows[rows["dr3"].isna() & rows["source_id"].notna()]
    info = {"n_rows": {k: int((rows["list"] == k).sum()) for k in HOST_LISTS},
            "n_dr3_from_nea": int(rows["dr3"].notna().sum()), "from_dr2_neighbourhood": []}
    if len(miss):
        ids = ", ".join(str(int(i)) for i in sorted(set(miss["source_id"])))
        t = ari_tap.query("SELECT n.dr2_source_id, n.dr3_source_id, n.angular_distance FROM gaiadr3.dr2_neighbourhood "
                          f"AS n WHERE n.dr2_source_id IN ({ids})", "dr3x_planet_hosts_neighbourhood").to_pandas()
        t.columns = [c.lower() for c in t.columns]
        best = t.sort_values("angular_distance").groupby("dr2_source_id").first()["dr3_source_id"]
        for i in miss.index:
            d2 = int(rows.at[i, "source_id"])
            rows.at[i, "dr3"] = int(best[d2]) if d2 in best.index else pd.NA
            info["from_dr2_neighbourhood"].append({"hostname": rows.at[i, "hostname"], "found": d2 in best.index})
    ids = set(rows["dr3"].dropna().astype(np.int64))
    info["n_unique_dr3_ids"] = len(ids)
    return ids, info


# ----------------------------------------------------------------------------- Method B field
def process_leaf(path, cfg, exclude):
    df = pd.read_csv(path)
    n = {"rows": len(df)}

    def keep(name, m):
        nonlocal df
        df = df[np.asarray(m, dtype=bool)].reset_index(drop=True)
        n[name] = len(df)

    keep("rv present", df["radial_velocity"].notna())
    keep("astrometric_params_solved = 31", df["astrometric_params_solved"] == 31)
    keep("cuts 1,3-10 (8 at 0.015)", cuts.apply(df, cfg["cuts"], FIELD_CUTS)[0])
    keep("BJ21 r_med_geo present", df["r_med_geo"].notna())
    keep("not a planet host (D17)", ~df["source_id"].isin(exclude))
    df["dist_pc"] = df["r_med_geo"].astype(float)
    extinction.deredden(df, cfg, ebv=extinction.ebv_los(df["l"].to_numpy(), df["b"].to_numpy(),
                                                        df["dist_pc"].to_numpy(), cfg))
    keep("main sequence", main_sequence.on_main_sequence(df["M_G0"], df["bp_rp0"], cfg))
    df = kinematics.add_kinematics(df, cfg, rv_col="radial_velocity")
    return df[FIELD_KEEP], n


def build_field(cfg):
    leaves = sorted(FIELD_DIR.glob("field_*.csv.gz"))
    oks = sorted(FIELD_DIR.glob("field_*.csv.ok"))
    if len(leaves) != N_LEAVES or len(oks) != N_LEAVES:
        raise SystemExit(f"STOP: expected {N_LEAVES} verified leaves, found {len(leaves)} files / {len(oks)} .ok")
    exclude, info = planet_host_dr3_ids(cfg)
    log(f"planet-host DR3 ids: {info['n_unique_dr3_ids']} ({info})")
    parts, tot = [], {}
    for p in leaves:
        t0 = time.time()
        df, n = process_leaf(p, cfg, exclude)
        parts.append(df)
        for k, v in n.items():
            tot[k] = tot.get(k, 0) + v
        log(f"{p.name}: {n['rows']} -> {len(df)} MS ({time.time() - t0:.0f} s)")
    f = pd.concat(parts, ignore_index=True)
    assert f["source_id"].is_unique
    A.write_csv(FIELD_MS, f)
    A.write_json(A.OUT / "field_b_log.json", {"steps": tot, "planet_hosts": info, "n_ms": len(f),
                                              "leaves": [p.name for p in leaves]})
    log(f"field: {tot}")
    return f


def load_field_b():
    return pd.read_csv(FIELD_MS, comment="#")


# ----------------------------------------------------------------------------- test B (frozen)
def test_b(h, fld, cfg, label):
    """As wasp_amendment1.vs_field: Monte Carlo (eq. 1, z matching, 40,000) + percentile bootstrap."""
    n = cfg["monte_carlo"]["n_iter"]
    assert n == 40000 and cfg["monte_carlo"]["match_vars"] == ["z"] and cfg["monte_carlo"]["statistic"] == "mean_abs"
    s, _ = pipeline.monte_carlo(h, fld, cfg, label, n, "mean_abs", ckpt_dir=CKPT)
    b = ek.bootstrap(h, s["control_median"], cfg, n, label, "mean_abs")
    lo, hi = b["ci95"]
    c = wk.case(lo, hi)
    return {"n": len(h), "sigma": s["host_dispersion"], "control_median": s["control_median"],
            "control_sd": s["control_std"], "delta": b["delta_hat"], "ci95": b["ci95"], "bootstrap": b,
            "fraction_controls_le": s["p_value"], "z_gauss": s["z_gauss"], "n_iter": n,
            "iters_with_unmatched_host": s["iters_with_unmatched_host"],
            "unmatched_hosts": s["unmatched_hosts"], "hosts_in_empty_zbins": s["hosts_in_empty_zbins"],
            "n_usable_field": s["n_usable_field"],
            "six_case_for_comparison_only": {"case": c, "text": wk.CASES[c]}}


# ----------------------------------------------------------------------------- rows
def method_a_rows(base, nea_rv):
    fld, n_removed = A.load_field()
    out = {"field_rows": len(fld), "field_nea_hosts_removed": n_removed}
    for key, excl in (("main_with_companion_exclusions", True), ("without_companion_exclusions", False)):
        _, _, h = dr3x_hosts.method_a(base, nea_rv, exclude_companions=excl)
        h = h.copy()
        h["rv"] = h["dr3"].map(nea_rv).fillna(h["radial_velocity"])
        h = A.add_uvw(h)
        r = A.compare(h, fld, "dr3x_A_wasp")
        r.update(rv_source=h["rv_source"].value_counts().to_dict(), ext_source=h["ext_source"].value_counts().to_dict())
        out[key] = r
        log(f"Method A {key}: N {r['n_hosts']}, sigma_W {r['sigma_W']['host']:.2f} vs {r['sigma_W']['field_median']:.2f}")
    return out


def run(cfg):
    if not FIELD_MS.is_file():
        raise SystemExit(f"STOP: {FIELD_MS} missing; run `python -m src.dr3x_stage2 field` first")
    base, base_b, nea_rv = dr3x_hosts.bases(cfg)
    R = {"what": "Stage 2, EXPLORATORY rows (notes/dr3_exploratory_plan.md; DECISIONS E61)"}
    R["method_A"] = method_a_rows(base, nea_rv)

    _, hb = dr3x_hosts.method_b(base_b, cfg)
    hb = hb.copy()
    hb["hostname"], hb["dist_pc"] = hb["name"], hb["r_med_geo"].astype(float)
    hb["source_id"] = hb["dr3"].astype(np.int64)                # DR3 id; used in the checkpoint key
    hb = kinematics.add_kinematics(hb, cfg, rv_col="rv")
    assert len(hb) == 101 and hb[["U", "V", "W", "z", "bp_rp0"]].notna().all().all(), len(hb)
    fb = load_field_b()
    log(f"Method B DR3: {len(hb)} hosts, field {len(fb)}")
    R["method_B_dr3"] = test_b(hb, fb, cfg, "dr3x_B_wasp")
    R["method_B_dr3"]["field_n"] = len(fb)
    R["method_B_dr3"]["rv_source"] = hb["rv_source"].value_counts().to_dict()
    log(f"Method B DR3: Delta {R['method_B_dr3']['delta']:.2f} {np.round(R['method_B_dr3']['ci95'], 2).tolist()}")

    xm = pd.read_csv(A.CROSSMATCH, usecols=["archive_name_2019", "dr2_source_id"], dtype={"dr2_source_id": "Int64"})
    d2 = set(xm.loc[xm["archive_name_2019"].isin(hb["name"]), "dr2_source_id"].astype(np.int64))
    p1 = pd.read_csv(config.out_dir(cfg) / "hosts_hj_table1_final.csv", dtype={"source_id": "Int64"})
    h2 = p1[p1["source_id"].isin(d2)].reset_index(drop=True)
    assert len(d2) == 101 and len(h2) == 101, (len(d2), len(h2))
    f2 = pd.read_csv(config.out_dir(cfg) / "field_ms.csv.gz")
    R["method_B_dr2_same_hosts"] = test_b(h2, f2, cfg, "dr3x_B_dr2match")
    R["method_B_dr2_same_hosts"]["field_n"] = len(f2)
    log(f"Method B DR2 match: Delta {R['method_B_dr2_same_hosts']['delta']:.2f} "
        f"{np.round(R['method_B_dr2_same_hosts']['ci95'], 2).tolist()}")

    fz = json.loads(FROZEN_A1.read_text())["B_primary"]
    R["comparison"] = {"frozen_dr2_107_headline": {"n": fz["n"], "sigma": fz["sigma"],
                                                   "control_median": fz["control_median"],
                                                   "delta": fz["B"]["delta_hat"], "ci95": fz["B"]["ci95"],
                                                   "case": fz["case"], "source": str(FROZEN_A1.relative_to(config.ROOT))},
                       "hu2026_wasp_d600_published": HU_WASP_D600}
    A.write_json(A.OUT / "stage2_results.json", R)
    return R


def main(argv):
    cfg = config.load()
    if argv[:1] == ["field"]:
        build_field(cfg)
    elif argv[:1] == ["run"]:
        run(cfg)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
