"""POST-HOC DIAGNOSTIC of the EXPLORATORY Stage 2 (DECISIONS E62): why the Method B DR3 control median exceeds
DR2's. Changes nothing in the Stage 2 results (E61).

    python -m src.dr3x_diag     # -> data/processed/dr3_exploratory/diag_results.json

1. Replays chunk 0 (5,000 iterations) of dr3x_B_dr2match (DR2) and dr3x_B_wasp (DR3), checks the replay against
   the saved checkpoint, and summarises the picked control stars.
2. Method B DR3 with the field restricted to stars also in the DR2 Phase 1 field (dr2_neighbourhood, closest
   match), 40,000 + 40,000, label dr3x_B_diag_dr2field.
3. Method B DR3 with the D5 halo cut (v_tot <= 180 km/s), 5,000 + 5,000 as the frozen D5, label dr3x_B_diag_halo.
"""
import json

import numpy as np
import pandas as pd

from . import (ari_tap, config, cuts, dr3x_hosts, eblm_kinematics as ek, field, kinematics, matching, pipeline,
               seeds, wasp_kinematics as wk)
from . import dr3x_method_a as A
from . import dr3x_stage2 as S

DIAG = "POST-HOC DIAGNOSTIC (E62)"
QS = (5, 16, 50, 84, 95)
XM_NAME = "dr3x_diag_dr2_dr3_xmatch"


def hosts(cfg):
    _, base_b, _ = dr3x_hosts.bases(cfg)
    _, hb = dr3x_hosts.method_b(base_b, cfg)
    hb = hb.copy()
    hb["hostname"], hb["dist_pc"] = hb["name"], hb["r_med_geo"].astype(float)
    hb["source_id"] = hb["dr3"].astype(np.int64)
    hb = kinematics.add_kinematics(hb, cfg, rv_col="rv")
    xm = pd.read_csv(A.CROSSMATCH, usecols=["archive_name_2019", "dr2_source_id"], dtype={"dr2_source_id": "Int64"})
    d2 = set(xm.loc[xm["archive_name_2019"].isin(hb["name"]), "dr2_source_id"].astype(np.int64))
    p1 = pd.read_csv(config.out_dir(cfg) / "hosts_hj_table1_final.csv", dtype={"source_id": "Int64"})
    h2 = p1[p1["source_id"].isin(d2)].reset_index(drop=True)
    assert len(hb) == 101 and len(h2) == 101
    return hb, h2


def test_b_n(h, fld, cfg, label, n):
    s, _ = pipeline.monte_carlo(h, fld, cfg, label, n, "mean_abs", ckpt_dir=S.CKPT)
    b = ek.bootstrap(h, s["control_median"], cfg, n, label, "mean_abs")
    c = wk.case(*b["ci95"])
    return {"n": len(h), "field_n": len(fld), "sigma": s["host_dispersion"], "control_median": s["control_median"],
            "delta": b["delta_hat"], "ci95": b["ci95"], "n_iter": n, "fraction_controls_le": s["p_value"],
            "iters_with_unmatched_host": s["iters_with_unmatched_host"], "n_usable_field": s["n_usable_field"],
            "six_case_for_comparison_only": c}


# ----------------------------------------------------------------------------- 1. replay
def replay(h, fld, cfg, label, n=5000):
    smp = matching.ControlSampler(h, fld, cfg)
    rng = np.random.default_rng(np.random.SeedSequence([cfg["seed"], seeds.label_code(label), 0]))
    picks = np.stack([smp.control(rng) for _ in range(n)])
    assert (picks >= 0).all(), "unmatched host in replay"
    disp = kinematics.dispersion(smp.fU[picks], smp.fV[picks], smp.fW[picks], "mean_abs")
    ck = sorted(S.CKPT.glob(f"{label}_*/chunk0000_{n}.npz"))
    assert len(ck) == 1, ck
    saved = np.load(ck[0])["disp_mean_abs"]
    assert np.allclose(disp, saved, rtol=0, atol=1e-9), np.abs(disp - saved).max()
    rows = smp.field_idx[picks.ravel()]
    return rows, float(np.abs(disp - saved).max())


def summary(v):
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    return {"n": int(len(v)), "mean": float(v.mean()), **{f"p{q}": float(np.percentile(v, q)) for q in QS}}


def matched_controls(cfg, hb, h2, fb, f2):
    out = {}
    # DR2
    r2, e2 = replay(h2, f2, cfg, "dr3x_B_dr2match")
    p2 = f2.iloc[r2][["source_id", "dist_pc", "z", "U", "V", "W"]].reset_index(drop=True)
    raw2 = field.download(cfg)                                  # cached DR2 prefilter (Phase 1)
    raw2 = raw2[["source_id", "phot_g_mean_mag", "radial_velocity_error", "parallax", "parallax_error"]]
    raw2 = raw2[raw2["source_id"].isin(set(p2["source_id"]))]
    p2 = p2.merge(raw2, on="source_id", how="left", validate="many_to_one")
    p2["parallax_over_error"] = p2["parallax"] / p2["parallax_error"]
    p2["grvs_mag"] = np.nan
    # DR3
    r3, e3 = replay(hb, fb, cfg, "dr3x_B_wasp")
    p3 = fb.iloc[r3][["source_id", "dist_pc", "z", "U", "V", "W"]].reset_index(drop=True)
    want = set(p3["source_id"])
    raw3 = pd.concat([d[d["source_id"].isin(want)] for d in
                      (pd.read_csv(f, usecols=["source_id", "phot_g_mean_mag", "grvs_mag", "radial_velocity_error",
                                               "parallax_over_error"]) for f in sorted(S.FIELD_DIR.glob("field_*.csv.gz")))])
    p3 = p3.merge(raw3, on="source_id", how="left", validate="many_to_one")
    for tag, p, err in (("DR2", p2, e2), ("DR3", p3, e3)):
        p["abs_z"] = p["z"].abs()
        o = {"replay_max_abs_diff_vs_checkpoint": err, "n_picks": len(p), "n_unique_stars": int(p["source_id"].nunique())}
        for c in ("dist_pc", "phot_g_mean_mag", "grvs_mag", "radial_velocity_error", "parallax_over_error", "abs_z"):
            o[c] = summary(p[c]) if p[c].notna().any() else "not available"
        o["rms_UVW_of_picks"] = {c: float(np.sqrt(np.mean((p[c] - p[c].mean()) ** 2))) for c in "UVW"}
        o["mean_abs_eq1_of_all_picks"] = float(kinematics.dispersion(p["U"], p["V"], p["W"], "mean_abs"))
        o["frac_v_tot_gt_180"] = float((np.sqrt(p["U"] ** 2 + p["V"] ** 2 + p["W"] ** 2) > 180).mean())
        out[tag] = o
        S.log(f"replay {tag}: {len(p)} picks, checkpoint diff {err:.1e}")
    return out


# ----------------------------------------------------------------------------- 2. DR2-field stars only
def dr2_field_dr3_ids(cfg, f2):
    where = cuts.server_where(cfg["cuts"]).replace("{", "{{").replace("}", "}}")
    q = ("SELECT n.dr2_source_id, n.dr3_source_id, n.angular_distance\n"
         "FROM gaiadr2.gaia_source AS g JOIN gaiadr3.dr2_neighbourhood AS n ON n.dr2_source_id = g.source_id\n"
         f"WHERE {{range}}\n AND {where}")
    tbl, failed = ari_tap.query_ranges(q, XM_NAME, cfg["gaia"]["ari"]["range_level"])
    if failed:
        raise SystemExit(f"STOP: crossmatch ranges failed: {[f[0] for f in failed]}")
    t = tbl.to_pandas()
    t.columns = [c.lower() for c in t.columns]
    t = t[t["dr2_source_id"].isin(set(f2["source_id"]))]
    best = t.sort_values("angular_distance").groupby("dr2_source_id", as_index=False).first()
    info = {"dr2_field": len(f2), "dr2_field_with_dr3_candidate": len(best),
            "dr3_ids_unique": int(best["dr3_source_id"].nunique()),
            "candidates_per_dr2_star_max": int(t.groupby("dr2_source_id").size().max())}
    return set(best["dr3_source_id"].astype(np.int64)), info


def main():
    cfg = config.load()
    hb, h2 = hosts(cfg)
    fb = S.load_field_b()
    f2 = pd.read_csv(config.out_dir(cfg) / "field_ms.csv.gz")
    R = {"what": f"{DIAG}: Method B DR3 vs DR2 control medians", "label": DIAG}
    R["1_matched_controls"] = matched_controls(cfg, hb, h2, fb, f2)

    ids, info = dr2_field_dr3_ids(cfg, f2)
    fr = fb[fb["source_id"].isin(ids)].reset_index(drop=True)
    info["dr3_field_kept"] = len(fr)
    S.log(f"DR3 field restricted to DR2-field stars: {info}")
    R["2_dr3_field_restricted_to_dr2_field"] = {"crossmatch": info,
                                                **test_b_n(hb, fr, cfg, "dr3x_B_diag_dr2field", 40000)}

    vt = np.sqrt(fb["U"] ** 2 + fb["V"] ** 2 + fb["W"] ** 2)
    f5 = fb[vt <= 180].reset_index(drop=True)
    R["3_d5_halo_cut"] = {"field_removed": int(len(fb) - len(f5)), **test_b_n(hb, f5, cfg, "dr3x_B_diag_halo", 5000)}

    st = json.loads((A.OUT / "stage2_results.json").read_text())
    R["reference_stage2"] = {k: {"control_median": st[k]["control_median"], "delta": st[k]["delta"],
                                 "ci95": st[k]["ci95"], "field_n": st[k]["field_n"]}
                             for k in ("method_B_dr3", "method_B_dr2_same_hosts")}
    A.write_json(A.OUT / "diag_results.json", R)
    for k in ("2_dr3_field_restricted_to_dr2_field", "3_d5_halo_cut"):
        r = R[k]
        S.log(f"{k}: field {r['field_n']}, median {r['control_median']:.2f}, Delta {r['delta']:.2f} "
              f"{np.round(r['ci95'], 2).tolist()}")


if __name__ == "__main__":
    main()
