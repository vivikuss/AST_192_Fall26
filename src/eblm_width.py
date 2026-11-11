"""EXPLORATORY (not pre-registered; requested after STOP D2): why is the EBLM control
distribution wider (SD 2.84 km/s) than the HJ hosts' (1.46 km/s)?

    python -m src.eblm_width

(a) Per-star velocity spread in the matched control pools: residual speeds |v - v_pool| of
    the control stars actually picked, for the EBLM setup (79; z + distance) and the HJ
    setup (297; z only, Phase 1), each also with the other matching scheme.
(b) The HJ-host Monte Carlo on random subsamples of 79 HJ hosts (with their own matched
    controls): the control-distribution SD, with z-only and z + distance matching.
Sample-size expectation: SD scales as 1/sqrt(N), so 1.46 * sqrt(297/79) = 2.83 km/s.
Aggregates only -> data/processed/eblm/d2/width_exploratory.json.
"""
import json

import numpy as np
import pandas as pd

from . import config, kinematics, matching, seeds
from .config import ROOT

OUT = ROOT / "data" / "processed" / "eblm" / "d2"
N_POOL_ITER = 2000       # control samples drawn per setup for (a)
N_SUB = 200              # random 79-host HJ subsamples for (b)
N_ITER_SUB = 1000        # Monte Carlo iterations per subsample


def _sampler(h, fld, cfg, match_vars):
    return matching.ControlSampler(h, fld, cfg, match_vars=match_vars)


def pool_spread(h, fld, cfg, match_vars, label):
    """(a): residual speeds of the picked control stars about the pooled mean velocity."""
    s = _sampler(h, fld, cfg, match_vars)
    rng = np.random.default_rng(np.random.SeedSequence([cfg["seed"], seeds.label_code(label), 11]))
    st = cfg["monte_carlo"]["statistic"]
    picks, disp = [], []
    for _ in range(N_POOL_ITER):
        p = s.control(rng)
        p = p[p >= 0]
        picks.append(p)
        disp.append(kinematics.dispersion(s.fU[p], s.fV[p], s.fW[p], st))
    p = np.concatenate(picks)
    U, V, W = s.fU[p], s.fV[p], s.fW[p]
    speed = np.sqrt((U - U.mean()) ** 2 + (V - V.mean()) ** 2 + (W - W.mean()) ** 2)
    n = len(h)
    return {"match_vars": match_vars, "n_hosts": n, "n_picks": int(len(p)),
            "per_star_mean_kms": float(speed.mean()), "per_star_sd_kms": float(speed.std(ddof=1)),
            "per_star_p16_p84": [float(x) for x in np.percentile(speed, [16, 84])],
            "sd_UVW_kms": [float(np.std(x, ddof=1)) for x in (U, V, W)],
            "pool_median_dist_pc": float(np.median(fld["dist_pc"].to_numpy()[s.field_idx][p])),
            "pool_median_abs_z_pc": float(np.median(np.abs(fld["z"].to_numpy()[s.field_idx][p]))),
            "predicted_sd_per_star_sd_over_sqrtN": float(speed.std(ddof=1) / np.sqrt(n)),
            "control_sd_here": float(np.std(disp, ddof=1)), "control_median_here": float(np.median(disp))}


def subsample_widths(hj, fld, cfg, match_vars, label, n=79):
    """(b): control SD (and median) for N_SUB random n-host subsamples of the HJ hosts."""
    rng = np.random.default_rng(np.random.SeedSequence([cfg["seed"], seeds.label_code(label), 12]))
    st = cfg["monte_carlo"]["statistic"]
    sds, meds, hdisp = [], [], []
    for k in range(N_SUB):
        sub = hj.iloc[np.sort(rng.choice(len(hj), n, replace=False))].reset_index(drop=True)
        r = _sampler(sub, fld, cfg, match_vars).run(N_ITER_SUB, rng, st)
        sds.append(np.std(r["dispersion"], ddof=1))
        meds.append(np.median(r["dispersion"]))
        hdisp.append(r["host_dispersion_full"][st])
    sds, meds = np.array(sds), np.array(meds)
    return {"match_vars": match_vars, "n_sub": N_SUB, "n_per_sub": n, "n_iter_each": N_ITER_SUB,
            "control_sd_median": float(np.median(sds)), "control_sd_p16_p84": [float(x) for x in np.percentile(sds, [16, 84])],
            "control_median_median": float(np.median(meds)),
            "host_sigma_minus_control_median": {"median": float(np.median(np.array(hdisp) - meds)),
                                                "p16_p84": [float(x) for x in np.percentile(np.array(hdisp) - meds, [16, 84])]}}


def main():
    cfg = config.load()
    out = config.out_dir(cfg)
    fld = pd.read_csv(out / "field_ms.csv.gz")
    hj = pd.read_csv(out / "hosts_hj_table1_final.csv")
    e = pd.read_csv(ROOT / "data" / "private" / "eblm_d2" / "eblm_final_uvw.csv")
    ids = set(pd.read_csv(ROOT / "data" / "eblm" / "targets_eblm4.csv")["source_id"].dropna().astype(np.int64))
    fld_e = fld[~fld["source_id"].isin(ids)].reset_index(drop=True)        # the D2 field (D9)
    assert len(fld_e) == 376786, len(fld_e)
    summ = json.loads((out / "summary.json").read_text())["mc"]["hj_table1"]
    R = {"label": "EXPLORATORY, not pre-registered (requested after STOP D2)",
         "reference": {"eblm_control_sd_primary": json.loads((OUT / "results.json").read_text())["primary"]["control_sd"],
                       "eblm_control_sd_z_only_G1": json.loads((OUT / "results.json").read_text())["G"]["G1_z_only"]["control_sd"],
                       "hj297_control_sd_z_only": summ["control_std"],
                       "sqrtN_scaled_hj_sd_to_79": summ["control_std"] * np.sqrt(len(hj) / len(e))}}
    zd, z = ["z", "dist_pc"], ["z"]
    R["a_pool_spread"] = {"eblm_z_dist": pool_spread(e, fld_e, cfg, zd, "w_e_zd"),
                          "eblm_z_only": pool_spread(e, fld_e, cfg, z, "w_e_z"),
                          "hj_z_only": pool_spread(hj, fld, cfg, z, "w_h_z"),
                          "hj_z_dist": pool_spread(hj, fld, cfg, zd, "w_h_zd")}
    R["b_hj_subsamples_79"] = {"z_only": subsample_widths(hj, fld, cfg, z, "w_sub_z"),
                               "z_dist": subsample_widths(hj, fld, cfg, zd, "w_sub_zd")}
    (OUT / "width_exploratory.json").write_text(json.dumps(R, indent=1, default=float))
    print(json.dumps(R, indent=1, default=float))


if __name__ == "__main__":
    main()
