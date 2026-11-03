"""STOP D1: aggregate summary of the gamma fits and their cross-checks (no per-system values).

    python -m src.gamma_fit --summary

Reads the private per-system table and prints/saves aggregates only
(data/processed/eblm/gamma_d1_summary.json). The poor-fit criteria are descriptive
(they exclude nothing; pre-registration A excludes only uncomputable fits). They were
fixed before the summary was first run (DECISIONS D3).
"""
import json

import numpy as np
import pandas as pd

from .config import ROOT
from .gamma_fit import k_errors, systems

OUT = ROOT / "data" / "processed" / "eblm"
K_FRAC = 0.05          # |K_fit - K_pub| / K_pub above this: K disagrees
K_PULL = 5.0           # ... and more than this many published sigma
SIG_BIG = 0.5          # sigma_gamma above this [km/s]: poorly constrained
P_TOL = 0.01           # DR3 NSS period within 1 % of the EBLM IV period: same (inner) orbit


def _dist(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return {"n": 0}
    q = np.quantile(x, [0.1, 0.5, 0.9])
    return {"n": int(len(x)), "p10": round(q[0], 4), "median": round(q[1], 4), "p90": round(q[2], 4),
            "min": round(x.min(), 4), "max": round(x.max(), 4)}


def _robust_sd(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if len(x) else np.nan


def _compare(d, ref, ref_err, label):
    """gamma_fit - reference, its robust spread and pulls."""
    delta = d["gamma"] - d[ref]
    pull = delta / np.hypot(d["sigma_gamma"], d[ref_err].fillna(0))
    ok = delta.notna()
    return {label: {"n": int(ok.sum()), "delta_kms": _dist(delta), "robust_sd_kms": round(_robust_sd(delta), 4),
                    "abs_delta_gt_1kms": int((delta.abs() > 1).sum()), "abs_pull_gt_3": int((pull.abs() > 3).sum()),
                    "pull_robust_sd": round(_robust_sd(pull), 3)}}


def summary(cfg, od):
    f = pd.read_csv(od / "gamma_fits.csv")
    s = systems().set_index("name")
    s["K_err_pub"] = k_errors()
    in79 = set(pd.read_csv(OUT / "eblm4_sample.csv")["name"].str.replace("EBLM ", "", regex=False))
    prim = f[f["variant"] == "primary"].set_index("name")
    prim = prim.join(s[["K", "K_err_pub", "flag", "P_obs"]].rename(columns={"K": "K_pub", "flag": "eblm4_flag"}))
    prim["in79"] = prim.index.isin(in79)
    prim["dK_frac"] = (prim["K"] - prim["K_pub"]) / prim["K_pub"]
    prim["K_pull"] = (prim["K"] - prim["K_pub"]) / prim["K_err_pub"]
    prim["poor_K"] = (prim["dK_frac"].abs() > K_FRAC) & (prim["K_pull"].abs() > K_PULL)
    prim["poor_sig"] = prim["sigma_gamma"] > SIG_BIG
    out = {"n_systems": int(len(prim)), "flags": prim["flag"].value_counts().to_dict()}
    for lab, d in (("all118", prim), ("kin79", prim[prim["in79"]])):
        ok = d[d["flag"] == "ok"]
        out[lab] = {
            "n": int(len(d)), "n_ok": int(len(ok)),
            "sigma_gamma_kms": _dist(ok["sigma_gamma"]),
            "sigma_gamma_stat_kms": _dist(np.sqrt(np.clip(ok["sigma_gamma"] ** 2 - 0.020 ** 2, 0, None))),
            "n_sigma_gt_0.1": int((ok["sigma_gamma"] > 0.1).sum()),
            "n_sigma_gt_0.5": int((ok["sigma_gamma"] > 0.5).sum()),
            "n_sigma_gt_1": int((ok["sigma_gamma"] > 1).sum()),
            "gamma_minus_mean_kms": _dist(ok["gamma_minus_mean"]),
            "abs_gamma_minus_mean_kms": _dist(ok["gamma_minus_mean"].abs()),
            "dK_frac": _dist(ok["dK_frac"]),
            "abs_dK_frac": _dist(ok["dK_frac"].abs()),
            "n_poor_K": int(ok["poor_K"].sum()), "n_poor_sigma": int(ok["poor_sig"].sum()),
            "n_poor_any": int((ok["poor_K"] | ok["poor_sig"]).sum()),
            "poor_by_eblm4_flag": ok[ok["poor_K"] | ok["poor_sig"]]["eblm4_flag"].fillna("none").value_counts().to_dict(),
            "gamma_kms": _dist(ok["gamma"]),
        }
    # variants vs primary
    v = f.pivot(index="name", columns="variant", values="gamma")
    out["variants_minus_primary"] = {}
    for c in v.columns.drop("primary"):
        dd = v[c] - v["primary"]
        sub = {"all118": _dist(dd), "kin79": _dist(dd[dd.index.isin(in79)]),
               "n_abs_gt_0.1": int((dd.abs() > 0.1).sum()), "n_abs_gt_1": int((dd.abs() > 1).sum()),
               "n_abs_gt_0.1_kin79": int((dd[dd.index.isin(in79)].abs() > 0.1).sum())}
        out["variants_minus_primary"][c] = sub
    # Gaia DR3 two-body orbits (SB1 with a centre-of-mass velocity; inner period)
    nss = pd.read_csv(OUT / "eblm4_dr3_nss.csv")
    nss = nss[nss["center_of_mass_velocity"].notna()].copy()
    nss = nss.join(s["P_obs"], on="name")
    nss["same_P"] = (nss["period"] / nss["P_obs"] - 1).abs() < P_TOL
    out["dr3_nss"] = {"systems_with_solution": int(pd.read_csv(OUT / "eblm4_dr3_nss.csv")["name"].nunique()),
                      "with_com_velocity": int(nss["name"].nunique()), "same_period": int(nss["same_P"].sum())}
    n1 = nss[nss["same_P"]].drop_duplicates("name").set_index("name")
    d = prim[prim["flag"] == "ok"].join(n1[["center_of_mass_velocity", "center_of_mass_velocity_error",
                                           "semi_amplitude_primary"]], how="inner")
    out["dr3_nss"].update(_compare(d, "center_of_mass_velocity", "center_of_mass_velocity_error", "all"))
    out["dr3_nss"].update(_compare(d[d["in79"]], "center_of_mass_velocity", "center_of_mass_velocity_error", "kin79"))
    out["dr3_nss"]["K_fit_minus_K_nss"] = _dist(d["K"] - d["semi_amplitude_primary"])
    # published gamma values (C2/C9), if compiled
    pg = OUT / "published_gamma.csv"
    if pg.exists():
        p = pd.read_csv(pg)
        p = p[pd.to_numeric(p["gamma_kms"], errors="coerce").notna()]
        p["gamma_kms"] = p["gamma_kms"].astype(float)
        p["gamma_err_kms"] = pd.to_numeric(p["gamma_err_kms"], errors="coerce")
        # primary-star orbital gammas only: drop secondary-star CCF velocities and values the
        # compilation derived itself (gamma + offset), and rows marked EXCLUDED (false name matches);
        # then one row per system, CORALIE first
        n_rows = len(p)
        p = p[~p["notes"].astype(str).str.contains("SECONDARY|DERIVED|EXCLUDED", case=True)]
        inst = p["instrument"].astype(str)
        p["prio"] = np.select([inst.str.match(r"^CORALIE(07)?$"), inst.str.contains("CORALIE"),
                               inst.str.startswith("HARPS")], [0, 1, 2], 3)
        p1 = p.sort_values(["name", "prio"]).drop_duplicates("name").set_index("name")
        d = prim[prim["flag"] == "ok"].join(p1[["gamma_kms", "gamma_err_kms", "instrument", "prio"]], how="inner")
        d["in79"] = d["in79"].astype(bool)
        out["published"] = {"rows_compiled": int(n_rows), "rows_primary_star": int(len(p)),
                            "systems": int(len(d)), "systems_kin79": int(d["in79"].sum()),
                            "reference_instrument": d["instrument"].value_counts().to_dict()}
        out["published"].update(_compare(d, "gamma_kms", "gamma_err_kms", "one_per_system"))
        out["published"].update(_compare(d[d["prio"] <= 1], "gamma_kms", "gamma_err_kms", "coralie_only"))
        out["published"].update(_compare(d[d["prio"] > 1], "gamma_kms", "gamma_err_kms", "other_instruments"))
    # names only (no values) of the systems behind each count, for the terminal (not saved)
    names = {"poor_K": sorted(prim.index[prim["poor_K"] & (prim["flag"] == "ok")]),
             "poor_sigma": sorted(prim.index[prim["poor_sig"] & (prim["flag"] == "ok")]),
             "not_ok": sorted(prim.index[prim["flag"] != "ok"])}
    (OUT / "gamma_d1_summary.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, indent=1, default=str))
    print("names (terminal only):", names)
