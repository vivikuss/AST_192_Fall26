"""The frozen WASP-only analysis (PREREG_wasp_kinematics_2026-11-04.txt, sha256 fece9b66...).

    python -m src.wasp_kinematics

Section letters follow the pre-registration: B (primary), C (secondary contrast), D1-D6 (sensitivity),
E (exploratory sigma_W). Implementation details the text leaves open are in DECISIONS E30 (written
before the run). The hash is checked before anything runs. Aggregates only ->
data/processed/wasp/results.json; Monte Carlo chunks -> data/processed/wasp/mc_checkpoints/ (WASP hosts)
and data/private/eblm_d2/mc_checkpoints/ (EBLMs, as Stage D2).
"""
import hashlib
import json

import numpy as np
import pandas as pd

from . import config, eblm_kinematics as ek, hu2026_checks, kinematics, matching, pipeline, seeds
from .config import ROOT
from .gaia_tap import log

PREREG = ROOT / "PREREG_wasp_kinematics_2026-11-04.txt"
SHA = "fece9b6608624c8213d1fb5fa96a6223bc07f5328adf820bd133db4a3490ecc4"
OUT = ROOT / "data" / "processed" / "wasp"
CKPT = OUT / "mc_checkpoints"
DELTA_HS19 = 36.1825 - 42.8234          # section 2: -6.64 km/s (Phase 1 summary.json, mc.hj_table1)
RATIO_HS19 = 36.182533473514546 / 42.823367938740205      # section E: sigma_W HS19-sized ratio
ZD = ["z", "dist_pc"]


# ------------------------------------------------------------------ statistics
def total(U, V, W):
    return kinematics.dispersion(U, V, W, "mean_abs")


def sigma_w(U, V, W):
    W = np.asarray(W, float)
    return np.sqrt(((W - W.mean(-1, keepdims=True)) ** 2).mean(-1))


def boot_draws(d, cfg, n, label, fn):
    """Bootstrap statistics of a sample; same stream scheme as eblm_kinematics.bootstrap."""
    U, V, W = (d[c].to_numpy(float) for c in "UVW")
    N = len(d)
    rng = np.random.default_rng(np.random.SeedSequence([cfg["seed"], seeds.label_code(label), 99]))
    return np.concatenate([fn(*(x[rng.integers(0, N, (min(5000, n - s), N))] for x in (U, V, W)))
                           for s in range(0, n, 5000)])


def case(lo, hi, d=DELTA_HS19):
    """Section B decision rule ("contains x" means lo <= x <= hi)."""
    if lo > 0:
        return 6
    if hi < d:
        return 2
    c0, cd = lo <= 0 <= hi, lo <= d <= hi
    if c0 and cd:
        return 5
    if cd:
        return 1
    if c0:
        return 4
    return 3


CASES = {1: "Supports HS19: a deficit consistent with HS19's size.",
         2: "A cold offset larger than HS19's (supports HS19's direction; the size is not HS19's).",
         3: "A real but smaller deficit: partial support for HS19; heterogeneity explains part of the effect.",
         4: "Supports Hu et al.: no deficit, and an HS19-sized one is excluded.",
         5: "Inconclusive: no claim.",
         6: "Hotter than the field: predicted by neither paper; reported, no claim."}


# ------------------------------------------------------------------ Monte Carlo
def mc_total(h, fld, cfg, label, n_iter, match_vars):
    """HS19 Monte Carlo for the total dispersion, as Phase 1 (pipeline.monte_carlo, chunked streams)."""
    c = config.variant(cfg, **{"monte_carlo.match_vars": list(match_vars)})
    s, res = pipeline.monte_carlo(h, fld, c, label, n_iter, "mean_abs", ckpt_dir=CKPT)
    s.pop("unmatched_hosts", None)
    return s, res


def mc_sigma_w(h, fld, cfg, label, n_iter, match_vars):
    """Control sigma_W with the same matching (ControlSampler; chunk streams as D19; D13b pairing)."""
    lab = seeds.label_code(label)               # first: wasp_a1_* labels need the frozen amendment (E36)
    c = config.variant(cfg, **{"monte_carlo.match_vars": list(match_vars)})
    s = matching.ControlSampler(h, fld, c)
    every = cfg["monte_carlo"]["checkpoint_every"]
    ctrl, host = np.empty(n_iter), np.empty(n_iter)
    hW = h["W"].to_numpy(float)
    for ci, st in enumerate(range(0, n_iter, every)):
        rng = np.random.default_rng(np.random.SeedSequence([cfg["seed"], lab, ci]))
        for i in range(st, min(st + every, n_iter)):
            p = s.control(rng)
            ok = p >= 0
            ctrl[i] = sigma_w(None, None, s.fW[p[ok]])
            host[i] = sigma_w(None, None, hW[ok])
    full = float(sigma_w(None, None, hW))
    return {"host_sigma_w": full, "control_median": float(np.median(ctrl)), "control_sd": float(np.std(ctrl, ddof=1)),
            "fraction_controls_le": float(np.mean(ctrl <= host)), "n_iter": n_iter}


# ------------------------------------------------------------------ blocks
def vs_field(h, fld, cfg, label, n_mc, n_boot, match_vars):
    """B-style: Delta = sigma - median(control), percentile CI (eblm_kinematics.bootstrap), C-1 fraction."""
    s, res = mc_total(h, fld, cfg, label, n_mc, match_vars)
    b = ek.bootstrap(h, s["control_median"], cfg, n_boot, label, "mean_abs")
    lo, hi = b["ci95"]
    return {"n": len(h), "match_vars": list(match_vars), "sigma": s["host_dispersion"],
            "control_median": s["control_median"], "control_sd": s["control_std"], "B": b,
            "case": case(lo, hi), "case_text": CASES[case(lo, hi)],
            "fraction_controls_le": s["p_value"], "n_le": s["n_le"], "z_gauss": s["z_gauss"], "n_iter": n_mc,
            "iters_with_unmatched_host": s["iters_with_unmatched_host"], "hosts_in_empty_zbins": s["hosts_in_empty_zbins"]}, res


def contrast(w, e, ref_w, ref_e, cfg, n_boot, label, fn=total):
    """Delta_EBLM - Delta_W with each control median fixed; joint (independent) bootstrap."""
    dw = boot_draws(w, cfg, n_boot, label + "_W", fn) - ref_w
    de = boot_draws(e, cfg, n_boot, label + "_E", fn) - ref_e
    est = float((fn(*(e[c].to_numpy(float) for c in "UVW")) - ref_e) - (fn(*(w[c].to_numpy(float) for c in "UVW")) - ref_w))
    d = de - dw
    lo, hi = np.percentile(d, [2.5, 97.5])
    reading = ("whole CI > 0: supports an HJ-specific deficit among WASP hosts" if lo > 0 else
               "whole CI < 0: predicted by neither paper; no claim" if hi < 0 else
               "CI includes 0: inconclusive at this power")
    return {"contrast": est, "ci95": [float(lo), float(hi)], "boot_n": n_boot, "boot_sd": float(np.std(d, ddof=1)),
            "boot_bias": float(d.mean() - est), "n_wasp": len(w), "n_eblm": len(e), "reading": reading}


# ------------------------------------------------------------------ run
def main():
    got = hashlib.sha256(PREREG.read_bytes()).hexdigest()
    assert got == SHA, f"pre-registration hash mismatch: {got}"
    cfg = config.load()
    n_main = cfg["monte_carlo"]["n_iter"]
    n_sens = cfg["monte_carlo"]["sensitivity"]["n_iter"]
    assert n_main == 40000 and n_sens == 5000, (n_main, n_sens)
    OUT.mkdir(parents=True, exist_ok=True)

    # samples (A)
    base = config.out_dir(cfg)
    fin = pd.read_csv(base / "hosts_hj_table1_final.csv")
    # The frozen run (E31) used the hostname join, which misses 9 WASP hosts (DECISIONS E32). It is
    # pinned here so the reported results stay reproducible; the corrected membership (107) is for
    # the amendment's rerun, not this script.
    lab = hu2026_checks.wasp_hosts(cfg, join="hostname_legacy")["final (MS with RV; Phase 1 primary)"][
        ["hostname", "facility", "name_wasp"]]
    fin = fin.merge(lab, on="hostname", how="left", validate="one_to_one")
    w = fin[fin["facility"].notna()].reset_index(drop=True)
    w_name = fin[fin["facility"].notna() | fin["name_wasp"]].reset_index(drop=True)
    fld = pd.read_csv(base / "field_ms.csv.gz")                                   # Phase 1 field (D17)
    e = pd.read_csv(ROOT / "data" / "private" / "eblm_d2" / "eblm_final_uvw.csv")
    ids = set(pd.read_csv(ROOT / "data" / "eblm" / "targets_eblm4.csv")["source_id"].dropna().astype(np.int64))
    fld_e = fld[~fld["source_id"].isin(ids)].reset_index(drop=True)              # D2 field (D9)
    assert len(w) == 98 and len(w_name) == 100 and len(e) == 79 and len(fld_e) == 376786, \
        (len(w), len(w_name), len(e), len(fld_e))
    R = {"prereg": PREREG.name, "sha256": got, "delta_hs19": DELTA_HS19, "ratio_hs19_sigma_w": RATIO_HS19,
         "n": {"wasp": len(w), "wasp_by_name": len(w_name), "eblm": len(e), "field": len(fld), "field_eblm": len(fld_e)}}
    log(f"WASP {len(w)}  EBLM {len(e)}  field {len(fld)} / {len(fld_e)}")

    # B: primary
    B, res_b = vs_field(w, fld, cfg, "wasp_primary", n_main, n_main, ["z"])
    np.savez(OUT / "mc_wasp_primary.npz", disp=res_b["dispersion"], host=res_b["host_dispersion"])
    R["B_primary"] = B
    log(f"B: {B['B']['delta_hat']:.2f} {np.round(B['B']['ci95'], 2).tolist()} case {B['case']}")

    # C: secondary. EBLM half = the frozen EBLM primary Monte Carlo (same label, cached chunks)
    se, _ = ek.mc(e, fld_e, cfg, "eblm_primary", n_main, ZD)
    assert abs(se["control_median"] - 41.77) < 0.01, se["control_median"]
    R["C_secondary"] = contrast(w, e, B["control_median"], se["control_median"], cfg, n_main, "wasp_C")
    R["C_secondary"].update({"delta_wasp": B["B"]["delta_hat"], "delta_eblm": se["host_dispersion"] - se["control_median"],
                             "eblm_control_median": se["control_median"]})
    log(f"C: {R['C_secondary']['contrast']:.2f} {np.round(R['C_secondary']['ci95'], 2).tolist()}")

    # D: sensitivity (5,000 + 5,000). A half that a test does not change keeps its 40,000-sample reference (E30).
    D = {}
    ref_w, ref_e = B["control_median"], se["control_median"]

    def eblm_ref(sub, f, label):
        s, _ = ek.mc(sub, f, cfg, label, n_sens, ZD)
        return s["control_median"]

    b, _ = vs_field(w_name, fld, cfg, "wasp_D1", n_sens, n_sens, ["z"])
    D["D1_membership_by_name"] = {"B": b, "C": contrast(w_name, e, b["control_median"], ref_e, cfg, n_sens, "wasp_D1C")}
    b, _ = vs_field(w, fld, cfg, "wasp_D2", n_sens, n_sens, ZD)
    D["D2_wasp_z_dist"] = {"B": b, "C": contrast(w, e, b["control_median"], ref_e, cfg, n_sens, "wasp_D2C")}
    w3, e3 = w[w["dist_pc"] < 600].reset_index(drop=True), e[e["dist_pc"] < 600].reset_index(drop=True)
    b, _ = vs_field(w3, fld, cfg, "wasp_D3", n_sens, n_sens, ["z"])
    D["D3_d_lt_600"] = {"B": b, "C": contrast(w3, e3, b["control_median"], eblm_ref(e3, fld_e, "eblm_wD3"),
                                               cfg, n_sens, "wasp_D3C")}
    w4 = w[w["bp_rp0"].between(0.5, 1.5)].reset_index(drop=True)
    e4 = e[e["bp_rp0"].between(0.5, 1.5)].reset_index(drop=True)
    b, _ = vs_field(w4, fld, cfg, "wasp_D4", n_sens, n_sens, ["z"])
    D["D4_colour_0.5_1.5"] = {"B": b, "C": contrast(w4, e4, b["control_median"], eblm_ref(e4, fld_e, "eblm_wD4"),
                                                     cfg, n_sens, "wasp_D4C")}
    vt = lambda f: np.sqrt(f["U"] ** 2 + f["V"] ** 2 + f["W"] ** 2)
    f5, f5e = fld[vt(fld) <= 180].reset_index(drop=True), fld_e[vt(fld_e) <= 180].reset_index(drop=True)
    b, _ = vs_field(w, f5, cfg, "wasp_D5", n_sens, n_sens, ["z"])
    D["D5_halo_cut"] = {"B": b, "C": contrast(w, e, b["control_median"], eblm_ref(e, f5e, "eblm_wD5"),
                                               cfg, n_sens, "wasp_D5C"),
                        "field_removed": int(len(fld) - len(f5)), "field_eblm_removed": int(len(fld_e) - len(f5e))}
    e6 = e[e["short"] != "J1500-33"].reset_index(drop=True)
    D["D6_no_J1500-33_C_only"] = {"C": contrast(w, e6, ref_w, eblm_ref(e6, fld_e, "eblm_wD6"), cfg, n_sens, "wasp_D6C")}
    for k, v in D.items():
        if "B" in v:
            log(f"{k} B: n={v['B']['n']} {v['B']['B']['delta_hat']:.2f} {np.round(v['B']['B']['ci95'], 2).tolist()} case {v['B']['case']}")
        log(f"{k} C: {v['C']['contrast']:.2f} {np.round(v['C']['ci95'], 2).tolist()}")
    R["D_sensitivity"] = D

    # E: exploratory sigma_W (matching of A; 40,000 controls and resamples)
    sw_w = mc_sigma_w(w, fld, cfg, "wasp_E_W", n_main, ["z"])
    sw_e = mc_sigma_w(e, fld_e, cfg, "wasp_E_E", n_main, ZD)
    E = {"label": "EXPLORATORY (pre-registered section E); a null in sigma_W is uninformative"}
    dW = boot_draws(w, cfg, n_main, "wasp_Ei", sigma_w) - sw_w["control_median"]
    E["i_wasp_vs_field"] = {**sw_w, "delta": sw_w["host_sigma_w"] - sw_w["control_median"],
                            "ci95": [float(x) for x in np.percentile(dW, [2.5, 97.5])],
                            "hs19_sized_delta": -(1 - RATIO_HS19) * sw_w["control_median"]}
    E["ii_contrast"] = contrast(w, e, sw_w["control_median"], sw_e["control_median"], cfg, n_main, "wasp_Eii", sigma_w)
    bE = boot_draws(e, cfg, n_main, "wasp_Eiii", sigma_w)
    E["iii_eblm"] = {**sw_e, "ci95_sigma_w": [float(x) for x in np.percentile(bE, [2.5, 97.5])],
                     "delta": sw_e["host_sigma_w"] - sw_e["control_median"],
                     "ci95_delta": [float(x) for x in np.percentile(bE - sw_e["control_median"], [2.5, 97.5])]}
    for k in ("i_wasp_vs_field", "iii_eblm"):
        lo, hi = E[k]["ci95"] if k.startswith("i_") else E[k]["ci95_delta"]
        E[k]["ci_excludes_0"] = bool(lo > 0 or hi < 0)
    lo, hi = E["ii_contrast"]["ci95"]
    E["ii_contrast"]["ci_excludes_0"] = bool(lo > 0 or hi < 0)
    R["E_exploratory_sigma_w"] = E
    log(f"E: {json.dumps({k: round(v.get('delta', v.get('contrast', 0)), 2) for k, v in E.items() if isinstance(v, dict)})}")

    (OUT / "results.json").write_text(json.dumps(R, indent=1, default=pipeline._json))
    log("results written")


if __name__ == "__main__":
    main()
