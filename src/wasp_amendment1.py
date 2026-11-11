"""Rerun under amendment 1 (PREREG_wasp_kinematics_amendment1_2026-11-10.txt, sha256 533a6437...).

    python -m src.wasp_amendment1

Section 4 of the amendment: B, C, D1-D6, E(i)-(ii) with the corrected membership (source_id join, 107 hosts;
D1 by HS19 ArchiveName, 109), and the DIAGNOSTIC: the original 98-host list (hostname join; D1 100) through
the same computations under the same wasp_a1_* labels. E(iii) and every EBLM reference are taken from the
frozen run (E31): eblm_primary / eblm_wD3-wD6 control medians from their cached frozen chunks, wasp_E_E sigma_W
control median and E(iii) from data/processed/wasp/results.json.

The frozen functions of src/wasp_kinematics.py are reused unchanged. vs_field is repeated here only so that the
WASP-half checkpoints go to data/processed/wasp_amendment1/mc_checkpoints/ (frozen outputs untouched).
The guard (seeds.require_frozen_amendment) runs first; every wasp_a1_* stream checks it again (seeds.label_code).
Aggregates only -> data/processed/wasp_amendment1/results.json and diagnostic_98/results.json. DECISIONS E40.
"""
import hashlib
import json

import numpy as np
import pandas as pd

from . import config, eblm_kinematics as ek, hu2026_checks, pipeline, seeds
from . import wasp_kinematics as wk
from .config import ROOT
from .gaia_tap import log

AMEND = ROOT / "PREREG_wasp_kinematics_amendment1_2026-11-10.txt"
AMEND_SHA = "533a6437d3bd32348a71bf7436647113ffe761c5164c59076053cde456be4e21"
OUT = ROOT / "data" / "processed" / "wasp_amendment1"
CKPT = OUT / "mc_checkpoints"
FROZEN = ROOT / "data" / "processed" / "wasp" / "results.json"
EBLM_CKPT = ROOT / "data" / "private" / "eblm_d2" / "mc_checkpoints"
F = "final (MS with RV; Phase 1 primary)"
ZD = wk.ZD


def vs_field(h, fld, cfg, label, n_mc, n_boot, match_vars):
    """wk.vs_field with the WASP-half checkpoints in the amendment directory (logic unchanged)."""
    c = config.variant(cfg, **{"monte_carlo.match_vars": list(match_vars)})
    s, res = pipeline.monte_carlo(h, fld, c, label, n_mc, "mean_abs", ckpt_dir=CKPT)
    s.pop("unmatched_hosts", None)
    b = ek.bootstrap(h, s["control_median"], cfg, n_boot, label, "mean_abs")
    lo, hi = b["ci95"]
    return {"n": len(h), "match_vars": list(match_vars), "sigma": s["host_dispersion"],
            "control_median": s["control_median"], "control_sd": s["control_std"], "B": b,
            "case": wk.case(lo, hi), "case_text": wk.CASES[wk.case(lo, hi)],
            "fraction_controls_le": s["p_value"], "n_le": s["n_le"], "z_gauss": s["z_gauss"], "n_iter": n_mc,
            "iters_with_unmatched_host": s["iters_with_unmatched_host"],
            "hosts_in_empty_zbins": s["hosts_in_empty_zbins"]}


def eblm_ref(sub, f, cfg, label, n_iter):
    """Frozen EBLM control median from its cached chunks; refuses to recompute (amendment section 4)."""
    before = {p.name for p in EBLM_CKPT.iterdir()}
    s, _ = ek.mc(sub, f, cfg, label, n_iter, ZD)
    after = {p.name for p in EBLM_CKPT.iterdir()}
    assert before == after, f"{label}: no cached frozen chunks; the EBLM reference would be recomputed"
    return s["control_median"]


def samples(cfg, join):
    base = config.out_dir(cfg)
    fin = pd.read_csv(base / "hosts_hj_table1_final.csv")
    lab = hu2026_checks.wasp_hosts(cfg, join=join)[F][["source_id", "facility", "name_wasp"]]
    fin["source_id"] = fin["source_id"].astype("Int64")
    lab["source_id"] = lab["source_id"].astype("Int64")
    fin = fin.merge(lab, on="source_id", how="left", validate="one_to_one")
    w = fin[fin["facility"].notna()].reset_index(drop=True)
    w_name = fin[fin["facility"].notna() | fin["name_wasp"]].reset_index(drop=True)
    return w, w_name


def run(cfg, w, w_name, e, fld, fld_e, frozen, n_main, n_sens, diag=False):
    cs = (lambda x: "") if diag else (lambda x: f" case {x}")         # no case is read for the diagnostic
    R = {"n": {"wasp": len(w), "wasp_by_name": len(w_name), "eblm": len(e)}}
    ref_e = eblm_ref(e, fld_e, cfg, "eblm_primary", n_main)
    assert abs(ref_e - frozen["C_secondary"]["eblm_control_median"]) < 1e-9, ref_e
    # B
    B = vs_field(w, fld, cfg, "wasp_a1_primary", n_main, n_main, ["z"])
    R["B_primary"] = B
    log(f"  B: n={B['n']} {B['B']['delta_hat']:.2f} {np.round(B['B']['ci95'], 2).tolist()}{cs(B['case'])}")
    # C
    R["C_secondary"] = wk.contrast(w, e, B["control_median"], ref_e, cfg, n_main, "wasp_a1_C")
    R["C_secondary"]["eblm_control_median"] = ref_e
    log(f"  C: {R['C_secondary']['contrast']:.2f} {np.round(R['C_secondary']['ci95'], 2).tolist()}")
    # D (5,000 + 5,000); EBLM references from the frozen run
    D, ref_w = {}, B["control_median"]
    b = vs_field(w_name, fld, cfg, "wasp_a1_D1", n_sens, n_sens, ["z"])
    D["D1_membership_by_name"] = {"B": b, "C": wk.contrast(w_name, e, b["control_median"], ref_e, cfg, n_sens,
                                                            "wasp_a1_D1C")}
    b = vs_field(w, fld, cfg, "wasp_a1_D2", n_sens, n_sens, ZD)
    D["D2_wasp_z_dist"] = {"B": b, "C": wk.contrast(w, e, b["control_median"], ref_e, cfg, n_sens, "wasp_a1_D2C")}
    w3, e3 = w[w["dist_pc"] < 600].reset_index(drop=True), e[e["dist_pc"] < 600].reset_index(drop=True)
    b = vs_field(w3, fld, cfg, "wasp_a1_D3", n_sens, n_sens, ["z"])
    D["D3_d_lt_600"] = {"B": b, "C": wk.contrast(w3, e3, b["control_median"],
                                                  eblm_ref(e3, fld_e, cfg, "eblm_wD3", n_sens), cfg, n_sens,
                                                  "wasp_a1_D3C")}
    w4 = w[w["bp_rp0"].between(0.5, 1.5)].reset_index(drop=True)
    e4 = e[e["bp_rp0"].between(0.5, 1.5)].reset_index(drop=True)
    b = vs_field(w4, fld, cfg, "wasp_a1_D4", n_sens, n_sens, ["z"])
    D["D4_colour_0.5_1.5"] = {"B": b, "C": wk.contrast(w4, e4, b["control_median"],
                                                        eblm_ref(e4, fld_e, cfg, "eblm_wD4", n_sens), cfg, n_sens,
                                                        "wasp_a1_D4C")}
    vt = lambda f: np.sqrt(f["U"] ** 2 + f["V"] ** 2 + f["W"] ** 2)
    f5, f5e = fld[vt(fld) <= 180].reset_index(drop=True), fld_e[vt(fld_e) <= 180].reset_index(drop=True)
    b = vs_field(w, f5, cfg, "wasp_a1_D5", n_sens, n_sens, ["z"])
    D["D5_halo_cut"] = {"B": b, "C": wk.contrast(w, e, b["control_median"],
                                                  eblm_ref(e, f5e, cfg, "eblm_wD5", n_sens), cfg, n_sens,
                                                  "wasp_a1_D5C")}
    e6 = e[e["short"] != "J1500-33"].reset_index(drop=True)
    D["D6_no_J1500-33_C_only"] = {"C": wk.contrast(w, e6, ref_w, eblm_ref(e6, fld_e, cfg, "eblm_wD6", n_sens),
                                                    cfg, n_sens, "wasp_a1_D6C")}
    for k, v in D.items():
        if "B" in v:
            log(f"  {k} B: n={v['B']['n']} {v['B']['B']['delta_hat']:.2f} "
                f"{np.round(v['B']['B']['ci95'], 2).tolist()}{cs(v['B']['case'])}")
        log(f"  {k} C: {v['C']['contrast']:.2f} {np.round(v['C']['ci95'], 2).tolist()}")
    R["D_sensitivity"] = D
    # E(i), E(ii); EBLM sigma_W reference and E(iii) from the frozen run
    fe = frozen["E_exploratory_sigma_w"]
    sw_e_med = fe["iii_eblm"]["control_median"]
    sw_w = wk.mc_sigma_w(w, fld, cfg, "wasp_a1_E_W", n_main, ["z"])
    dW = wk.boot_draws(w, cfg, n_main, "wasp_a1_Ei", wk.sigma_w) - sw_w["control_median"]
    E = {"label": "EXPLORATORY (pre-registered section E); a null in sigma_W is uninformative"}
    E["i_wasp_vs_field"] = {**sw_w, "delta": sw_w["host_sigma_w"] - sw_w["control_median"],
                            "ci95": [float(x) for x in np.percentile(dW, [2.5, 97.5])],
                            "hs19_sized_delta": -(1 - wk.RATIO_HS19) * sw_w["control_median"]}
    E["ii_contrast"] = wk.contrast(w, e, sw_w["control_median"], sw_e_med, cfg, n_main, "wasp_a1_Eii", wk.sigma_w)
    E["ii_contrast"]["eblm_sigma_w_control_median_from_frozen_run"] = sw_e_med
    for k in ("i_wasp_vs_field", "ii_contrast"):
        lo, hi = E[k]["ci95"]
        E[k]["ci_excludes_0"] = bool(lo > 0 or hi < 0)
    E["iii_eblm"] = {**fe["iii_eblm"], "source": "frozen run (E31), reported unchanged (amendment section 4)"}
    R["E_exploratory_sigma_w"] = E
    log(f"  E: i {E['i_wasp_vs_field']['delta']:.2f} {np.round(E['i_wasp_vs_field']['ci95'], 2).tolist()}; "
        f"ii {E['ii_contrast']['contrast']:.2f} {np.round(E['ii_contrast']['ci95'], 2).tolist()}")
    return R


def _strip_readings(d):
    """DIAGNOSTIC: no decision rule or reading applies (amendment section 4)."""
    if isinstance(d, dict):
        return {k: _strip_readings(v) for k, v in d.items()
                if k not in ("case", "case_text", "reading", "interpretation")}
    return d


def main():
    amend = seeds.require_frozen_amendment()                      # guard; refuses unless frozen and on origin
    got = hashlib.sha256(amend.read_bytes()).hexdigest()
    assert amend == AMEND and got == AMEND_SHA, (amend, got)
    assert hashlib.sha256(wk.PREREG.read_bytes()).hexdigest() == wk.SHA
    cfg = config.load()
    n_main, n_sens = cfg["monte_carlo"]["n_iter"], cfg["monte_carlo"]["sensitivity"]["n_iter"]
    assert n_main == 40000 and n_sens == 5000, (n_main, n_sens)
    frozen = json.loads(FROZEN.read_text())
    assert frozen["sha256"] == wk.SHA
    OUT.mkdir(parents=True, exist_ok=True)
    base = config.out_dir(cfg)
    fld = pd.read_csv(base / "field_ms.csv.gz")
    e = pd.read_csv(ROOT / "data" / "private" / "eblm_d2" / "eblm_final_uvw.csv")
    ids = set(pd.read_csv(ROOT / "data" / "eblm" / "targets_eblm4.csv")["source_id"].dropna().astype(np.int64))
    fld_e = fld[~fld["source_id"].isin(ids)].reset_index(drop=True)
    assert len(e) == 79 and len(fld_e) == 376786 and len(fld) == 376791

    head = {"amendment": AMEND.name, "amendment_sha256": got, "prereg": wk.PREREG.name, "prereg_sha256": wk.SHA,
            "delta_hs19": wk.DELTA_HS19, "ratio_hs19_sigma_w": wk.RATIO_HS19,
            "reused_from_frozen_run": "E(iii); EBLM control medians eblm_primary, eblm_wD3-wD6 (cached frozen "
                                      "chunks); wasp_E_E sigma_W control median (results.json)"}
    w, w_name = samples(cfg, "source_id")
    assert len(w) == 107 and len(w_name) == 109, (len(w), len(w_name))
    log(f"HEADLINE: WASP {len(w)} (name {len(w_name)})  EBLM {len(e)}")
    R = {**head, "run": "headline (corrected membership, source_id join)",
         **run(cfg, w, w_name, e, fld, fld_e, frozen, n_main, n_sens)}
    (OUT / "results.json").write_text(json.dumps(R, indent=1, default=pipeline._json))

    w, w_name = samples(cfg, "hostname_legacy")
    assert len(w) == 98 and len(w_name) == 100, (len(w), len(w_name))
    log(f"DIAGNOSTIC: WASP {len(w)} (name {len(w_name)})  EBLM {len(e)}")
    Rd = {**head, "run": "DIAGNOSTIC: original 98-host list under the wasp_a1_* labels; not a result; "
                         "no decision rule or reading applies (amendment section 4)",
          **_strip_readings(run(cfg, w, w_name, e, fld, fld_e, frozen, n_main, n_sens, diag=True))}
    (OUT / "diagnostic_98").mkdir(exist_ok=True)
    (OUT / "diagnostic_98" / "results.json").write_text(json.dumps(Rd, indent=1, default=pipeline._json))
    log("results written")


if __name__ == "__main__":
    main()
