"""Stage D2: the pre-registered EBLM kinematics (PREREG_eblm_kinematics_2026-11-02.txt, frozen).

    python -m src.eblm_kinematics              # B, C-1, C-2, G1-G11, I; resumable Monte Carlo chunks
    python -m src.eblm_kinematics --figures    # Figs. 3 and 4 from the saved results

Section letters follow the pre-registration. G11 (without J1500-33) was added by the user after
STOP D1 and before any kinematics (DECISIONS D8); it is labelled as such everywhere.
Per-system tables (gamma-derived UVW) stay in data/private/eblm_d2/ (C7); aggregates go to
data/processed/eblm/d2/results.json.
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

from . import config, figures, kinematics, main_sequence, matching, pipeline
from .config import ROOT
from .gaia_tap import log

PRIV = ROOT / "data" / "private" / "eblm_d2"
GAMMA = ROOT / "data" / "private" / "eblm_gamma" / "gamma_fits.csv"
OUT = ROOT / "data" / "processed" / "eblm" / "d2"
FIG = ROOT / "figures" / "eblm"
N_SENS = None            # from config monte_carlo.sensitivity.n_iter


# ------------------------------------------------------------------ samples
def target_list(kind):
    """Private target CSV with gamma (primary fit) or the simple average (G10) as the RV."""
    t = pd.read_csv(ROOT / "data" / "eblm" / "targets_eblm4.csv", dtype={"source_id": "Int64"})
    g = pd.read_csv(GAMMA)
    g = g[g["variant"] == "primary"].set_index("name")
    key = t["name"].str.replace("EBLM ", "", regex=False)
    ok = g.loc[key, "flag"].to_numpy() == "ok"
    gamma = g.loc[key, "gamma"].to_numpy()
    rv = gamma if kind == "gamma" else gamma - g.loc[key, "gamma_minus_mean"].to_numpy()
    t["rv"] = np.where(ok, rv, np.nan)
    t["rv_err"] = np.where(ok, g.loc[key, "sigma_gamma"].to_numpy(), np.nan)
    t["rv_kind"] = "gamma"            # both are systemic values from the amended fit (G10 = method check)
    t["rv_ref"] = "07 amendment fit, drift t_ref mid" if kind == "gamma" else "simple average (G10)"
    PRIV.mkdir(parents=True, exist_ok=True)
    p = PRIV / f"targets_eblm4_{kind}.csv"
    t.to_csv(p, index=False)
    return p


def eblm_samples(cfg, kind="gamma"):
    final, lg, after_cuts = pipeline.target_sample(cfg, target_list(kind), "eblm4")   # same Gaia rows as Stage C
    sysf = pd.read_csv(ROOT / "data" / "processed" / "eblm" / "eblm4_systems.csv").set_index("name")
    for d in (final, after_cuts):
        d["short"] = d["name"].str.replace("EBLM ", "", regex=False)
        d["eblm4_flag"] = sysf.loc[d["short"], "flag"].fillna("").to_numpy()
        d["abs_z"] = d["z"].abs()
    return final, lg, after_cuts


def field(cfg, eblm_ids):
    """Phase 1 field (same cuts and exclusions), with the EBLM IV primaries also removed (D9)."""
    out = config.out_dir(cfg)
    excl = set(int(i) for i in eblm_ids)
    for k in ("hj_table1", "hj", "longp"):
        excl |= set(pd.read_csv(out / f"hosts_{k}_after_cuts.csv")["source_id"].dropna().astype(np.int64))
    fld, lg, _ = pipeline.field_sample(cfg, excl)
    fld["abs_z"] = fld["z"].abs()
    return fld, lg


# ------------------------------------------------------------------ statistics
def _disp_boot(U, V, W, idx, stat):
    return kinematics.dispersion(U[idx], V[idx], W[idx], stat)


def bootstrap(e, ref, cfg, n, label, stat):
    """B: percentile CI of sigma(EBLM resample) - ref; bias; BCa when |bias| > 0.2 half-width."""
    U, V, W = (e[c].to_numpy(float) for c in "UVW")
    N = len(e)
    rng = np.random.default_rng(np.random.SeedSequence([cfg["seed"], sum(map(ord, label)), 99]))
    d = np.concatenate([_disp_boot(U, V, W, rng.integers(0, N, (min(5000, n - s), N)), stat)
                        for s in range(0, n, 5000)]) - ref
    est = float(kinematics.dispersion(U, V, W, stat) - ref)
    lo, hi = np.percentile(d, [2.5, 97.5])
    bias = float(d.mean() - est)
    out = {"delta_hat": est, "ci95": [float(lo), float(hi)], "boot_n": n, "boot_bias": bias,
           "half_width": float((hi - lo) / 2), "boot_sd": float(np.std(d, ddof=1))}
    out["bca_needed"] = bool(abs(bias) > 0.2 * out["half_width"])
    if out["bca_needed"]:
        z0 = stats.norm.ppf(np.mean(d < est))
        jack = np.array([kinematics.dispersion(np.delete(U, i), np.delete(V, i), np.delete(W, i), stat)
                         for i in range(N)])
        dj = jack.mean() - jack
        a = float(np.sum(dj ** 3) / (6 * np.sum(dj ** 2) ** 1.5))
        za = stats.norm.ppf([0.025, 0.975])
        q = stats.norm.cdf(z0 + (z0 + za) / (1 - a * (z0 + za)))
        out["bca_ci95"] = [float(x) for x in np.percentile(d, 100 * q)]
        out["bca_z0_a"] = [float(z0), a]
    lo, hi = out["ci95"]
    out["interpretation"] = ("colder than the field" if hi < 0 else "hotter than the field" if lo > 0
                             else "consistent with the field")
    return out


def mc(e, fld, cfg, label, n_iter, match_vars, stat=None, ckpt=True):
    """HS19 Monte Carlo (pipeline.monte_carlo) with the given step-(a) matching variables."""
    c = config.variant(cfg, **{"monte_carlo.match_vars": list(match_vars)})
    if stat:
        c["monte_carlo"]["statistic"] = stat
    s, res = pipeline.monte_carlo(e, fld, c, label, n_iter, stat,
                                  ckpt_dir=PRIV / "mc_checkpoints" if ckpt else None)
    s.pop("unmatched_hosts", None)                # names + per-host counts: kept out of the aggregates
    s["n_hosts_ever_unmatched"] = int((np.asarray(res["host_unmatched_count"]) > 0).sum())
    return s, res


def primary_block(e, fld, cfg, label, n_iter, n_boot, match_vars, stat=None):
    stat = stat or cfg["monte_carlo"]["statistic"]
    s, res = mc(e, fld, cfg, label, n_iter, match_vars, stat)
    b = bootstrap(e, s["control_median"], cfg, n_boot, label, stat)
    return {"n": len(e), "match_vars": list(match_vars), "statistic": stat, "sigma_eblm": s["host_dispersion"],
            "control_median": s["control_median"], "control_sd": s["control_std"],
            "B": b, "C1": {"fraction_controls_le": s["p_value"], "n_le": s["n_le"], "n_iter": n_iter,
                           "z_gauss": s["z_gauss"]}, "mc": s}, res


def c2(e, hj, cfg, n_iter, label, tol=0.05, stat=None):
    """C-2: EBLM bootstrap vs colour (|d| <= tol) and EBLM-distance-quintile matched HJ hosts."""
    stat = stat or cfg["monte_carlo"]["statistic"]
    edges = matching.bin_edges(e["dist_pc"], {"bins": "quantile", "n_quantile_bins": 5})
    eb = np.searchsorted(edges, e["dist_pc"].to_numpy(), side="right") - 1
    hb = np.searchsorted(edges, hj["dist_pc"].to_numpy(), side="right") - 1
    hb = np.where((hb >= 0) & (hb < len(edges) - 1), hb, -1)
    ec, hc = e["bp_rp0"].to_numpy(), hj["bp_rp0"].to_numpy()
    cands = [np.flatnonzero((np.abs(hc - ec[i]) <= tol) & (hb == eb[i])) for i in range(len(e))]
    nc = np.array([len(c) for c in cands])
    off = np.r_[0, np.cumsum(nc)[:-1]]
    flat = np.concatenate(cands) if nc.sum() else np.array([], int)
    EU, EV, EW = (e[c].to_numpy(float) for c in "UVW")
    HU, HV, HW = (hj[c].to_numpy(float) for c in "UVW")
    rng = np.random.default_rng(np.random.SeedSequence([cfg["seed"], sum(map(ord, label)), 2]))
    out = np.empty(n_iter)
    n_drop = np.empty(n_iter, int)
    for k in range(n_iter):
        idx = rng.integers(0, len(e), len(e))
        idx = idx[nc[idx] > 0]
        n_drop[k] = len(e) - len(idx)
        pick = flat[off[idx] + (rng.random(len(idx)) * nc[idx]).astype(int)]
        out[k] = (kinematics.dispersion(EU[idx], EV[idx], EW[idx], stat)
                  - kinematics.dispersion(HU[pick], HV[pick], HW[pick], stat))
    return {"n_eblm": len(e), "n_hj": len(hj), "n_eblm_without_hj_candidate": int((nc == 0).sum()),
            "median_hj_candidates": float(np.median(nc)), "mean_dropped_per_iter": float(n_drop.mean()),
            "delta_median": float(np.median(out)), "ci95": [float(x) for x in np.percentile(out, [2.5, 97.5])],
            "n_iter": n_iter, "tol": tol, "distance_edges": [float(x) for x in edges]}


# ------------------------------------------------------------------ I (exploratory)
def evolved_fraction(E, H, cfg, bw, n_boot, label):
    """I: colour-reweighted evolved fraction, EBLM minus HJ hosts, with a joint bootstrap."""
    lo, hi = E["bp_rp0"].min(), E["bp_rp0"].max()
    H = H[(H["bp_rp0"] >= lo) & (H["bp_rp0"] <= hi)]
    edges = np.arange(lo, hi + bw, bw)
    if edges[-1] <= hi:
        edges = np.append(edges, edges[-1] + bw)
    ebin = np.clip(np.searchsorted(edges, E["bp_rp0"].to_numpy(), "right") - 1, 0, len(edges) - 2)
    hbin = np.clip(np.searchsorted(edges, H["bp_rp0"].to_numpy(), "right") - 1, 0, len(edges) - 2)
    ev_e, ev_h = E["evolved"].to_numpy(bool), H["evolved"].to_numpy(bool)
    nb = len(edges) - 1

    def one(ie, ih):
        ne = np.bincount(ebin[ie], minlength=nb)
        nh = np.bincount(hbin[ih], minlength=nb)
        keep = (ne > 0) & (nh > 0)
        me = keep[ebin[ie]]
        fe = ev_e[ie][me].mean()
        w = np.where(keep, ne / np.maximum(nh, 1), 0)[hbin[ih]]
        fh = np.sum(w * ev_h[ih]) / np.sum(w)
        return fe, fh, int((~me).sum()), int((keep == 0)[np.unique(ebin[ie])].sum())

    fe, fh, n_e_dropped, n_bins_dropped = one(np.arange(len(E)), np.arange(len(H)))
    rng = np.random.default_rng(np.random.SeedSequence([cfg["seed"], sum(map(ord, label)), 3]))
    bd = np.array([np.subtract(*one(rng.integers(0, len(E), len(E)), rng.integers(0, len(H), len(H)))[:2])
                   for _ in range(n_boot)])
    return {"bin_width": bw, "n_eblm": len(E), "n_hj_in_colour_range": len(H),
            "raw_f_eblm": float(ev_e.mean()), "raw_f_hj_in_range": float(ev_h.mean()),
            "f_eblm": float(fe), "f_hj_weighted": float(fh), "delta_f": float(fe - fh),
            "ci95": [float(x) for x in np.percentile(bd, [2.5, 97.5])], "n_boot": n_boot,
            "eblm_dropped_no_hj_bin": n_e_dropped, "colour_bins_dropped": n_bins_dropped,
            "colour_range": [float(lo), float(hi)]}


def evolved_inputs(cfg, e_after_cuts, delta_mag=None, hj_table=None):
    c = config.variant(cfg, **({"main_sequence.delta_mag": delta_mag} if delta_mag else {}))
    E = e_after_cuts[e_after_cuts["dist_pc"].notna() & np.isfinite(e_after_cuts["bp_rp0"])
                     & np.isfinite(e_after_cuts["M_G0"])].copy()
    H = (hj_table if hj_table is not None
         else pd.read_csv(config.out_dir(cfg) / "hosts_hj_after_cuts.csv"))
    H = H[H["dist_pc"].notna() & np.isfinite(H["bp_rp0"]) & np.isfinite(H["M_G0"])].copy()
    for d in (E, H):
        d["evolved"] = ~main_sequence.on_main_sequence(d["M_G0"], d["bp_rp0"], c)
    return E, H


# ------------------------------------------------------------------ run
def main():
    cfg = config.load()
    OUT.mkdir(parents=True, exist_ok=True)
    PRIV.mkdir(parents=True, exist_ok=True)
    mcc = cfg["monte_carlo"]
    n_main, n_sens = mcc["n_iter"], mcc["sensitivity"]["n_iter"]
    R = {"prereg": "PREREG_eblm_kinematics_2026-11-02.txt", "statistic": mcc["statistic"]}

    e, e_log, e_cut = eblm_samples(cfg, "gamma")
    assert len(e) == 79, len(e)
    R["sample"] = {"n_final": len(e), "fgk": pipeline.counts_fgk(e), "log": {k: v for k, v in e_log.items()
                                                                           if k != "cuts"}}
    e.to_csv(PRIV / "eblm_final_uvw.csv", index=False)
    fld, f_log = field(cfg, pd.read_csv(ROOT / "data" / "eblm" / "targets_eblm4.csv")["source_id"].dropna())
    R["field"] = {"n_ms_with_uvw": len(fld), "n_after_cuts": f_log["n_after_cuts"]}
    hj = pd.read_csv(config.out_dir(cfg) / "hosts_hj_table1_final.csv")
    log(f"EBLM {len(e)}  field {len(fld)}  HJ {len(hj)}")

    # B + C-1 (primary: z + distance)
    R["primary"], res = primary_block(e, fld, cfg, "eblm_primary", n_main, n_main, ["z", "dist_pc"])
    np.savez(OUT / "mc_eblm_primary.npz", **{f"disp_{k}": v for k, v in res["dispersion_all"].items()},
             host_disp=res["host_dispersion"])
    log(f"PRIMARY: {json.dumps(R['primary']['B'])}")
    # C-2
    R["C2"] = c2(e, hj, cfg, n_main, "eblm_c2")
    log(f"C-2: {R['C2']['delta_median']:.2f} {R['C2']['ci95']}")

    # G sensitivity (5,000 iterations; bootstrap 5,000)
    G = {}
    zd = ["z", "dist_pc"]
    G["G1_z_only"] = primary_block(e, fld, cfg, "eblm_G1", n_sens, n_sens, ["z"])[0]
    cq = config.variant(cfg, **{"monte_carlo.z_match.bins": "quantile", "monte_carlo.z_match.n_quantile_bins": 5})
    G["G2_z_dist_quantile5x5"] = primary_block(e, fld, cq, "eblm_G2", n_sens, n_sens, zd)[0]
    G["G3_absz_dist"] = primary_block(e, fld, cfg, "eblm_G3", n_sens, n_sens, ["abs_z", "dist_pc"])[0]
    G["G4_absz_only"] = primary_block(e, fld, cfg, "eblm_G4", n_sens, n_sens, ["abs_z"])[0]
    G["G5_field_dec_lt_0"] = primary_block(e, fld[fld["dec"] < 0].reset_index(drop=True), cfg, "eblm_G5",
                                           n_sens, n_sens, zd)[0]
    sub = e[~e["eblm4_flag"].isin(["drift", "triple"])].reset_index(drop=True)
    G["G6_no_drift_triple"] = primary_block(sub, fld, cfg, "eblm_G6", n_sens, n_sens, zd)[0]
    sub = e[e["eblm4_flag"] != "active"].reset_index(drop=True)
    G["G7_no_active"] = primary_block(sub, fld, cfg, "eblm_G7", n_sens, n_sens, zd)[0]
    G["G8_rms"] = primary_block(e, fld, cfg, "eblm_G8", n_sens, n_sens, zd, stat="rms")[0]
    eF = e[e["sptype"] == "F"].reset_index(drop=True)
    G["G9_F_only"] = primary_block(eF, fld, cfg, "eblm_G9", n_sens, n_sens, zd)[0]
    G["G9_F_only"]["C2_vs_F_HJ"] = c2(eF, hj[hj["sptype"] == "F"].reset_index(drop=True), cfg, n_sens, "eblm_G9c2")
    e_avg, _, _ = eblm_samples(cfg, "average")
    G["G10_simple_average"] = primary_block(e_avg, fld, cfg, "eblm_G10", n_sens, n_sens, zd)[0]
    sub = e[e["short"] != "J1500-33"].reset_index(drop=True)
    G["G11_no_J1500-33_POST_D1"] = primary_block(sub, fld, cfg, "eblm_G11", n_sens, n_sens, zd)[0]
    for k, v in G.items():
        log(f"{k}: n={v['n']} dhat={v['B']['delta_hat']:.2f} CI={np.round(v['B']['ci95'], 2).tolist()} "
            f"C1={v['C1']['fraction_controls_le']:.4f}")
    R["G"] = G

    # I (exploratory, pre-registered)
    I = {}
    E, H = evolved_inputs(cfg, e_cut)
    I["primary_bw0.1"] = evolved_fraction(E, H, cfg, 0.1, 10000, "I_main")
    I["bw0.05"] = evolved_fraction(E, H, cfg, 0.05, 10000, "I_bw005")
    I["bw0.2"] = evolved_fraction(E, H, cfg, 0.2, 10000, "I_bw02")
    E5, H5 = evolved_inputs(cfg, e_cut, delta_mag=0.5)
    I["threshold_0.5mag"] = evolved_fraction(E5, H5, cfg, 0.1, 10000, "I_d05")
    try:
        cn = config.variant(cfg, **{"hosts.disc_year_max": None})
        _, _, h_now = pipeline.host_sample(cn, "hj")
        _, Hn = evolved_inputs(cfg, e_cut, hj_table=h_now)
        I["hj_current_archive"] = evolved_fraction(E, Hn, cfg, 0.1, 10000, "I_now")
    except Exception as ex:                       # network; recorded, not silently skipped
        I["hj_current_archive"] = {"error": repr(ex)}
    R["I_exploratory"] = I
    (OUT / "results.json").write_text(json.dumps(R, indent=1, default=pipeline._json))
    log("results written")
    make_figures(cfg, e, fld, hj)


def make_figures(cfg, e=None, fld=None, hj=None):
    FIG.mkdir(parents=True, exist_ok=True)
    R = json.loads((OUT / "results.json").read_text())
    if e is None:
        e = pd.read_csv(PRIV / "eblm_final_uvw.csv")
        fld, _ = field(cfg, pd.read_csv(ROOT / "data" / "eblm" / "targets_eblm4.csv")["source_id"].dropna())
    figures.fig_uvw(fld, e, FIG / "fig3_uvw_eblm.png", label="EBLM primaries", points=False)
    st = cfg["monte_carlo"]["statistic"]
    d = np.load(OUT / "mc_eblm_primary.npz")[f"disp_{st}"]
    s = dict(R["primary"]["mc"])
    hjs = json.loads((config.out_dir(cfg) / "summary.json").read_text())["mc"]["hj_table1"]
    dh = np.load(config.out_dir(cfg) / "mc_hj_table1.npz")[f"disp_{st}"]
    figures.fig_mc({"eblm": (s, d), "hj_table1": (hjs, dh)}, FIG / "fig4_mc_eblm.png",
                   {"eblm": "EBLM primaries", "hj_table1": "HJ hosts (HS19 Table 1)"}, shared_x=True)


if __name__ == "__main__":
    if "--figures" in sys.argv:
        make_figures(config.load())
    else:
        main()
