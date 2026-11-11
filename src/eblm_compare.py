"""Stage C (no RVs): EBLM IV primaries vs HJ hosts, and the expected EBLM precision.

    python -m src.eblm_compare

Writes data/processed/eblm/{compare_stats.csv, precision.json, eblm4_sample.csv} and
figures/eblm/{fig_eblm_vs_hj_distributions.png, fig_eblm_vs_hj_sky.png}.
Uses positions, photometry and astrometry only; no RV of any EBLM is read.
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

from . import config, kinematics, matching, pipeline, seeds, targets
from .config import ROOT

OUT = ROOT / "data" / "processed" / "eblm"
FIG = ROOT / "figures" / "eblm"
VARS = [("dist_pc", "Distance [pc]"), ("phot_g_mean_mag", "G [mag]"), ("bp_rp0", "(BP-RP)$_0$ [mag]"),
        ("z", "z [pc]"), ("abs_z", "|z| [pc]"), ("dec", "Dec [deg]"), ("b", "b [deg]")]


def samples(cfg):
    e_final, e_log, e = pipeline.target_sample(cfg, targets.default_path("eblm4"), "eblm4")
    e = e[e["on_ms"]].reset_index(drop=True)                   # after cuts 1-4 and the MS cut; gamma pending
    out = config.out_dir(cfg)
    hj = pd.read_csv(out / "hosts_hj_table1_final.csv")        # Phase 1 primary HJ sample (297)
    for d in (e, hj):
        d["abs_z"] = d["z"].abs()
    return e, hj, e_log


def compare(e, hj):
    rows = []
    for v, _ in VARS:
        a, b = e[v].dropna().to_numpy(), hj[v].dropna().to_numpy()
        ks = stats.ks_2samp(a, b)
        ad = stats.anderson_ksamp([a, b])
        rows.append({"variable": v, "eblm_n": len(a), "eblm_median": np.median(a),
                     "eblm_p16": np.percentile(a, 16), "eblm_p84": np.percentile(a, 84),
                     "hj_n": len(b), "hj_median": np.median(b), "hj_p16": np.percentile(b, 16),
                     "hj_p84": np.percentile(b, 84), "ks_D": ks.statistic, "ks_p": ks.pvalue,
                     "ad_p_capped": ad.significance_level})
    df = pd.DataFrame(rows)
    df.loc[len(df)] = {"variable": "frac_south(dec<0)", "eblm_n": len(e), "eblm_median": float((e["dec"] < 0).mean()),
                       "hj_n": len(hj), "hj_median": float((hj["dec"] < 0).mean())}
    return df


def figures(e, hj):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .figures import ACCENT, HOST, INK, MUTED
    FIG.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    for ax, (v, lab) in zip(axes.flat, VARS[:5] + [("b", "b [deg]")]):
        a, b = e[v].dropna(), hj[v].dropna()
        lo, hi = np.nanpercentile(np.r_[a, b], [0.5, 99.5])
        bins = np.linspace(lo, hi, 25)
        ax.hist(b, bins=bins, density=True, color=HOST, alpha=0.55, label=f"HJ hosts (N={len(b)})")
        ax.hist(a, bins=bins, density=True, histtype="step", lw=2, color=ACCENT, label=f"EBLM primaries (N={len(a)})")
        ax.axvline(np.median(b), color=HOST, ls="--", lw=1)
        ax.axvline(np.median(a), color=ACCENT, ls="--", lw=1)
        ax.set(xlabel=lab, ylabel="density")
    axes.flat[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_eblm_vs_hj_distributions.png", dpi=150)
    plt.close(fig)

    fig = plt.figure(figsize=(12, 4.6))
    for k, (x, y, title) in enumerate([("ra", "dec", "Equatorial (RA, Dec)"), ("l", "b", "Galactic (l, b)")]):
        ax = fig.add_subplot(1, 2, k + 1, projection="mollweide")
        for d, col, mk, lab in ((hj, HOST, "o", "HJ hosts"), (e, ACCENT, "^", "EBLM primaries")):
            lon = np.radians(-(((d[x] + 180) % 360) - 180))    # astronomical convention: east to the left
            ax.scatter(lon, np.radians(d[y]), s=12, color=col, marker=mk, alpha=0.8, lw=0, label=lab)
        ax.set_title(title, color=INK, fontsize=10, pad=12)
        ax.set_xticklabels([])
        ax.grid(color=MUTED, lw=0.4, alpha=0.5)
    fig.axes[0].legend(frameon=False, fontsize=8, loc="lower left", bbox_to_anchor=(0, -0.15))
    fig.savefig(FIG / "fig_eblm_vs_hj_sky.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def control_width(e, fld, cfg, match_vars, n_iter=5000, seed_label="eblm_precision"):
    """Control-dispersion distribution for a host sample with e's colours, z (and distances),
    using field velocities only (host velocities are not needed and not used)."""
    h = e.copy()
    h["U"] = h["V"] = h["W"] = 0.0                              # placeholders; only controls are used
    s = matching.ControlSampler(h, fld, cfg, match_vars=match_vars)
    rng = np.random.default_rng(np.random.SeedSequence([cfg["seed"], seeds.label_code(seed_label)]))
    st = cfg["monte_carlo"]["statistic"]
    d = np.array([kinematics.dispersion(s.fU[p[p >= 0]], s.fV[p[p >= 0]], s.fW[p[p >= 0]], st)
                  for p in (s.control(rng) for _ in range(n_iter))])
    return {"match_vars": match_vars, "n_iter": n_iter, "median": float(np.median(d)),
            "sd": float(np.std(d, ddof=1)), "p16": float(np.percentile(d, 16)), "p84": float(np.percentile(d, 84))}


def astrometric_error(d, cfg, rv_col, n=200, seed=1):
    """Per-star 3D velocity error [km/s] from DR2 parallax/proper-motion errors (MC, BJ18-like
    distance = 1/parallax draw scaled to the adopted distance). RV error not included here."""
    rng = np.random.default_rng(seed)
    errs = []
    for _, r in d.iterrows():
        plx = r["parallax"] + rng.normal(0, r["parallax_error"], n)
        dist = r["dist_pc"] * r["parallax"] / np.where(plx > 0, plx, np.nan)
        df = pd.DataFrame({"ra": r["ra"], "dec": r["dec"], "dist_pc": dist,
                           "pmra": r["pmra"] + rng.normal(0, r["pmra_error"], n),
                           "pmdec": r["pmdec"] + rng.normal(0, r["pmdec_error"], n),
                           "rv": r[rv_col] if np.isfinite(r[rv_col]) else 0.0}).dropna()
        U, V, W = kinematics.uvw(kinematics.skycoord(df), cfg["kinematics"]["solar_motion_kms"])
        errs.append(float(np.sqrt(np.var(U) + np.var(V) + np.var(W))))
    return np.array(errs)


def gamma_error_estimate(e_names):
    """Rough sigma_gamma from EBLM IV's PUBLISHED summary columns only (tab:obs):
    sqrt(sigma_median^2 + sigma_add^2) / sqrt(n_rv). Not from the RVs."""
    s = pd.read_csv(OUT / "eblm4_systems.csv")
    s = s[s["name"].isin(e_names)]
    sig = np.sqrt(s["sig_median"] ** 2 + s["sig_add"].clip(lower=0) ** 2) / np.sqrt(s["n_rv"]) / 1000.0  # km/s
    return sig.to_numpy()


def inflation(sigma_true, err):
    """Inflation of an RMS-like dispersion by independent errors: sqrt(s^2 + <e^2>) - s."""
    return float(np.sqrt(sigma_true ** 2 + np.mean(np.asarray(err) ** 2)) - sigma_true)


def main():
    cfg = config.load()
    OUT.mkdir(parents=True, exist_ok=True)
    e, hj, e_log = samples(cfg)
    e.to_csv(OUT / "eblm4_sample.csv", index=False)
    st = compare(e, hj)
    st.to_csv(OUT / "compare_stats.csv", index=False)
    figures(e, hj)
    fld = pd.read_csv(config.out_dir(cfg) / "field_ms.csv.gz")
    prec = {"n_eblm_after_cuts_ms": len(e), "n_hj": len(hj), "eblm_cuts": e_log["cuts"],
            "eblm_fgk_after_ms": pipeline.counts_fgk(e)}
    prec["control_z_colour"] = control_width(e, fld, cfg, ["z"])
    prec["control_z_dist_colour"] = control_width(e, fld, cfg, ["z", "dist_pc"])
    prec["control_hj_z_colour"] = control_width(hj, fld, cfg, ["z"], seed_label="hj_precision")
    # bootstrap SD of the HJ hosts' own dispersion (for the direct EBLM-vs-HJ comparison)
    rng = np.random.default_rng(np.random.SeedSequence([cfg["seed"], 7]))
    boot = [kinematics.dispersion(*(hj[c].to_numpy()[i] for c in "UVW"), cfg["monte_carlo"]["statistic"])
            for i in (rng.integers(0, len(hj), len(hj)) for _ in range(5000))]
    prec["hj_bootstrap_sd"] = float(np.std(boot, ddof=1))
    sd_e = prec["control_z_dist_colour"]["sd"]
    prec["mdd_vs_field_2sigma_80pct"] = (1.96 + 0.84) * sd_e     # minimum detectable difference, two-sided 5%, 80% power
    sd_diff = float(np.hypot(sd_e, prec["hj_bootstrap_sd"]))
    prec["sd_diff_eblm_minus_hj"] = sd_diff
    prec["mdd_vs_hj_2sigma_80pct"] = (1.96 + 0.84) * sd_diff
    prec["hj_minus_field_observed"] = 42.82 - 36.18             # Phase 1 primary, km/s
    # error budgets (km/s): astrometric (both samples), gamma (EBLM, from published summaries)
    ea = astrometric_error(e, cfg, "rv_used")
    ha = astrometric_error(hj, cfg, "rv_used")
    eg = gamma_error_estimate(e["name"].str.replace("EBLM ", "", regex=False))
    sig = prec["control_z_dist_colour"]["median"]
    prec["errors"] = {
        "eblm_astrometric_median": float(np.median(ea)), "eblm_astrometric_rms": float(np.sqrt(np.mean(ea ** 2))),
        "hj_astrometric_median": float(np.median(ha)), "hj_astrometric_rms": float(np.sqrt(np.mean(ha ** 2))),
        "eblm_gamma_est_median": float(np.median(eg)), "eblm_gamma_est_max": float(np.max(eg)),
        "inflation_eblm_astrometric": inflation(sig, ea), "inflation_hj_astrometric": inflation(sig, ha),
        "inflation_eblm_gamma": inflation(sig, eg),
        "note": "inflation = sqrt(sigma^2 + <err^2>) - sigma at sigma = EBLM control median; RMS approximation"}
    (OUT / "precision.json").write_text(json.dumps(prec, indent=1, default=float))
    pd.set_option("display.width", 200)
    print(st.round(3).to_string(index=False))
    print(json.dumps({k: v for k, v in prec.items() if k != "eblm_cuts"}, indent=1, default=float))


if __name__ == "__main__":
    main()
