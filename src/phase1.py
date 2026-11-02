"""Phase 1 driver: faithful DR2 run, sensitivity variants, figures.

    python -m src.phase1            # everything
    python -m src.phase1 --primary  # primary run + figures only
    python -m src.phase1 --test 1000  # short Monte Carlo (Stage A), separate output folder
    python -m src.phase1 --figures    # redraw Figs. 2-4 from saved outputs
"""
import json
import sys
import time

import numpy as np
import pandas as pd

from . import config, figures, pipeline
from .config import ROOT

LABELS = {"hj_table1": "HJ hosts (HS19 Table 1)", "hj": "HJ hosts (archive, disc. <= 2019)",
          "longp": "Long-period giant hosts"}

# Variants that change the samples (full re-selection), run at the sensitivity n_iter.
SAMPLE_VARIANTS = {
    "cut4_literal": {"cuts.uwe.form": "literal"},
    "cut8_0.0015": {"cuts.bp_rp_excess.lo_b": 0.0015},
    "longp_without_hj_hosts": {"hosts.exclude_hj_hosts_from_longp": True},
    "z_heliocentric": {"kinematics.z_definition": "heliocentric"},
    "disc_year_le_2018": {"hosts.disc_year_max": 2018},
    "host_rv_fill_dr3": {"hosts.rv_fill_dr3": True},
    "drop_mass_upper_limits": {"hosts.drop_mass_upper_limits": True},
}
# Main-sequence threshold shifts: keep M_G0 > MS - (1 + shift) (user request, Stage B)
for _d in (-0.3, -0.2, -0.1, 0.1, 0.2, 0.3):
    SAMPLE_VARIANTS[f"ms_shift_{_d:+.1f}"] = {"main_sequence.delta_mag": round(1.0 + _d, 2)}


def matching_variants(cfg):
    """Variants that only change step (a) of the Monte Carlo (samples reused)."""
    s = cfg["monte_carlo"]["sensitivity"]
    out = {}
    for w in s["bin_widths_pc"]:
        out[f"zbin_{w:g}pc"] = {"monte_carlo.z_match.bins": "fixed", "monte_carlo.z_match.bin_width_pc": w}
    for q in s["quantile_bins"]:
        out[f"zbin_q{q}"] = {"monte_carlo.z_match.bins": "quantile", "monte_carlo.z_match.n_quantile_bins": q}
    out["weighted_resample"] = {"monte_carlo.z_match.method": "weighted_resample"}
    out["no_z_matching"] = {"monte_carlo.z_match.bins": "none"}
    out["control_with_replacement"] = {"monte_carlo.control_without_replacement": False}
    return out


def row(step, kind, res, mc):
    return {"variant": step, "sample": kind, "N_host": mc["n_host"], "N_field": res,
            "host_disp": round(mc["host_dispersion"], 2),
            "ctrl_median": round(mc["control_median"], 2),
            "ctrl_16": round(mc["control_pct"][16], 2), "ctrl_84": round(mc["control_pct"][84], 2),
            "ctrl_2": round(mc["control_pct"][2], 2), "ctrl_98": round(mc["control_pct"][98], 2),
            "ctrl_std": round(mc["control_std"], 2), "z_gauss": round(mc["z_gauss"], 2),
            "p": mc["p_value"], "n_le": mc["n_le"], "n_iter": mc["n_iter"],
            "host_disp_rms": round(mc["host_dispersion_all"]["rms"], 2),
            "p_rms": mc["p_value_all"]["rms"],
            "iters_unmatched": mc["iters_with_unmatched_host"]}


def make_figures(res, cfg):
    fd = config.out_dir(cfg, "figures")
    out = config.out_dir(cfg)
    fld, hj = res["field"], res["samples"]["hj_table1"]
    figures.fig_cmd(fld, hj, cfg, fd / "fig2_cmd_hj_table1.png", "HS19 Table 1 hosts, Gaia DR2")
    if "hj" in res["samples"]:
        figures.fig_cmd(fld, res["samples"]["hj"], cfg, fd / "fig2_cmd_hj_archive.png", "Archive HJ hosts (disc. <= 2019)")
    figures.fig_uvw(fld, hj, fd / "fig3_uvw_hj_table1.png")
    mcres = {}
    for kind in ("hj_table1", "longp"):
        d = np.load(out / f"mc_{kind}.npz")[f"disp_{cfg['monte_carlo']['statistic']}"]
        mcres[kind] = (res["mc"][kind], d)
    figures.fig_mc(mcres, fd / "fig4_mc.png", LABELS)
    if "hj" in res["mc"]:
        d = np.load(out / "mc_hj.npz")[f"disp_{cfg['monte_carlo']['statistic']}"]
        figures.fig_mc({"hj": (res["mc"]["hj"], d)}, fd / "fig4_mc_hj_archive.png", LABELS)


def main(primary_only=False):
    cfg = config.load()
    out = config.out_dir(cfg)
    t0 = time.time()
    res = pipeline.run(cfg)
    make_figures(res, cfg)
    rows = [row("primary", k, len(res["field"]), m) for k, m in res["mc"].items()]
    pd.DataFrame(rows).to_csv(out / "results_primary.csv", index=False)
    print(f"primary done in {time.time() - t0:.0f} s", flush=True)
    if primary_only:
        return
    n_sens = cfg["monte_carlo"]["sensitivity"]["n_iter"]
    for name, ov in matching_variants(cfg).items():
        v = config.variant(cfg, **ov)
        for kind in ("hj_table1", "longp"):
            mc, _ = pipeline.monte_carlo(res["samples"][kind], res["field"], v, kind, n_sens)
            rows.append(row(name, kind, len(res["field"]), mc))
            print(name, kind, rows[-1]["host_disp"], rows[-1]["ctrl_median"], rows[-1]["p"], flush=True)
    # HJ hosts limited to Table 1 hosts that are also in the archive list (334 before cuts)
    arch = set(res["full"]["hj"]["source_id"].dropna().astype(np.int64)) | set(
        pipeline.hosts.build(cfg, "hj")["source_id"].dropna().astype(np.int64))
    t1 = res["samples"]["hj_table1"]
    sub = t1[t1["source_id"].astype("Int64").isin(arch)].reset_index(drop=True)
    mc, _ = pipeline.monte_carlo(sub, res["field"], cfg, "hj_table1", n_sens)
    rows.append(row("table1_in_archive", "hj_table1", len(res["field"]), mc))
    print("table1_in_archive", len(sub), rows[-1]["host_disp"], rows[-1]["p"], flush=True)
    variants = dict(SAMPLE_VARIANTS)
    refit = out / "pleiades_refit.json"
    if refit.exists():                       # MS from our own Pleiades refit instead of HS19's coefficients
        variants["ms_refit_poly"] = {"main_sequence.coeffs": json.load(open(refit))["refit_coeffs"]}
    for name, ov in variants.items():
        v = config.variant(cfg, run_name=f"{cfg['run_name']}_{name}", **ov)
        r = pipeline.run(v, n_iter=n_sens, save=True)
        for kind, mc in r["mc"].items():
            rows.append(row(name, kind, len(r["field"]), mc))
        json.dump({k: r["logs"][k] for k in r["logs"]}, open(config.out_dir(v) / "logs.json", "w"),
                  indent=1, default=pipeline._json)
    pd.DataFrame(rows).to_csv(out / "results_sensitivity.csv", index=False)
    print(f"all done in {time.time() - t0:.0f} s")


def figures_from_saved():
    """Redraw Figs. 2-4 from the saved primary outputs (no recomputation)."""
    cfg = config.load()
    out = config.out_dir(cfg)
    summ = json.load(open(out / "summary.json"))
    res = {"field": pd.read_csv(out / "field_ms.csv.gz"),
           "samples": {k: pd.read_csv(out / f"hosts_{k}_final.csv") for k in summ["mc"]},
           "mc": summ["mc"]}
    make_figures(res, cfg)


def test(n_iter, **ov):
    """Short Monte Carlo for all host samples; summary + per-iteration runtime."""
    cfg = config.load()
    tag = "_excl" if ov else ""
    v = config.variant(cfg, run_name=f"{cfg['run_name']}_test{n_iter}{tag}", **ov)
    res = pipeline.run(v, n_iter=n_iter, save=True)
    rows = []
    for kind, mc in res["mc"].items():
        r = row(f"test{n_iter}", kind, len(res["field"]), mc)
        r.update({"host_disp_iter_mean": round(mc["host_dispersion_iter_mean"], 2),
                  "mean_unmatched_per_iter": round(mc["mean_unmatched_per_iter"], 2),
                  "n_excluded": len(mc["excluded_hosts"]),
                  "seconds": mc["seconds"], "ms_per_iter": round(1000 * mc["seconds"] / n_iter, 2)})
        rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(config.out_dir(v) / "test_mc.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    if "--figures" in sys.argv:
        figures_from_saved()
    elif "--test" in sys.argv:
        ov = {"monte_carlo.exclude_unmatchable_hosts": True} if "--exclude-unmatchable" in sys.argv else {}
        test(int(sys.argv[sys.argv.index("--test") + 1]), **ov)
    else:
        main(primary_only="--primary" in sys.argv)
