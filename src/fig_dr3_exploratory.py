"""EXPLORATORY figures for the DR3 Stage 2 rows and the E62 POST-HOC DIAGNOSTIC (DECISIONS E63). No host names.

    python -m src.fig_dr3_exploratory 1 2 3 4     # -> notes/figures/dr3_exploratory/fig{1..4}_*.png

Values come from the saved results files, not retyped:
  - stage2_results.json (E61), diag_results.json (E62), data/processed/wasp_amendment1/results.json (via stage2).
  - Figs 1-2: the Method A main control draws are not stored, so they are regenerated with the same registered
    label (dr3x_A_wasp) and the same RNG consumption as src.dr3x_method_a.draws; the script asserts that every
    host value, median and 1/2-sigma percentile reproduces stage2_results.json to 1e-9 before plotting.
    sigma_U and sigma_V come from these same draws and are written to fig2_method_a_components.json (new
    exploratory calculation, E63).
  - Fig 3: the Method B matched control stars are the chunk-0 picks (5,000 iterations) replayed by
    src.dr3x_diag.replay, which asserts the replay against the saved MC checkpoint.
"""
import json
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

from . import config, dr3x_hosts, field, seeds
from . import dr3x_diag as D
from . import dr3x_method_a as A
from . import dr3x_stage2 as S
from .figures import ACCENT, HOST, INK, MUTED, SEQ
from .fig_wasp_amendment1 import DELTA_HS19
from .wasp_kinematics import DELTA_HS19 as DELTA_HS19_SRC

OUT = config.ROOT / "notes" / "figures" / "dr3_exploratory"
STAGE2 = A.OUT / "stage2_results.json"
DIAG_JSON = A.OUT / "diag_results.json"
COMP_JSON = A.OUT / "fig2_method_a_components.json"
LABEL_A = "dr3x_A_wasp"
STAT_NAMES = {"sigma_W": r"$\sigma_W$", "sigma_tot": r"$\sigma_{tot}$", "hs19_eq1": "HS19 Eq. 1"}
DR2_COL, DR3_COL = MUTED, HOST
assert round(DELTA_HS19_SRC, 2) == DELTA_HS19


def stamp(fig):
    fig.text(0.5, -0.01, A.HEADER, ha="center", va="top", fontsize=10, color=ACCENT, weight="bold")


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    stamp(fig)
    fig.savefig(OUT / name)
    plt.close(fig)
    print(f"wrote {OUT / name}")


# ----------------------------------------------------------------------------- Method A (figs 1, 2)
def method_a_draws():
    """Method A main (97): hosts, pools, field and draws incl. sigma_U, sigma_V; checked against stage2."""
    cfg = config.load()
    base, _, nea_rv = dr3x_hosts.bases(cfg)
    _, _, h = dr3x_hosts.method_a(base, nea_rv, exclude_companions=True)
    h = h.copy()
    h["rv"] = h["dr3"].map(nea_rv).fillna(h["radial_velocity"])
    h = A.add_uvw(h)
    fld, _ = A.load_field()
    pl = A.pools(h, fld)
    assert all(len(p) for p in pl)
    uvw = fld[["U", "V", "W"]].to_numpy()
    del fld
    # same RNG stream and consumption as A.draws (n = N_DRAWS, block 500)
    rng = np.random.default_rng(np.random.SeedSequence(list(seeds.stream_entropy(LABEL_A))))
    cnt = np.array([len(p) for p in pl])
    flat, off = np.concatenate(pl), np.concatenate([[0], np.cumsum(cnt)[:-1]])
    stats = dict(A.STATS, sigma_U=lambda v: np.std(v[..., 0], axis=-1), sigma_V=lambda v: np.std(v[..., 1], axis=-1))
    mc = {k: np.empty(A.N_DRAWS) for k in stats}
    for i0 in range(0, A.N_DRAWS, 500):
        k = min(500, A.N_DRAWS - i0)
        v = uvw[flat[off + (rng.random((k, len(pl))) * cnt).astype(np.int64)]]
        for s, fn in stats.items():
            mc[s][i0:i0 + k] = fn(v)
    hv = h[["U", "V", "W"]].to_numpy()
    host = {s: float(fn(hv)) for s, fn in stats.items()}
    ref = json.loads(STAGE2.read_text())["method_A"]["main_with_companion_exclusions"]
    assert ref["n_hosts"] == len(h) == 97
    for s in A.STATS:
        q = np.percentile(mc[s], [2.5, 16, 50, 84, 97.5])
        got = [host[s], q[2], q[1], q[3], q[0], q[4]]
        want = [ref[s]["host"], ref[s]["field_median"], *ref[s]["one_sigma"], *ref[s]["two_sigma"]]
        assert np.allclose(got, want, rtol=0, atol=1e-9), (s, got, want)
    # expected draw weight of each field star: sum over hosts of 1/pool size
    w = np.zeros(len(uvw))
    for p in pl:
        np.add.at(w, p, 1.0 / len(p))
    return h, mc, host, ref, uvw, w


def band(ax, lo2, lo1, hi1, hi2):
    ax.axvspan(lo2, hi2, color=HOST, alpha=0.12, lw=0, label="2σ (2.5–97.5 %)")
    ax.axvspan(lo1, hi1, color=HOST, alpha=0.25, lw=0, label="1σ (16–84 %)")


def fig1(mc, ref):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, s in zip(axes, A.STATS):
        r = ref[s]
        band(ax, r["two_sigma"][0], *r["one_sigma"], r["two_sigma"][1])
        ax.hist(mc[s], bins=60, histtype="stepfilled", color=MUTED, alpha=0.55, lw=0, label="control draws")
        ax.axvline(r["field_median"], color=INK, ls="--", lw=1.5, label="field median")
        ax.axvline(r["host"], color=ACCENT, lw=2, label="hosts")
        ax.text(0.97, 0.95, f"{r['one_sided_sigma']:+.2f}σ\n(one-sided)", transform=ax.transAxes,
                ha="right", va="top", fontsize=10, color=INK)
        ax.set(xlabel=f"{STAT_NAMES[s]} (km/s)", ylabel="draws" if s == "sigma_W" else None)
        ax.set_title(STAT_NAMES[s], fontsize=10, color=INK)
    hd, lb = axes[0].get_legend_handles_labels()
    fig.legend(hd, lb, frameon=False, fontsize=9, loc="lower center", ncol=5, bbox_to_anchor=(0.5, 0.93))
    fig.suptitle(f"Method A main: {ref['n_hosts']} WASP hosts vs {len(mc['sigma_W']):,} matched control draws "
                 "(Gaia DR3)", fontsize=11, color=INK, y=1.08)
    save(fig, "fig1_method_a_main_draws.png")


def fig2(h, mc, host, uvw, w):
    comp = {}
    for c in "UVW":
        s = f"sigma_{c}"
        q = np.percentile(mc[s], [2.5, 16, 50, 84, 97.5])
        p, sig = A.one_sided_sigma(host[s], mc[s])
        comp[s] = {"host": host[s], "field_median": q[2], "one_sigma": [q[1], q[3]], "two_sigma": [q[0], q[4]],
                   "tail_fraction": p, "one_sided_sigma": sig}
    A.write_json(COMP_JSON, {"header": A.HEADER, "what": "Method A main (97 hosts): per-component dispersions "
                             "(np.std) on the dr3x_A_wasp control draws; new exploratory calculation (E63)",
                             "n_hosts": len(h), "n_draws": len(mc["sigma_W"]), **comp})
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1, 1.25]})
    for i, c in enumerate("UVW"):
        r = comp[f"sigma_{c}"]
        a.plot(r["two_sigma"], [i, i], color=HOST, alpha=0.35, lw=10, solid_capstyle="butt",
               label="control 2σ" if i == 0 else None)
        a.plot(r["one_sigma"], [i, i], color=HOST, lw=10, solid_capstyle="butt", label="control 1σ" if i == 0 else None)
        a.plot(r["field_median"], i, "|", color=INK, ms=18, mew=2, label="control median" if i == 0 else None)
        a.plot(r["host"], i, "o", color=ACCENT, ms=9, mec="white", mew=1.5, label="hosts" if i == 0 else None)
        a.text(r["two_sigma"][1] + 0.6, i, f"{r['one_sided_sigma']:+.2f}σ", va="center", fontsize=9, color=INK)
    a.set_yticks(range(3), [r"$\sigma_U$", r"$\sigma_V$", r"$\sigma_W$"])
    a.set_ylim(-0.6, 2.6)
    a.invert_yaxis()
    a.set_xlabel("dispersion (km/s)")
    a.set_title("(a) per-component dispersion", fontsize=10, color=INK)
    a.legend(frameon=False, fontsize=8, loc="lower right")
    m = w > 0
    hb = b.hist2d(uvw[m, 0], uvw[m, 1], weights=w[m], bins=[np.linspace(-150, 150, 61), np.linspace(-150, 100, 51)],
                  cmap=SEQ, norm=LogNorm(), cmin=1e-6)
    fig.colorbar(hb[3], ax=b, label="expected draws per bin", pad=0.02)
    b.scatter(h["U"], h["V"], s=22, facecolor="white", edgecolor=INK, lw=0.9, zorder=3, label=f"hosts (N={len(h)})")
    b.set(xlabel="U (km/s)", ylabel="V (km/s)", xlim=(-150, 150), ylim=(-150, 100))
    b.set_title("(b) U–V: matched field (weighted by draw probability) and hosts", fontsize=10, color=INK)
    b.legend(frameon=False, fontsize=8, loc="lower left")
    fig.suptitle("Method A main, same draws as Fig. 1", fontsize=11, color=INK, y=1.04)
    save(fig, "fig2_method_a_components_uv.png")
    out_of_range = int(((np.abs(h["U"]) > 150) | (h["V"] < -150) | (h["V"] > 100)).sum())
    return comp, out_of_range


# ----------------------------------------------------------------------------- Method B (fig 3)
def picks():
    """DR2 and DR3 matched control picks (chunk 0, checkpoint-verified) with their properties, and the hosts."""
    cfg = config.load()
    hb, h2 = D.hosts(cfg)
    f2 = pd.read_csv(config.out_dir(cfg) / "field_ms.csv.gz")
    r2, _ = D.replay(h2, f2, cfg, "dr3x_B_dr2match")
    p2 = f2.iloc[r2][["source_id", "dist_pc", "z"]].reset_index(drop=True)
    del f2
    raw2 = field.download(cfg)[["source_id", "phot_g_mean_mag", "radial_velocity_error"]]
    p2 = p2.merge(raw2[raw2["source_id"].isin(set(p2["source_id"]))], on="source_id", how="left",
                  validate="many_to_one")
    del raw2
    fb = S.load_field_b()
    r3, _ = D.replay(hb, fb, cfg, "dr3x_B_wasp")
    p3 = fb.iloc[r3][["source_id", "dist_pc", "z"]].reset_index(drop=True)
    del fb
    want = set(p3["source_id"])
    raw3 = pd.concat([d[d["source_id"].isin(want)] for d in
                      (pd.read_csv(f, usecols=["source_id", "phot_g_mean_mag", "radial_velocity_error"])
                       for f in sorted(S.FIELD_DIR.glob("field_*.csv.gz")))])
    p3 = p3.merge(raw3, on="source_id", how="left", validate="many_to_one")
    for p in (p2, p3, hb, h2):
        p["abs_z"] = p["z"].abs()
    return p2, p3, hb, h2


def fig3(p2, p3, hb):
    panels = [("dist_pc", "distance (pc)", np.linspace(0, 1800, 73), False),
              ("phot_g_mean_mag", "G (mag)", np.linspace(7, 15.5, 69), False),
              ("radial_velocity_error", "Gaia radial_velocity_error (km/s)", np.logspace(-1.5, 1.3, 57), True),
              ("abs_z", "|z| (pc)", np.linspace(0, 600, 61), False)]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, (c, xl, bins, logx) in zip(axes.ravel(), panels):
        for tag, d in (("DR2", p2), ("DR3", p3), ("hosts", hb)):
            if c in d:
                v = d[c].dropna()
                print(f"{c:22s} {tag:5s} outside plotted range: {((v < bins[0]) | (v > bins[-1])).mean():.4%}")
        ax.hist(p2[c].dropna(), bins=bins, density=True, histtype="step", lw=2, color=DR2_COL,
                label=f"DR2 matched controls (n={p2[c].notna().sum():,})")
        ax.hist(p3[c].dropna(), bins=bins, density=True, histtype="step", lw=2, color=DR3_COL,
                label=f"DR3 matched controls (n={p3[c].notna().sum():,})")
        if c in hb:
            ax.hist(hb[c].dropna(), bins=bins, density=True, histtype="step", lw=1.5, color=ACCENT, ls="--",
                    label=f"hosts, DR3 (n={hb[c].notna().sum()})")
        else:
            ax.text(0.97, 0.6, "hosts not shown:\nMethod B uses archive RVs", transform=ax.transAxes, ha="right",
                    fontsize=8, color=MUTED)
        if logx:
            ax.set_xscale("log")
        ax.set(xlabel=xl, ylabel="normalized density")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("Method B, 101 hosts: matched control stars (chunk 0, 5,000 iterations, picks with repetition), "
                 "DR2 vs DR3", fontsize=11, color=INK, y=1.02)
    fig.tight_layout()
    save(fig, "fig3_method_b_controls_dr2_vs_dr3.png")


# ----------------------------------------------------------------------------- forest (fig 4)
def fig4():
    st, dg = json.loads(STAGE2.read_text()), json.loads(DIAG_JSON.read_text())
    fz, b2, b3 = st["comparison"]["frozen_dr2_107_headline"], st["method_B_dr2_same_hosts"], st["method_B_dr3"]
    d2, d5 = dg["2_dr3_field_restricted_to_dr2_field"], dg["3_d5_halo_cut"]
    rows = [("DR2 frozen headline (pre-registered)", fz["n"], fz["delta"], fz["ci95"], INK),
            ("DR2, same 101 hosts", b2["n"], b2["delta"], b2["ci95"], DR2_COL),
            ("DR3", b3["n"], b3["delta"], b3["ci95"], DR3_COL),
            ("DR3, field = DR2-field stars (diagnostic)", d2["n"], d2["delta"], d2["ci95"], DR3_COL),
            ("DR3, D5 halo cut (diagnostic)", d5["n"], d5["delta"], d5["ci95"], DR3_COL)]
    fig, ax = plt.subplots(figsize=(8.5, 3.8))
    for i, (lab, n, d, (lo, hi), col) in enumerate(rows):
        y = len(rows) - 1 - i
        ax.plot([lo, hi], [y, y], color=col, lw=2, solid_capstyle="round")
        ax.plot(d, y, "o", ms=8, color=col, mfc=col if "diagnostic" not in lab else "white", mew=2)
        ax.text(1.02, y, f"{d:.2f} [{lo:.2f}, {hi:.2f}]", va="center", fontsize=8, color=INK,
                transform=ax.get_yaxis_transform())
    ax.set_yticks(range(len(rows)), [f"{r[0]} (N={r[1]})" for r in rows][::-1])
    ax.axvline(0, color=INK, lw=1)
    ax.axvline(DELTA_HS19, color=ACCENT, lw=1.2, ls="--")
    ax.text(DELTA_HS19, len(rows) - 0.45, f"δ_HS19 = {DELTA_HS19:.2f}", ha="center", va="bottom", fontsize=8,
            color=INK)
    ax.set_xlabel("Δ = host − control median, HS19 Eq. 1 (km/s), 95 % bootstrap interval")
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_title("Method B: Δ by data release and diagnostic variant (open = post-hoc diagnostic)",
                 fontsize=10, color=INK, pad=18)
    save(fig, "fig4_method_b_forest.png")


def main(argv):
    which = set(argv) or {"1", "2", "3", "4"}
    if which & {"1", "2"}:
        h, mc, host, ref, uvw, w = method_a_draws()
        if "1" in which:
            fig1(mc, ref)
        if "2" in which:
            comp, oor = fig2(h, mc, host, uvw, w)
            print(json.dumps(comp, indent=1), f"\nhosts outside the U-V plot range: {oor}")
        del uvw, w
    if "3" in which:
        p2, p3, hb, _ = picks()
        fig3(p2, p3, hb)
    if "4" in which:
        fig4()


if __name__ == "__main__":
    main(sys.argv[1:])
