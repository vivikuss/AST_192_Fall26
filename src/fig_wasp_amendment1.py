"""Forest plot for the amendment 1 correction: 107-host rerun vs the frozen 98-host run (presentation only).

    python -m src.fig_wasp_amendment1    # -> notes/figures/fig_wasp_amendment1_forest.png

Every value is read from data/processed/wasp_amendment1/results.json (107 hosts) and the frozen
data/processed/wasp/results.json (98 hosts); nothing is retyped or recomputed. The diagnostic is not shown.
Left: WASP hosts vs matched field (B, D1-D5 B, E(i) sigma_W). Right: contrasts Delta_EBLM - Delta_W
(C, D1-D6 C, E(ii) sigma_W). 95 % percentile intervals. sigma_W rows are exploratory. DECISIONS E44.
"""
import json

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from . import figures                      # project palette and rcParams
from .config import ROOT
from .figures import HOST, INK, MUTED

NEW = ROOT / "data" / "processed" / "wasp_amendment1" / "results.json"
OLD = ROOT / "data" / "processed" / "wasp" / "results.json"
OUT = ROOT / "notes" / "figures" / "fig_wasp_amendment1_forest.png"
DELTA_HS19 = -6.64                         # section B of the frozen pre-registration (delta_hs19 checked below)
D_NAMES = {"D1_membership_by_name": "D1 WASP by name", "D2_wasp_z_dist": "D2 matched on z + distance",
           "D3_d_lt_600": "D3 d < 600 pc", "D4_colour_0.5_1.5": "D4 0.5 ≤ (BP−RP)₀ ≤ 1.5",
           "D5_halo_cut": "D5 field halo cut", "D6_no_J1500-33_C_only": "D6 without J1500-33"}


def rows_field(R):
    """(label, n, estimate, [lo, hi], exploratory) for B, D1-D5 B, E(i)."""
    out = [("B primary", R["B_primary"]["n"], R["B_primary"]["B"]["delta_hat"], R["B_primary"]["B"]["ci95"], False)]
    for k, v in R["D_sensitivity"].items():
        if "B" in v:
            out.append((D_NAMES[k], v["B"]["n"], v["B"]["B"]["delta_hat"], v["B"]["B"]["ci95"], False))
    e = R["E_exploratory_sigma_w"]["i_wasp_vs_field"]
    out.append(("E(i) σ_W (exploratory)", R["B_primary"]["n"], e["delta"], e["ci95"], True))
    return out


def rows_contrast(R):
    """(label, 'n_wasp/n_eblm', estimate, [lo, hi], exploratory) for C, D1-D6 C, E(ii)."""
    c = R["C_secondary"]
    out = [("C contrast", f"{c['n_wasp']}/{c['n_eblm']}", c["contrast"], c["ci95"], False)]
    for k, v in R["D_sensitivity"].items():
        out.append((D_NAMES[k], f"{v['C']['n_wasp']}/{v['C']['n_eblm']}", v["C"]["contrast"], v["C"]["ci95"], False))
    e = R["E_exploratory_sigma_w"]["ii_contrast"]
    out.append(("E(ii) σ_W (exploratory)", f"{e['n_wasp']}/{e['n_eblm']}", e["contrast"], e["ci95"], True))
    return out


def panel(ax, new, old, title, ref_lines):
    n = len(new)
    assert [r[0] for r in new] == [r[0] for r in old]
    for i, (rn, ro) in enumerate(zip(new, old)):
        y = n - 1 - i
        for r, dy, filled, col in ((rn, +0.17, True, HOST), (ro, -0.17, False, MUTED)):
            lo, hi = r[3]
            ax.plot([lo, hi], [y + dy, y + dy], color=col, lw=2, solid_capstyle="round", zorder=2)
            ax.plot(r[2], y + dy, "o", ms=8, mec=col, mew=2, mfc=col if filled else "white", zorder=3)
    ax.set_yticks(range(n))
    ax.set_yticklabels([f"{rn[0]}  ({rn[1]} vs {ro[1]})" for rn, ro in zip(new, old)][::-1])
    for t, r in zip(ax.get_yticklabels()[::-1], new):
        if r[4]:
            t.set_color(MUTED)
            t.set_style("italic")
    if new[-1][4]:                                       # separate the exploratory row
        ax.axhline(0.5, color=MUTED, lw=0.8, ls=":", zorder=1)
    for x, ls, lab in ref_lines:
        ax.axvline(x, color=INK, lw=1.2, ls=ls, zorder=1)
        if lab:                                          # label above the panel, clear of the line
            ax.text(x, 1.01, lab, ha="center", va="bottom", fontsize=10, color=INK,
                    transform=ax.get_xaxis_transform(), clip_on=False)
    ax.set_ylim(-0.6, n - 0.2)
    ax.set_xlabel("Δ (km/s)")
    ax.set_title(title, loc="left", fontsize=13, pad=22)
    ax.grid(axis="x", color="#e6e5e1", lw=0.8, zorder=0)
    ax.set_axisbelow(True)


def main():
    new, old = json.loads(NEW.read_text()), json.loads(OLD.read_text())
    assert abs(new["delta_hs19"] - DELTA_HS19) < 0.005 and abs(old["delta_hs19"] - DELTA_HS19) < 0.005
    assert new["prereg_sha256"] == old["sha256"]
    plt.rcParams.update({"font.size": 11})
    fig, (a, b) = plt.subplots(1, 2, figsize=(15, 6.4))
    panel(a, rows_field(new), rows_field(old), "WASP hosts vs matched field (Δ_W)",
          [(0, "-", "0"), (DELTA_HS19, "--", "δ_HS19 = −6.64")])
    panel(b, rows_contrast(new), rows_contrast(old), "Contrast Δ_EBLM − Δ_W",
          [(0, "-", "0")])
    handles = [Line2D([], [], color=HOST, lw=2, marker="o", ms=8, mfc=HOST, mec=HOST,
                      label="Corrected membership (107 WASP hosts)"),
               Line2D([], [], color=MUTED, lw=2, marker="o", ms=8, mfc="white", mec=MUTED,
                      label="Original run (98 hosts; membership join fault)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("WASP kinematics, corrected membership (107) vs original (98)", x=0.01, ha="left",
                 fontsize=15, fontweight="bold")
    fig.text(0.01, 0.925, "Points: estimate; bars: 95 % percentile bootstrap interval. Labels give N (107 run vs 98 run); "
             "contrasts give N_WASP/N_EBLM. σ_W rows are exploratory.", fontsize=10, color=MUTED, ha="left")
    fig.tight_layout(rect=(0, 0.05, 1, 0.92))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print(OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
