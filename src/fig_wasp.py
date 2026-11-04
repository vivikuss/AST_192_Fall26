"""Figure of the frozen WASP-only results (presentation only; no new statistics).

    python -m src.fig_wasp      # -> figures/fig_wasp_primary.png

(a) the primary test's 40,000 matched-field total dispersions (data/processed/wasp/mc_wasp_primary.npz)
    with the WASP hosts' value and the field median; styled like figures/phase1_dr2/fig4_mc.png.
(b) Delta_W and its 95 % percentile interval for B and D1-D5, from data/processed/wasp/results.json.
No sigma_W result and no EBLM value is shown.
"""
import json

import numpy as np
import matplotlib.pyplot as plt

from . import figures                      # project palette and rcParams
from .config import ROOT
from .figures import ACCENT, HOST, INK, MUTED

WASP = ROOT / "data" / "processed" / "wasp"
OUT = ROOT / "figures" / "fig_wasp_primary.png"
ROWS = [("B_primary", "Primary (98 WASP hosts)"),
        ("D1_membership_by_name", "D1  WASP by name (100)"),
        ("D2_wasp_z_dist", "D2  matched on z + distance"),
        ("D3_d_lt_600", "D3  d < 600 pc (97)"),
        ("D4_colour_0.5_1.5", "D4  0.5 ≤ (BP−RP)₀ ≤ 1.5 (95)"),
        ("D5_halo_cut", "D5  field halo cut")]


def main():
    R = json.loads((WASP / "results.json").read_text())
    z = np.load(WASP / "mc_wasp_primary.npz")
    d, B = z["disp"], R["B_primary"]
    # the saved controls must be the ones results.json summarises
    assert np.median(d) == B["control_median"] and np.std(d, ddof=1) == B["control_sd"]
    assert np.mean(d <= z["host"]) == B["fraction_controls_le"]

    plt.rcParams.update({"font.size": 13})
    fig, (a, b) = plt.subplots(1, 2, figsize=(14, 5.6), gridspec_kw={"width_ratios": [1.05, 1]})

    # (a)
    host, med = B["sigma"], B["control_median"]
    lo, hi = np.floor(min(d.min(), host) - 0.5), np.ceil(max(d.max(), host) + 0.5)
    a.hist(d, bins=np.arange(lo, hi + 0.1, 0.1), color=HOST, lw=0, label="Matched field (40,000 samples)")
    for q, ls in ((2, ":"), (16, "--"), (84, "--"), (98, ":")):
        a.axvline(np.percentile(d, q), color=INK, ls=ls, lw=1.2)
    a.axvline(med, color=INK, lw=2, label=f"Field median: {med:.2f} km/s")
    a.axvline(host, color=ACCENT, lw=3.5, label=f"WASP hosts: {host:.2f} km/s")
    a.text(0.98, 0.66, f"{B['n_le']} of {B['n_iter']:,} controls\nas cold as the hosts\n"
           f"(Gaussian-equivalent {B['z_gauss']:.2f}σ)\n\ndashed: 16th/84th pct\ndotted: 2nd/98th pct",
           transform=a.transAxes, ha="right", va="top", fontsize=11, color=INK)
    a.set(xlim=(lo, hi), xlabel="Total velocity dispersion (mean residual speed) [km/s]", ylabel="Control samples")
    a.set_ylim(0, a.get_ylim()[1] * 1.45)
    a.legend(loc="upper left", fontsize=12, frameon=True, framealpha=0.95, edgecolor="none")
    a.set_title("(a) Primary test: WASP hosts vs matched field", loc="left", fontsize=14)

    # (b)
    y = np.arange(len(ROWS))[::-1]
    dh = R["delta_hs19"]
    b.axvline(0, color=MUTED, ls="--", lw=1.5, label="No offset")
    b.axvline(dh, color=INK, ls=":", lw=2, label=f"HS19-sized offset ({dh:.2f} km/s)")
    for yi, (key, lab) in zip(y, ROWS):
        r = B if key == "B_primary" else R["D_sensitivity"][key]["B"]
        est, (l, u) = r["B"]["delta_hat"], r["B"]["ci95"]
        assert u < 0 and l <= dh <= u and r["case"] == 1, (key, l, u)       # backs the footer text
        primary = key == "B_primary"
        col = ACCENT if primary else HOST
        b.plot([l, u], [yi, yi], color=col, lw=4 if primary else 2.5, solid_capstyle="round")
        b.plot(est, yi, "o", color=col, ms=11 if primary else 8, mec="white", mew=1.5, zorder=3)
        b.text(1.5, yi, f"{est:+.2f}  [{l:+.2f}, {u:+.2f}]", va="center", fontsize=12, color=INK)
    b.set_yticks(y, [lab for _, lab in ROWS])
    b.set_xlim(-12, 7.5)
    b.set_ylim(-1.9, len(ROWS) - 0.3)
    b.set_xlabel("WASP hosts minus matched-field median [km/s]\n(point and 95 % percentile interval)")
    b.legend(frameon=False, loc="lower left", fontsize=11, ncol=2)
    b.set_title("(b) Primary and sensitivity tests D1–D5", loc="left", fontsize=14)
    for ax in (a, b):
        ax.spines[["top", "right"]].set_visible(False)

    fig.text(0.01, 0.005, "Frozen pre-registration PREREG_wasp_kinematics_2026-11-04 (sha256 fece9b66…). "
             "Gaia DR2; HS19 eq. 1 statistic. All intervals exclude 0 and contain −6.64 km/s (case 1).",
             fontsize=10, color=MUTED)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(OUT, dpi=200)
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
