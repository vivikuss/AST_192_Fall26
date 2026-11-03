"""Step 11: paper figures (Fig. 1 Pleiades MS, Fig. 2 CMD, Fig. 3 UVW ellipses, Fig. 4 MC)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from matplotlib.patches import Ellipse

from .main_sequence import ms_poly

HOST = "#2a78d6"        # categorical slot 1 (hosts)
ACCENT = "#eb6834"      # categorical slot 2 (MS line / host value)
INK = "#3a3a38"
MUTED = "#8a8985"
SEQ = LinearSegmentedColormap.from_list(
    "blue_seq", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#104281", "#0d366b"])

plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 150, "savefig.bbox": "tight"})


def fig_pleiades(c0, M0, coeffs_pub, coeffs_refit=None, spline=None, path=None, title=None):
    """Fig. 1: Pleiades members, published polynomial, and (optionally) our refit."""
    fig, ax = plt.subplots(figsize=(4.5, 5))
    ax.scatter(c0, M0, s=6, color=HOST, alpha=0.6, lw=0, label="Pleiades member")
    g = np.linspace(-0.15, 3.7, 400)                       # polynomial over HS19's colour range
    ax.plot(g, ms_poly(g, coeffs_pub), color=INK, lw=2, label="HS19 polynomial")
    g = np.linspace(np.nanmin(c0), np.nanmax(c0), 400)     # refit only where members constrain it
    if coeffs_refit is not None:
        ax.plot(g, ms_poly(g, coeffs_refit), color=ACCENT, lw=2, ls="--", label="refit polynomial")
    if spline is not None:
        ax.plot(g, spline(g), color=ACCENT, lw=1, ls=":", label="refit spline")
    ax.set(xlabel=r"$(G_{BP}-G_{RP})_0$", ylabel=r"$M_{G,0}$", xlim=(-0.3, 3.9), ylim=(15, -1))
    if title:
        ax.set_title(title, fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=8)
    _save(fig, path)


def fig_cmd(field, host, cfg, path=None, title=None):
    """Fig. 2: field density + hosts + MS line and the 1-mag-above boundary."""
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    h = ax.hist2d(field["bp_rp0"], field["M_G0"], bins=[np.linspace(-0.1, 2.2, 116), np.linspace(0, 11, 111)],
                  cmap=SEQ, norm=LogNorm(), cmin=1)
    fig.colorbar(h[3], ax=ax, label="Control sample stars", pad=0.02)
    g = np.linspace(-0.1, 2.2, 300)
    ms = ms_poly(g, cfg["main_sequence"]["coeffs"])
    ax.plot(g, ms, color=ACCENT, lw=2, label="Pleiades main sequence")
    ax.plot(g, ms - cfg["main_sequence"]["delta_mag"], color=ACCENT, lw=1, ls="--",
            label="1 mag above (selection edge)")
    ax.scatter(host["bp_rp0"], host["M_G0"], s=18, facecolor="white", edgecolor=INK, lw=0.8,
               zorder=3, label=f"HJ host (N={len(host)})")
    ax.set(xlim=(-0.1, 2.2), ylim=(11, 0), xlabel=r"$(G_{BP}-G_{RP})_0$", ylabel=r"$M_{G,0}$")
    if title:
        ax.set_title(title, fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    _save(fig, path)


def _cov_ellipse(ax, x, y, n_sigma, **kw):
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    ang = np.degrees(np.arctan2(vecs[1, 1], vecs[0, 1]))
    w, h = 2 * n_sigma * np.sqrt(vals[::-1])
    ax.add_patch(Ellipse((np.mean(x), np.mean(y)), w, h, angle=ang, **kw))


def fig_uvw(field, host, path=None, label="HJ hosts", points=True):
    """Fig. 3: 1/2/3-sigma Gaussian ellipses for field (grey fill) and hosts (blue lines), by F/G/K.

    points=False draws the ellipses only (no per-star velocities, e.g. for private-RV samples)."""
    pairs = [("U", "V"), ("U", "W"), ("V", "W")]
    fig, axes = plt.subplots(3, 3, figsize=(10, 10), sharex=True, sharey=True)
    for i, sp in enumerate(("F", "G", "K")):
        f, h = field[field["sptype"] == sp], host[host["sptype"] == sp]
        for j, (a, b) in enumerate(pairs):
            ax = axes[i, j]
            for k, (alpha, ns) in enumerate([(0.18, 3), (0.28, 2), (0.4, 1)]):
                _cov_ellipse(ax, f[a], f[b], ns, facecolor=MUTED, alpha=alpha, lw=0)
            for ns in (1, 2, 3):
                _cov_ellipse(ax, h[a], h[b], ns, facecolor="none", edgecolor=HOST, lw=2 if ns == 1 else 1.2,
                             alpha=1.0 if ns == 1 else 0.6)
            if points:
                ax.scatter(h[a], h[b], s=5, color=HOST, lw=0)
            ax.text(0.04, 0.94, f"{sp} stars", transform=ax.transAxes, fontsize=9, color=INK)
            ax.text(0.96, 0.94, f"N={len(h)}", transform=ax.transAxes, ha="right", fontsize=9, color=INK)
            ax.set(xlim=(-100, 100), ylim=(-100, 100))
            ax.grid(color="#e5e4df", lw=0.6)
            ax.set_ylabel(f"{b} [km/s]")
            ax.yaxis.set_tick_params(labelleft=True)
            if i == 2:
                ax.set_xlabel(f"{a} [km/s]")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(facecolor=MUTED, alpha=0.4, label="Field MS stars in the same colour bin: 1, 2, 3$\\sigma$"),
                        Line2D([], [], color=HOST, lw=2, label=f"{label}: 1, 2, 3$\\sigma$")],
               loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    _save(fig, path)


def fig_mc(results, path=None, labels=None, shared_x=False):
    """Fig. 4: control-sample dispersion histograms with host values and percentiles.

    shared_x=True puts every panel on one horizontal range (and the same 0.1 km/s bins)."""
    n = len(results)
    fig, axes = plt.subplots(n, 1, figsize=(6.5, 3.0 * n), squeeze=False)
    if shared_x:
        glo = min(min(d.min(), s["host_dispersion"]) for s, d in results.values()) - 0.5
        ghi = max(max(d.max(), s["host_dispersion"]) for s, d in results.values()) + 0.5
    for ax, (key, (summary, disp)) in zip(axes[:, 0], results.items()):
        lo = min(disp.min(), summary["host_dispersion"]) - 0.5
        hi = max(disp.max(), summary["host_dispersion"]) + 0.5
        if shared_x:
            lo, hi = glo, ghi
        ax.hist(disp, bins=np.arange(np.floor(lo), np.ceil(hi) + 0.1, 0.1), color=HOST, lw=0)
        for q, ls in ((2, ":"), (16, "--"), (50, "-"), (84, "--"), (98, ":")):
            ax.axvline(np.percentile(disp, q), color=INK, ls=ls, lw=1.2)
        ax.axvline(summary["host_dispersion"], color=ACCENT, lw=3,
                   label=f"{(labels or {}).get(key, key)}: {summary['host_dispersion']:.2f} km/s")
        p = summary["p_value"]
        ptxt = (f"0 of {summary['n_iter']} controls as cold" if summary["n_le"] == 0
                else f"{summary['n_le']} of {summary['n_iter']} as cold (p = {p:.3f})")
        ztxt = f"Gaussian-equivalent {summary['z_gauss']:.1f}$\\sigma$"
        ax.text(0.98, 0.95, f"N = {summary['n_host']}\n{ptxt}\n{ztxt}",
                transform=ax.transAxes, ha="right", va="top", fontsize=9, color=INK)
        ax.set(xlabel="Sample velocity dispersion (mean residual speed) [km/s]", ylabel="Control samples")
        ax.set_ylim(0, ax.get_ylim()[1] * 1.35)          # headroom for legend and text
        if shared_x:
            ax.set_xlim(np.floor(lo), np.ceil(hi))
        ax.legend(frameon=False, loc="upper left", fontsize=9)
    fig.tight_layout()
    _save(fig, path)


def _save(fig, path):
    if path is not None:
        fig.savefig(path)
    plt.close(fig)
