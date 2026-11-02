"""HS19 appendix quality cuts (Lindegren et al. 2018; Marchetti et al. 2018).

Each cut_N(df, cfg) returns a boolean numpy array (True = passes). Missing
values fail. Cut numbering follows the HS19 appendix.
"""
import numpy as np

# Gaia DR2 columns each cut needs (DR3 mapping lives in Phase 2).
CUT_COLUMNS = {
    1: ["parallax_over_error"],
    2: ["mean_varpi_factor_al"],
    3: ["visibility_periods_used"],
    4: ["astrometric_chi2_al", "astrometric_n_good_obs_al", "phot_g_mean_mag"],
    5: ["astrometric_gof_al"],
    6: ["astrometric_excess_noise_sig"],
    7: ["rv_nb_transits"],
    8: ["phot_bp_rp_excess_factor", "bp_rp"],
    9: ["phot_bp_mean_flux_over_error"],
    10: ["phot_rp_mean_flux_over_error"],
}


def _a(df, col):
    return np.asarray(df[col], dtype=float)


def unit_weight_error(df):
    """u = sqrt(astrometric_chi2_al / (astrometric_n_good_obs_al - 5))."""
    return np.sqrt(_a(df, "astrometric_chi2_al") / (_a(df, "astrometric_n_good_obs_al") - 5.0))


def uwe_threshold(g, a=1.2, b=-0.2, g0=19.5, form="lindegren"):
    """Cut-4 threshold on u.

    lindegren: a * max(1, exp(b * (G - g0)))   -- Lindegren+2018 eq. C.1
    literal:   a * max(1, exp(b * G - g0))     -- as printed in HS19 (no parentheses)
    """
    g = np.asarray(g, dtype=float)
    if form == "lindegren":
        e = np.exp(b * (g - g0))
    elif form == "literal":
        e = np.exp(b * g - g0)
    else:
        raise ValueError(f"unknown cut-4 form {form!r}")
    return a * np.maximum(1.0, e)


def cut_1(df, c):
    return _a(df, "parallax_over_error") > c["parallax_over_error_min"]


def cut_2(df, c):
    lo, hi = c["mean_varpi_factor_al_range"]
    x = _a(df, "mean_varpi_factor_al")
    return (x > lo) & (x < hi)


def cut_3(df, c):
    return _a(df, "visibility_periods_used") > c["visibility_periods_used_min"]


def cut_4(df, c):
    u = c["uwe"]
    thr = uwe_threshold(_a(df, "phot_g_mean_mag"), u["a"], u["b"], u["g0"], u["form"])
    return unit_weight_error(df) < thr


def cut_5(df, c):
    return _a(df, "astrometric_gof_al") < c["astrometric_gof_al_max"]


def cut_6(df, c):
    return _a(df, "astrometric_excess_noise_sig") < c["astrometric_excess_noise_sig_max"]


def cut_7(df, c):
    return _a(df, "rv_nb_transits") > c["rv_nb_transits_min"]


def bp_rp_excess_bounds(bp_rp, lo_a=1.0, lo_b=0.015, hi_a=1.3, hi_b=0.06):
    c2 = np.asarray(bp_rp, dtype=float) ** 2
    return lo_a + lo_b * c2, hi_a + hi_b * c2


def cut_8(df, c):
    e = c["bp_rp_excess"]
    lo, hi = bp_rp_excess_bounds(_a(df, "bp_rp"), e["lo_a"], e["lo_b"], e["hi_a"], e["hi_b"])
    x = _a(df, "phot_bp_rp_excess_factor")
    return (x > lo) & (x < hi)


def cut_9(df, c):
    return _a(df, "phot_bp_mean_flux_over_error") > c["bp_flux_over_error_min"]


def cut_10(df, c):
    return _a(df, "phot_rp_mean_flux_over_error") > c["rp_flux_over_error_min"]


CUTS = {i: globals()[f"cut_{i}"] for i in range(1, 11)}


def apply(df, cuts_cfg, which):
    """Return (mask passing all cuts in `which`, dict cut -> individual mask).

    NaN comparisons are False, so missing values fail every cut.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        masks = {i: CUTS[i](df, cuts_cfg) for i in which}
    mask = np.ones(len(df), dtype=bool)
    for m in masks.values():
        mask &= m
    return mask, masks


def cut_table(df, cuts_cfg, which):
    """Sequential pass counts, for logging."""
    rows, mask = [], np.ones(len(df), dtype=bool)
    _, masks = apply(df, cuts_cfg, which)
    for i in which:
        mask &= masks[i]
        rows.append({"cut": i, "pass_alone": int(masks[i].sum()), "remaining": int(mask.sum())})
    return rows


def server_where(cuts_cfg, alias="g"):
    """ADQL WHERE clause for the field prefilter.

    Applies cuts 1-3, 5-7, 9, 10 exactly. Cut 8 uses the looser lower bound
    (server_lo_b) and cut 4 is omitted, so both readings of each can be applied
    locally. The Gaia archive ADQL has no CASE/GREATEST, which cut 4 needs.
    """
    c, e, a = cuts_cfg, cuts_cfg["bp_rp_excess"], alias
    lo, hi = c["mean_varpi_factor_al_range"]
    lob = min(e["lo_b"], c["server_lo_b"])
    return "\n AND ".join([
        f"{a}.radial_velocity IS NOT NULL",
        f"{a}.parallax IS NOT NULL",                                   # 5-parameter solution
        f"{a}.parallax_over_error > {c['parallax_over_error_min']}",
        f"{a}.mean_varpi_factor_al > {lo} AND {a}.mean_varpi_factor_al < {hi}",
        f"{a}.visibility_periods_used > {c['visibility_periods_used_min']}",
        f"{a}.astrometric_gof_al < {c['astrometric_gof_al_max']}",
        f"{a}.astrometric_excess_noise_sig < {c['astrometric_excess_noise_sig_max']}",
        f"{a}.rv_nb_transits > {c['rv_nb_transits_min']}",
        f"{a}.phot_bp_rp_excess_factor > {e['lo_a']} + {lob}*{a}.bp_rp*{a}.bp_rp",
        f"{a}.phot_bp_rp_excess_factor < {e['hi_a']} + {e['hi_b']}*{a}.bp_rp*{a}.bp_rp",
        f"{a}.phot_bp_mean_flux_over_error > {c['bp_flux_over_error_min']}",
        f"{a}.phot_rp_mean_flux_over_error > {c['rp_flux_over_error_min']}",
    ])
