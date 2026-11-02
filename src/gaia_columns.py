"""Gaia DR2 columns selected for every sample (astrometry, RV, photometry, cut inputs)."""
from .cuts import CUT_COLUMNS

ASTROMETRY = ["source_id", "ra", "dec", "l", "b", "parallax", "parallax_error",
              "pmra", "pmra_error", "pmdec", "pmdec_error",
              "radial_velocity", "radial_velocity_error",
              "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag"]
CUT_INPUTS = sorted({c for cols in CUT_COLUMNS.values() for c in cols} - set(ASTROMETRY))
GAIA_COLUMNS = ASTROMETRY + [c for c in CUT_INPUTS if c not in ASTROMETRY]
DIST_COLUMNS = ["r_est", "r_lo", "r_hi", "result_flag"]


def select_list(gaia_alias="g", dist_alias="d"):
    cols = [f"{gaia_alias}.{c}" for c in GAIA_COLUMNS]
    cols += [f"{dist_alias}.{c}" for c in DIST_COLUMNS]
    return ", ".join(cols)
