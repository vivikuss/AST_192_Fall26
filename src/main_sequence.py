"""Step 7: Pleiades zero-age main sequence and the 'less than 1 mag above' cut."""
import numpy as np
from scipy.interpolate import UnivariateSpline


def ms_poly(color0, coeffs):
    """M_G0 of the main sequence at (BP-RP)_0; coeffs in increasing order."""
    return np.polynomial.polynomial.polyval(np.asarray(color0, dtype=float), coeffs)


def on_main_sequence(M_G0, color0, cfg, relation=None):
    """True where the star is less than delta_mag above (brighter than) the MS.

    Magnitudes increase downward, so 'above' = smaller M_G0. Keep M_G0 > MS(c) - delta.
    Stars below the MS are kept (no lower limit in HS19).
    """
    ms = cfg["main_sequence"]
    f = relation if relation is not None else (lambda c: ms_poly(c, ms["coeffs"]))
    c = np.asarray(color0, dtype=float)
    keep = np.asarray(M_G0, dtype=float) > f(c) - ms["delta_mag"]
    if ms["color_range"] is not None:
        lo, hi = ms["color_range"]
        keep &= (c >= lo) & (c <= hi)
    return keep


def refit(color0, M_G0, order=8, s=None, n_dense=1000):
    """HS19 procedure: smoothing spline through the Pleiades, then an
    order-`order` polynomial fit to the densely sampled spline.

    Returns (poly coeffs increasing order, spline, dense colour grid).
    """
    c = np.asarray(color0, dtype=float)
    m = np.asarray(M_G0, dtype=float)
    o = np.argsort(c)
    c, m = c[o], m[o]
    # UnivariateSpline needs strictly increasing x: average exact duplicates
    uc, inv = np.unique(c, return_inverse=True)
    um = np.bincount(inv, weights=m) / np.bincount(inv)
    w = np.sqrt(np.bincount(inv))
    spl = UnivariateSpline(uc, um, w=w, k=3, s=s)
    grid = np.linspace(uc[0], uc[-1], n_dense)
    coeffs = np.polynomial.polynomial.polyfit(grid, spl(grid), order)
    return coeffs, spl, grid
