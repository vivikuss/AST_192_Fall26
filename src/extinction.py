"""Step 6: line-of-sight E(B-V) through the Capitanio et al. (2017) STILISM cube.

The cube holds differential reddening dE(B-V)/ds [mag/pc] on a 5 pc grid,
x, y in [-1997.5, 1997.5] pc, z in [-297.5, 297.5] pc, Sun at index
(399.5, 399.5, 59.5). Axis convention (x -> l=0, y -> l=90, z -> b=90) is
checked in tests/test_extinction.py against known clouds.
"""
from functools import lru_cache

import h5py
import numpy as np
from scipy.ndimage import map_coordinates

from .config import ROOT


@lru_cache(maxsize=2)
def load_cube(path):
    with h5py.File(path, "r") as f:
        d = f["stilism/cube_datas"]
        cube = d[...].astype(np.float32)
        step = np.asarray(d.attrs["gridstep_values"], dtype=float)
        sun = np.asarray(d.attrs["sun_position"], dtype=float)
    return cube, step, sun


def ebv_los(l_deg, b_deg, dist_pc, cfg, chunk=20000):
    """Integrated E(B-V) from the Sun to dist_pc along (l, b).

    Midpoint rule with n = ceil(d/step) equal sub-steps; trilinear interpolation;
    points outside the cube contribute zero (cfg beyond_cube = 'zero').
    """
    ec = cfg["extinction"]
    if ec["beyond_cube"] != "zero":
        raise NotImplementedError(ec["beyond_cube"])
    cube, grid, sun = load_cube(str(ROOT / ec["cube_file"]))
    l = np.radians(np.asarray(l_deg, dtype=float))
    b = np.radians(np.asarray(b_deg, dtype=float))
    d = np.asarray(dist_pc, dtype=float)
    out = np.full(d.shape, np.nan)
    good = np.isfinite(l) & np.isfinite(b) & np.isfinite(d) & (d > 0)
    idx = np.flatnonzero(good)
    unit = np.stack([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], axis=-1)
    for s in range(0, len(idx), chunk):
        ii = idx[s:s + chunk]
        n = np.maximum(1, np.ceil(d[ii] / ec["step_pc"]).astype(int))
        h = d[ii] / n                                   # actual step per star
        star = np.repeat(np.arange(len(ii)), n)
        k = np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)    # 0..n-1 per star
        dist = (k + 0.5) * h[star]
        xyz = unit[ii][star] * dist[:, None]            # pc, heliocentric Galactic
        coords = (xyz / grid + sun).T                   # fractional cube indices
        val = map_coordinates(cube, coords, order=1, mode="constant", cval=0.0)
        out[ii] = np.bincount(star, weights=val * h[star], minlength=len(ii))
    return out


def gaia_extinction(ebv, cfg):
    """A_G and E(BP-RP) from E(B-V) with Casagrande & VandenBerg (2018) mean R."""
    ec = cfg["extinction"]
    ebv = np.asarray(ebv, dtype=float)
    return ec["R_G"] * ebv, (ec["R_BP"] - ec["R_RP"]) * ebv


def deredden(df, cfg, ebv=None, dist_col="dist_pc"):
    """Add ebv, A_G, E_BPRP, bp_rp0, M_G0 columns (in place) and return df."""
    if ebv is None:
        ebv = ebv_los(df["l"], df["b"], df[dist_col], cfg)
    df["ebv"] = ebv
    df["A_G"], df["E_BPRP"] = gaia_extinction(ebv, cfg)
    df["bp_rp0"] = df["bp_rp"] - df["E_BPRP"]
    df["M_G0"] = df["phot_g_mean_mag"] - 5 * np.log10(df[dist_col] / 10.0) - df["A_G"]
    return df
