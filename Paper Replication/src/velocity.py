"""Per-star velocity likelihoods L_v(v) and L_u(u) (Morton & Winn 2014 Sec. 4.2)."""

import numpy as np
from scipy import stats

from .inclination import gaussian_likelihood

KM_S_PER_RSUN_PER_DAY = 2 * np.pi * 695_700.0 / 86_400.0  # 2πR/P with R in R_sun, P in d


def sample_split_normal(mu, err_hi, err_lo, n, rng):
    """Asymmetric errors as a two-piece normal: σ = err_hi above mu, err_lo below."""
    z = rng.standard_normal(n)
    return mu + np.where(z > 0, z * err_hi, z * err_lo)


def sample_equatorial_velocity(row, n, rng, *, alpha, lat_mean_deg, lat_sd_deg):
    """Monte Carlo over stellar radius, measured period and spot latitude:
        P_rot(ℓ) = P_eq / (1 - α sin²ℓ)   =>   P_eq = P_rot (1 - α sin²ℓ)
        v = 2π R / P_eq
    Draws with R <= 0 or P <= 0 are rejected."""
    R = sample_split_normal(row.R_star, row.R_star_err_hi, row.R_star_err_lo, n, rng)
    P = rng.normal(row.P_rot, row.P_rot_err, n)
    lat = np.deg2rad(rng.normal(lat_mean_deg, lat_sd_deg, n))
    ok = (R > 0) & (P > 0)
    P_eq = P[ok] * (1.0 - alpha * np.sin(lat[ok]) ** 2)
    return KM_S_PER_RSUN_PER_DAY * R[ok] / P_eq


class GridLikelihood:
    """A density tabulated on a grid, callable anywhere (zero outside the grid)."""

    def __init__(self, grid, values):
        self.grid, self.values = np.asarray(grid), np.asarray(values)

    def __call__(self, x):
        return np.interp(x, self.grid, self.values, left=0.0, right=0.0)


def velocity_likelihood_from_samples(v_samples, n_grid=1500, pad_bw=6.0):
    """Smooth the Monte Carlo v distribution with a Gaussian KDE and tabulate it.
    Returns (L_v, v_grid), where v_grid spans the support with padding and is the
    grid the cos i integrals should run over."""
    sub = v_samples[: min(len(v_samples), 40_000)]
    kde = stats.gaussian_kde(sub)
    bw = np.sqrt(kde.covariance[0, 0])
    lo = max(0.0, v_samples.min() - pad_bw * bw)
    hi = v_samples.max() + pad_bw * bw
    grid = np.linspace(lo, hi, n_grid)
    return GridLikelihood(grid, kde(grid)), grid


def vsini_likelihood(row):
    """Gaussian L_u from the table value. The table value is already corrected
    for differential rotation (MW14 Table 1 note f); it is used unchanged."""
    assert row.vsini_diffrot_corrected, "v sin i must be the once-corrected table value"
    return gaussian_likelihood(row.vsini, row.vsini_err)


def vsini_excess_sigma(row, rng, *, alpha, lat_mean_deg, lat_sd_deg, n=200_000):
    """(v sin i − v) / σ, with v = 2πR/P_eq from the same Monte Carlo as L_v
    (spot-latitude corrected). v is the median draw, σ_v half the 16–84% range,
    combined in quadrature with the v sin i error. Values > 2 are unphysical
    at 2σ (v sin i cannot exceed v)."""
    v = sample_equatorial_velocity(row, n, rng, alpha=alpha,
                                   lat_mean_deg=lat_mean_deg, lat_sd_deg=lat_sd_deg)
    q16, q50, q84 = np.percentile(v, [16, 50, 84])
    sigma_v = 0.5 * (q84 - q16)
    return dict(v=q50, sigma_v=sigma_v,
                excess_sigma=(row.vsini - q50) / np.hypot(row.vsini_err, sigma_v))
