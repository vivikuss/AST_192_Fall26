"""Population inference for the Fisher concentration κ (Morton & Winn 2014 Sec. 2).

    L(κ) ≈ ∏_n (1/K) Σ_k f_cosI(c_nk | κ) / π_0(c_nk),   π_0 = 1 (uniform cos i)
    p(κ | data) ∝ L(κ) (1 + κ²)^(-3/4)
"""

from multiprocessing import Pool
from pathlib import Path

import numpy as np

from .inclination import posterior_cosi
from .obliquity import f_cosi, log_kappa_prior
from .velocity import (sample_equatorial_velocity, velocity_likelihood_from_samples,
                       vsini_likelihood)

COSI_GRID = np.linspace(0.0, 1.0, 4001)
# κ grid: fine where the posteriors live, logarithmic tail out to 1e4
KAPPA_GRID = np.unique(np.concatenate([np.linspace(0, 100, 401),
                                       np.geomspace(100, 1e4, 121)]))
# z grid for tabulating f_cosI: uniform plus dense near 0 (width ~ 1/sqrt(κ))
Z_GRID = np.unique(np.concatenate([np.linspace(0, 1, 801), np.geomspace(1e-4, 0.2, 400)]))


# ------------------------------------------------------------- per star
def star_cosi_posterior(row, *, method, rng, n_mc=200_000, alpha=0.23,
                        lat_mean_deg=20.0, lat_sd_deg=20.0, vsini_err_scale=1.0,
                        cosi_grid=COSI_GRID):
    """cos i posterior for one star on cosi_grid under `method` ('correct'|'mw14')."""
    v = sample_equatorial_velocity(row, n_mc, rng, alpha=alpha,
                                   lat_mean_deg=lat_mean_deg, lat_sd_deg=lat_sd_deg)
    L_v, v_grid = velocity_likelihood_from_samples(v)
    L_u = vsini_likelihood(row, err_scale=vsini_err_scale)
    return posterior_cosi(cosi_grid, L_v, L_u, method=method, v_grid=v_grid)


def sample_from_grid(grid, density, k, rng):
    """K inverse-CDF draws from a density tabulated on a grid."""
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(grid))])
    cdf /= cdf[-1]
    return np.interp(rng.uniform(size=k), cdf, grid)


def grid_quantile(grid, density, q):
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(grid))])
    return np.interp(q, cdf / cdf[-1], grid)


def inclination_upper_limit_deg(grid, density, level=0.95):
    """i such that P(i' < i) = level, i.e. cos i at the (1 - level) quantile."""
    return np.rad2deg(np.arccos(grid_quantile(grid, density, 1.0 - level)))


# ------------------------------------------------------------- f_cosI table
def _row(kappa):
    return f_cosi(Z_GRID, kappa)


def fcosi_table(cache=None, kappa_grid=KAPPA_GRID, processes=8):
    """f_cosI(z|κ) on Z_GRID × kappa_grid, cached to an .npz file."""
    if cache is not None and Path(cache).exists():
        d = np.load(cache)
        if np.array_equal(d["kappa"], kappa_grid) and np.array_equal(d["z"], Z_GRID):
            return d["table"]
    with Pool(processes) as pool:
        table = np.array(pool.map(_row, list(kappa_grid)))
    if cache is not None:
        np.savez_compressed(cache, kappa=kappa_grid, z=Z_GRID, table=table)
    return table


def per_star_log_likelihood(samples, table):
    """Each star's term log[(1/K) Σ_k f_κ(c_nk) / π_0(c_nk)], shape (N_stars, N_κ).
    samples: (N_stars, K) cos i draws taken under the uniform prior, π_0 = 1."""
    samples = np.atleast_2d(samples)
    out = np.zeros((samples.shape[0], table.shape[0]))
    for j, f_row in enumerate(table):
        f = np.interp(samples, Z_GRID, f_row)                # f_κ(c_nk)
        out[:, j] = np.log(np.mean(f / 1.0, axis=1))         # π_0(c_nk) = 1
    return out


def log_likelihood_kappa(samples, table):
    """Hogg et al. (2010) importance-sampled log L(κ) for each κ in the table."""
    return per_star_log_likelihood(samples, table).sum(axis=0)


def kappa_posterior(loglike, kappa_grid=KAPPA_GRID):
    logp = loglike + log_kappa_prior(kappa_grid)
    p = np.exp(logp - logp.max())
    return p / np.trapz(p, kappa_grid)


def summarize_kappa(p, kappa_grid=KAPPA_GRID):
    """Mode, median, and two 68% intervals (around the mode by highest density,
    and 16th–84th percentiles)."""
    mode = kappa_grid[np.argmax(p)]
    q16, q50, q84 = (grid_quantile(kappa_grid, p, q) for q in (0.16, 0.5, 0.84))
    # highest-density 68% interval on a fine resampling
    fine = np.linspace(kappa_grid[0], kappa_grid[-1], 400_001)
    pf = np.interp(fine, kappa_grid, p)
    order = np.argsort(pf)[::-1]
    mass = np.cumsum(pf[order]) * (fine[1] - fine[0])
    inside = order[: np.searchsorted(mass, 0.68) + 1]
    return dict(mode=mode, hpd68_lo=fine[inside].min(), hpd68_hi=fine[inside].max(),
                median=q50, q16=q16, q84=q84)


# ------------------------------------------------- comparing κ posteriors
def tvd_mw14(p_a, p_b):
    """MW14 Eq. 15 as printed: max |p_A(κ) − p_B(κ)| (a sup-norm on densities)."""
    return np.max(np.abs(p_a - p_b))


def hellinger_sq(p_a, p_b, kappa_grid=KAPPA_GRID):
    """MW14 Eq. 16: 1 − ∫ sqrt(p_A p_B) dκ."""
    return 1.0 - np.trapz(np.sqrt(p_a * p_b), kappa_grid)


def permutation_test(ll_star, n_a, n_trials, rng, observed_mask, kappa_grid=KAPPA_GRID):
    """MW14 Sec. 4.3: random splits into groups of n_a and N − n_a stars.
    ll_star: (N, N_κ) per-star log-likelihoods on kappa_grid. observed_mask: the real group A.
    Returns observed metrics, trial metrics, and the fraction of trials at least
    as different as observed."""
    n = ll_star.shape[0]

    def metrics(mask):
        p_a = kappa_posterior(ll_star[mask].sum(0), kappa_grid)
        p_b = kappa_posterior(ll_star[~mask].sum(0), kappa_grid)
        return tvd_mw14(p_a, p_b), hellinger_sq(p_a, p_b, kappa_grid)

    obs = np.array(metrics(observed_mask))
    trials = np.empty((n_trials, 2))
    for t in range(n_trials):
        mask = np.zeros(n, bool)
        mask[rng.choice(n, n_a, replace=False)] = True
        trials[t] = metrics(mask)
    return obs, trials, (trials >= obs).mean(axis=0)


# ------------------------------------------------------- full-sample run
def run_sample(df, table, *, method, vsini_err_scale=1.0, k_samples=1000,
               sample_seed=100, kappa_grid=KAPPA_GRID):
    """Per-star posteriors (seeded by KOI number), K samples each (one shared
    sample seed), and per-star log-likelihoods on kappa_grid. Every run with the
    same seeds sees the same Monte Carlo draws, so runs differ only in the
    method or the inputs."""
    post = {int(r.koi): star_cosi_posterior(r, method=method, vsini_err_scale=vsini_err_scale,
                                            rng=np.random.default_rng(int(r.koi)))
            for _, r in df.iterrows()}
    rng = np.random.default_rng(sample_seed)
    samples = np.array([sample_from_grid(COSI_GRID, post[int(k)], k_samples, rng) for k in df.koi])
    return dict(post=post, samples=samples, ll_star=per_star_log_likelihood(samples, table))


def kl_divergence(p, q, grid=COSI_GRID):
    """KL(p || q) for densities tabulated on the same grid (nats)."""
    ok = p > 0
    return np.trapz(np.where(ok, p * np.log(np.where(ok, p, 1.0) / np.maximum(q, 1e-300)), 0.0), grid)
