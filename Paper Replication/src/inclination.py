"""Per-star inference of cos(i) from v = 2*pi*R/P_rot and u = v sin i.

Notation follows Masuda & Winn (2020, AJ 159, 81; "MW20"):
    L_v(v)  likelihood for the equatorial velocity v
    L_u(u)  likelihood for the projected velocity u = v sin i
    P_v(v)  prior on v (uniform here)
A uniform prior on cos i over [0, 1] is assumed everywhere and applied exactly
once: since it is constant, the posterior on a cos i grid is simply the
normalized marginal likelihood.

Two marginal likelihoods are central to this project and are named explicitly:
    marginal_likelihood_correct  -- MW20 Eq. 9
    marginal_likelihood_mw14     -- Morton & Winn (2014) Eq. 13 = MW20 Eq. 18
The other "incorrect_*" functions exist only to reproduce MW20 Figure 1.
(MW20 Eqs. 20 and 21 are deliberately not implemented: their form could not
be read reliably from the PDF.)
"""

import numpy as np
from scipy import stats

DEFAULT_V_GRID = np.linspace(0.0, 100.0, 20001)  # km/s, support of uniform P_v
METHODS = ("correct", "mw14")


# ---------------------------------------------------------------- likelihoods
def gaussian_likelihood(mu, sigma):
    """Gaussian likelihood in a velocity, as a vectorized callable."""
    return lambda x: stats.norm.pdf(x, loc=mu, scale=sigma)


def flat_upper_limit_likelihood(upper):
    """MW20 Fig. 1 toy upper limit: constant for 0 <= u < upper, zero above."""
    return lambda x: ((np.asarray(x) >= 0) & (np.asarray(x) < upper)).astype(float)


def _sin_from_cos(cosi):
    return np.sqrt(np.clip(1.0 - np.asarray(cosi, dtype=float) ** 2, 0.0, 1.0))


def _normalize(cosi_grid, density):
    density = np.asarray(density, dtype=float)
    return density / np.trapz(density, cosi_grid)


# ------------------------------------------------------ marginal likelihoods
def marginal_likelihood_correct(cosi, L_v, L_u, v_grid=DEFAULT_V_GRID, P_v=None):
    """MW20 Eq. 9:  p(D | cos i) ∝ ∫ L_v(v) L_u(v sqrt(1-cos^2 i)) P_v(v) dv.

    The same v enters both likelihoods. P_v defaults to uniform on v_grid.
    """
    sini = _sin_from_cos(np.atleast_1d(cosi))
    prior = np.ones_like(v_grid) if P_v is None else P_v(v_grid)
    weight = L_v(v_grid) * prior                                  # (Nv,)
    integrand = weight[None, :] * L_u(sini[:, None] * v_grid[None, :])
    return np.trapz(integrand, v_grid, axis=1)


def marginal_likelihood_mw14(cosi, L_v, L_u, v_grid=DEFAULT_V_GRID):
    """Morton & Winn (2014) Eq. 13 = MW20 Eq. 18 (INCORRECT):
        p(D | cos i) ∝ ∫ L_v(u / sqrt(1-cos^2 i)) L_u(u) du.
    Missing the 1/sin i Jacobian and the P_v weighting relative to Eq. 9.

    Evaluated after the exact substitution u = v sin i (du = sin i dv):
        = sin i ∫ L_v(v) L_u(v sin i) dv,
    which keeps the integrand resolved on the v grid as sin i -> 0. The direct
    u-space form is marginal_likelihood_mw14_uspace (used to test this one).
    """
    sini = _sin_from_cos(np.atleast_1d(cosi))
    integrand = L_v(v_grid)[None, :] * L_u(sini[:, None] * v_grid[None, :])
    return sini * np.trapz(integrand, v_grid, axis=1)


def marginal_likelihood_mw14_uspace(cosi, L_v, L_u, u_grid=DEFAULT_V_GRID):
    """MW14 Eq. 13 integrated literally over u. Reference for tests only."""
    sini = np.maximum(_sin_from_cos(np.atleast_1d(cosi)), 1e-300)
    integrand = L_v(u_grid[None, :] / sini[:, None]) * L_u(u_grid)[None, :]
    return np.trapz(integrand, u_grid, axis=1)


def posterior_cosi(cosi_grid, L_v, L_u, *, method, v_grid=DEFAULT_V_GRID):
    """Normalized posterior for cos i on a grid. `method` is required and must be
    'correct' (MW20 Eq. 9) or 'mw14' (MW14 Eq. 13); there is no default.
    P_v is uniform on v_grid in both cases."""
    if method == "correct":
        like = marginal_likelihood_correct(cosi_grid, L_v, L_u, v_grid=v_grid)
    elif method == "mw14":
        like = marginal_likelihood_mw14(cosi_grid, L_v, L_u, v_grid=v_grid)
    else:
        raise ValueError(f"method must be one of {METHODS}, got {method!r}")
    return _normalize(cosi_grid, like)


# ------------------------------------- other incorrect procedures (MW20 Fig. 1)
def incorrect_product_rule(cosi_grid, L_v, L_u, v_grid=DEFAULT_V_GRID):
    """MW20 Sec. 4.2 / Eq. 12: treat v and u as independent, z = sin i = u/v,
    p_z(z) = ∫ p_v(v) p_u(z v) v dv, then change variables to cos i.
    Mass at sin i > 1 is discarded (renormalized over 0 <= sin i <= 1).
    The density diverges as 1/sin i at cos i = 1 whenever p_u(0) > 0; that
    point is returned as nan rather than allowed to swamp the normalization."""
    p_v = L_v(v_grid) / np.trapz(L_v(v_grid), v_grid)
    u_norm = np.trapz(L_u(v_grid), v_grid)

    def p_z(z):
        z = np.atleast_1d(z)
        return np.trapz(p_v[None, :] * L_u(z[:, None] * v_grid[None, :]) / u_norm
                        * v_grid[None, :], v_grid, axis=1)

    z_fine = np.linspace(0.0, 1.0, 4001)
    mass = np.trapz(p_z(z_fine), z_fine)
    cosi_grid = np.asarray(cosi_grid, dtype=float)
    sini = _sin_from_cos(cosi_grid)
    out = np.full(cosi_grid.shape, np.nan)
    ok = sini > 0
    out[ok] = p_z(sini[ok]) * cosi_grid[ok] / sini[ok] / mass
    return out


# ------------------------------------------------------- Monte Carlo methods
def incorrect_independent_mc(v_samples, u_samples):
    """MW20 Sec. 4.1: sample v and u independently, sin i = u/v.
    Samples with sin i outside [0, 1] are discarded. Returns cos i samples."""
    sini = np.asarray(u_samples) / np.asarray(v_samples)
    sini = sini[(sini >= 0) & (sini <= 1)]
    return np.sqrt(1.0 - sini ** 2)


def incorrect_divide_and_clip(v_samples, u_samples):
    """MW14 Sec. 2 / Fig. 2 'intuitive' procedure: sin i = u/v, clipped to <= 1,
    which piles samples up at cos i = 0. Returns cos i samples."""
    sini = np.clip(np.asarray(u_samples) / np.asarray(v_samples), 0.0, 1.0)
    return np.sqrt(1.0 - sini ** 2)


def mcmc_posterior_cosi(L_v, L_u, *, v_max=DEFAULT_V_GRID[-1], v_init=10.0,
                        nwalkers=32, nsteps=6000, burn=1000, seed=0):
    """Independent check of Eq. 10: MCMC over (v, cos i) with uniform priors on
    v in (0, v_max) and cos i in (0, 1), likelihood L_v(v) L_u(v sin i).
    Returns flattened cos i samples."""
    import emcee

    def log_prob(theta):
        v, c = theta[:, 0], theta[:, 1]
        ok = (v > 0) & (v < v_max) & (c >= 0) & (c <= 1)
        out = np.full(len(v), -np.inf)
        with np.errstate(divide="ignore"):
            out[ok] = (np.log(L_v(v[ok]))
                       + np.log(L_u(v[ok] * np.sqrt(1 - c[ok] ** 2))))
        return np.where(np.isnan(out), -np.inf, out)

    rng = np.random.default_rng(seed)
    p0 = np.column_stack([v_init + 0.1 * rng.standard_normal(nwalkers),
                          rng.uniform(0.05, 0.95, nwalkers)])
    # start every walker at finite probability (matters for hard upper limits)
    bad = ~np.isfinite(log_prob(p0))
    while bad.any():
        p0[bad, 1] = rng.uniform(0.05, 0.999, bad.sum())
        bad = ~np.isfinite(log_prob(p0))
    sampler = emcee.EnsembleSampler(nwalkers, 2, log_prob, vectorize=True)
    sampler.random_state = np.random.RandomState(seed).get_state()
    sampler.run_mcmc(p0, nsteps, progress=False)
    return sampler.get_chain(discard=burn, flat=True)[:, 1]
