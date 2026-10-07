"""Fisher obliquity distribution and the induced distribution of cos i.

Morton & Winn (2014, ApJ 796, 47; "MW14") Section 2.
    θ  obliquity (star spin vs. planet orbit normal), Fisher(κ)
    i  stellar inclination to the line of sight;  cos i = sin θ cos φ  (Eq. 5)
"""

import numpy as np
from scipy import integrate


def log_kappa_over_sinh(kappa):
    """log(κ / sinh κ), stable for κ -> 0 and for κ >~ 710 where sinh overflows."""
    kappa = np.asarray(kappa, dtype=float)
    small = kappa < 1e-6
    k = np.where(small, 1.0, kappa)
    big = np.log(k) - (k + np.log1p(-np.exp(-2.0 * k)) - np.log(2.0))
    return np.where(small, -kappa ** 2 / 6.0, big)


def log_kappa_prior(kappa):
    """Fabrycky & Winn (2009) Eq. 3, unnormalized: π(κ) ∝ (1 + κ²)^(-3/4)."""
    return -0.75 * np.log1p(np.asarray(kappa, dtype=float) ** 2)


def fisher_theta_pdf(theta, kappa):
    """MW14 Eq. 1: f(θ|κ) = κ/(2 sinh κ) exp(κ cos θ) sin θ, in log-safe form."""
    theta = np.asarray(theta, dtype=float)
    return 0.5 * np.exp(log_kappa_over_sinh(kappa) + kappa * np.cos(theta)) * np.sin(theta)


def sample_fisher_theta(kappa, n, rng):
    """Exact inverse-CDF draws of θ from the Fisher distribution."""
    u = rng.uniform(size=n)
    if kappa < 1e-8:
        cos_t = 2.0 * u - 1.0
    else:
        cos_t = 1.0 + np.log(u + (1.0 - u) * np.exp(-2.0 * kappa)) / kappa
    return np.arccos(np.clip(cos_t, -1.0, 1.0))


def sample_cosi_given_kappa(kappa, n, rng):
    """Direct Monte Carlo of cos i = |sin θ cos φ| with θ ~ Fisher(κ), φ ~ U(0, 2π)."""
    theta = sample_fisher_theta(kappa, n, rng)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n)
    return np.abs(np.sin(theta) * np.cos(phi))


def _f_cosi_scalar(z, kappa):
    log_pref = log_kappa_over_sinh(kappa)
    if z >= 1.0:
        # limit z -> 1: the integral tends to π g(1), giving f(1|κ) = κ / sinh κ
        return float(np.exp(log_pref))

    # cosh(κ w) κ/sinh κ with w = sqrt(1-y^2), evaluated in log space
    def cosh_term(y):
        w = np.sqrt(max(1.0 - y * y, 0.0))
        return np.exp(log_pref + kappa * w + np.log1p(np.exp(-2.0 * kappa * w)) - np.log(2.0))

    # Integrand of MW14 Eq. 11 written as  g(y) (y-z)^(-1/2) (1-y)^(-1/2)  with
    # g smooth; the two endpoint singularities go into QUADPACK's QAWS weights.
    def g(y):
        if y + z <= 0.0:          # z = 0 and y = 0: y / sqrt(y) -> 0
            return 0.0
        return cosh_term(y) * y / (np.sqrt(1.0 + y) * np.sqrt(y + z))

    # Split so each piece carries one singular endpoint; for large κ the mass
    # sits within ~1/sqrt(κ) of y = 0, so keep the first piece narrow then.
    m = z + min(0.5 * (1.0 - z), 5.0 / np.sqrt(max(kappa, 1.0)))
    opts = dict(weight="alg", limit=200, epsabs=0.0, epsrel=1e-10)
    a, _ = integrate.quad(lambda y: g(y) / np.sqrt(1.0 - y), z, m, wvar=(-0.5, 0.0), **opts)
    b, _ = integrate.quad(lambda y: g(y) / np.sqrt(y - z), m, 1.0, wvar=(0.0, -0.5), **opts)
    return 2.0 / np.pi * (a + b)


def f_cosi(z, kappa):
    """MW14 Eq. 11: probability density of cos i on [0, 1] given κ.
        f(z|κ) = 2κ/(π sinh κ) ∫_z^1 cosh(κ√(1−y²))/√(1−y²) · 1/√(1−(z/y)²) dy
    """
    z = np.asarray(z, dtype=float)
    out = np.array([_f_cosi_scalar(zz, kappa) for zz in z.ravel()])
    return out.reshape(z.shape)
