"""Stage 0 validation tests. Runs under pytest, or without it via
`python tests/run_tests.py` (pytest is not in the ast192 environment)."""

import sys
from pathlib import Path

import numpy as np
from scipy import integrate

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.inclination import (flat_upper_limit_likelihood, gaussian_likelihood,
                             incorrect_divide_and_clip, incorrect_product_rule,
                             mcmc_posterior_cosi, posterior_cosi)
from src.obliquity import f_cosi, log_kappa_over_sinh, sample_cosi_given_kappa

COSI = np.linspace(0.0, 1.0, 1001)
L_V = gaussian_likelihood(10.0, 2.0)
MW20_PANELS = {  # Masuda & Winn (2020) Figure 1
    "vsini=9.8+-2.0": gaussian_likelihood(9.8, 2.0),
    "vsini=4.0+-2.0": gaussian_likelihood(4.0, 2.0),
    "vsini<5": flat_upper_limit_likelihood(5.0),
    "vsini<30": flat_upper_limit_likelihood(30.0),
}


def _max_cdf_distance(samples, density):
    cdf = integrate.cumulative_trapezoid(density, COSI, initial=0.0)
    ecdf = np.searchsorted(np.sort(samples), COSI, side="right") / len(samples)
    return np.max(np.abs(ecdf - cdf))


# ---------------------------------------------------------------- f_cosI
def test_f_cosi_normalized():
    for kappa in [0.0, 1e-3, 0.1, 1, 5, 20, 100, 500, 1000, 5000]:
        pts = [min(0.99, s / np.sqrt(max(kappa, 1.0))) for s in (1, 3, 10)]
        norm, _ = integrate.quad(lambda z: f_cosi(z, kappa), 0, 1, points=pts, limit=200)
        assert abs(norm - 1.0) < 1e-6, (kappa, norm)


def test_f_cosi_isotropic_limit():
    assert np.allclose(f_cosi(np.linspace(0, 0.99, 50), 0.0), 1.0, atol=1e-8)


def test_f_cosi_endpoint_limit():
    for kappa in [0.0, 1.0, 5.0, 20.0]:
        expected = 1.0 if kappa == 0 else kappa / np.sinh(kappa)
        assert np.isclose(f_cosi(1.0, kappa), expected)
        assert np.isclose(f_cosi(1.0 - 1e-7, kappa), expected, rtol=1e-3)


def test_f_cosi_matches_monte_carlo():
    """Residuals vs. direct Monte Carlo are consistent with Poisson noise
    (chi^2/dof near 1, no mean offset) when compared to bin-integrated f."""
    edges = np.linspace(0, 1, 41)
    n = 1_000_000
    for kappa in [1, 5, 20, 100]:
        counts, _ = np.histogram(sample_cosi_given_kappa(kappa, n, np.random.default_rng(7)), edges)
        expect = n * np.array([integrate.quad(lambda z: f_cosi(z, kappa), a, b)[0]
                               for a, b in zip(edges[:-1], edges[1:])])
        use = expect > 20
        z = (counts[use] - expect[use]) / np.sqrt(expect[use])
        assert np.mean(z ** 2) < 1.6 and abs(np.mean(z)) < 3 / np.sqrt(use.sum()), kappa


def test_kappa_over_sinh_does_not_overflow():
    vals = log_kappa_over_sinh(np.array([0.0, 1e-9, 1.0, 700.0, 710.0, 1e4, 1e6]))
    assert np.all(np.isfinite(vals))
    assert np.isclose(vals[2], np.log(1.0 / np.sinh(1.0)))


# ------------------------------------------------------ cos i posteriors
def test_flat_prior_recovered_for_uninformative_vsini():
    """MW20 Fig. 1 bottom right: the single best test of the implementation."""
    post = posterior_cosi(COSI, L_V, MW20_PANELS["vsini<30"], method="correct")
    assert np.allclose(post, 1.0, atol=1e-3)
    samples = mcmc_posterior_cosi(L_V, MW20_PANELS["vsini<30"], seed=3)
    assert _max_cdf_distance(samples, np.ones_like(COSI)) < 0.02


def test_mcmc_agrees_with_analytic_all_panels():
    for seed, (name, L_u) in enumerate(MW20_PANELS.items()):
        post = posterior_cosi(COSI, L_V, L_u, method="correct")
        samples = mcmc_posterior_cosi(L_V, L_u, seed=seed)
        d = _max_cdf_distance(samples, post)
        assert d < 0.02, (name, d)


def test_method_must_be_explicit():
    try:
        posterior_cosi(COSI, L_V, L_V)
    except TypeError:
        pass
    else:
        raise AssertionError("posterior_cosi accepted a call with no method")
    try:
        posterior_cosi(COSI, L_V, L_V, method="Correct")
    except ValueError:
        pass
    else:
        raise AssertionError("posterior_cosi accepted an unknown method")


def test_mw14_bias_direction():
    """MW20 Sec. 5: Eq. 18 over-weights small cos i and under-weights large cos i."""
    for L_u in (MW20_PANELS["vsini=9.8+-2.0"], MW20_PANELS["vsini<30"]):
        good = posterior_cosi(COSI, L_V, L_u, method="correct")
        bad = posterior_cosi(COSI, L_V, L_u, method="mw14")
        assert bad[0] > good[0] and bad[900] < good[900]


def test_values_read_off_mw20_figure1_top_left():
    """Densities at cos i = 0 read by eye from the published top-left panel
    (so only ±0.1): correct ~1.35, Eq.18 ~1.5, product-rule peak ~1.95
    near cos i ~ 0.65."""
    L_u = MW20_PANELS["vsini=9.8+-2.0"]
    assert abs(posterior_cosi(COSI, L_V, L_u, method="correct")[0] - 1.35) < 0.1
    assert abs(posterior_cosi(COSI, L_V, L_u, method="mw14")[0] - 1.5) < 0.1
    prod = incorrect_product_rule(COSI, L_V, L_u)
    assert abs(np.nanmax(prod) - 1.95) < 0.1 and 0.6 < COSI[np.nanargmax(prod)] < 0.72


def test_mw14_figure2():
    """MW14 Fig. 2: V sin I = V = 10 ± 0.5. Their Eq. 13 curve is flat at ~2.6
    for cos I < 0.2; divide-and-clip puts ~50% of samples at cos I = 0."""
    L = gaussian_likelihood(10.0, 0.5)
    mw14 = posterior_cosi(COSI, L, L, method="mw14")
    assert abs(mw14[0] - 2.6) < 0.15 and abs(mw14[200] - mw14[0]) < 0.15
    rng = np.random.default_rng(0)
    c = incorrect_divide_and_clip(rng.normal(10, 0.5, 200_000), rng.normal(10, 0.5, 200_000))
    assert 0.45 < np.mean(c == 0.0) < 0.55
