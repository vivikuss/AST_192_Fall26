"""gamma fit (Stage D): simulated data only."""
import numpy as np

from src.gamma_fit import PERMITTED, drift, fit_gamma, rv_shape, tperi_from_t0, true_anomaly, _err


def _sim(e, omega, n=25, seed=1, gamma=12.345, K=18.0, P=5.3, t0=2456000.2):
    rng = np.random.default_rng(seed)
    t = np.sort(2455000 + rng.uniform(0, 2000, n))
    tp = tperi_from_t0(t0, P, e, omega)
    rv = gamma + K * rv_shape(t, P, tp, e, omega) + rng.normal(0, 0.02, n)
    return t, rv, np.full(n, 0.02), rv_shape(t, P, tp, e, omega)


def test_circular_reduces_to_minus_sine():
    t = np.linspace(0, 20, 50)
    P, t0 = 3.0, 1.2
    tp = tperi_from_t0(t0, P, 0.0, 37.0)
    assert np.allclose(rv_shape(t, P, tp, 0.0, 37.0), -np.sin(2 * np.pi * (t - t0) / P), atol=1e-9)


def test_eclipse_at_conjunction_eccentric():
    P, t0, e, w = 7.0, 3.0, 0.4, 120.0
    nu = true_anomaly(np.array([t0]), P, tperi_from_t0(t0, P, e, w), e)
    assert np.isclose(np.mod(nu[0] + np.radians(w), 2 * np.pi), np.pi / 2, atol=1e-9)


def test_fit_recovers_gamma_and_K():
    for e, w in [(0.0, 0.0), (0.3, 75.0)]:
        t, rv, err, shape = _sim(e, w)
        r = fit_gamma(t, rv, err, shape, np.zeros_like(t), 0.0, 0.020)
        assert r["flag"] == "ok"
        assert abs(r["gamma"] - 12.345) < 0.02 and abs(r["K"] - 18.0) < 0.03
        assert 0.020 < r["sigma_gamma"] < 0.03


def test_fit_returns_only_permitted_outputs():
    t, rv, err, shape = _sim(0.1, 10.0)
    assert set(fit_gamma(t, rv, err, shape, np.zeros_like(t), 0.0, 0.02)) == set(PERMITTED)
    assert fit_gamma(t[:2], rv[:2], err[:2], shape[:2], np.zeros(2), 0.0, 0.02)["flag"] == "too_few_rv"
    assert fit_gamma(t, rv, err, np.zeros_like(t), np.zeros_like(t), 0.0, 0.02)["flag"] == "singular"


def test_fixed_drift_subtracted():
    t, rv, err, shape = _sim(0.0, 0.0)
    tm = 0.5 * (t.min() + t.max())
    d = drift(t, 150.0, -20.0, 3.0, tm)
    r = fit_gamma(t, rv + d, err, shape, d, 0.0, 0.02)
    assert abs(r["gamma"] - 12.345) < 0.02


def test_published_error_parsing():
    assert np.isclose(_err("16.6460(60)"), 0.0060)
    assert np.isclose(_err("23.4(2.2)"), 2.2)
    assert np.isclose(_err("3062(91)"), 91.0)
