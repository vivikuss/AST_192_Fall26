"""Hu et al. (2026) checks: sigma_W definition and the power calculation."""
import numpy as np

from src.hu2026_checks import _stats, power


def test_sigma_w_is_rms_with_1_over_n():
    W = np.array([[1.0, 2.0, 3.0, 6.0]])
    _, sw = _stats(np.zeros_like(W), np.zeros_like(W), W)
    assert np.isclose(sw[0], np.sqrt(np.mean((W - W.mean()) ** 2)))
    assert not np.isclose(sw[0], np.std(W, ddof=1))


def test_power_null_is_alpha_and_grows_with_offset():
    X = np.random.default_rng(1).normal(40.0, 2.0, 200000)
    assert abs(power(X, 1.0) - 0.05) < 0.005
    assert power(X, 0.9) < power(X, 0.85) < power(X, 0.8)
