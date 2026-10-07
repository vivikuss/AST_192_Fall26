"""Tests for the Stage 2 helpers: error scaling, KL divergence, permutation test."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.data import load_mw14_table
from src.hierarchical import (COSI_GRID, KAPPA_GRID, hellinger_sq, kl_divergence,
                              permutation_test, tvd_mw14)
from src.velocity import vsini_likelihood


def test_vsini_error_scale():
    row = load_mw14_table(p_rot_err_if_zero=0.005).iloc[0]
    u = np.linspace(-5, 10, 30001)
    for f in (1.0, 1.35):
        L = vsini_likelihood(row, err_scale=f)(u)
        sd = np.sqrt(np.trapz(L * (u - row.vsini) ** 2, u) / np.trapz(L, u))
        assert abs(sd - f * row.vsini_err) < 1e-3


def test_kl_divergence():
    p = np.ones_like(COSI_GRID)
    q = 2 * (1 + COSI_GRID) / 3
    assert abs(kl_divergence(p, p)) < 1e-12
    # KL(U || 2(1+x)/3) = log(3/2) - (2 log 2 - 1) analytically
    assert abs(kl_divergence(p, q) - (np.log(1.5) - 2 * np.log(2) + 1)) < 1e-6


def test_distance_metrics_identical_and_disjoint():
    a = np.exp(-0.5 * (KAPPA_GRID - 5) ** 2); a /= np.trapz(a, KAPPA_GRID)
    b = np.exp(-0.5 * (KAPPA_GRID - 50) ** 2); b /= np.trapz(b, KAPPA_GRID)
    assert tvd_mw14(a, a) == 0 and abs(hellinger_sq(a, a)) < 1e-9
    assert hellinger_sq(a, b) > 0.999


def test_permutation_null_is_uniform():
    """With exchangeable stars, the observed split is just another split."""
    rng = np.random.default_rng(1)
    ll = rng.normal(0, 0.3, (30, len(KAPPA_GRID))).cumsum(axis=1) * 0.01
    pvals = []
    for seed in range(20):
        mask = np.zeros(30, bool); mask[np.random.default_rng(seed).choice(30, 18, replace=False)] = True
        pvals.append(permutation_test(ll, 18, 100, np.random.default_rng(100 + seed), mask)[2][1])
    assert 0.25 < np.mean(pvals) < 0.75
