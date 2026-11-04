"""WASP-only pre-registration: decision rule (section B) and the sigma_W statistic."""
import numpy as np

from src.wasp_kinematics import DELTA_HS19, case, sigma_w


def test_delta_hs19_is_the_frozen_value():
    assert round(DELTA_HS19, 2) == -6.64


def test_each_case_of_the_decision_rule():
    d = DELTA_HS19
    assert case(-9.0, -3.0) == 1          # below 0, contains delta
    assert case(-12.0, -7.0) == 2         # entirely below delta
    assert case(-5.0, -1.0) == 3          # between delta and 0, touching neither
    assert case(-4.0, 3.0) == 4           # contains 0, excludes delta
    assert case(-8.0, 2.0) == 5           # contains both
    assert case(0.5, 4.0) == 6            # entirely above 0
    assert case(d, -1.0) == 1 and case(-3.0, 0.0) == 4 and case(d, 0.0) == 5   # closed ends


def test_rule_is_exhaustive_and_exclusive():
    rng = np.random.default_rng(0)
    for lo, w in zip(rng.uniform(-15, 5, 2000), rng.uniform(0, 12, 2000)):
        assert case(lo, lo + w) in range(1, 7)


def test_sigma_w_rms_1_over_n():
    W = np.array([1.0, 2.0, 3.0, 6.0])
    assert np.isclose(sigma_w(None, None, W), np.sqrt(np.mean((W - W.mean()) ** 2)))
