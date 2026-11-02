import numpy as np

from src import config, main_sequence as msq

CFG = config.load()
A = CFG["main_sequence"]["coeffs"]


def test_published_coefficients():
    assert len(A) == 9                                  # 8th order
    assert np.isclose(msq.ms_poly(0.0, A), 1.07857)
    # sum of coefficients = value at c = 1
    assert np.isclose(msq.ms_poly(1.0, A), sum(A))


def test_ms_is_monotonic_over_fgk_range():
    c = np.linspace(0.3, 2.5, 500)
    assert np.all(np.diff(msq.ms_poly(c, A)) > 0)       # redder -> fainter


def test_sign_convention():
    c = np.array([0.8, 0.8, 0.8, 0.8])
    ms = msq.ms_poly(0.8, A)
    M = np.array([ms - 0.5, ms - 1.5, ms + 0.7, ms - 0.999])   # brighter = smaller M
    keep = msq.on_main_sequence(M, c, CFG)
    assert keep.tolist() == [True, False, True, True]


def test_refit_recovers_known_polynomial():
    rng = np.random.default_rng(1)
    c = np.sort(rng.uniform(0.0, 3.5, 600))
    M = msq.ms_poly(c, A) + rng.normal(0, 0.05, c.size)
    coeffs, spl, grid = msq.refit(c, M, s=c.size * 0.05**2)
    g = np.linspace(0.2, 3.3, 200)
    assert np.max(np.abs(msq.ms_poly(g, coeffs) - msq.ms_poly(g, A))) < 0.1
