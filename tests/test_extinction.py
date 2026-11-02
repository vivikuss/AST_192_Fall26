import numpy as np

from src import config, extinction as ex

CFG = config.load()


def _ebv(l, b, d, cfg=CFG):
    return ex.ebv_los(np.atleast_1d(l), np.atleast_1d(b), np.atleast_1d(d), cfg)[0]


def test_cube_orientation_known_clouds():
    # Ophiuchus is north of the plane, Taurus south; mirrored directions are clean
    assert _ebv(353, 17, 200) > 5 * _ebv(353, -17, 200)
    assert _ebv(172, -15, 250) > 3 * _ebv(172, 15, 250)
    assert _ebv(172, -15, 250) > 3 * _ebv(278, -15, 250)     # x<->y swap check


def test_pleiades_reddening():
    # Gaia Collab. 2018a adopt E(B-V) = 0.045 for the Pleiades
    assert abs(_ebv(166.57, -23.52, 136.0) - 0.045) < 0.015


def test_monotonic_in_distance_and_zero_at_origin():
    d = np.array([1.0, 50, 100, 200, 400, 800, 1600])
    e = ex.ebv_los(np.full(7, 30.0), np.full(7, 2.0), d, CFG)
    assert np.all(np.diff(e) >= -1e-9) and e[0] < 1e-3


def test_step_convergence():
    l, b, d = np.array([30.0, 120, 250]), np.array([2.0, -5, 10]), np.array([800.0, 500, 300])
    fine = ex.ebv_los(l, b, d, config.variant(CFG, **{"extinction.step_pc": 0.5}))
    coarse = ex.ebv_los(l, b, d, CFG)
    assert np.allclose(fine, coarse, atol=0.005)


def test_casagrande_vandenberg_coefficients():
    ag, ebprp = ex.gaia_extinction(1.0, CFG)
    assert np.isclose(ag, 2.740) and np.isclose(ebprp, 3.374 - 2.035)
