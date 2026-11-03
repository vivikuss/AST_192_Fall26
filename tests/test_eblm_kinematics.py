"""Stage D2 statistics on synthetic samples."""
import numpy as np
import pandas as pd

from src.eblm_kinematics import bootstrap, c2, evolved_fraction

CFG = {"seed": 1, "monte_carlo": {"statistic": "mean_abs"}}


def _stars(n, sig, seed, color=(0.5, 1.2)):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({"U": rng.normal(0, sig, n), "V": rng.normal(0, sig, n), "W": rng.normal(0, sig, n),
                         "bp_rp0": rng.uniform(*color, n), "dist_pc": rng.uniform(100, 600, n)})


def test_bootstrap_ci_covers_estimate_and_sign():
    e = _stars(80, 20.0, 1)
    b = bootstrap(e, ref=100.0, cfg=CFG, n=2000, label="t", stat="mean_abs")
    assert b["ci95"][0] < b["delta_hat"] < b["ci95"][1]
    assert b["interpretation"] == "colder than the field"


def test_c2_same_population_is_centred_on_zero():
    e, h = _stars(80, 25.0, 2), _stars(400, 25.0, 3)
    r = c2(e, h, CFG, 2000, "t")
    assert r["ci95"][0] < 0 < r["ci95"][1]
    cold = h.assign(U=h["U"] * 0.5, V=h["V"] * 0.5, W=h["W"] * 0.5)
    assert c2(e, cold, CFG, 500, "t")["delta_median"] > 0


def test_evolved_fraction_reweights_colour():
    # same evolved fraction per colour bin, different colour mix -> delta ~ 0 after reweighting
    rng = np.random.default_rng(4)
    E = pd.DataFrame({"bp_rp0": rng.uniform(0.5, 0.7, 300)})
    H = pd.DataFrame({"bp_rp0": np.r_[rng.uniform(0.5, 0.7, 100), rng.uniform(0.7, 1.0, 900)]})
    for d in (E, H):
        d["evolved"] = rng.random(len(d)) < np.where(d["bp_rp0"] < 0.7, 0.3, 0.0)
    r = evolved_fraction(E, H, CFG, 0.1, 300, "t")
    assert r["raw_f_hj_in_range"] > 0.2 and abs(r["delta_f"]) < 0.12
    assert r["ci95"][0] < 0 < r["ci95"][1]
