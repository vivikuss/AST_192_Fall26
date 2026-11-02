import numpy as np
import pandas as pd

from src import config, matching

CFG = config.load()


def _mock(n_field=20000, n_host=150, seed=0):
    rng = np.random.default_rng(seed)
    field = pd.DataFrame(dict(z=rng.normal(0, 200, n_field), bp_rp0=rng.uniform(0.4, 2.0, n_field),
                              U=rng.normal(0, 40, n_field), V=rng.normal(-20, 25, n_field),
                              W=rng.normal(0, 20, n_field)))
    host = pd.DataFrame(dict(z=rng.normal(50, 80, n_host), bp_rp0=rng.uniform(0.5, 1.2, n_host),
                             U=rng.normal(0, 30, n_host), V=rng.normal(-15, 20, n_host),
                             W=rng.normal(0, 15, n_host)))
    return host, field


def test_histogram_weights_reproduce_host_distribution():
    host, field = _mock()
    edges = [matching.bin_edges(host["z"], {"bins": "fixed", "bin_width": 50.0, "n_quantile_bins": 0})]
    w, info = matching.histogram_weights(host[["z"]].to_numpy(), field[["z"]].to_numpy(), edges)
    hh = np.histogram(host["z"], edges[0])[0] / len(host)
    fw = np.histogram(field["z"], edges[0], weights=w)[0]
    assert np.allclose(fw / fw.sum(), hh, atol=1e-12)
    assert (w[(field["z"] < edges[0][0]) | (field["z"] > edges[0][-1])] == 0).all()


def test_controls_are_colour_matched_and_unique():
    host, field = _mock()
    s = matching.ControlSampler(host, field, CFG)
    rng = np.random.default_rng(0)
    for _ in range(20):
        pick = s.control(rng)
        assert (pick >= 0).all()
        assert np.all(np.abs(s.fc[pick] - s.hc) <= CFG["monte_carlo"]["color_tol"] + 1e-12)
        assert len(np.unique(pick)) == len(pick)                  # without replacement


def test_zmatched_controls_follow_host_z():
    host, field = _mock(n_field=60000, n_host=300)
    fz = field["z"].to_numpy()
    s = matching.ControlSampler(host, field, CFG)
    rng = np.random.default_rng(1)
    zc = np.concatenate([fz[s.field_idx[s.control(rng)]] for _ in range(200)])
    # control z is pulled from the field's N(0,200) to the hosts' N(50,80)
    assert abs(np.median(zc) - np.median(host["z"])) < 15
    assert abs(np.std(zc) - np.std(host["z"])) < 15


def test_null_hosts_drawn_from_field_give_uniform_p():
    # if hosts are a random field subsample, p-values should not be extreme
    _, field = _mock(n_field=30000)
    rng = np.random.default_rng(5)
    host = field.sample(150, random_state=5)
    fld = field.drop(host.index)
    s = matching.ControlSampler(host, fld, CFG)
    res = s.run(400, rng)
    from src.kinematics import dispersion
    p, _ = matching.p_value(dispersion(host.U, host.V, host.W), res["dispersion"])
    assert 0.01 < p < 0.99


def test_weighted_resample_method_runs():
    host, field = _mock()
    cfg = config.variant(CFG, **{"monte_carlo.z_match.method": "weighted_resample"})
    s = matching.ControlSampler(host, field, cfg)
    res = s.run(50, np.random.default_rng(2))
    assert np.isfinite(res["dispersion"]).all()


def test_unmatched_host_dropped_from_both_sides():
    """A host with no field star of its colour is excluded from the host dispersion
    of every iteration in which it is unmatched (paired comparison, D13b)."""
    host, field = _mock()
    host.loc[0, ["bp_rp0", "U", "V", "W"]] = [3.0, 500.0, 500.0, 500.0]   # unmatchable outlier
    s = matching.ControlSampler(host, field, CFG)
    res = s.run(50, np.random.default_rng(2), "mean_abs", chunk=25)
    assert (res["n_unmatched"] >= 1).all()
    assert res["host_unmatched_count"][0] == 50
    expect = matching.dispersion(host.U[1:], host.V[1:], host.W[1:], "mean_abs")
    assert np.allclose(res["host_dispersion"], expect)
    assert res["host_dispersion_full"]["mean_abs"] > expect                # outlier only in the full value
    assert (matching.candidates(host, field, CFG["monte_carlo"]["color_tol"])[0] == 0)
