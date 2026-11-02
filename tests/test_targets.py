import numpy as np
import pandas as pd

from src import targets


def _t(**kw):
    base = dict(is_binary=[True, True, True, False, False], rv=[5.0, 7.0, np.nan, 3.0, np.nan],
                rv_kind=["gamma", "systemic", "", "", ""], radial_velocity=[50.0, 50.0, 50.0, 9.0, 11.0])
    base.update(kw)
    return pd.DataFrame(base)


def test_binary_rv_only_from_gamma():
    t = targets.assign_rv(_t())
    rv = t["rv_used_input"].to_numpy()
    assert rv[0] == 5.0                           # binary with orbit-fit gamma
    assert np.isnan(rv[1]) and np.isnan(rv[2])    # binary: non-gamma value or nothing -> no RV
    assert list(t["rv_source"][:3]) == ["gamma", "rejected_not_gamma", "none_binary_no_gamma"]
    assert list(t["flag_binary_no_gamma"]) == [False, True, True, False, False]


def test_binary_never_uses_gaia_rv():
    t = targets.assign_rv(_t(rv=[np.nan] * 5, rv_kind=[""] * 5))
    assert np.isnan(t["rv_used_input"][:3]).all()       # Gaia RV present but ignored for binaries


def test_single_star_override_then_gaia():
    t = targets.assign_rv(_t())
    assert t["rv_used_input"][3] == 3.0 and t["rv_source"][3] == "override"
    assert t["rv_used_input"][4] == 11.0 and t["rv_source"][4] == "gaia_dr2"
