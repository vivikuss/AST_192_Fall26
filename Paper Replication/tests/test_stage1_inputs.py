"""Guards on the Stage 1 inputs: Table 1 bookkeeping and the MW14 likelihood."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.data import check_vsini_not_double_corrected, load_mw14_table
from src.inclination import (gaussian_likelihood, marginal_likelihood_correct,
                             marginal_likelihood_mw14, marginal_likelihood_mw14_uspace)
from src.velocity import vsini_likelihood


def test_multiplicity_mapping():
    df = load_mw14_table(p_rot_err_if_zero=0.005)
    assert (df.loc[df.n_cands == 1, "group"] == "single").all()
    assert (df.loc[df.n_cands >= 2, "group"] == "multi").all()
    assert df.group.value_counts().to_dict() == {"single": 45, "multi": 25}


def test_vsini_used_once_corrected_koi975():
    df = load_mw14_table(p_rot_err_if_zero=0.005)
    row = df.loc[df.koi == 975].iloc[0]
    L_u = vsini_likelihood(row)
    u = np.linspace(0, 20, 20001)
    assert abs(u[np.argmax(L_u(u))] - 8.47) < 1e-3  # not 9.57, not 7.5


def test_double_correction_is_caught():
    df = load_mw14_table(p_rot_err_if_zero=0.005)
    df.loc[:, "vsini"] = df.vsini / (1 - 0.23 / 2)
    try:
        check_vsini_not_double_corrected(df)
    except AssertionError:
        return
    raise AssertionError("a second (1 - α/2) division went undetected")


def test_p_rot_floor_only_touches_koi63():
    a = load_mw14_table(p_rot_err_if_zero=0.005)
    b = load_mw14_table(p_rot_err_if_zero=0.05)
    diff = a.koi[a.P_rot_err != b.P_rot_err].tolist()
    assert diff == [63]
    assert a.loc[a.koi == 63, "P_rot_err"].item() == 0.005


def test_mw14_vspace_equals_uspace():
    """The v-space evaluation of MW14 Eq. 13 is an exact change of variables."""
    L_v, L_u = gaussian_likelihood(10, 2), gaussian_likelihood(4, 2)
    c = np.linspace(0, 0.95, 40)
    v_grid = np.linspace(0, 40, 40001)
    a = marginal_likelihood_mw14(c, L_v, L_u, v_grid=v_grid)
    b = marginal_likelihood_mw14_uspace(c, L_v, L_u, u_grid=v_grid)
    assert np.allclose(a, b, rtol=1e-4)


def test_mw14_is_sini_times_correct():
    """With uniform P_v, Eq. 13 = sin i × Eq. 9: the missing Jacobian (MW20)."""
    L_v, L_u = gaussian_likelihood(10, 2), gaussian_likelihood(9.8, 2)
    c = np.linspace(0, 1, 101)
    a = marginal_likelihood_mw14(c, L_v, L_u)
    b = marginal_likelihood_correct(c, L_v, L_u) * np.sqrt(1 - c ** 2)
    assert np.allclose(a, b)
