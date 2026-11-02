import numpy as np
import pandas as pd

from src import config, cuts

CFG = config.load()["cuts"]


def _row(**kw):
    base = dict(parallax_over_error=20.0, mean_varpi_factor_al=0.0, visibility_periods_used=12,
                astrometric_chi2_al=100.0, astrometric_n_good_obs_al=105, phot_g_mean_mag=12.0,
                astrometric_gof_al=0.0, astrometric_excess_noise_sig=0.0, rv_nb_transits=10,
                phot_bp_rp_excess_factor=1.2, bp_rp=1.0, phot_bp_mean_flux_over_error=100.0,
                phot_rp_mean_flux_over_error=100.0)
    base.update(kw)
    return pd.DataFrame([base])


def test_good_star_passes_all():
    m, _ = cuts.apply(_row(), CFG, range(1, 11))
    assert m.all()


def test_strict_inequalities_at_boundaries():
    # every threshold in HS19 is strict
    for kw, i in [(dict(parallax_over_error=10.0), 1), (dict(mean_varpi_factor_al=-0.23), 2),
                  (dict(mean_varpi_factor_al=0.36), 2), (dict(visibility_periods_used=8), 3),
                  (dict(astrometric_gof_al=3.0), 5), (dict(astrometric_excess_noise_sig=2.0), 6),
                  (dict(rv_nb_transits=5), 7), (dict(phot_bp_mean_flux_over_error=10.0), 9),
                  (dict(phot_rp_mean_flux_over_error=10.0), 10)]:
        assert not cuts.CUTS[i](_row(**kw), CFG)[0], (kw, i)


def test_nan_fails():
    m, masks = cuts.apply(_row(parallax_over_error=np.nan, bp_rp=np.nan), CFG, range(1, 11))
    assert not masks[1][0] and not masks[8][0] and not m[0]


def test_uwe_threshold_forms():
    # Lindegren C.1: 1.2 * max(1, exp(-0.2 (G - 19.5)))
    assert np.isclose(cuts.uwe_threshold(19.5), 1.2)
    assert np.isclose(cuts.uwe_threshold(20.5), 1.2)                   # faint end floors at 1.2
    assert np.isclose(cuts.uwe_threshold(14.5), 1.2 * np.e)            # exp(1)
    assert np.isclose(cuts.uwe_threshold(9.5), 1.2 * np.exp(2.0))
    # literal (as printed, no parentheses): exp(-0.2 G - 19.5) << 1 -> constant 1.2
    for g in (5.0, 12.0, 19.0):
        assert np.isclose(cuts.uwe_threshold(g, form="literal"), 1.2)


def test_cut4_uses_unit_weight_error():
    # u = sqrt(chi2 / (n - 5)); chi2=400, n=105 -> u = 2
    r = _row(astrometric_chi2_al=400.0, astrometric_n_good_obs_al=105, phot_g_mean_mag=12.0)
    assert np.isclose(cuts.unit_weight_error(r)[0], 2.0)
    lind = dict(CFG, uwe=dict(CFG["uwe"], form="lindegren"))
    lit = dict(CFG, uwe=dict(CFG["uwe"], form="literal"))
    assert cuts.cut_4(r, lind)[0]          # threshold 1.2*exp(1.5) = 5.4
    assert not cuts.cut_4(r, lit)[0]       # threshold 1.2


def test_cut8_bounds():
    lo, hi = cuts.bp_rp_excess_bounds(2.0)          # 1 + 0.015*4, 1.3 + 0.06*4
    assert np.isclose(lo, 1.06) and np.isclose(hi, 1.54)
    lo_p, _ = cuts.bp_rp_excess_bounds(2.0, lo_b=0.0015)
    assert np.isclose(lo_p, 1.006)
    assert not cuts.cut_8(_row(bp_rp=2.0, phot_bp_rp_excess_factor=1.05), CFG)[0]
    assert cuts.cut_8(_row(bp_rp=2.0, phot_bp_rp_excess_factor=1.07), CFG)[0]
    assert not cuts.cut_8(_row(bp_rp=2.0, phot_bp_rp_excess_factor=1.55), CFG)[0]


def test_server_where_is_superset_of_local():
    # the server prefilter must never remove a star the local cut set keeps
    w = cuts.server_where(CFG)
    assert "0.0015" in w and "astrometric_chi2_al" not in w
