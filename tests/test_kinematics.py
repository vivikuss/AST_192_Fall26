import numpy as np
import pandas as pd

from src import kinematics as k

# HD 6755, worked example in the IDL Astronomy Library gal_uvw.pro (Landsman):
# ra=01:09:42.3 dec=+61:32:49.5, pm=(628.42, 76.65) mas/yr, d=139 pc, vr=-321.4 km/s
#   /lsr ===> u=141.2 v=-491.7 w=93.9   (u positive toward the ANTI-centre)
# The printed result predates the 2011 solar-motion update: it equals the
# heliocentric velocity plus the Dehnen & Binney (1998) solar motion
# (10.00, 5.25, 7.17) km/s, i.e. heliocentric (U_GC, V, W) = (-151.2, -496.95, 86.73).
HD6755 = pd.DataFrame(dict(ra=[(1 + 9 / 60 + 42.3 / 3600) * 15], dec=[61 + 32 / 60 + 49.5 / 3600],
                           dist_pc=[139.0], pmra=[628.42], pmdec=[76.65], rv=[-321.4]))


def test_uvw_hd6755_against_published():
    c = k.skycoord(HD6755)
    U, V, W = k.uvw(c, solar_motion_kms=(10.00, 5.25, 7.17))
    assert abs(-U[0] - 141.2) < 0.15            # sign: published U is anticentre-positive
    assert abs(V[0] - (-491.7)) < 0.15
    assert abs(W[0] - 93.9) < 0.15


def _js_matrix_uvw(ra, dec, plx, pmra, pmdec, vr):
    """Independent implementation: literal port of IDL gal_uvw.pro (Johnson &
    Soderblom 1987 with the Hipparcos J2000 matrix). Returns U anticentre-positive."""
    AG = np.array([[0.0548755604, +0.4941094279, -0.8676661490],     # IDL A_G[i, j]: i = column
                   [0.8734370902, -0.4448296300, -0.1980763734],
                   [0.4838350155, 0.7469822445, +0.4559837762]])
    cd, sd, ca, sa = np.cos(np.radians(dec)), np.sin(np.radians(dec)), np.cos(np.radians(ra)), np.sin(np.radians(ra))
    v1, v2, v3 = vr, 4.74047 * pmra / plx, 4.74047 * pmdec / plx
    out = []
    for r in range(3):                     # IDL A_G[r, c] -> AG[c, r] in row-major numpy
        a0, a1, a2 = AG[0, r], AG[1, r], AG[2, r]
        out.append((a0 * ca * cd + a1 * sa * cd + a2 * sd) * v1 + (-a0 * sa + a1 * ca) * v2
                   + (-a0 * ca * sd - a1 * sa * sd + a2 * cd) * v3)
    return out


def test_uvw_matches_independent_matrix_method():
    rng = np.random.default_rng(3)
    n = 20
    df = pd.DataFrame(dict(ra=rng.uniform(0, 360, n), dec=np.degrees(np.arcsin(rng.uniform(-1, 1, n))),
                           dist_pc=rng.uniform(20, 1000, n), pmra=rng.normal(0, 50, n),
                           pmdec=rng.normal(0, 50, n), rv=rng.normal(0, 40, n)))
    U, V, W = k.uvw_heliocentric(k.skycoord(df))
    for i, r in df.iterrows():
        u, v, w = _js_matrix_uvw(r.ra, r.dec, 1e3 / r.dist_pc, r.pmra, r.pmdec, r.rv)
        # IDL U is anticentre-positive; ours points to the Galactic centre
        assert np.allclose([-u, v, w], [U[i], V[i], W[i]], atol=0.02), i


def test_lsr_shift_is_additive():
    c = k.skycoord(HD6755)
    h = np.array(k.uvw_heliocentric(c)).ravel()
    l = np.array(k.uvw(c, (11.1, 12.24, 7.25))).ravel()
    assert np.allclose(l - h, [11.1, 12.24, 7.25], atol=1e-6)


def test_dispersion_formula():
    U = np.array([1.0, -1.0, 3.0, -3.0]); V = np.zeros(4); W = np.zeros(4)
    assert np.isclose(k.dispersion(U, V, W, "rms"), np.sqrt(5.0))
    assert np.isclose(k.dispersion(U, V, W, "mean_abs"), 2.0)
    # constant offsets (e.g. solar motion) do not change it
    assert np.isclose(k.dispersion(U + 10, V - 4, W + 2), np.sqrt(5.0))
    # 2-D input: one value per row
    assert np.allclose(k.dispersion(np.vstack([U, 2 * U]), np.zeros((2, 4)), np.zeros((2, 4))),
                       [np.sqrt(5.0), 2 * np.sqrt(5.0)])


def test_height_z_heliocentric_vs_galactocentric():
    df = pd.DataFrame(dict(ra=[0.0], dec=[90.0], dist_pc=[100.0], pmra=[0.0], pmdec=[0.0], rv=[0.0]))
    c = k.skycoord(df)
    zh = k.height_z(c, {"z_definition": "heliocentric"})[0]
    b = c.galactic.b.radian[0]
    assert np.isclose(zh, 100 * np.sin(b))
    zg = k.height_z(c, {"z_definition": "galactocentric", "galcen_distance_kpc": 8.3, "z_sun_pc": 27.0})[0]
    assert abs(zg - (zh + 27.0)) < 1.0                  # tilt is tiny near the Sun
