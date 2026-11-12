"""DR3 EXPLORATORY run (DECISIONS E52, E53): Method A building blocks on synthetic data only."""
import pathlib
import tempfile

import numpy as np
import pandas as pd

from src import dr3x_method_a as A
from src import seeds


def test_dr3x_labels_registered_and_share_no_code():
    codes = {}
    for label in seeds.REGISTRY:
        codes.setdefault(sum(map(ord, label)), set()).add(label)
    for label in seeds.DR3X:
        assert label in seeds.REGISTRY and label.startswith(seeds.DR3X_PREFIX)
        assert codes[sum(map(ord, label))] == {label}, codes[sum(map(ord, label))]


def test_healpix_nested_hierarchy_and_equal_area():
    rng = np.random.default_rng(0)
    n = 400_000
    lon = rng.uniform(0, 360, n)
    lat = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))
    for nside in (2, 4, 8):
        p = A.ang2pix_nest(nside, lon, lat)
        assert p.min() == 0 and p.max() == 12 * nside ** 2 - 1
        c = np.bincount(p, minlength=12 * nside ** 2)
        e = n / (12 * nside ** 2)
        assert np.all(np.abs(c - e) < 6 * np.sqrt(e)), (nside, c.min(), c.max(), e)   # equal-area pixels
        assert np.array_equal(A.ang2pix_nest(nside // 2, lon, lat), p // 4)          # NESTED parent = child // 4
    # base pixels: north pole in face 0-3, south pole in 8-11, equator on the ecliptic of faces 4-7 or 0-11
    assert A.ang2pix_nest(1, [45.0], [89.0])[0] == 0 and A.ang2pix_nest(1, [45.0], [-89.0])[0] == 8
    assert A.ang2pix_nest(1, [0.0], [0.0])[0] == 4


def _fake_bayestar(path, missing_pix=(8, 9, 10, 11)):
    import h5py
    nside = 1
    pix = np.array([p for p in range(12) if p not in missing_pix])
    info = np.zeros(len(pix), dtype=[("nside", "<u4"), ("healpix_index", "<i8")])
    info["nside"], info["healpix_index"] = nside, pix
    dm = np.linspace(4.0, 19.0, 4)                                  # 4 distance bins
    prof = np.array([0.0, 0.1, 0.2, 0.3])
    samples = np.repeat(np.tile(prof, (len(pix), 1))[:, None, :], 5, axis=1)   # 5 identical samples
    samples[:, 0, :] += 10.0                                        # one outlier sample: median ignores it
    with h5py.File(path, "w") as f:
        d = f.create_dataset("pixel_info", data=info)
        d.attrs["DM_bin_edges"] = dm
        f.create_dataset("samples", data=samples)
    return dm, prof


def test_bayestar_stops_when_missing_and_follows_the_rules():
    with tempfile.TemporaryDirectory() as t:
        missing = pathlib.Path(t) / "none.h5"
        try:
            A.bayestar_ebv(np.array([10.0]), np.array([10.0]), np.array([100.0]), missing)
        except RuntimeError as e:
            assert "no zero fallback" in str(e)
        else:
            raise AssertionError("missing Bayestar19 file did not stop")
        p = pathlib.Path(t) / "b.h5"
        dm, prof = _fake_bayestar(p)
        # a point in base pixel 0 (b ~ +60 deg in Galactic coordinates): pick it in ICRS via astropy
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        north = SkyCoord(l=45 * u.deg, b=60 * u.deg, frame="galactic").icrs
        south = SkyCoord(l=45 * u.deg, b=-60 * u.deg, frame="galactic").icrs
        ra = np.array([north.ra.deg] * 4 + [south.ra.deg])
        dec = np.array([north.dec.deg] * 4 + [south.dec.deg])
        d = 10 ** ((np.array([3.0, 9.0, 19.0, 25.0, 9.0]) + 5) / 5)   # DM 3 (too close), 9, 19 (last), 25, 9
        ebv, why = A.bayestar_ebv(ra, dec, d, p)
        assert list(why) == ["zero: closer than first bin", "bayestar19", "bayestar19", "bayestar19",
                             "zero: outside map"]
        exp9 = np.interp(9.0, dm, prof) * A.BAYESTAR_TO_EBV
        assert np.allclose(ebv, [0.0, exp9, 0.3 * 0.884, 0.3 * 0.884, 0.0])   # held at the last bin


def test_photometry_prefers_gspphot_and_never_zero_fills_silently():
    h = pd.DataFrame({"ra": [10.0, 10.0], "dec": [10.0, 10.0], "parallax": [10.0, 10.0],
                      "phot_g_mean_mag": [10.0, 10.0], "bp_rp": [1.0, 1.0],
                      "ag_gspphot": [0.1, np.nan], "ebpminrp_gspphot": [0.05, np.nan]})
    try:
        A.photometry(h, bayestar_path="/nonexistent/bayestar2019.h5")
    except RuntimeError:
        pass
    else:
        raise AssertionError("a host without GSP-Phot was handled without Bayestar19")
    out = A.photometry(h.iloc[:1], bayestar_path="/nonexistent/bayestar2019.h5")   # GSP-Phot only: no map needed
    assert out["ext_source"].iat[0] == "gspphot" and np.isclose(out["color0"].iat[0], 0.95)
    assert np.isclose(out["d"].iat[0], 100.0)


def test_pools_are_inclusive_and_draws_reproducible():
    field = pd.DataFrame({"ra": [0.0, 0.0, 0.0, 0.0, 20.0], "dec": [0.0] * 5,
                          "d": [100.0, 150.0, 151.0, 0.0, 100.0], "phot_g_mean_mag": [10.0, 9.0, 10.1, 10.0, 10.0],
                          "color0": [1.2, 0.8, 1.0, 1.0, 1.0], "U": [1.0, 2, 3, 4, 5], "V": [0.0] * 5,
                          "W": [1.0, -1, 2, -2, 3]})
    hosts = pd.DataFrame({"name": ["h"], "ra": [0.0], "dec": [0.0], "d": [100.0], "phot_g_mean_mag": [10.0],
                          "color0": [1.0]})
    pl = A.pools(hosts, field)
    # star0: dcol 0.2 (edge, in); star1: dd +50 (edge), dg -1 (edge), dcol 0.2 (in); star2: dd +51 out;
    # star3: dd -100 (edge, in); star4: 20 deg away (out)
    assert list(pl[0]) == [0, 1, 3], pl
    uvw = field[["U", "V", "W"]].to_numpy()
    a = A.draws(pl, uvw, "dr3x_A_val_tess", n=50)
    b = A.draws(pl, uvw, "dr3x_A_val_tess", n=50)
    assert all(np.array_equal(a[k], b[k]) for k in a)
    try:
        A.draws(pl, uvw, "dr3x_unregistered", n=5)
    except AssertionError:
        pass
    else:
        raise AssertionError("unregistered label accepted")


def test_statistics_and_one_sided_sigma():
    v = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [2.0, 0.0, 1.0], [-2.0, 0.0, -1.0]])
    assert np.isclose(A.sigma_w(v), 1.0)
    assert np.isclose(A.sigma_tot(v), np.sqrt(2.0 + 1.0))
    assert np.isclose(A.eq1(v), np.mean(np.linalg.norm(v, axis=1)))
    mc = np.arange(1, 101, dtype=float)
    p, s = A.one_sided_sigma(10.0, mc)
    assert np.isclose(p, 0.10) and s < 0
    p, s = A.one_sided_sigma(1000.0, mc)
    assert p == 0.0 and s > 0 and np.isclose(s, 2.326, atol=1e-3)          # floor 1/N


def test_wasp_run_refuses_without_stage2():
    try:
        A.main(["wasp"])
    except SystemExit as e:
        assert "Stage 2" in str(e)
    else:
        raise AssertionError


def test_add_uvw_solar_motion_and_axes():
    # a star with zero heliocentric motion has UVW = SUN_UVW; RV toward the Galactic centre is +U
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    gc = SkyCoord(l=0 * u.deg, b=0 * u.deg, frame="galactic").icrs
    h = pd.DataFrame({"ra": [gc.ra.deg] * 2, "dec": [gc.dec.deg] * 2, "d": [100.0] * 2,
                      "pmra": [0.0] * 2, "pmdec": [0.0] * 2, "rv": [0.0, 10.0]})
    v = A.add_uvw(h)[["U", "V", "W"]].to_numpy()
    assert np.allclose(v[0], A.SUN_UVW, atol=1e-6)
    assert np.allclose(v[1], np.array(A.SUN_UVW) + [10.0, 0, 0], atol=1e-6)


def test_nea_radv_by_name_rule():
    with tempfile.TemporaryDirectory() as t:
        p = pathlib.Path(t) / "nea.csv"
        pd.DataFrame({"pl_name": ["X c", "X b", "Y b", "Y c", "Z b"], "hostname": ["X", "X", "Y", "Y", "Z"],
                      "pl_orbper": [3.0, 2.0, 50.0, 3.0, 80.0], "pl_bmassj": [1.0, 1.0, 1.0, 1.0, 1.0],
                      "disc_year": [2010, 2010, 2001, 2012, 2005], "st_radv": [5.0, np.nan, 7.0, 9.0, 4.0]}
                     ).to_csv(p, index=False)
        r = A.nea_radv_by_name(["X", "Y", "Z"], p)
        # X: HJ rows sorted (2010, "X b") first but null -> next non-null 5.0; Y: only the HJ row counts -> 9.0;
        # Z: no HJ row -> all rows -> 4.0
        assert r.to_dict() == {"X": 5.0, "Y": 9.0, "Z": 4.0}


def test_write_csv_header_and_gzip():
    import gzip
    with tempfile.TemporaryDirectory() as t:
        for name in ("x.csv", "x.csv.gz"):
            p = pathlib.Path(t) / name
            A.write_csv(p, pd.DataFrame({"a": [1, 2]}))
            raw = p.read_bytes()
            assert raw.startswith(b"\x1f\x8b") == name.endswith(".gz")
            text = gzip.decompress(raw).decode() if name.endswith(".gz") else raw.decode()
            assert text.splitlines()[0] == f"# {A.HEADER}"
            assert pd.read_csv(p, comment="#")["a"].tolist() == [1, 2]
