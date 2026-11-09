"""WASP membership join (DECISIONS E32): facilities keyed by Gaia DR2 source_id, not host name."""
import pandas as pd

from src.hu2026_checks import facility_label, facility_map


def _archive():
    planets = pd.DataFrame({
        "pl_name": ["WASP-161 b", "HAT-P-1 b", "WASP-53 b"],
        "hostname": ["WASP-161", "HAT-P-1", "WASP-53"],
        "gaia_dr2_id": ["Gaia DR2 5751177095877066496", "Gaia DR2 1", None],
    })
    disc = pd.DataFrame({"disc_facility": ["SuperWASP-South", "HATNet", "Multiple Facilities"]},
                        index=pd.Index(["WASP-161 b", "HAT-P-1 b", "WASP-53 b"], name="pl_name"))
    return planets, disc


def test_source_id_join_labels_host_with_catalogue_simbad_name():
    planets, disc = _archive()
    hosts = pd.DataFrame({"hostname": ["TYC 5435-1109-1"],          # HS19 SimbadName of WASP-161
                          "source_id": pd.array([5751177095877066496], dtype="Int64")})
    fac = facility_map(planets, disc, "source_id")
    assert hosts["source_id"].map(fac).fillna("").map(facility_label).tolist() == ["WASP-South"]
    # the legacy hostname join misses it
    fac_old = facility_map(planets, disc, "hostname")
    assert hosts["hostname"].map(fac_old).fillna("").map(facility_label).isna().all()


def test_source_ids_beyond_float_precision_and_missing_ids():
    planets, disc = _archive()
    fac = facility_map(planets, disc, "source_id")
    assert 5751177095877066496 in fac.index           # exact 19-digit id, no float rounding
    assert len(fac) == 2                               # the planet without a DR2 id is dropped


def test_facility_label():
    assert facility_label("SuperWASP-North") == "SuperWASP-North"
    assert facility_label("SuperWASP-South; WASP-South") == "WASP-South"
    assert facility_label("SuperWASP") == "SuperWASP (hemisphere not given)"
    assert facility_label("HATNet") is None
