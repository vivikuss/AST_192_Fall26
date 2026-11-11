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


def test_protected_stage_e_outputs_cannot_be_overwritten():
    import sys
    from src import hu2026_checks as hc
    before = {p: p.read_bytes() for p in hc.PROTECTED if p.exists()}
    for p in hc.PROTECTED:
        try:
            hc._write(p, "x")
        except RuntimeError:
            pass
        else:
            raise AssertionError(f"{p} was writable")
    argv = sys.argv
    for cmd in ("wasp", "power", "contrast"):
        sys.argv = ["hu2026_checks", cmd]
        try:
            hc.main()
        except SystemExit:
            pass
        else:
            raise AssertionError(cmd)
        finally:
            sys.argv = argv
    assert all(p.read_bytes() == b for p, b in before.items())


def test_hemisphere_csv_readers_use_comment_hash():
    """notes/wasp_hemisphere_labels.csv has a leading '#' header (E35). Read without comment='#', pandas
    silently returns 99 rows x 1 column. Every reader in the repository must pass comment='#' (E36)."""
    import re
    from src.config import ROOT
    path = ROOT / "notes" / "wasp_hemisphere_labels.csv"
    assert pd.read_csv(path, comment="#").shape[0] == 98
    for f in list((ROOT / "src").rglob("*.py")) + list((ROOT / "tests").rglob("*.py")) + \
            list((ROOT / "notebooks").rglob("*.ipynb")):
        txt = f.read_text()
        for m in re.finditer(r"read_csv\([^)]*wasp_hemisphere_labels[^)]*\)", txt):
            assert "comment=" in m.group(0), (f, m.group(0))
