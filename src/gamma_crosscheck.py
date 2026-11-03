"""Stage D cross-checks for the EBLM gamma fits (C9: checks only, never a source of gamma).

    python -m src.gamma_crosscheck nss     # Gaia DR3 two-body orbits of the EBLM IV primaries (ARI)

Public inputs only. Output: data/processed/eblm/eblm4_dr3_nss.csv
"""
import sys

import pandas as pd
from astropy.table import Table

from . import ari_tap, config
from .config import ROOT

OUT = ROOT / "data" / "processed" / "eblm"


def nss(cfg):
    ids = pd.read_csv(OUT / "eblm4_dr3_crosscheck.csv")[["name", "dr3_source_id"]].dropna()
    up = Table.from_pandas(ids.rename(columns={"dr3_source_id": "source_id"}).astype({"source_id": "int64"}))
    q = ("SELECT t.name, n.source_id, n.nss_solution_type, n.period, n.period_error, "
         "n.center_of_mass_velocity, n.center_of_mass_velocity_error, "
         "n.semi_amplitude_primary, n.semi_amplitude_primary_error, n.eccentricity, n.eccentricity_error, "
         "n.goodness_of_fit, n.significance, n.flags\n"
         "FROM TAP_UPLOAD.t AS t JOIN gaiadr3.nss_two_body_orbit AS n ON n.source_id = t.source_id")
    m = ari_tap.query(q, "eblm/eblm4_dr3_nss", refresh=cfg["refresh"], uploads={"t": up}).to_pandas()
    m.to_csv(OUT / "eblm4_dr3_nss.csv", index=False)
    print(f"DR3 NSS two-body solutions: {len(m)} rows for {m['name'].nunique()} of {len(ids)} systems")
    print(m["nss_solution_type"].value_counts().to_dict())


def main():
    cfg = config.load()
    ari_tap.install_signal_handlers()
    if sys.argv[1:2] == ["nss"]:
        nss(cfg)


if __name__ == "__main__":
    main()
