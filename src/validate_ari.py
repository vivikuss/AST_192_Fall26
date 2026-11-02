"""A6 validation: ARI Heidelberg vs ESA for the same DR2 field chunks and BJ18 distances.

    python -u -m src.validate_ari
Writes logs/ari_validation.json. Exit code 1 on any difference.
"""
import json
import sys

import numpy as np
import pandas as pd
from astropy.table import Table

from . import ari_tap, config, field, gaia_tap, hosts
from .config import RAW, ROOT
from .gaia_tap import log

REGIONS = [0, 109, 110, 111]          # L2 regions whose ESA download matched ESA COUNT(*)
N_HOSTS = 10


def _frame(tbl):
    df = tbl.to_pandas()
    df.columns = [c.lower() for c in df.columns]
    return df.set_index("source_id").sort_index()


def compare(esa, ari, label):
    """Same source_ids, and every column identical (NaN == NaN), compared in float64."""
    res = {"label": label, "n_esa": len(esa), "n_ari": len(ari),
           "only_esa": len(esa.index.difference(ari.index)), "only_ari": len(ari.index.difference(esa.index)),
           "columns": {}}
    common = esa.index.intersection(ari.index)
    missing = sorted(set(esa.columns) ^ set(ari.columns))
    res["column_set_diff"] = missing
    for c in sorted(set(esa.columns) & set(ari.columns)):
        a = pd.to_numeric(esa.loc[common, c], errors="coerce").to_numpy(float)
        b = pd.to_numeric(ari.loc[common, c], errors="coerce").to_numpy(float)
        same = (a == b) | (np.isnan(a) & np.isnan(b))
        d = np.abs(a - b)
        res["columns"][c] = {"n_diff": int((~same).sum()),
                             "max_abs_diff": float(np.nanmax(d)) if (~same).any() and np.isfinite(d).any() else 0.0,
                             "dtype_esa": str(esa[c].dtype), "dtype_ari": str(ari[c].dtype)}
    res["identical"] = (res["only_esa"] == res["only_ari"] == 0 and not missing
                        and all(v["n_diff"] == 0 for v in res["columns"].values()))
    return res


def main():
    cfg = config.load()
    out = []
    esa_t, ari_t = field.field_adql(cfg, "esa"), field.field_adql(cfg, "ari")
    name = field.chunk_name(cfg)
    for p in REGIONS:
        esa = _frame(Table.read(RAW / "gaia" / name / f"L2_{p:05d}.fits"))
        ari = _frame(ari_tap.query(gaia_tap.chunk_query(ari_t, 2, p), f"validate/{name}_L2_{p:05d}"))
        out.append(compare(esa, ari, f"field L2 region {p}"))

    ids = hosts.build(cfg, "hj_table1")["source_id"].dropna().astype("int64")
    ids = sorted(ids.sample(N_HOSTS, random_state=cfg["seed"]).tolist())
    cols = "source_id, r_est, r_lo, r_hi, result_flag"
    idlist = ", ".join(map(str, ids))
    esa_d = gaia_tap.query(f"SELECT {cols} FROM {cfg['gaia']['distance_table']} WHERE source_id IN ({idlist})",
                           "validate/bj18_hosts_esa")
    ari_d = ari_tap.query(f"SELECT {cols} FROM {cfg['gaia']['ari']['distance_table']} WHERE source_id IN ({idlist})",
                          "validate/bj18_hosts_ari")
    r = compare(_frame(esa_d), _frame(ari_d), f"BJ18 distances, {N_HOSTS} HJ hosts")
    r["source_ids"] = ids
    out.append(r)

    ok = True
    for r in out:
        bad = {c: v for c, v in r["columns"].items() if v["n_diff"]}
        log(f"{r['label']}: ESA {r['n_esa']} ARI {r['n_ari']} only-ESA {r['only_esa']} only-ARI {r['only_ari']} "
            f"columns compared {len(r['columns'])} differing {len(bad)} "
            f"{'IDENTICAL' if r['identical'] else 'DIFFERENT'}")
        for c, v in bad.items():
            log(f"    {c}: {v}")
        if r["column_set_diff"]:
            log(f"    column set differs: {r['column_set_diff']}")
        ok &= r["identical"]
    (ROOT / "logs" / "ari_validation.json").write_text(json.dumps(out, indent=1))
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
