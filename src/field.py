"""Step 3: Gaia DR2 field sample (5-parameter astrometry + RV), with BJ18 distances."""
import numpy as np

from . import ari_tap, cuts, gaia_tap
from .gaia_columns import select_list


def tables(cfg, service=None):
    """(source_table, distance_table) on the given service (default gaia.field_service)."""
    g = cfg["gaia"]
    s = g if (service or g["field_service"]) == "esa" else g[service or g["field_service"]]
    return s["source_table"], s["distance_table"]


def field_adql(cfg, service="esa"):
    """ADQL template; `{range}` is a source_id range (ESA chunks) or `1=1` (one job)."""
    src, dist = tables(cfg, service)
    where = cuts.server_where(cfg["cuts"]).replace("{", "{{").replace("}", "}}")
    return (f"SELECT {select_list()}\n"
            f"FROM {src} AS g\n"
            f"LEFT OUTER JOIN {dist} AS d ON g.source_id = d.source_id\n"
            f"WHERE {{range}}\n AND {where}")


def count_adql(cfg, service="esa"):
    """COUNT(*) with the same WHERE as field_adql (the distance join is 1:1, LEFT)."""
    src, _ = tables(cfg, service)
    where = cuts.server_where(cfg["cuts"]).replace("{", "{{").replace("}", "}}")
    return f"SELECT COUNT(*) AS n\nFROM {src} AS g\nWHERE {{range}}\n AND {where}"


def chunk_name(cfg):
    return f"field_{cfg['gaia']['release']}_prefilter"


def download(cfg, allow_missing=False):
    """Field prefilter (cuts 1-3,5-7,9,10 exact; 8 loose; 4 local).

    gaia.field_service = ari: one ARI Heidelberg job per HEALPix source_id range
    (A6, A7); esa: ESA archive chunks split on failure (A4). Cached per range.
    Raises if any range failed, unless allow_missing (then returns (df, failed)).
    """
    if cfg["gaia"]["field_service"] == "ari":
        tbl, failed = ari_tap.query_ranges(field_adql(cfg, "ari"), chunk_name(cfg),
                                           cfg["gaia"]["ari"]["range_level"], refresh=cfg["refresh"])
    else:
        tbl, failed = gaia_tap.query_chunked(field_adql(cfg), chunk_name(cfg), level=2, refresh=cfg["refresh"])
    if allow_missing:
        return tbl.to_pandas(), failed
    if failed:
        raise RuntimeError(f"field download incomplete, failed ranges: "
                           f"{[f[:-1] for f in failed]}; rerun python -u -m src.download_field")
    return tbl.to_pandas()


def select(df, cfg):
    """Apply the field cut set locally. Returns (subset, sequential cut table)."""
    which = cfg["cuts"]["apply"]["field"]
    mask, _ = cuts.apply(df, cfg["cuts"], which)
    return df[mask].reset_index(drop=True), cuts.cut_table(df, cfg["cuts"], which)
