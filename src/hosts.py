"""Steps 1-2: hot-Jupiter and longer-period giant host samples.

NASA Exoplanet Archive (pscomppars) -> planet selection -> unique hosts ->
Gaia DR2 source_id (SIMBAD, fallback archive) -> Gaia DR2 astrometry + BJ18
distances -> RV (archive systemic RV, else Gaia DR2).
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from . import ari_tap, gaia_tap
from .config import RAW, ROOT
from .gaia_columns import select_list

ARCHIVE_COLS = ("pl_name,hostname,hd_name,hip_name,disc_year,discoverymethod,pl_orbper,"
                "pl_orbperlim,pl_bmassj,pl_bmassjlim,pl_bmassprov,pl_msinij,pl_msinijlim,"
                "st_radv,st_radverr1,st_radverr2,st_radvlim,st_radv_reflink,gaia_dr2_id,"
                "gaia_dr3_id,ra,dec,sy_dist,sy_pnum,pl_pubdate")
MEASURED_MASS = ["Mass", "Msini", "Msin(i)/sin(i)"]     # excludes "M-R relationship" (D4)


# ------------------------------------------------------------------ archive
def archive_table(cfg, snapshot="2026-10-28"):
    """Cached full pscomppars download (all planets; selection happens later)."""
    out = RAW / "exoarchive" / f"pscomppars_{snapshot}.csv"
    if not out.exists() or cfg["refresh"]:
        out.parent.mkdir(parents=True, exist_ok=True)
        r = requests.get(cfg["hosts"]["archive_tap"], timeout=600, params={
            "query": f"select {ARCHIVE_COLS} from {cfg['hosts']['archive_table']}",
            "format": "csv"})
        r.raise_for_status()
        if r.text.lstrip().startswith("<?xml"):
            raise RuntimeError(r.text[:500])
        out.write_text(r.text)
    return pd.read_csv(out)


def select_planets(pl, cfg, kind):
    """Planets meeting the HJ ('hj') or longer-period ('longp') definition.

    Mass: pl_bmassj with a measured provenance (Mass / Msini); upper limits
    kept unless hosts.drop_mass_upper_limits (D4). Period: pl_orbper, non-null.
    """
    h = cfg["hosts"]
    m = pl[h["mass_column"]]
    ok = pl["pl_bmassprov"].isin(MEASURED_MASS)
    if h["drop_mass_upper_limits"]:
        ok &= pl["pl_bmassjlim"].fillna(0) != 1
    if kind == "hj":
        sel = ok & (pl["pl_orbper"] < h["hj"]["period_max"]) & (m > h["hj"]["mass_min"])
    elif kind == "longp":
        sel = ok & (pl["pl_orbper"] > h["longp"]["period_min"]) & (m > h["longp"]["mass_min"])
    else:
        raise ValueError(kind)
    if h["disc_year_max"] is not None:
        sel &= pl["disc_year"] <= h["disc_year_max"]
    return pl[sel]


def host_list(pl, cfg, kind):
    """One row per host for sample `kind`, with archive RV columns."""
    sel = select_planets(pl, cfg, kind)
    if kind == "longp" and cfg["hosts"]["exclude_hj_hosts_from_longp"]:
        hj_hosts = set(select_planets(pl, cfg, "hj")["hostname"])
        sel = sel[~sel["hostname"].isin(hj_hosts)]
    keep = ["hostname", "hd_name", "hip_name", "ra", "dec", "sy_dist", "st_radv", "st_radverr1",
            "st_radverr2", "st_radvlim", "st_radv_reflink", "gaia_dr2_id", "gaia_dr3_id"]
    hosts = sel.sort_values("disc_year").groupby("hostname", as_index=False).agg(
        **{k: (k, "first") for k in keep if k != "hostname"},
        n_planets_sel=("pl_name", "size"), planets=("pl_name", lambda s: ";".join(s)),
        disc_year_min=("disc_year", "min"))
    return hosts


# ------------------------------------------------------------------ SIMBAD
def _name_variants(name):
    v = [name, name.replace("-", " "), re.sub(r"\s+", " ", name)]
    v.append(re.sub(r"^(Qatar|KELT|XO|TrES|CoRoT|Kepler|K2|HAT-P|HATS|WASP|NGTS|TOI)[ -]?", r"\1 ", name))
    return list(dict.fromkeys(v))


def simbad_gaia_ids(names, cache_name, refresh=False):
    """Map host name -> 'Gaia DR2 <id>' via SIMBAD sim-script (Harvard mirror).

    Each query line returns one record; we verify that the queried name is in
    the returned identifier list before accepting it. Cached as CSV.
    """
    out = RAW / "simbad" / f"{cache_name}.csv"
    if out.exists() and not refresh:
        return pd.read_csv(out, dtype=str)
    out.parent.mkdir(parents=True, exist_ok=True)
    url = "https://simbad.harvard.edu/simbad/sim-script"
    rows = []
    pending = {n: _name_variants(n) for n in names}
    for attempt in range(4):
        todo = [(n, v[attempt]) for n, v in pending.items() if attempt < len(v)]
        if not todo:
            break
        for i in range(0, len(todo), 200):
            chunk = todo[i:i + 200]
            script = ['output console=off script=off',
                      'format object "@@REC %IDLIST(1) | %IDLIST(Gaia DR2) | %IDLIST(Gaia DR3) | %IDLIST"']
            script += [f"query id {q}" for _, q in chunk]
            r = requests.post(url, data={"script": "\n".join(script)}, timeout=300)
            r.raise_for_status()
            txt = r.text
            err_lines = {int(x) for x in re.findall(r"^\[(\d+)\]", txt.split("::data::")[0], re.M)}
            body = txt.split("::data::")[-1]          # no ::data:: header when no errors
            data = [rec.strip() for rec in body.split("@@REC ")[1:]]     # one record per object
            ok = [c for j, c in enumerate(chunk) if (j + 3) not in err_lines]   # script line numbers
            if len(ok) != len(data):
                raise RuntimeError(f"SIMBAD response mismatch ({len(ok)} vs {len(data)})")
            for (name, q), line in zip(ok, data):
                main, dr2, dr3, ids = [x.strip() for x in line.split("|", 3)]
                norm = lambda s: re.sub(r"[^0-9a-z+]", "", s.lower())
                if norm(q) not in {norm(x) for x in ids.splitlines()}:
                    continue                                # mismatched record, try next variant
                rows.append({"hostname": name, "simbad_query": q, "simbad_main_id": main,
                             "simbad_dr2": dr2, "simbad_dr3": dr3})
                pending.pop(name, None)
    res = pd.DataFrame(rows, columns=["hostname", "simbad_query", "simbad_main_id",
                                      "simbad_dr2", "simbad_dr3"])
    res.to_csv(out, index=False)
    return res.astype(str)


def _dr2_ids(series):
    """'Gaia DR2 <id>' strings -> nullable Int64 (source_ids exceed float64 precision)."""
    return series.astype(str).str.extract(r"Gaia DR2 (\d+)")[0].astype("Int64")


def attach_gaia_ids(hosts, cfg, cache_name):
    """Add `source_id` (Gaia DR2) and `id_source` (simbad|archive|none)."""
    sim = simbad_gaia_ids(hosts["hostname"].tolist(), cache_name, cfg["refresh"])
    h = hosts.merge(sim, on="hostname", how="left")
    sid_sim = _dr2_ids(h["simbad_dr2"])
    sid_arc = _dr2_ids(h["gaia_dr2_id"])
    primary, other = (sid_sim, sid_arc) if cfg["hosts"]["gaia_id_source"] == "simbad" else (sid_arc, sid_sim)
    h["source_id"] = primary.where(primary.notna(), other)
    h["id_source"] = np.where(primary.notna(), cfg["hosts"]["gaia_id_source"],
                              np.where(other.notna(), "fallback", "none"))
    h["id_disagree"] = sid_sim.notna() & sid_arc.notna() & (sid_sim != sid_arc)
    return h


# ------------------------------------------------------------------ Gaia
def gaia_rows(source_ids, cfg, name):
    """Gaia DR2 rows + BJ18 distances for a list of source_ids (cached)."""
    ids = sorted({int(s) for s in source_ids if pd.notna(s)})
    g = cfg["gaia"]
    if g["host_service"] == "ari":              # one job per list; no 60 s cap (A8)
        src, dist = g["ari"]["source_table"], g["ari"]["distance_table"]
        q = (f"SELECT {select_list()}\nFROM {src} AS g\n"
             f"LEFT OUTER JOIN {dist} AS d ON g.source_id = d.source_id\n"
             f"WHERE g.source_id IN ({', '.join(map(str, ids))})")
        return ari_tap.query(q, f"{name}_ari", refresh=cfg["refresh"]).to_pandas()
    parts = []
    for i in range(0, len(ids), 50):            # small batches stay under the 60 s anonymous limit
        chunk = ids[i:i + 50]
        q = (f"SELECT {select_list()}\nFROM {g['source_table']} AS g\n"
             f"LEFT OUTER JOIN {g['distance_table']} AS d ON g.source_id = d.source_id\n"
             f"WHERE g.source_id IN ({', '.join(map(str, chunk))})")
        parts.append(gaia_tap.query(q, f"{name}/batch{i // 50:03d}", refresh=cfg["refresh"]).to_pandas())
    return pd.concat(parts, ignore_index=True)


def assign_rv(h, cfg):
    """rv = archive st_radv unless missing/limit/rejected reference, else Gaia DR2 RV."""
    ref = h["st_radv_reflink"].fillna("").str.replace(r"<[^>]*>", "", regex=True)
    pat = "|".join(re.escape(p) for p in cfg["hosts"]["archive_rv_reject_ref_patterns"])
    arc_ok = h["st_radv"].notna() & (h["st_radvlim"].fillna(0) == 0) & ~ref.str.contains(pat)
    h["rv_archive_ref"] = ref.str.strip()
    h["rv"] = np.where(arc_ok, h["st_radv"], h["radial_velocity"])
    h["rv_source"] = np.where(arc_ok, "archive",
                              np.where(h["radial_velocity"].notna(), "gaia_dr2", "none"))
    return h


def table1_hosts(pl, cfg):
    """HS19 Table 1 hosts with archive host-level RV columns (matched on ArchiveName).

    SIMBAD is queried with the paper's own SimbadName.
    """
    t1 = load_table1(cfg)
    hosts_all = (pl.sort_values("disc_year").groupby("hostname", as_index=False)
                 .first()[["hostname", "hd_name", "hip_name", "ra", "dec", "sy_dist", "st_radv",
                           "st_radverr1", "st_radverr2", "st_radvlim", "st_radv_reflink",
                           "gaia_dr2_id", "gaia_dr3_id"]])
    h = t1.rename(columns={"ArchiveName": "archive_name_2019", "SimbadName": "simbad_name_2019"})
    h["hostname"] = h["simbad_name_2019"]                 # SIMBAD lookup key
    h = h.merge(hosts_all.rename(columns={"hostname": "archive_name_2019"}),
                on="archive_name_2019", how="left")
    h["in_archive_now"] = h["gaia_dr2_id"].notna() | h["ra"].notna()
    return h


def dr3_rv(dr2_ids, cfg, name):
    """Gaia DR3 radial_velocity for DR2 source_ids via gaiadr3.dr2_neighbourhood
    (closest DR3 match by angular distance; magnitude difference kept for checks)."""
    ids = sorted({int(s) for s in dr2_ids if pd.notna(s)})
    parts = []
    for i in range(0, len(ids), 50):
        q = ("SELECT n.dr2_source_id, n.dr3_source_id, n.angular_distance, n.magnitude_difference, "
             "g.radial_velocity AS rv_dr3, g.radial_velocity_error AS rv_dr3_error "
             "FROM gaiadr3.dr2_neighbourhood AS n JOIN gaiadr3.gaia_source AS g "
             "ON n.dr3_source_id = g.source_id "
             f"WHERE n.dr2_source_id IN ({', '.join(map(str, ids[i:i + 50]))})")
        parts.append(gaia_tap.query(q, f"{name}/batch{i // 50:03d}", refresh=cfg["refresh"]).to_pandas())
    df = pd.concat(parts, ignore_index=True).sort_values("angular_distance")
    df = df.drop_duplicates("dr2_source_id").rename(columns={"dr2_source_id": "source_id"})
    df["source_id"] = df["source_id"].astype("Int64")
    return df


def build(cfg, kind):
    """Full host table for `kind` ('hj' | 'longp' | 'hj_table1') before quality cuts."""
    pl = archive_table(cfg)
    if kind == "hj_table1":
        hosts = table1_hosts(pl, cfg)
        tag = "hj_table1"
    else:
        hosts = host_list(pl, cfg, kind)
        tag = f"{kind}_dy{cfg['hosts']['disc_year_max'] or 'all'}"
    h = attach_gaia_ids(hosts, cfg, f"simbad_{tag}")
    gaia = gaia_rows(h["source_id"], cfg, f"hosts_{tag}_{cfg['gaia']['release']}")
    gaia["source_id"] = gaia["source_id"].astype("Int64")
    h = h.merge(gaia, on="source_id", how="left", suffixes=("_archive", ""))
    h["in_gaia"] = h["parallax"].notna() | h["phot_g_mean_mag"].notna()
    h = assign_rv(h, cfg)
    if cfg["hosts"].get("rv_fill_dr3"):
        # Variant (not in HS19): hosts with neither an archive nor a DR2 RV get the DR3 RV
        miss = h["rv"].isna() & h["source_id"].notna()
        if miss.any():
            d3 = dr3_rv(h.loc[miss, "source_id"], cfg, f"hosts_{tag}_dr3rv")
            h = h.merge(d3[["source_id", "rv_dr3"]], on="source_id", how="left")
            fill = h["rv"].isna() & h["rv_dr3"].notna()
            h.loc[fill, "rv"] = h.loc[fill, "rv_dr3"]
            h.loc[fill, "rv_source"] = "gaia_dr3_fill"
    return h


def load_table1(cfg):
    """HS19 machine-readable Table 1, if present (else None)."""
    p = ROOT / cfg["hosts"]["table1_file"]
    if not p.exists():
        return None
    from astropy.io import ascii
    return ascii.read(p).to_pandas()
