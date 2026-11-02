"""Stage C: EBLM IV primaries (Triaud et al. 2017) -> Gaia DR2 targets.

    python -m src.eblm      # parse EBLM IV tables, match Gaia DR2, write the target list

Public inputs only (arXiv source tables). No RVs are read here (C6).
Outputs: data/processed/eblm/eblm4_systems.csv, eblm4_dr2_match.csv,
data/eblm/targets_eblm4.csv (the target-list format of src/targets.py).
"""
import re
import tarfile

import numpy as np
import pandas as pd
from astropy.table import Table

from . import ari_tap, config
from .config import RAW, ROOT

SRC = RAW / "arxiv" / "1707.07521_source.tar.gz"
TEX = "EBLM_V12-nobold.tex"
OUT = ROOT / "data" / "processed" / "eblm"


def _tex():
    with tarfile.open(SRC) as t:
        return t.extractfile(TEX).read().decode("utf-8", "replace")


def _num(s):
    """'1.60(11)' -> 1.60; '$<$0.015' -> nan; '$--$' -> nan."""
    s = s.strip().replace("$", "")
    m = re.match(r"^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    return float(m.group(0)) if m and not s.startswith("<") else np.nan


def _upper(s):
    s = s.strip().replace("$", "")
    return float(s[1:]) if s.startswith("<") else np.nan


def _rows(tex, label, ncol):
    i = tex.index(rf"\label{{{label}}}")
    j = tex.index(r"\end{longtable}", i)
    out = []
    for line in tex[i:j].splitlines():
        line = line.strip()
        if line.startswith("EBLM J"):
            cells = [c.strip() for c in line.rstrip("\\").split("&")]
            if len(cells) >= ncol:
                out.append(cells[:ncol])
    return out


def systems():
    """One row per EBLM IV binary: positions, ephemeris, primary properties, orbit."""
    tex = _tex()
    obs = pd.DataFrame(_rows(tex, "tab:obs", 12),
                       columns=["name", "ra_hms", "dec_dms", "P_obs", "T0_pri", "T0_sec", "sig_1800s",
                                "sig_median", "sig_add", "span_yr", "n_rv", "flag"])
    pri = pd.DataFrame(_rows(tex, "tab:primary", 9),
                       columns=["name", "m_A", "R_A", "Teff", "ebv_eblm4", "Vmag", "Rmag", "Jmag", "SpT"])
    par = pd.DataFrame(_rows(tex, "tab:params", 13),
                       columns=["name", "P", "a", "K", "e", "omega", "T_peri", "f_m", "m_A_par", "m_B",
                                "drift_lin", "drift_quad", "drift_cubic"])
    for d in (obs, pri, par):
        d["name"] = d["name"].str.replace("EBLM ", "", regex=False).str.strip()
    ra = obs["ra_hms"].str.split(expand=True).astype(float)
    de = obs["dec_dms"].str.split(expand=True)
    sign = np.where(de[0].str.startswith("-"), -1.0, 1.0)
    obs["ra_deg"] = 15 * (ra[0] + ra[1] / 60 + ra[2] / 3600)
    obs["dec_deg"] = sign * (de[0].astype(float).abs() + de[1].astype(float) / 60 + de[2].astype(float) / 3600)
    for c in ("P_obs", "T0_pri", "T0_sec", "sig_1800s", "sig_median", "sig_add", "span_yr", "n_rv"):
        obs[c] = obs[c].map(_num)
    for c in ("m_A", "R_A", "Teff", "ebv_eblm4", "Vmag", "Rmag", "Jmag"):
        pri[c] = pri[c].map(_num)
    par["e_upper"] = par["e"].map(_upper)
    par["drift_lin_upper"] = par["drift_lin"].map(_upper)
    for c in ("P", "a", "K", "e", "omega", "T_peri", "f_m", "m_A_par", "m_B", "drift_lin", "drift_quad",
              "drift_cubic"):
        par[c] = par[c].map(_num)
    df = obs.merge(pri, on="name", how="left").merge(par, on="name", how="left")
    df["T0_pri_bjd"] = df["T0_pri"] + 2455000.0
    return df


def dr2_match(sys_df, cfg, radius_arcsec=10.0):
    """All DR2 sources within radius of each EBLM IV J2000 position (one ARI job, uploaded list).

    Gaia positions (epoch 2015.5) are propagated back to J2000 with the DR2 proper
    motion before the separation is computed. The match is the nearest source at
    J2000; its G - V and the next-nearest candidate are kept for checks.
    """
    up = Table.from_pandas(sys_df[["name", "ra_deg", "dec_deg"]].rename(
        columns={"ra_deg": "ra", "dec_deg": "dec"}))
    a = cfg["gaia"]["ari"]
    q = (f"SELECT t.name, g.source_id, g.ra, g.dec, g.pmra, g.pmdec, g.parallax, g.phot_g_mean_mag, "
         f"g.bp_rp\n"                          # no Gaia RVs: never used for binaries (spec)
         f"FROM TAP_UPLOAD.t AS t JOIN {a['source_table']} AS g\n"
         f"ON 1 = CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', t.ra, t.dec, {radius_arcsec / 3600}))")
    m = ari_tap.query(q, "eblm/eblm4_dr2_cone", refresh=cfg["refresh"], uploads={"t": up}).to_pandas()
    m["name"] = m["name"].astype(str)
    s = sys_df.set_index("name")
    dt = 2000.0 - 2015.5
    pmra = m["pmra"].fillna(0) / 3.6e6
    pmde = m["pmdec"].fillna(0) / 3.6e6
    ra0 = m["ra"] + pmra * dt / np.cos(np.radians(m["dec"]))
    de0 = m["dec"] + pmde * dt
    tra, tde = s.loc[m["name"], "ra_deg"].to_numpy(), s.loc[m["name"], "dec_deg"].to_numpy()
    m["sep_j2000"] = 3600 * np.hypot((ra0 - tra) * np.cos(np.radians(tde)), de0 - tde)
    m["G_minus_V"] = m["phot_g_mean_mag"] - s.loc[m["name"], "Vmag"].to_numpy()
    m = m.sort_values(["name", "sep_j2000"])
    first = m.groupby("name").head(1).set_index("name")
    second = m.groupby("name").nth(1).set_index("name")
    first["n_within_radius"] = m.groupby("name").size()
    first["sep2_arcsec"] = second["sep_j2000"]
    first["dG_second"] = second["phot_g_mean_mag"] - first.loc[second.index, "phot_g_mean_mag"]
    first = first.reindex(sys_df["name"])
    first.index.name = "name"
    return first.reset_index()


def match_flags(m, cfg):
    """Quality flag per match (logged; the cut values are in config eblm.match)."""
    c = cfg["eblm"]["match"]
    ok_sep = m["sep_j2000"] < c["max_sep_arcsec"]
    ok_gv = m["G_minus_V"].between(*c["G_minus_V_range"])
    # a second source within radius that is bright enough to be confused with the primary
    confused = m["dG_second"] < c["min_dG_second"]
    flag = np.where(m["source_id"].isna(), "no_source",
                    np.where(~ok_sep, "far", np.where(~ok_gv, "G_minus_V", np.where(confused, "confused", "ok"))))
    return flag


def main():
    cfg = config.load()
    OUT.mkdir(parents=True, exist_ok=True)
    s = systems()
    s.to_csv(OUT / "eblm4_systems.csv", index=False)
    m = dr2_match(s, cfg)
    m["match_flag"] = match_flags(m, cfg)
    m.to_csv(OUT / "eblm4_dr2_match.csv", index=False)
    # target list (src/targets.py format); gamma to be filled by the pre-registered fit (C6)
    tl = pd.DataFrame({"source_id": m["source_id"].astype("Int64"), "name": "EBLM " + m["name"],
                       "is_binary": True, "rv": np.nan, "rv_err": np.nan, "rv_kind": "",
                       "rv_ref": "", "match_flag": m["match_flag"]})
    (ROOT / "data" / "eblm").mkdir(parents=True, exist_ok=True)
    tl.to_csv(ROOT / "data" / "eblm" / "targets_eblm4.csv", index=False)
    print(f"EBLM IV systems parsed: {len(s)} (tab:obs), with V: {s['Vmag'].notna().sum()}, "
          f"with K: {s['K'].notna().sum()}")
    print("match flags:", m["match_flag"].value_counts().to_dict())
    print("sep_j2000 [arcsec] median/max (ok):",
          round(float(m.loc[m.match_flag == 'ok', 'sep_j2000'].median()), 2),
          round(float(m.loc[m.match_flag == 'ok', 'sep_j2000'].max()), 2))
    print(m.loc[m.match_flag != "ok", ["name", "source_id", "sep_j2000", "G_minus_V", "n_within_radius",
                                       "sep2_arcsec", "dG_second", "match_flag"]].to_string(index=False))
    return s, m


if __name__ == "__main__":
    main()
