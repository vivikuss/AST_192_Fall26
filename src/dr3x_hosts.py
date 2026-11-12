"""EXPLORATORY (not pre-registered; DECISIONS E52-E54). Stage 1 host step counts for the 107 WASP hosts,
Method A and Method B, by name. No kinematics: RVs are used only for presence and source, never their values.

    python -m src.dr3x_hosts      # -> data/processed/dr3_exploratory/host_steps_wasp107.json

Gaia DR3 rows of the 107 (DR3 IDs from E48) come from one ARI query joined to the Bailer-Jones et al. (2021)
distances (ARI gaiaedr3_complements.distances), cached by src.ari_tap.

Method A (Kayla's steps, E54): close-companion exclusion (main result; the count without it is also given),
no DR3 row, App. A cuts 1-3, RV (NEA st_radv, else DR3 RV), extinction (GSP-Phot, else Bayestar19 web query),
MS + 0.5 <= (BP-RP)_0 <= 1.5, d = 1000/parallax < 600 pc. Hosts without GSP-Phot take Bayestar19 from the local
bayestar2019.h5 in Kayla's order (median of the samples per tabulated distance, then interpolate; E59, E60);
zero only where the map has no value (no coverage, or closer than its first distance bin).

Method B (frozen DR2 rules ported, E53/E54): no DR3 row, astrometric_params_solved = 31, host cuts 1, 3, 4
(cut 2 dropped), RV (archive st_radv under the frozen D6 rule, else Gaia DR3 RV), BJ21 r_med_geo present, MS
(STILISM v4.1 along r_med_geo, published polynomial, no colour limit). No exclusions, no distance cut.
"""
import numpy as np
import pandas as pd

from . import ari_tap, config, cuts, extinction, hosts, main_sequence
from . import dr3x_method_a as A

COLS = ["source_id", "ra", "dec", "l", "b", "parallax", "pmra", "pmdec", "radial_velocity", "phot_g_mean_mag",
        "bp_rp", "ag_gspphot", "ebpminrp_gspphot", "parallax_over_error", "visibility_periods_used",
        "astrometric_chi2_al", "astrometric_n_good_obs_al", "ruwe", "astrometric_params_solved"]
CLOSE_COMPANION = ["KELT-2 A", "KELT-4 A", "WASP-70 A", "WASP-94 A", "WASP-94 B", "NGTS-3 A", "WASP-173 A"]


def gaia_rows(dr3_ids, name):
    ids = ", ".join(str(int(i)) for i in sorted(set(dr3_ids)))
    q = (f"SELECT {', '.join('g.' + c for c in COLS)}, d.r_med_geo, d.flag AS bj21_flag "
         f"FROM gaiadr3.gaia_source AS g LEFT JOIN gaiaedr3_complements.distances AS d "
         f"ON d.source_id = g.source_id WHERE g.source_id IN ({ids})")
    t = ari_tap.query(q, name).to_pandas()
    t.columns = [c.lower() for c in t.columns]
    return t.rename(columns={"source_id": "dr3"})


class Steps:
    def __init__(self, h):
        self.h, self.rows = h.reset_index(drop=True), [{"step": "start", "n_after": len(h), "removed": []}]

    def cut(self, step, keep, note=None):
        keep = np.asarray(keep, dtype=bool)
        r = {"step": step, "removed": self.h.loc[~keep, "name"].tolist(), "n_after": int(keep.sum())}
        if note:
            r["note"] = note
        self.rows.append(r)
        self.h = self.h[keep].reset_index(drop=True)


def method_a(base, nea_rv, exclude_companions):
    s = Steps(base)
    if exclude_companions:
        s.cut("close stellar companion (Kayla's list)", ~s.h["name"].isin(CLOSE_COMPANION))
    s.cut("no Gaia DR3 row", s.h["parallax"].notna())
    s.cut("App. A astrometric cuts 1-3", A.astrometry_ok(s.h))
    has = s.h["dr3"].map(nea_rv).notna() | s.h["radial_velocity"].notna()
    s.h["rv_source"] = np.where(s.h["dr3"].map(nea_rv).notna(), "NEA st_radv",
                                np.where(s.h["radial_velocity"].notna(), "Gaia DR3", "none"))
    s.cut("no RV (NEA or Gaia DR3)", has)
    s.h = A.photometry(s.h)                                    # GSP-Phot, else Bayestar19 (E60)
    ext = {n: e for n, e in zip(s.h["name"], s.h["ext_source"]) if e != "gspphot"}
    s.cut("main sequence + 0.5<=(BP-RP)_0<=1.5",
          ((s.h["zams_offset"] > -A.MS_MAX_ABOVE) & s.h["color0"].between(*A.COLOUR_WINDOW)).to_numpy(),
          note=f"extinction not from GSP-Phot: {ext}")
    s.cut("d = 1000/parallax < 600 pc", (1000.0 / s.h["parallax"]) < A.D_MAX)
    return s.rows, ext, s.h


def method_b(base, cfg):
    s = Steps(base)
    s.cut("no Gaia DR3 row", s.h["parallax"].notna())
    s.cut("astrometric_params_solved = 31 (5-parameter)", s.h["astrometric_params_solved"] == 31)
    m, _ = cuts.apply(s.h, cfg["cuts"], [1, 3, 4])
    s.cut("host cuts 1, 3, 4 (cut 2 dropped: no DR3 column)", m)
    r = hosts.assign_rv(s.h.copy(), cfg)                       # frozen D6 rule, DR3 RV as the fallback
    s.h["rv_source"] = np.where(r["rv_source"] == "gaia_dr2", "Gaia DR3", r["rv_source"])
    s.cut("no RV (archive st_radv under D6, else Gaia DR3)", r["rv_source"].ne("none").to_numpy())
    s.cut("no BJ21 r_med_geo", s.h["r_med_geo"].notna())
    h = s.h.copy()
    ebv = extinction.ebv_los(h["l"].to_numpy(), h["b"].to_numpy(), h["r_med_geo"].to_numpy(), cfg)
    extinction.deredden(h, cfg, ebv=ebv, dist_col="r_med_geo")
    s.cut("main sequence (STILISM, published polynomial, no colour limit)",
          main_sequence.on_main_sequence(h["M_G0"], h["bp_rp0"], cfg))
    return s.rows


def main():
    cfg = config.load()
    w = A.wasp_hosts()
    g = gaia_rows(w["dr3"], "dr3x_hosts_wasp107")
    base = w.merge(g, on="dr3", how="left", validate="one_to_one")
    # archive RV columns for Method B (our 2026-10-28 snapshot, as frozen), joined on the DR2 id (E48)
    xm = pd.read_csv(A.CROSSMATCH, usecols=["archive_name_2019", "dr2_source_id"], dtype={"dr2_source_id": "Int64"})
    arc = pd.read_csv(config.out_dir(cfg) / "hosts_hj_table1_final.csv",
                      usecols=["source_id", "st_radv", "st_radvlim", "st_radv_reflink"], dtype={"source_id": "Int64"})
    xm = xm.merge(arc, left_on="dr2_source_id", right_on="source_id", how="left", validate="one_to_one")
    base_b = base.merge(xm.rename(columns={"archive_name_2019": "name"})[["name", "st_radv", "st_radvlim",
                                                                            "st_radv_reflink"]],
                        on="name", how="left", validate="one_to_one")
    nea_rv = A.nea_radv(w["dr3"])
    a_main, a_ext, a_left = method_a(base, nea_rv, exclude_companions=True)
    a_extra, _, _ = method_a(base, nea_rv, exclude_companions=False)
    b = method_b(base_b, cfg)
    out = {"what": "Stage 1 host steps, 107 WASP hosts (no kinematics)",
           "method_A_main_with_companion_exclusions": a_main,
           "method_A_extra_without_companion_exclusions": a_extra,
           "method_A_extinction_not_gspphot": a_ext,
           "method_B": b,
           "gspphot_available_of_107": int((base["ag_gspphot"].notna() & base["ebpminrp_gspphot"].notna()).sum())}
    A.write_json(A.OUT / "host_steps_wasp107.json", out)
    for title, rows in (("METHOD A (main: with close-companion exclusions)", a_main),
                        ("METHOD A (extra: without exclusions)", a_extra), ("METHOD B", b)):
        print(f"\n{title}")
        for r in rows:
            print(f"  {r['step']}: {r['n_after']}" + (f"  removed {r['removed']}" if r["removed"] else "")
                  + (f"\n    {r['note']}" if r.get("note") else ""))


if __name__ == "__main__":
    main()
