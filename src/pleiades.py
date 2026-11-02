"""Step 7 check: Pleiades members (Gaia Collab. 2018a) and the main-sequence refit."""
import numpy as np
import pandas as pd

from . import cuts, gaia_tap, main_sequence
from .config import RAW
from .gaia_columns import GAIA_COLUMNS

VIZIER_SERVERS = ["vizier.cds.unistra.fr", "vizier.cfa.harvard.edu", "vizier.iucaa.in", "vizier.nao.ac.jp"]


def members(cfg):
    """Pleiades rows of J/A+A/616/A10/tablea1a (source_id, optimised plx). Cached CSV."""
    out = RAW / "vizier" / "J_A+A_616_A10" / "tablea1a_pleiades.csv"
    if out.exists() and not cfg["refresh"]:
        return pd.read_csv(out)
    tsv = out.with_suffix(".tsv")                     # ASU-TSV download (CDS outage fallback)
    if tsv.exists():
        lines = [l for l in tsv.read_text().splitlines() if l and not l.startswith("#")]
        rows = [l.split("\t") for l in lines[3:]]     # header, units, dashes
        df = pd.DataFrame(rows, columns=[c.strip() for c in lines[0].split("\t")])
        df = df.rename(columns={"Source": "source_id", "plx": "plx_opt", "e_plx": "e_plx_opt"})
        df["source_id"] = df["source_id"].str.strip().astype("int64")
        for c in ("RAdeg", "DEdeg", "Gmag", "plx_opt", "e_plx_opt"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df["Cluster"] = df["Cluster"].str.strip()
        df.to_csv(out, index=False)
        return df
    from astroquery.vizier import Vizier
    last = None
    for server in VIZIER_SERVERS:
        try:
            v = Vizier(columns=["Source", "Cluster", "RAdeg", "DEdeg", "Gmag", "plx", "e_plx"],
                       row_limit=-1, vizier_server=server)
            v.TIMEOUT = 120
            t = v.query_constraints(catalog=f"{cfg['pleiades']['vizier_catalog']}/tablea1a",
                                    Cluster="Pleiades")[0]
            df = t.to_pandas().rename(columns={"Source": "source_id", "plx": "plx_opt",
                                               "e_plx": "e_plx_opt"})
            out.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(out, index=False)
            return df
        except Exception as err:          # CDS outage: try the next mirror
            last = err
    raise RuntimeError(f"VizieR unreachable on all mirrors: {last}")


def gaia_data(ids, cfg):
    g = cfg["gaia"]
    parts = []
    ids = sorted(int(i) for i in ids)
    for i in range(0, len(ids), 300):
        q = (f"SELECT {', '.join('g.' + c for c in GAIA_COLUMNS)} FROM {g['source_table']} AS g "
             f"WHERE g.source_id IN ({', '.join(map(str, ids[i:i + 300]))})")
        parts.append(gaia_tap.query(q, f"pleiades_{g['release']}/batch{i // 300:03d}",
                                    refresh=cfg["refresh"]).to_pandas())
    return pd.concat(parts, ignore_index=True)


def sample(cfg):
    """Members passing the Pleiades cut set, dereddened with the single E(B-V)."""
    m = members(cfg)
    df = m.merge(gaia_data(m["source_id"], cfg), on="source_id", how="left")
    p = cfg["pleiades"]
    mask, _ = cuts.apply(df, cfg["cuts"], cfg["cuts"]["apply"]["pleiades"])
    log = {"n_members": len(m), "n_after_cuts": int(mask.sum())}
    df = df[mask & ~df["source_id"].isin(p["exclude_source_ids"])].copy()
    log["n_after_wd_removal"] = len(df)
    ec = cfg["extinction"]
    df["bp_rp0"] = df["bp_rp"] - (ec["R_BP"] - ec["R_RP"]) * p["ebv"]
    df["M_G0"] = df["phot_g_mean_mag"] + 5 * np.log10(df["plx_opt"] / 100.0) - ec["R_G"] * p["ebv"]
    return df, log


def refit(cfg):
    df, log = sample(cfg)
    p = cfg["pleiades"]
    s = p["spline_smoothing"]
    if s == "auto":                      # residual sum ~ N sigma^2, sigma = robust scatter about the MS
        r = df["M_G0"] - main_sequence.ms_poly(df["bp_rp0"], cfg["main_sequence"]["coeffs"])
        sigma = 1.4826 * np.median(np.abs(r - np.median(r)))
        s = len(df) * sigma**2
        log["spline_sigma"] = float(sigma)
    coeffs, spl, grid = main_sequence.refit(df["bp_rp0"], df["M_G0"], order=p["poly_order"], s=s)
    g = np.linspace(0.44, 2.02, 200)
    pub = main_sequence.ms_poly(g, cfg["main_sequence"]["coeffs"])
    diff = main_sequence.ms_poly(g, coeffs) - pub
    log.update({"refit_coeffs": coeffs.tolist(), "max_abs_diff_FGK": float(np.max(np.abs(diff))),
                "rms_diff_FGK": float(np.sqrt(np.mean(diff**2)))})
    return df, coeffs, spl, log


PRINTED_CUTS = [1, 5, 8, 9, 10]          # HS19 Sec. 2 as printed (D10b)


def main():
    """Refit with the cut set used (D10b) and with the printed one; Fig. 1 for each.

        python -m src.pleiades
    """
    import json
    from . import config, figures
    cfg = config.load()
    out, fd = config.out_dir(cfg), config.out_dir(cfg, "figures")
    for tag, which in (("", cfg["cuts"]["apply"]["pleiades"]), ("_printed_cuts", PRINTED_CUTS)):
        v = config.variant(cfg, **{"cuts.apply.pleiades": list(which)})
        df, coeffs, spl, log = refit(v)
        log["cuts"] = list(which)
        log["color_range"] = [float(df["bp_rp0"].min()), float(df["bp_rp0"].max())]
        (out / f"pleiades_refit{tag}.json").write_text(json.dumps(log, indent=1))
        figures.fig_pleiades(df["bp_rp0"], df["M_G0"], cfg["main_sequence"]["coeffs"], coeffs, spl,
                             fd / f"fig1_pleiades{tag}.png",
                             title=f"Pleiades, cuts {', '.join(map(str, which))}: N = {len(df)}")
        print(tag or "used", {k: log[k] for k in ("n_members", "n_after_cuts", "n_after_wd_removal",
                                                  "max_abs_diff_FGK", "rms_diff_FGK", "color_range")})


if __name__ == "__main__":
    main()
