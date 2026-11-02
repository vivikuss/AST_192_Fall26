"""Stage A: step-by-step survivor counts for the field, HJ-host and long-period-host samples.

    python -u -m src.counts      # writes data/processed/<run>/count_table.csv (+ .md)

Counts are sequential: each row is the number left after that step and every
step above it. The final rows equal what pipeline.run feeds the Monte Carlo.
"""
import numpy as np
import pandas as pd

from . import config, cuts, hosts, main_sequence, pipeline
from .gaia_tap import log

PAPER = {"field": 385036, "hj_table1": 338, "longp": 367, "F": 119, "G": 146, "K": 66}
CUT_NAMES = {
    1: "cut 1: parallax_over_error > 10",
    2: "cut 2: 5-parameter solution / mean_varpi_factor_al range",
    3: "cut 3: visibility_periods_used > 8",
    4: "cut 4: unit weight error (Lindegren C.1)",
    5: "cut 5: astrometric_gof_al < 3",
    6: "cut 6: astrometric_excess_noise_sig < 2",
    7: "cut 7: rv_nb_transits > 5",
    8: "cut 8: BP/RP excess factor (Lindegren C.2)",
    9: "cut 9: phot_bp_mean_flux_over_error > 10",
    10: "cut 10: phot_rp_mean_flux_over_error > 10",
}


def _tail(df, cfg, rows, rv_step):
    """Distance -> extinction -> MS -> (RV/UVW) -> F/G/K, appending (step, n)."""
    m = df["dist_pc"].notna() & (df["dist_pc"] > 0)
    rows.append(("BJ18 distance available", int(m.sum())))
    m &= np.isfinite(df["M_G0"]) & np.isfinite(df["bp_rp0"])
    rows.append(("extinction/reddening corrected (finite M_G0, (BP-RP)_0)", int(m.sum())))
    m &= df["on_ms"].astype(bool)
    rows.append(("main-sequence cut (M_G0 > MS - 1 mag)", int(m.sum())))
    if rv_step:
        m &= df["rv_used"].notna()
        rows.append(("RV available (UVW computable)", int(m.sum())))
    m &= df["U"].notna()
    rows.append(("final sample (UVW computed)", int(m.sum())))
    sp = df.loc[m, "sptype"]
    for k in ("F", "G", "K"):
        rows.append((f"  {k} bin", int((sp == k).sum())))
    rows.append(("  outside F/G/K bins", int((sp == "").sum())))
    return rows


def _sequential_cuts(df, cfg, which, rows, prefix_mask=None):
    m = np.ones(len(df), bool) if prefix_mask is None else prefix_mask.copy()
    _, masks = cuts.apply(df, cfg["cuts"], which)
    for i in sorted(masks):
        m &= masks[i]
        rows.append((CUT_NAMES[i], int(m.sum())))
    return m


def host_counts(cfg, kind):
    h = hosts.build(cfg, kind)
    rows = [("listed", len(h))]
    rows.append(("Gaia DR2 source_id found", int(h["source_id"].notna().sum())))
    h = h[h["in_gaia"]].reset_index(drop=True)
    rows.append(("row returned from gaiadr2.gaia_source", len(h)))
    m = _sequential_cuts(h, cfg, cfg["cuts"]["apply"]["hosts"], rows)
    h = pipeline.prepare(h[m].copy(), cfg, f"hosts_{kind}", "rv")
    return _tail(h, cfg, rows, rv_step=True), set(h["source_id"].dropna().astype(np.int64))


def field_counts(cfg, exclude_ids):
    raw = pipeline.field_prepared(cfg)
    raw["on_ms"] = main_sequence.on_main_sequence(raw["M_G0"], raw["bp_rp0"], cfg)
    rows = [("pre-filtered download (server: RV, plx, cuts 1-3, 5-7, 9, 10; loose 8)", len(raw))]
    m = _sequential_cuts(raw, cfg, cfg["cuts"]["apply"]["field"], rows)
    m &= ~raw["source_id"].isin(exclude_ids).to_numpy()
    rows.append(("excluding planet hosts", int(m.sum())))
    return _tail(raw[m], cfg, rows, rv_step=False)


def table(cfg):
    cols, excl = {}, set()
    for kind in ("hj_table1", "hj", "longp"):
        cols[kind], ids = host_counts(cfg, kind)
        excl |= ids
        log(f"{kind}: final {dict(cols[kind])['final sample (UVW computed)']}")
    cols["field"] = field_counts(cfg, excl)
    log(f"field: final {dict(cols['field'])['final sample (UVW computed)']}")
    pre = ["pre-filtered download (server: RV, plx, cuts 1-3, 5-7, 9, 10; loose 8)", "listed",
           "Gaia DR2 source_id found", "row returned from gaiadr2.gaia_source"]
    cut_rows = [CUT_NAMES[i] for i in sorted(CUT_NAMES)]
    tail = ["excluding planet hosts", "BJ18 distance available",
            "extinction/reddening corrected (finite M_G0, (BP-RP)_0)", "main-sequence cut (M_G0 > MS - 1 mag)",
            "RV available (UVW computable)", "final sample (UVW computed)",
            "  F bin", "  G bin", "  K bin", "  outside F/G/K bins"]
    order = pre + cut_rows + tail
    df = pd.DataFrame({"step": order})
    for c in ("field", "hj_table1", "hj", "longp"):
        d = dict(cols[c])
        df[c] = [d.get(s, pd.NA) for s in order]
    df = df.rename(columns={"hj_table1": "hj_hosts_table1", "hj": "hj_hosts_archive_le2019",
                            "longp": "longp_hosts"})
    paper = {"final sample (UVW computed)": (PAPER["field"], PAPER["hj_table1"], PAPER["longp"]),
             "  F bin": (pd.NA, PAPER["F"], pd.NA), "  G bin": (pd.NA, PAPER["G"], pd.NA),
             "  K bin": (pd.NA, PAPER["K"], pd.NA)}
    df["paper_field"] = [paper.get(s, (pd.NA,) * 3)[0] for s in order]
    df["paper_hj"] = [paper.get(s, (pd.NA,) * 3)[1] for s in order]
    df["paper_longp"] = [paper.get(s, (pd.NA,) * 3)[2] for s in order]
    return df


def markdown(df):
    cells = df.astype(object).where(df.notna(), "")
    lines = ["| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)]
    lines += ["| " + " | ".join(str(v) for v in r) + " |" for r in cells.itertuples(index=False)]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    cfg = config.load()
    df = table(cfg)
    out = config.out_dir(cfg)
    df.to_csv(out / "count_table.csv", index=False)
    (out / "count_table.md").write_text(markdown(df))
    print(df.to_string(index=False))
