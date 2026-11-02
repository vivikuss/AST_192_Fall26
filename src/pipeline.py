"""End-to-end pipeline: samples -> cuts -> distances -> extinction -> MS -> UVW -> Monte Carlo.

    python -m src.pipeline                 # config.yaml as is
    python -m src.pipeline cuts.uwe.form=literal run_name=phase1_c4literal
"""
import hashlib
import json
import sys
import time

import numpy as np
import pandas as pd

from . import config, cuts, extinction, field, hosts, kinematics, main_sequence, matching
from .config import ROOT
from .gaia_tap import log

HOST_KINDS = ("hj_table1", "hj", "longp")


# ------------------------------------------------------------------ helpers
def _cache_ebv(df, cfg, tag):
    """E(B-V) is expensive for ~1e6 stars: cache by source_id per (release, distance, map, step)."""
    ec = cfg["extinction"]
    p = ROOT / "data" / "processed" / "cache"
    p.mkdir(parents=True, exist_ok=True)
    f = p / f"ebv_{tag}_{cfg['gaia']['release']}_{cfg['gaia']['distance_column']}_{ec['map']}_step{ec['step_pc']}.npz"
    sid = np.asarray(df["source_id"], dtype=np.int64)
    if f.exists() and not cfg["refresh"]:
        z = np.load(f)
        lut = pd.Series(z["ebv"], index=z["source_id"])
        if np.isin(sid, z["source_id"]).all():
            return lut.loc[sid].to_numpy()
    ebv = extinction.ebv_los(df["l"], df["b"], df["dist_pc"], cfg)
    np.savez(f, source_id=sid, ebv=ebv)
    return ebv


def prepare(df, cfg, tag, rv_col):
    """Distances, extinction, MS flag, kinematics, spectral bin (adds columns)."""
    df = df.copy()
    df["dist_pc"] = df[cfg["gaia"]["distance_column"]].astype(float)
    df["rv_used"] = df[rv_col].astype(float)
    extinction.deredden(df, cfg, ebv=_cache_ebv(df, cfg, tag))
    df["on_ms"] = main_sequence.on_main_sequence(df["M_G0"], df["bp_rp0"], cfg)
    ok = df["dist_pc"].notna() & df["rv_used"].notna()
    for c in ("U", "V", "W", "z"):
        df[c] = np.nan
    if ok.any():
        sub = kinematics.add_kinematics(df.loc[ok].copy(), cfg, rv_col="rv_used")
        for c in ("U", "V", "W", "z"):
            df.loc[ok, c] = sub[c]
    pos = df["dist_pc"].notna() & ~ok                 # z needs only a distance (e.g. binaries awaiting gamma)
    if pos.any():
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        c = SkyCoord(ra=df.loc[pos, "ra"].to_numpy() * u.deg, dec=df.loc[pos, "dec"].to_numpy() * u.deg,
                     distance=df.loc[pos, "dist_pc"].to_numpy() * u.pc)
        df.loc[pos, "z"] = kinematics.height_z(c, cfg["kinematics"])
    df["sptype"] = spectral_bin(df["bp_rp0"], cfg)
    return df


def spectral_bin(color0, cfg):
    c = np.asarray(color0, dtype=float)
    out = np.full(len(c), "", dtype=object)
    for name, (lo, hi) in cfg["spectral_bins"].items():
        out[(c > lo) & (c < hi)] = name          # strict, as printed in HS19 Sec. 3
    return out


def counts_fgk(df):
    return {k: int((df["sptype"] == k).sum()) for k in ("F", "G", "K")}


# ------------------------------------------------------------------ samples
_FIELD_CACHE = {}


def field_prepared(cfg):
    """Prefiltered field with extinction/MS/UVW for every star (cut-independent), memoised."""
    k = cfg["extinction"]
    key = (cfg["gaia"]["release"], cfg["gaia"]["distance_column"], k["map"], k["step_pc"],
           str(cfg["kinematics"]))                    # MS flag is recomputed per cfg in field_sample
    if key not in _FIELD_CACHE:
        _FIELD_CACHE.clear()
        _FIELD_CACHE[key] = prepare(field.download(cfg), cfg, "field", "radial_velocity")
    return _FIELD_CACHE[key]


def field_sample(cfg, exclude_ids=()):
    raw = field_prepared(cfg)
    raw["on_ms"] = main_sequence.on_main_sequence(raw["M_G0"], raw["bp_rp0"], cfg)
    df, cut_rows = field.select(raw, cfg)
    df = df[~df["source_id"].isin(set(exclude_ids))]
    log = {"n_prefilter": len(raw), "cuts": cut_rows, "n_after_cuts": len(df),
           "n_no_distance": int(df["dist_pc"].isna().sum()), "n_ms": int(df["on_ms"].sum())}
    ms = df[df["on_ms"] & df["U"].notna()].reset_index(drop=True)
    log["n_ms_with_uvw"] = len(ms)
    return ms, log, df


def host_sample(cfg, kind):
    return table_sample(cfg, hosts.build(cfg, kind), f"hosts_{kind}")


def target_sample(cfg, path, tag):
    """Any target list (src/targets.py CSV) through the same steps as the hosts."""
    from . import targets
    return table_sample(cfg, targets.build(cfg, path, tag), f"targets_{tag}")


def table_sample(cfg, h, tag):
    log = {"n_hosts": len(h), "n_in_gaia": int(h["in_gaia"].sum())}
    if "id_source" in h:
        log["id_source"] = h["id_source"].value_counts().to_dict()
    h = h[h["in_gaia"]].copy()
    mask, masks = cuts.apply(h, cfg["cuts"], cfg["cuts"]["apply"]["hosts"])
    log["cuts"] = cuts.cut_table(h, cfg["cuts"], cfg["cuts"]["apply"]["hosts"])
    h = h[mask].copy()
    h = prepare(h, cfg, tag, "rv")
    log["n_after_cuts"] = len(h)
    log["n_no_distance"] = int(h["dist_pc"].isna().sum())
    log["n_no_rv"] = int(h["rv_used"].isna().sum())
    log["rv_source_after_cuts"] = h["rv_source"].value_counts().to_dict()
    log["n_ms"] = int(h["on_ms"].sum())
    final = h[h["on_ms"] & h["U"].notna()].reset_index(drop=True)
    log["n_final"] = len(final)
    log["fgk"] = counts_fgk(final)
    log["rv_source_final"] = final["rv_source"].value_counts().to_dict()
    return final, log, h


# ------------------------------------------------------------------ Monte Carlo
def _run_chunked(sampler, host, fld, cfg, label, n_iter, statistic, ckpt_dir=None):
    """Run the Monte Carlo in chunks of monte_carlo.checkpoint_every iterations.

    Chunk c uses its own stream SeedSequence([seed, label code, c]), so results do
    not depend on interruptions. With ckpt_dir, each finished chunk is saved and
    reused on restart; the directory name hashes the config and the samples, so
    a changed input never reuses stale chunks (D19).
    """
    mc = cfg["monte_carlo"]
    every = mc["checkpoint_every"]
    lab = sum(map(ord, label))
    d = None
    if ckpt_dir is not None:
        key = hashlib.sha1(json.dumps([cfg["seed"], cfg["monte_carlo"], cfg["main_sequence"], label,
                                       np.asarray(host["source_id"], dtype=np.int64).tolist(), len(fld),
                                       float(np.nansum(fld["U"]))], default=_json).encode()).hexdigest()[:12]
        d = ckpt_dir / f"{label}_{key}"
        d.mkdir(parents=True, exist_ok=True)
    parts = []
    for c, s in enumerate(range(0, n_iter, every)):
        m = min(every, n_iter - s)
        f = d / f"chunk{c:04d}_{m}.npz" if d is not None else None
        if f is not None and f.exists():
            z = np.load(f)
            parts.append({k: z[k] for k in z.files})
            continue
        rng = np.random.default_rng(np.random.SeedSequence([cfg["seed"], lab, c]))
        r = sampler.run(m, rng, statistic)
        part = {"disp_rms": r["dispersion_all"]["rms"], "disp_mean_abs": r["dispersion_all"]["mean_abs"],
                "host_rms": r["host_dispersion_all"]["rms"], "host_mean_abs": r["host_dispersion_all"]["mean_abs"],
                "n_unmatched": r["n_unmatched"], "host_unmatched_count": r["host_unmatched_count"]}
        if f is not None:
            np.savez(f, **part)
        parts.append(part)
        log(f"[{label}] chunk {c + 1}/{-(-n_iter // every)} done ({s + m}/{n_iter} iterations)")
    cat = {k: np.concatenate([p[k] for p in parts]) for k in parts[0] if k != "host_unmatched_count"}
    hfull = {k: float(kinematics.dispersion(host["U"], host["V"], host["W"], k)) for k in ("rms", "mean_abs")}
    return {"dispersion": cat[f"disp_{statistic}"],
            "dispersion_all": {k: cat[f"disp_{k}"] for k in ("rms", "mean_abs")},
            "host_dispersion": cat[f"host_{statistic}"],
            "host_dispersion_all": {k: cat[f"host_{k}"] for k in ("rms", "mean_abs")},
            "host_dispersion_full": hfull, "n_unmatched": cat["n_unmatched"],
            "host_unmatched_count": sum(p["host_unmatched_count"] for p in parts)}


def monte_carlo(host, fld, cfg, label, n_iter=None, statistic=None, ckpt_dir=None):
    mc = cfg["monte_carlo"]
    n_iter = n_iter or mc["n_iter"]
    statistic = statistic or mc["statistic"]
    t0 = time.time()
    ncand = matching.candidates(host, fld, mc["color_tol"])
    excluded = []
    if mc["exclude_unmatchable_hosts"]:          # sensitivity (D13b): hosts with no field star of their colour
        drop = ncand < mc["min_field_candidates"]
        excluded = host.loc[drop, ["hostname", "bp_rp0"]].to_dict("records") if "hostname" in host else []
        host, ncand = host[~drop].reset_index(drop=True), ncand[~drop]
    sampler = matching.ControlSampler(host, fld, cfg)
    res = _run_chunked(sampler, host, fld, cfg, label, n_iter, statistic, ckpt_dir)
    d = res["dispersion"]
    p, k = matching.p_value(res["host_dispersion"], d)              # paired per iteration
    pa = {s: matching.p_value(res["host_dispersion_all"][s], res["dispersion_all"][s]) for s in ("rms", "mean_abs")}
    names = host["hostname"].astype(str).to_numpy() if "hostname" in host else np.arange(len(host)).astype(str)
    hu = res["host_unmatched_count"]
    unmatched = [{"host": names[i], "bp_rp0": round(float(host["bp_rp0"].iloc[i]), 3),
                  "z_pc": round(float(host["z"].iloc[i])), "field_candidates": int(ncand[i]),
                  "frac_iters_unmatched": round(hu[i] / n_iter, 4)}
                 for i in np.argsort(-hu) if hu[i] > 0]
    summary = {
        "label": label, "n_host": len(host), "n_iter": n_iter, "statistic": statistic,
        "host_dispersion": res["host_dispersion_full"][statistic],
        "host_dispersion_all": res["host_dispersion_full"],
        "host_dispersion_iter_mean": float(np.mean(res["host_dispersion"])),
        "control_median": float(np.median(d)),
        "control_std": float(np.std(d, ddof=1)),
        # Gaussian-equivalent significance: (control median - host) / control std; assumes Gaussian tails
        "z_gauss": float((np.median(d) - res["host_dispersion_full"][statistic]) / np.std(d, ddof=1)),
        "control_pct": {q: float(np.percentile(d, q)) for q in (2, 16, 50, 84, 98)},
        "control_pct_all": {s: {q: float(np.percentile(res["dispersion_all"][s], q)) for q in (2, 16, 50, 84, 98)}
                            for s in ("rms", "mean_abs")},
        "p_value": p, "n_le": k, "p_value_all": {s: v[0] for s, v in pa.items()},
        "n_usable_field": sampler.n_usable, "mean_zmatched_subsample": sampler.mean_subsample,
        "hosts_in_empty_zbins": sampler.bin_info["n_hosts_in_empty_bins"],
        "iters_with_unmatched_host": int((res["n_unmatched"] > 0).sum()),
        "mean_unmatched_per_iter": float(res["n_unmatched"].mean()),
        "unmatched_hosts": unmatched, "excluded_hosts": excluded,
        "z_edges": [e.tolist() for e in sampler.edges], "seconds": round(time.time() - t0, 1),
    }
    return summary, res


# ------------------------------------------------------------------ run
def run(cfg, kinds=HOST_KINDS, n_iter=None, save=True):
    out = config.out_dir(cfg)
    t0 = time.time()
    samples, logs, full = {}, {}, {}
    for kind in kinds:
        samples[kind], logs[kind], full[kind] = host_sample(cfg, kind)
        print(f"[{kind}] final {logs[kind]['n_final']}  FGK {logs[kind]['fgk']}", flush=True)
    exclude = set()
    for h in full.values():
        exclude |= set(h["source_id"].dropna().astype(np.int64))
    fld, logs["field"], fld_all = field_sample(cfg, exclude)
    print(f"[field] after cuts {logs['field']['n_after_cuts']}  MS {logs['field']['n_ms_with_uvw']}", flush=True)
    logs["field"]["fgk"] = counts_fgk(fld)
    mc = {}
    for kind in kinds:
        mc[kind], res = monte_carlo(samples[kind], fld, cfg, kind, n_iter,
                                    ckpt_dir=out / "mc_checkpoints" if save else None)
        print(f"[{kind}] sigma_host={mc[kind]['host_dispersion']:.2f} control median="
              f"{mc[kind]['control_median']:.2f} p={mc[kind]['p_value']:.2e} ({mc[kind]['seconds']} s)", flush=True)
        if save:
            np.savez(out / f"mc_{kind}.npz", **{f"disp_{s}": v for s, v in res["dispersion_all"].items()},
                     **{f"host_disp_{s}": v for s, v in res["host_dispersion_all"].items()},
                     n_unmatched=res["n_unmatched"], host_unmatched_count=res["host_unmatched_count"])
    if save:
        for kind, h in samples.items():
            h.to_csv(out / f"hosts_{kind}_final.csv", index=False)
            full[kind].to_csv(out / f"hosts_{kind}_after_cuts.csv", index=False)
        fld[["source_id", "bp_rp0", "M_G0", "U", "V", "W", "z", "dist_pc", "ebv", "sptype"]].to_csv(
            out / "field_ms.csv.gz", index=False)
        with open(out / "summary.json", "w") as f:
            json.dump({"config": cfg, "logs": logs, "mc": mc,
                       "runtime_s": round(time.time() - t0, 1)}, f, indent=1, default=_json)
    return {"samples": samples, "field": fld, "field_all": fld_all, "logs": logs, "mc": mc,
            "full": full}


def _json(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def _parse_overrides(argv):
    ov = {}
    for a in argv:
        k, v = a.split("=", 1)
        try:
            import yaml
            v = yaml.safe_load(v)
        except Exception:
            pass
        ov[k] = v
    return ov


if __name__ == "__main__":
    run(config.load(**_parse_overrides(sys.argv[1:])))
