"""EXPLORATORY (not pre-registered; DECISIONS E71). Metallicity-matched Method B (HS19-style).

    python -m src.dr3x_feh_b field    # rebuild the E62 diagnostic field (124,142) -> data/processed/dr3_exploratory/
    python -m src.dr3x_feh_b fetch    # slow resumable GSP-Spec fetch for field stars not fetched in E68
    python -m src.dr3x_feh_b run      # -> data/processed/dr3_exploratory/feh_b_results.json (aggregates)

Field: DR3 Method B field restricted to DR2-field stars (E62 item 2). Hosts: the 101 Method B hosts
(dr3x_diag.hosts). [M/H]: GSP-Spec, E66 calibration and medium quality, hosts and field alike.
Runs on the hosts with a GSP-Spec [M/H]: dr3x_B_feh_none (Method B as the diagnostic), dr3x_B_feh_02 (control
within +-0.2 dex of its host, main) and dr3x_B_feh_01 (+-0.1 dex, check). 40,000 controls + 40,000 bootstrap.
Host names dropped for an empty pool go to data/private/feh_b/ only.
"""
import json
import sys

import numpy as np
import pandas as pd

from . import config, kinematics, matching, eblm_kinematics as ek
from . import dr3x_diag as D
from . import dr3x_method_a as A
from . import dr3x_stage2 as S
from . import feh_field as F
from . import seeds
from .feh_compare import calibrate, quality

HEADER = "EXPLORATORY — not pre-registered (DECISIONS E71): metallicity-matched Method B"
FIELD = A.OUT / "feh_b_field.csv.gz"                       # E62 diagnostic field (public Gaia-derived; gitignored)
PREFIX = "feh/fieldB_dr2_gspspec"
STATUS = F.PROC / "fetch_status_fieldB.json"
GSP_OUT = F.PROC / "field_gspspec_fieldB.csv.gz"
LOG = config.ROOT / "logs" / "feh_fieldB_fetch.log"
GIVE_UP_IF_NO_SUCCESS_BY = "2026-11-19T09:00:00"
OUT = A.OUT / "feh_b_results.json"
PRIV = config.ROOT / "data" / "private" / "feh_b"
N = 40000
RUNS = (("dr3x_B_feh_none", None), ("dr3x_B_feh_02", 0.2), ("dr3x_B_feh_01", 0.1))


# ----------------------------------------------------------------------------- field and fetch
def build_field():
    cfg = config.load()
    fb = S.load_field_b()
    f2 = pd.read_csv(config.out_dir(cfg) / "field_ms.csv.gz")
    ids, info = D.dr2_field_dr3_ids(cfg, f2)                # cached ARI crossmatch (E62)
    fr = fb[fb["source_id"].isin(ids)].reset_index(drop=True)
    assert len(fr) == 124142, len(fr)
    fr.to_csv(FIELD, index=False, compression="gzip")
    print(len(fr), info)


def fetch_chunks():
    ids = pd.read_csv(FIELD, usecols=["source_id"])["source_id"].to_numpy(np.int64)
    have = set(np.concatenate(F.pool_chunks()).tolist())
    new = np.sort(np.array([i for i in ids if i not in have], dtype=np.int64))
    print(f"field {len(ids)}: {len(ids) - len(new)} already fetched in E68, {len(new)} to fetch", flush=True)
    return [new[k:k + F.POOL_CHUNK] for k in range(0, len(new), F.POOL_CHUNK)]


def fetch():
    F.fetch_ids(fetch_chunks(), PREFIX, STATUS, GSP_OUT, LOG, GIVE_UP_IF_NO_SUCCESS_BY)


# ----------------------------------------------------------------------------- [M/H]
def mh_table(g):
    g = g.copy()
    med = g["flags_gspspec"].map(lambda f: quality(f)[0])
    g["mh"] = calibrate(g["mh_gspspec"], g["logg_gspspec"]).where(g["mh_gspspec"].notna() & med)
    return g


def field_mh(fld):
    g = pd.concat([pd.read_csv(f, dtype={"flags_gspspec": str}) for f in
                   (F.PROC / "field_gspspec_pool.csv.gz", GSP_OUT)], ignore_index=True)
    g = mh_table(g.drop_duplicates("source_id"))
    m = fld[["source_id"]].merge(g[["source_id", "mh_gspspec", "mh"]], on="source_id", how="left", validate="1:1")
    cov = {"field": len(fld), "with_ap_row": int(fld["source_id"].isin(g["source_id"]).sum()),
           "with_mh_gspspec": int(m["mh_gspspec"].notna().sum()), "with_medium_mh": int(m["mh"].notna().sum())}
    return m["mh"].to_numpy(float), cov


def host_mh(hb):
    m = pd.read_csv(config.ROOT / "data" / "private" / "feh" / "feh_gspspec_merged.csv",
                    usecols=["gaia_dr3_source_id", "gsp"], dtype={"gaia_dr3_source_id": str})
    m = m.set_index("gaia_dr3_source_id")["gsp"]
    return hb["dr3"].astype(np.int64).astype(str).map(m).to_numpy(float)


# ----------------------------------------------------------------------------- sampler
class FehSampler(matching.ControlSampler):
    """ControlSampler with an extra step-(b) condition |[M/H]_control - [M/H]_host| <= window (None = no window).
    Step (a) is unchanged. With window=None it consumes the random stream exactly as ControlSampler."""

    def __init__(self, host, field, cfg, host_mh=None, field_mh=None, window=None):
        super().__init__(host, field, cfg)
        assert self.method == "thinning"
        lo = np.searchsorted(self.fc, self.hc - self.tol, side="left")
        hi = np.searchsorted(self.fc, self.hc + self.tol, side="right")
        fmh = None if window is None else np.asarray(field_mh, dtype=float)[self.field_idx]
        self.cand = []
        for h in range(self.n_host):
            pos = np.arange(lo[h], hi[h])
            if window is not None:
                pos = pos[np.abs(fmh[pos] - host_mh[h]) <= window + 1e-9]       # NaN never passes
            self.cand.append(pos)
        self.pool_sizes = np.array([len(c) for c in self.cand])

    def control(self, rng):
        keep = rng.random(self.n_usable) < self.p_keep                         # = ControlSampler.subsample
        cs = [c[keep[c]] for c in self.cand]
        n = np.array([len(c) for c in cs])
        pick = np.full(self.n_host, -1, dtype=np.int64)
        has = np.flatnonzero(n > 0)
        k = (rng.random(len(has)) * n[has]).astype(np.int64)
        for j, h in enumerate(has):
            pick[h] = cs[h][k[j]]
        if self.without_replacement:
            pick = self._dedupe_c(pick, cs, rng)
        return pick

    def _dedupe_c(self, pick, cs, rng):
        order = rng.permutation(self.n_host)
        p = pick[order]
        _, first = np.unique(p, return_index=True)
        dup = np.ones(len(p), dtype=bool)
        dup[first] = False
        dup &= p >= 0
        if not dup.any():
            return pick
        used = set(p[~dup & (p >= 0)].tolist())
        for j in np.flatnonzero(dup):
            h = order[j]
            cand = [c for c in cs[h].tolist() if c not in used]
            if cand:
                c = cand[int(rng.random() * len(cand))]
                pick[h] = c
                used.add(c)
            else:
                pick[h] = -1
        return pick


def same_stream_check(h, fld, cfg, n=200):
    """FehSampler(window=None) must reproduce ControlSampler's picks exactly."""
    a, b = matching.ControlSampler(h, fld, cfg), FehSampler(h, fld, cfg)
    ra, rb = (np.random.default_rng(np.random.SeedSequence([cfg["seed"], 12345, 0])) for _ in range(2))
    return all(np.array_equal(a.control(ra), b.control(rb)) for _ in range(n))


# ----------------------------------------------------------------------------- statistics
def _stats(U, V, W):
    return {"hs19_eq1": kinematics.dispersion(U, V, W, "mean_abs"), "sigma_tot": kinematics.dispersion(U, V, W, "rms"),
            "sigma_U": np.std(U, axis=-1), "sigma_V": np.std(V, axis=-1), "sigma_W": np.std(W, axis=-1)}


def mc(smp, cfg, label, n=N):
    every = cfg["monte_carlo"]["checkpoint_every"]
    lab = seeds.label_code(label)
    picks = []
    for c, s in enumerate(range(0, n, every)):
        rng = np.random.default_rng(np.random.SeedSequence([cfg["seed"], lab, c]))     # as pipeline._run_chunked
        picks.append(np.stack([smp.control(rng) for _ in range(min(every, n - s))]))
        S.log(f"[{label}] chunk {c + 1}/{-(-n // every)}")
    P = np.concatenate(picks)
    miss = P < 0
    ok = ~miss.any(axis=1)
    ctl = {k: np.empty(n) for k in ("hs19_eq1", "sigma_tot", "sigma_U", "sigma_V", "sigma_W")}
    hst = {k: np.empty(n) for k in ctl}
    full = _stats(smp.hU, smp.hV, smp.hW)
    pc = P[ok]
    for k, v in _stats(smp.fU[pc], smp.fV[pc], smp.fW[pc]).items():
        ctl[k][ok] = v
        hst[k][ok] = full[k]
    for i in np.flatnonzero(~ok):                       # D13b: matched hosts only, for control and host alike
        m = ~miss[i]
        p = P[i][m]
        for k, v in _stats(smp.fU[p], smp.fV[p], smp.fW[p]).items():
            ctl[k][i] = v
        for k, v in _stats(smp.hU[m], smp.hV[m], smp.hW[m]).items():
            hst[k][i] = v
    return P, ctl, hst, {k: float(v) for k, v in full.items()}


def run_one(h, fld, cfg, label, window, hmh, fmh):
    smp = FehSampler(h, fld, cfg, hmh, fmh, window)
    empty = np.flatnonzero(smp.pool_sizes == 0)
    if len(empty):                                      # drop never-matchable hosts; z weights from the rest
        keep = smp.pool_sizes > 0
        h, hmh = h[keep].reset_index(drop=True), hmh[keep]
        smp = FehSampler(h, fld, cfg, hmh, fmh, window)
        assert (smp.pool_sizes > 0).all()
    P, ctl, hst, full = mc(smp, cfg, label)
    r = {"label": label, "window_dex": window, "n_hosts": len(h), "n_dropped_empty_pool": int(len(empty)),
         "pool_size_median": float(np.median(smp.pool_sizes)), "pool_size_min": int(smp.pool_sizes.min()),
         "n_usable_field": smp.n_usable, "hosts_in_empty_zbins": smp.bin_info["n_hosts_in_empty_bins"],
         "iters_with_unmatched_host": int((P < 0).any(axis=1).sum()),
         "mean_unmatched_per_iter": float((P < 0).sum(axis=1).mean()), "n_iter": N, "bootstrap_n": N}
    for k in ctl:
        r[k] = {"host": full[k], "control_median": float(np.median(ctl[k])),
                "control_16_84": [float(x) for x in np.percentile(ctl[k], [16, 84])],
                "fraction_controls_le_host": float(np.mean(ctl[k] <= hst[k]))}
    for k, stat in (("hs19_eq1", "mean_abs"), ("sigma_tot", "rms")):
        b = ek.bootstrap(h, r[k]["control_median"], cfg, N, label, stat)
        r[k].update({"delta": b["delta_hat"], "ci95": b["ci95"], "bca_needed": b["bca_needed"],
                     **({"bca_ci95": b["bca_ci95"]} if "bca_ci95" in b else {})})
    return r, empty


def run():
    cfg = config.load()
    hb, _ = D.hosts(cfg)
    fld = pd.read_csv(FIELD)
    assert len(fld) == 124142
    _, err = D.replay(hb, fld, cfg, "dr3x_B_diag_dr2field")            # guard: E62 chunk 0 reproduces exactly
    S.log(f"guard: E62 dr3x_B_diag_dr2field chunk 0 replayed, max diff {err:.1e}")
    fmh, cov = field_mh(fld)
    hmh_all = host_mh(hb)
    has = np.isfinite(hmh_all)
    h = hb[has].reset_index(drop=True)
    hmh = hmh_all[has]
    assert same_stream_check(h, fld, cfg), "FehSampler(window=None) differs from ControlSampler"
    lo, hi = np.nanmin(hmh), np.nanmax(hmh)
    cov["with_medium_mh_in_host_range"] = int(((fmh >= lo) & (fmh <= hi)).sum())
    R = {"header": HEADER, "guard_replay_max_abs_diff": err, "same_stream_check": True,
         "hosts_method_b": len(hb), "hosts_with_gsp_medium": int(has.sum()), "field_coverage": cov}
    dropped = {}
    for label, w in RUNS:
        R[label], e = run_one(h, fld, cfg, label, w, hmh, fmh)
        dropped[label] = h["name"].iloc[e].tolist()
        r = R[label]
        S.log(f"{label}: N {r['n_hosts']}, pool median {r['pool_size_median']:.0f} min {r['pool_size_min']}, "
              f"eq1 {r['hs19_eq1']['host']:.2f} vs {r['hs19_eq1']['control_median']:.2f}, "
              f"Delta {r['hs19_eq1']['delta']:.2f} {np.round(r['hs19_eq1']['ci95'], 2).tolist()}")
    A.write_json(OUT, R)
    PRIV.mkdir(parents=True, exist_ok=True)
    (PRIV / "dropped_empty_pool.json").write_text(json.dumps(dropped, indent=1))
    return R


if __name__ == "__main__":
    {"field": build_field, "fetch": fetch, "run": run}[sys.argv[1]]()
