"""Download the Gaia field prefilter (ARI source_id ranges, or ESA HEALPix chunks), resuming
from the cache, then verify it.

    python -u -m src.download_field             # status, download, status
    python -u -m src.download_field --status    # ranges complete / remaining only
    python -u -m src.download_field --time L P  # time one ARI range (HEALPix level L, pixel P)
    python -u -m src.download_field --verify    # per-range COUNT(*) vs cached rows + duplicate check
    python -u -m src.download_field --repair    # verify, delete mismatched ranges; rerun to refetch
"""
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from astropy.table import Table

from . import ari_tap, config, field, gaia_tap
from .config import RAW, ROOT
from .gaia_tap import log


def status(cfg):
    if cfg["gaia"]["field_service"] == "ari":
        lev = cfg["gaia"]["ari"]["range_level"]
        done, todo = ari_tap.status(field.field_adql(cfg, "ari"), field.chunk_name(cfg), lev)
        log(f"ARI level-{lev} ranges: {len(done)} complete, {len(todo)} remaining of {len(done) + len(todo)}")
    else:
        done, todo = gaia_tap.chunk_status(field.field_adql(cfg), field.chunk_name(cfg), level=2)
        log(f"ESA level-2 regions: {len(done)} complete, {len(todo)} remaining of {len(done) + len(todo)}")
    if todo:
        log(f"remaining: {todo}")
    return done, todo


def download(cfg):
    t0 = time.time()
    df, failed = field.download(cfg, allow_missing=True)
    log(f"{len(df)} rows in {time.time() - t0:.0f} s")
    out = ROOT / "logs" / f"{field.chunk_name(cfg)}_failed.json"
    out.write_text(json.dumps(failed, indent=1))
    if failed:
        log(f"FAILED ranges ({len(failed)}); rerun to retry them:")
        for f in failed:
            log(f"  {f}")
    else:
        log("no failed ranges")
    return df, failed


def time_range(cfg, level, pix):
    """Run (and cache) one ARI range, reporting wall time and rows."""
    t0 = time.time()
    name = f"{field.chunk_name(cfg)}/L{level}_{pix:05d}"
    tbl = ari_tap.query(gaia_tap.chunk_query(field.field_adql(cfg, "ari"), level, pix), name, refresh=cfg["refresh"])
    log(f"TIMING L{level} pixel {pix}: {len(tbl)} rows in {time.time() - t0:.0f} s")


def verify_ari(cfg, delete_bad=False):
    """COUNT(*) on ARI for every cached range (same WHERE) vs cached rows; duplicate source_ids."""
    name, lev = field.chunk_name(cfg), cfg["gaia"]["ari"]["range_level"]
    done, todo = ari_tap.status(field.field_adql(cfg, "ari"), name, lev)
    if todo:
        log(f"{len(todo)} ranges not downloaded yet; verify after the download finishes")
        return False
    ctmpl = field.count_adql(cfg, "ari")

    def one(p):
        n_srv = int(ari_tap.query(gaia_tap.chunk_query(ctmpl, lev, p), f"{name}_count/L{lev}_{p:05d}",
                                  refresh=cfg["refresh"])["n"][0])
        return p, n_srv, len(Table.read(RAW / "ari" / name / f"L{lev}_{p:05d}.fits"))

    with ThreadPoolExecutor(cfg["gaia"]["ari"]["workers"]) as ex:
        res = list(ex.map(one, done))
    bad = [r for r in res if r[1] != r[2]]
    n_srv, n_loc = sum(r[1] for r in res), sum(r[2] for r in res)
    log(f"{len(res)} ARI ranges: server COUNT(*) = {n_srv}, cached rows = {n_loc}, mismatched = {len(bad)}")
    ids = np.concatenate([np.asarray(Table.read(RAW / "ari" / name / f"L{lev}_{p:05d}.fits")["source_id"],
                                     dtype=np.int64) for p in done])
    n_dup = len(ids) - len(np.unique(ids))
    log(f"duplicate source_ids: {n_dup}")
    for p, a, b in bad:
        log(f"  L{lev} pixel {p}: server {a} vs cached {b}")
        if delete_bad:
            for ext in ("fits", "adql"):
                (RAW / "ari" / name / f"L{lev}_{p:05d}.{ext}").unlink()
    summary = {"service": cfg["gaia"]["ari"]["tap_url"], "range_level": lev, "ranges": len(res),
               "server_count": n_srv, "cached_rows": n_loc, "per_range": [list(r) for r in res],
               "mismatched": [list(r) for r in bad], "duplicate_source_ids": int(n_dup)}
    (ROOT / "logs" / f"{name}_verify.json").write_text(json.dumps(summary, indent=1))
    return not bad and n_dup == 0


def verify(cfg, delete_bad=False):
    """COUNT(*) on the server for exactly the chunks in the cache; duplicates across chunks."""
    if cfg["gaia"]["field_service"] == "ari":
        return verify_ari(cfg, delete_bad)
    tmpl, name = field.field_adql(cfg), field.chunk_name(cfg)
    files = []
    for p in range(12 * 4**2):
        f = gaia_tap.chunk_files(tmpl, name, 2, p, gaia_tap.net()["max_level"])
        if f is None:
            log(f"L2 pixel {p} incomplete; verify after the download finishes")
            return False
        files += f
    ctmpl = field.count_adql(cfg)

    def one(lp):
        lev, pix = lp
        n_srv = int(gaia_tap.query(gaia_tap.chunk_query(ctmpl, lev, pix), f"{name}_count/L{lev}_{pix:05d}",
                                   refresh=cfg["refresh"])["n"][0])
        n_loc = len(Table.read(RAW / "gaia" / name / f"L{lev}_{pix:05d}.fits"))
        return lev, pix, n_srv, n_loc

    with ThreadPoolExecutor(gaia_tap.net()["workers"]) as ex:
        res = list(ex.map(one, files))
    bad = [r for r in res if r[2] != r[3]]
    n_srv, n_loc = sum(r[2] for r in res), sum(r[3] for r in res)
    log(f"{len(files)} chunks: server COUNT(*) = {n_srv}, cached rows = {n_loc}, mismatched chunks = {len(bad)}")
    ids = np.concatenate([np.asarray(Table.read(RAW / "gaia" / name / f"L{l}_{p:05d}.fits")["source_id"],
                                     dtype=np.int64) for l, p in files])
    n_dup = len(ids) - len(np.unique(ids))
    log(f"duplicate source_ids: {n_dup}")
    for r in bad:
        log(f"  L{r[0]} pixel {r[1]}: server {r[2]} vs cached {r[3]}")
        if delete_bad:
            for ext in ("fits", "adql"):
                (RAW / "gaia" / name / f"L{r[0]}_{r[1]:05d}.{ext}").unlink()
    if delete_bad and bad:
        log(f"deleted {len(bad)} mismatched chunks; rerun python -u -m src.download_field")
    summary = {"chunks": len(files), "server_count": n_srv, "cached_rows": n_loc,
               "mismatched": [list(r) for r in bad], "duplicate_source_ids": int(n_dup)}
    (ROOT / "logs" / f"{name}_verify.json").write_text(json.dumps(summary, indent=1))
    return not bad and n_dup == 0


if __name__ == "__main__":
    cfg = config.load()
    ari_tap.install_signal_handlers()          # SIGTERM / Ctrl-C abort running ARI jobs
    if "--status" in sys.argv:
        status(cfg)
    elif "--time" in sys.argv:
        i = sys.argv.index("--time")
        time_range(cfg, int(sys.argv[i + 1]), int(sys.argv[i + 2]))
    elif "--verify" in sys.argv or "--repair" in sys.argv:
        sys.exit(0 if verify(cfg, delete_bad="--repair" in sys.argv) else 1)
    else:
        status(cfg)
        download(cfg)
        status(cfg)
