"""Cached ADQL async queries against the Gaia archive (astroquery).

Every query has a socket timeout and a wall-clock job timeout, and is retried
with increasing waits on timeouts, dropped connections and HTTP 500 (A5).
"""
import hashlib
import socket
import time
from datetime import datetime

from astropy.table import Table
from astroquery.gaia import GaiaClass

from . import config
from .config import RAW

_NET = None


def net():
    """Network settings (config.yaml gaia.network); also sets the socket timeout."""
    global _NET
    if _NET is None:
        _NET = config.load()["gaia"]["network"]
        socket.setdefaulttimeout(_NET["socket_timeout_s"])   # astroquery's http.client has none
    return _NET


def log(msg):
    print(f"{datetime.now():%H:%M:%S} {msg}", flush=True)


class JobTimeout(TimeoutError):
    """Async job did not finish within job_timeout_s."""


class JobFailed(RuntimeError):
    """Job ended in ERROR/ABORTED, or the server answered HTTP 500 (usually the 60 s cap)."""


def _is_server_error(err):
    return isinstance(err, JobFailed) or "HTTPError" in type(err).__name__


def _run_job(adql, kw):
    n = net()
    gaia = GaiaClass(show_server_messages=False)     # one client per call (thread-safe)
    job = gaia.launch_job_async(adql, verbose=False, background=True, **kw)
    t0 = time.time()
    while True:
        phase = job.get_phase(update=True).upper().strip()
        if phase not in ("PENDING", "QUEUED", "EXECUTING"):
            break
        if time.time() - t0 > n["job_timeout_s"]:
            try:
                job.abort()
            except Exception:
                pass
            raise JobTimeout(f"job {job.jobid} still {phase} after {n['job_timeout_s']} s")
        time.sleep(n["poll_s"])
    if phase != "COMPLETED":                          # ABORTED jobs can carry partial results
        raise JobFailed(f"job {job.jobid} ended {phase}")
    tbl = job.get_results()
    if len(tbl) >= n["max_rows"]:                     # anonymous output cap: result truncated
        raise JobFailed(f"job {job.jobid} returned {len(tbl)} rows = output cap")
    return tbl


def cached(name, adql):
    out, qfile = RAW / "gaia" / f"{name}.fits", RAW / "gaia" / f"{name}.adql"
    digest = hashlib.sha1(adql.encode()).hexdigest()
    return out.exists() and qfile.exists() and qfile.read_text().splitlines()[0] == f"-- sha1 {digest}"


def query(adql, name, refresh=False, upload=None, upload_name=None, retries=None, retries_500=None):
    """Run `adql` as an async job and cache the result as data/raw/gaia/<name>.fits.

    The ADQL text is stored next to the result; a changed query invalidates the cache.
    Timeouts / dropped connections are retried `retries` times and HTTP 500 /
    failed jobs `retries_500` times, waiting backoff_s[i] before attempt i+1.
    """
    n = net()
    retries = n["retries"] if retries is None else retries
    retries_500 = n["retries_500_max_level"] if retries_500 is None else retries_500
    out, qfile = RAW / "gaia" / f"{name}.fits", RAW / "gaia" / f"{name}.adql"
    out.parent.mkdir(parents=True, exist_ok=True)
    if cached(name, adql) and not refresh:
        return Table.read(out)
    t0 = time.time()
    kw = {}
    if upload is not None:
        kw = dict(upload_resource=str(upload), upload_table_name=upload_name)
    n_net = n_500 = 0
    while True:
        try:
            tbl = _run_job(adql, kw)
            break
        except Exception as err:
            if _is_server_error(err):
                n_500 += 1
                if n_500 > retries_500:
                    raise
            else:                                   # JobTimeout, socket timeout, dropped connection
                n_net += 1
                if n_net > retries:
                    raise
            wait = n["backoff_s"][min(n_net + n_500 - 1, len(n["backoff_s"]) - 1)]
            log(f"  {name}: {type(err).__name__} {str(err)[:80]!r}; retry in {wait} s")
            time.sleep(wait)
    for c in tbl.colnames:                      # astroquery may upper/lower-case names
        if c != c.lower():
            tbl.rename_column(c, c.lower())
    tbl.meta = {"QUERYSEC": round(time.time() - t0, 1)}
    tbl.write(out, overwrite=True)
    qfile.write_text(f"-- sha1 {hashlib.sha1(adql.encode()).hexdigest()}\n{adql}\n")
    return tbl


# source_id // 2**35 is the HEALPix level-12 (nested) index, so contiguous
# source_id ranges are sky pixels; level-L pixel p spans 4**(12-L) level-12 pixels.
def _healpix_range(level, pix):
    span = 2**35 * 4 ** (12 - level)
    return pix * span, (pix + 1) * span - 1


def chunk_query(adql_template, lev, pix):
    lo, hi = _healpix_range(lev, pix)
    return adql_template.format(range=f"g.source_id BETWEEN {lo} AND {hi}")


def chunk_files(adql_template, name, lev, pix, max_level):
    """Cached files covering pixel (lev, pix), or None if any part is missing."""
    f = f"{name}/L{lev}_{pix:05d}"
    if cached(f, chunk_query(adql_template, lev, pix)):
        return [(lev, pix)]
    if lev >= max_level:
        return None
    kids = [chunk_files(adql_template, name, lev + 1, 4 * pix + k, max_level) for k in range(4)]
    return None if any(k is None for k in kids) else [c for k in kids for c in k]


def any_cached(name, lev, pix, max_level):
    """True if pixel (lev, pix) or any of its descendants has a cache file."""
    if (RAW / "gaia" / name / f"L{lev}_{pix:05d}.fits").exists():
        return True
    return lev < max_level and any(any_cached(name, lev + 1, 4 * pix + k, max_level) for k in range(4))


def chunk_status(adql_template, name, level=2, max_level=None):
    """(complete, remaining) lists of level-`level` pixels."""
    max_level = max_level or net()["max_level"]
    done = [p for p in range(12 * 4**level) if chunk_files(adql_template, name, level, p, max_level)]
    return done, sorted(set(range(12 * 4**level)) - set(done))


def query_chunked(adql_template, name, level=1, refresh=False, max_level=None, workers=None):
    """Run `adql_template` (containing `{range}`, filled with a source_id BETWEEN
    clause) once per HEALPix pixel, so each job stays under the anonymous 60 s
    limit. A pixel that keeps failing is split into its 4 children; at
    max_level it is logged and skipped.
    Returns (vstacked Table, failed list of (level, pixel, error)).
    """
    from astropy.table import vstack
    from concurrent.futures import ThreadPoolExecutor

    n = net()
    max_level = max_level or n["max_level"]
    workers = workers or n["workers"]
    failed = []

    def run(lev, pix):
        deepest = lev >= max_level
        q = chunk_query(adql_template, lev, pix)
        if (not deepest and not refresh and not cached(f"{name}/L{lev}_{pix:05d}", q)
                and any(any_cached(name, lev + 1, 4 * pix + k, max_level) for k in range(4))):
            return [t for k in range(4) for t in run(lev + 1, 4 * pix + k)]   # split on an earlier run
        try:
            return [query(q, f"{name}/L{lev}_{pix:05d}", refresh=refresh,
                          retries_500=n["retries_500_max_level"] if deepest else n["retries_500"])]
        except Exception as err:
            if deepest:
                log(f"  {name} L{lev} pixel {pix} FAILED ({type(err).__name__}: {str(err)[:80]}); skipped")
                failed.append((lev, pix, f"{type(err).__name__}: {err}"[:200]))
                return []
            log(f"  {name} L{lev} pixel {pix} failed ({type(err).__name__}); splitting")
            return [t for k in range(4) for t in run(lev + 1, 4 * pix + k)]

    def top(pix):
        out = run(level, pix)
        log(f"{name} L{level} pixel {pix} done ({sum(len(t) for t in out)} rows)")
        return out

    (RAW / "gaia" / name).mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(workers) as ex:
        parts = [t for ts in ex.map(top, range(12 * 4**level)) for t in ts]
    parts = [p for p in parts if len(p)]
    for p in parts:
        p.meta = {}
    return vstack(parts, metadata_conflicts="silent"), failed
