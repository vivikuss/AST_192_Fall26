"""Cached ADQL async queries against the ARI Heidelberg Gaia TAP service (pyvo).

Same cache layout as gaia_tap (data/raw/ari/<name>.fits + .adql with the query's
sha1). The field query runs as one job per source_id range (HEALPix pixel), each
saved only if it reaches COMPLETED, so a restart resumes from finished ranges (A6, A7).
"""
import hashlib
import signal
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import pyvo
from astropy.table import Table, vstack

from . import config
from .config import RAW
from .gaia_tap import chunk_query, log, net


def settings():
    return config.load()["gaia"]["ari"]


def _paths(name):
    d = RAW / "ari"
    return d / f"{name}.fits", d / f"{name}.adql"


def _key(adql, uploads=None):
    """Cache key: the ADQL plus the content of any uploaded tables."""
    h = hashlib.sha1(adql.encode())
    for k, t in sorted((uploads or {}).items()):
        h.update(k.encode())
        for c in t.colnames:
            h.update(c.encode() + repr(list(t[c])).encode())
    return h.hexdigest()


def cached(name, adql, uploads=None):
    out, qfile = _paths(name)
    digest = _key(adql, uploads)
    return out.exists() and qfile.exists() and qfile.read_text().splitlines()[0] == f"-- sha1 {digest}"


class JobFailed(RuntimeError):
    """Job ended ABORTED / ERROR, or exceeded max_wait_s."""


ACTIVE = {}                       # job_id -> job, for abort on SIGTERM / SIGINT
_LOCK = threading.Lock()


def abort_active(*_):
    """Signal handler: abort every running ARI job of this process, then exit."""
    with _LOCK:
        jobs = list(ACTIVE.items())
    for jid, job in jobs:
        try:
            job.abort()
            log(f"  signal: aborted ARI job {jid}")
        except Exception:
            log(f"  signal: could not abort ARI job {jid}; abort it by hand")
    sys.exit(1)


def install_signal_handlers():
    signal.signal(signal.SIGTERM, abort_active)
    signal.signal(signal.SIGINT, abort_active)


def _run_job(adql, name, uploads=None, max_wait_s=None):
    """Submit, then poll the phase (non-blocking GET) every poll_s; log every status_log_s.

    pyvo's job.wait() long-polls with WAIT=-1 and has no total deadline, which
    is how a 4 h ABORTED job went unnoticed (A7); this loop does not use it.
    """
    a = settings()
    max_wait_s = a["max_wait_s"] if max_wait_s is None else max_wait_s
    svc = pyvo.dal.TAPService(a["tap_url"])
    job = svc.submit_job(adql, maxrec=a["maxrec"], uploads=uploads)
    job.execution_duration = a["execution_duration_s"]
    job.run()
    with _LOCK:
        ACTIVE[job.job_id] = job
    t0 = last_log = time.time()
    log(f"  {name}: ARI job {job.job_id} submitted")
    n, poll_fail = net(), 0
    try:
        while True:
            try:
                phase = job.phase
                poll_fail = 0
            except (ConnectionError, TimeoutError, OSError, pyvo.dal.DALServiceError) as err:
                poll_fail += 1                           # keep the same job; never resubmit here
                if poll_fail > n["retries"]:
                    raise
                wait = n["backoff_s"][min(poll_fail - 1, len(n["backoff_s"]) - 1)]
                log(f"  {name}: poll of job {job.job_id} failed ({type(err).__name__}); retry in {wait} s")
                time.sleep(wait)
                continue
            el = time.time() - t0
            if phase == "COMPLETED":
                break
            if phase in ("ABORTED", "ERROR"):
                raise JobFailed(f"ARI job {job.job_id} {phase} after {el:.0f} s: {job.job.errorsummary}")
            if el > max_wait_s:
                raise JobFailed(f"ARI job {job.job_id} still {phase} after {el:.0f} s (max_wait_s)")
            if time.time() - last_log >= a["status_log_s"]:
                log(f"  {name}: ARI job {job.job_id} {phase}, {el / 60:.1f} min")
                last_log = time.time()
            time.sleep(a["poll_s"])
    except BaseException:                                # incl. KeyboardInterrupt: leave no orphan job
        try:
            job.abort()
            log(f"  {name}: ARI job {job.job_id} aborted")
        except Exception:
            log(f"  {name}: could not abort ARI job {job.job_id}; abort it by hand")
        raise
    finally:
        with _LOCK:
            ACTIVE.pop(job.job_id, None)
    tbl = job.fetch_result(max_retries=net()["retries"]).to_table()   # same job, no resubmit
    try:
        job.delete()
    except Exception:
        pass
    log(f"  {name}: ARI job {job.job_id} COMPLETED in {time.time() - t0:.0f} s, {len(tbl)} rows")
    return tbl


def query(adql, name, refresh=False, retries=None, uploads=None, max_wait_s=None):
    """Run `adql` as one async job; cache the result only if the job COMPLETED.

    Dropped connections are retried with increasing waits; a job that ends
    ABORTED / ERROR is not retried (raises JobFailed). A result with maxrec rows
    is treated as truncated. `max_wait_s` overrides gaia.ari.max_wait_s for this call.
    """
    a, n = settings(), net()
    retries = n["retries"] if retries is None else retries
    out, qfile = _paths(name)
    out.parent.mkdir(parents=True, exist_ok=True)
    if cached(name, adql, uploads) and not refresh:
        return Table.read(out)
    t0 = time.time()
    for attempt in range(retries + 1):
        try:
            tbl = _run_job(adql, name, uploads, max_wait_s)
            break
        except JobFailed:
            raise
        except (ConnectionError, TimeoutError, OSError, pyvo.dal.DALServiceError) as err:
            if attempt == retries:
                raise
            wait = n["backoff_s"][min(attempt, len(n["backoff_s"]) - 1)]
            log(f"  {name}: {type(err).__name__} {str(err)[:80]!r}; retry in {wait} s")
            time.sleep(wait)
    if len(tbl) >= a["maxrec"]:
        raise JobFailed(f"{name}: {len(tbl)} rows = maxrec, result truncated")
    for c in tbl.colnames:
        if c != c.lower():
            tbl.rename_column(c, c.lower())
    for c in tbl.colnames:                      # variable-length strings cannot be written to FITS
        if tbl[c].dtype.kind == "O":
            tbl[c] = tbl[c].astype(str)
    tbl.meta = {"QUERYSEC": round(time.time() - t0, 1), "SERVICE": a["tap_url"]}
    tbl.write(out, overwrite=True)
    qfile.write_text(f"-- sha1 {_key(adql, uploads)}\n{adql}\n")
    return tbl


def range_names(name, level):
    return [(p, f"{name}/L{level}_{p:05d}") for p in range(12 * 4**level)]


def status(adql_template, name, level):
    """(done, todo) range indices at this level."""
    done = [p for p, nm in range_names(name, level) if cached(nm, chunk_query(adql_template, level, p))]
    return done, sorted(set(range(12 * 4**level)) - set(done))


def query_ranges(adql_template, name, level, refresh=False, workers=None):
    """One job per HEALPix level-`level` source_id range. Completed ranges are
    cached and skipped on restart; failed ranges are logged and listed.
    Returns (vstacked Table of completed ranges, failed [(pixel, error)]).
    """
    workers = workers or settings()["workers"]
    failed = []

    def run(item):
        p, nm = item
        try:
            return query(chunk_query(adql_template, level, p), nm, refresh=refresh)
        except Exception as err:
            log(f"  {nm} FAILED ({type(err).__name__}: {str(err)[:120]})")
            failed.append((p, f"{type(err).__name__}: {err}"[:300]))
            return None

    with ThreadPoolExecutor(workers) as ex:
        parts = [t for t in ex.map(run, range_names(name, level)) if t is not None and len(t)]
    for t in parts:
        t.meta = {}
    return (vstack(parts, metadata_conflicts="silent") if parts else Table()), sorted(failed)
