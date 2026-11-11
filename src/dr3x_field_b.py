"""EXPLORATORY (not pre-registered; DECISIONS E52, E53). Method B field: Gaia DR3 field download for the frozen
DR2 primary pipeline ported to DR3.

    python -m src.dr3x_field_b plan          # print the queries and the cut table (no network)
    python -m src.dr3x_field_b download      # fetch every source_id range (resumable)

Chunking follows Kayla's approach (read as text; DECISIONS E50): 48 equal source_id ranges up to the Gaia
source_id limit 12*4^12*2^35; a range whose job times out, or whose result fills the row cap, is split in half
(names <range>a / <range>b); a split range never keeps a partial file. Results are gzipped CSV, one per leaf,
with the ADQL stored next to it (a changed query invalidates the leaf).

Server-side WHERE (only cuts whose DR3 form is exact, or the looser side of an open choice, so that every
pending choice in notes/dr3_exploratory_plan.md can still be applied locally):
  radial_velocity IS NOT NULL; cut 1 parallax_over_error > 10; cut 3 visibility_periods_used > 8; cut 5
  astrometric_gof_al < 3; cut 6 astrometric_excess_noise_sig < 2; cut 7 rv_nb_transits > 5; cut 8 with the
  LOOSER lower coefficient 0.0015 (0.015 is applied locally, as the frozen DR2 run did: config cuts.server_lo_b);
  cuts 9-10 flux_over_error > 10.
Not on the server: cut 4 (UWE; applied locally, as in DR2), cut 2 (no DR3 column), the 5-parameter restriction,
the Katz et al. (2023) RV correction, and the choice of distance column; their inputs are downloaded.
Distances: LEFT JOIN the Bailer-Jones et al. (2021) EDR3 distances (ARI gaiaedr3_complements.distances =
ESA external.gaiaedr3_distance): r_med_geo, r_lo_geo, r_hi_geo, r_med_photogeo, flag.
Service: ARI Heidelberg (ESA anonymous jobs stop at 60 s; E53).
Robustness (E57): network errors (dropped connections, DNS failures, timeouts, HTTP 429/5xx) are retried
with waits of 1, 2, 4, 8, 16, 30... min (10 tries per range), each retry logged; a bad query stops the run.
A range counts as complete only if its whole gzip reads, its rows are below the row limit and its last
line parses; anything else is deleted and rerun from scratch. Run under `caffeinate -i`.
"""
import gzip
import hashlib
import shutil
import sys
import time
from pathlib import Path

ROOT_DIR = Path.home() / "data" / "dr3x_field_b"          # outside every repository (large)
# ARI Heidelberg, as the frozen DR2 field download (A6): anonymous jobs up to 4 h and 10M rows. ESA aborts
# anonymous jobs after 60 s (smoke test 2026-11-11, E53), too short even for 1/1.6M of the source_id range.
TAP = "https://gaia.ari.uni-heidelberg.de/tap"
DIST_TABLE = "gaiaedr3_complements.distances"       # ESA: external.gaiaedr3_distance
MAXREC = 10_000_000
JOB_TIMEOUT_S = 14_400
MIN_FREE_GB = 4.0
N_TOP = 48
MAX_ID = 12 * 4 ** 12 * 2 ** 35

COLUMNS = ["g.source_id", "g.ra", "g.dec", "g.l", "g.b", "g.parallax", "g.parallax_error", "g.parallax_over_error",
           "g.pmra", "g.pmdec", "g.radial_velocity", "g.radial_velocity_error", "g.rv_nb_transits",
           "g.rv_template_teff", "g.grvs_mag", "g.phot_g_mean_mag", "g.bp_rp", "g.phot_bp_rp_excess_factor",
           "g.phot_bp_mean_flux_over_error", "g.phot_rp_mean_flux_over_error", "g.visibility_periods_used",
           "g.astrometric_chi2_al", "g.astrometric_n_good_obs_al", "g.astrometric_gof_al",
           "g.astrometric_excess_noise_sig", "g.astrometric_params_solved", "g.ruwe",
           "d.r_med_geo", "d.r_lo_geo", "d.r_hi_geo", "d.r_med_photogeo", "d.flag AS bj21_flag"]
WHERE = ["g.radial_velocity IS NOT NULL",
         "g.parallax_over_error > 10",                                   # 1
         "g.visibility_periods_used > 8",                                # 3
         "g.astrometric_gof_al < 3",                                     # 5
         "g.astrometric_excess_noise_sig < 2",                           # 6
         "g.rv_nb_transits > 5",                                         # 7
         "g.phot_bp_rp_excess_factor > 1.0 + 0.0015 * POWER(g.bp_rp, 2)",   # 8, looser lower bound
         "g.phot_bp_rp_excess_factor < 1.3 + 0.06 * POWER(g.bp_rp, 2)",
         "g.phot_bp_mean_flux_over_error > 10",                          # 9
         "g.phot_rp_mean_flux_over_error > 10"]                          # 10


def query(lo, hi):
    return (f"SELECT {', '.join(COLUMNS)} FROM gaiadr3.gaia_source AS g "
            f"LEFT JOIN {DIST_TABLE} AS d ON d.source_id = g.source_id "
            f"WHERE {' AND '.join(WHERE)} AND g.source_id >= {lo} AND g.source_id < {hi}")


def top_ranges():
    e = [MAX_ID * i // N_TOP for i in range(N_TOP + 1)]
    return [(f"{i:02d}", e[i], e[i + 1]) for i in range(N_TOP)]


class Split(Exception):
    """The job timed out or filled the row cap: split the range."""


class BadQuery(RuntimeError):
    """The service rejected the query or the job ended in ERROR for a reason other than time: never retried."""


class Incomplete(RuntimeError):
    """A downloaded file failed the completeness check (truncated / unparseable / at the row limit)."""


def _free_gb(p):
    return shutil.disk_usage(p).free / 1e9


# ----------------------------------------------------------------------------- retries (E57)
RETRY_TRIES = 10
RETRY_WAITS_MIN = (1, 2, 4, 8, 16, 30, 30, 30, 30)        # between the 10 tries; doubling, capped at 30 min


def _chain(e):
    seen = []
    while e is not None and e not in seen:
        seen.append(e)
        e = e.__cause__ or e.__context__
    return seen


def is_transient(e):
    """Dropped connections, DNS failures, timeouts, HTTP 429 / 5xx: worth retrying. Never Split / BadQuery."""
    import http.client
    import socket
    import requests
    import urllib3
    transient = (requests.exceptions.ConnectionError, requests.exceptions.Timeout,
                 requests.exceptions.ChunkedEncodingError, urllib3.exceptions.ProtocolError,
                 urllib3.exceptions.NewConnectionError, urllib3.exceptions.NameResolutionError,
                 http.client.RemoteDisconnected, http.client.IncompleteRead, socket.gaierror, ConnectionError,
                 TimeoutError)
    for x in _chain(e):
        if isinstance(x, (Split, BadQuery)):
            return False
        if isinstance(x, transient):
            return True
        name, msg = type(x).__name__, str(x)
        if name == "DALRateLimitError" or ("DALServiceError" in name and any(
                f"{c} " in msg or f"HTTP {c}" in msg or f"{c} Server Error" in msg for c in (500, 502, 503, 504))):
            return True
    return False


def with_retries(fn, what, tries=RETRY_TRIES, waits_min=RETRY_WAITS_MIN, sleep=time.sleep, log=print):
    """fn() with retries on transient network errors only; every retry is logged. Other errors propagate."""
    for attempt in range(1, tries + 1):
        try:
            return fn()
        except Exception as e:                                       # noqa: BLE001 (classified below)
            if not is_transient(e) or attempt == tries:
                if is_transient(e):
                    log(f"{what}: network error on try {attempt}/{tries}, giving up: {type(e).__name__}: "
                        f"{str(e)[:200]}")
                raise
            wait = waits_min[min(attempt - 1, len(waits_min) - 1)]
            log(f"{what}: network error on try {attempt}/{tries} ({type(e).__name__}: {str(e)[:160]}); "
                f"retrying in {wait} min")
            sleep(60 * wait)


# ----------------------------------------------------------------------------- one job
def job_table(adql):
    """Run one ARI async job; returns an astropy Table. Raises Split (time limit / row cap), BadQuery, or the
    network error itself (for with_retries)."""
    import pyvo
    svc = pyvo.dal.TAPService(TAP)
    try:
        job = svc.submit_job(adql, maxrec=MAXREC)
    except Exception as e:                                           # noqa: BLE001
        if is_transient(e):
            raise
        raise BadQuery(f"submission rejected: {type(e).__name__}: {str(e)[:300]}") from e
    try:
        job.execution_duration = JOB_TIMEOUT_S
        job.run()
        t0 = time.time()
        while True:
            ph = job.phase
            if ph in ("COMPLETED", "ERROR", "ABORTED"):
                break
            if time.time() - t0 > JOB_TIMEOUT_S:
                job.abort()
                raise Split(f"timeout after {JOB_TIMEOUT_S} s")
            time.sleep(15)
        if ph != "COMPLETED":
            msg = str((job.job.errorsummary.message.content if job.job.errorsummary else "") or ph)
            if "time" in msg.lower() or "abort" in msg.lower() or ph == "ABORTED":
                raise Split(f"{ph}: {msg[:200]}")
            raise BadQuery(f"job {ph}: {msg[:300]}")
        tbl = job.fetch_result().to_table()
    finally:
        try:
            job.delete()
        except Exception:                                            # noqa: BLE001
            pass
    if len(tbl) >= MAXREC:
        raise Split(f"{len(tbl)} rows = row cap")
    return tbl


# ----------------------------------------------------------------------------- completeness (E57)
def verify(path, max_rows=MAXREC):
    """Read the whole gzip (a truncated stream raises), count data rows, and parse the last line with the
    header's number of fields. Returns the row count; raises Incomplete."""
    import csv
    import zlib
    try:
        with gzip.open(path, "rt", newline="") as fh:
            header = fh.readline()
            if not header.strip():
                raise Incomplete(f"{path.name}: empty")
            n_fields = len(next(csv.reader([header])))
            rows, last = 0, None
            for line in fh:
                rows += 1
                last = line
    except (EOFError, OSError, zlib.error) as e:
        raise Incomplete(f"{path.name}: unreadable gzip ({type(e).__name__}: {e})") from e
    if rows >= max_rows:
        raise Incomplete(f"{path.name}: {rows} rows = row limit")
    if last is not None:
        if not last.endswith("\n"):
            raise Incomplete(f"{path.name}: last line has no newline")
        if len(next(csv.reader([last]))) != n_fields:
            raise Incomplete(f"{path.name}: last line has the wrong number of fields")
    return rows


def write_verified(tbl, out):
    tmp = out.with_suffix(".tmp")
    with gzip.open(tmp, "wt") as fh:
        tbl.to_pandas().to_csv(fh, index=False)
    n = verify(tmp)
    if n != len(tbl):
        tmp.unlink()
        raise Incomplete(f"{out.name}: wrote {n} rows, expected {len(tbl)}")
    tmp.rename(out)
    _mark(out, n)
    return n


def _mark(out, n):
    st = out.stat()
    out.with_suffix(".ok").write_text(f"rows={n} bytes={st.st_size}\n")


def complete(out, qf, q, log=print):
    """A leaf counts as complete only if its ADQL matches and the file passes verify(). A marker records a
    verified file's size, so unchanged files are not re-read on every restart. A failed file is deleted."""
    if not (out.exists() and qf.exists() and qf.read_text() == q):
        return False
    ok = out.with_suffix(".ok")
    if ok.exists() and ok.read_text().split()[1] == f"bytes={out.stat().st_size}":
        return True
    try:
        n = verify(out)
    except Incomplete as e:
        log(f"{out.name}: incomplete ({e}); deleting and rerunning")
        out.unlink()
        ok.unlink(missing_ok=True)
        return False
    _mark(out, n)
    log(f"{out.name}: verified {n} rows")
    return True


def fetch(name, lo, hi, root=ROOT_DIR, log=print, job=None, sleep=time.sleep):
    """Leaf files for [lo, hi): <root>/field_<name>.csv.gz (+ .adql, .ok). Splits on timeout / row cap,
    retries network errors with increasing waits, stops on a bad query."""
    job = job or job_table
    root.mkdir(parents=True, exist_ok=True)
    out, qf = root / f"field_{name}.csv.gz", root / f"field_{name}.adql"
    if (root / f"field_{name}a.adql").exists():                       # split on an earlier run
        mid = (lo + hi) // 2
        return (fetch(name + "a", lo, mid, root, log, job, sleep)
                + fetch(name + "b", mid, hi, root, log, job, sleep))
    q = query(lo, hi)
    if complete(out, qf, q, log):
        return [out]
    out.with_suffix(".tmp").unlink(missing_ok=True)                   # an interrupted write: start over
    if _free_gb(root) < MIN_FREE_GB:
        raise RuntimeError(f"less than {MIN_FREE_GB} GB free on {root}; stopping before range {name}")
    qf.write_text(q)
    t0 = time.time()
    try:
        tbl = with_retries(lambda: job(q), f"range {name}", sleep=sleep, log=log)
    except Split as e:
        qf.unlink()
        if hi - lo < 2:
            raise
        log(f"{name}: {e}; splitting")
        mid = (lo + hi) // 2
        (root / f"field_{name}a.adql").write_text(query(lo, mid))     # marks the split before fetching
        return (fetch(name + "a", lo, mid, root, log, job, sleep)
                + fetch(name + "b", mid, hi, root, log, job, sleep))
    n = write_verified(tbl, out)
    log(f"{name}: {n} rows in {time.time() - t0:.0f} s (verified); free {_free_gb(root):.1f} GB")
    return [out]


def main(argv):
    if not argv or argv[0] not in ("plan", "download"):
        raise SystemExit(__doc__)
    if argv[0] == "plan":
        for name, lo, hi in top_ranges()[:1]:
            print(query(lo, hi))
        print(f"{N_TOP} top-level ranges; output under {ROOT_DIR}")
        return
    files = []
    for name, lo, hi in top_ranges():
        files += fetch(name, lo, hi, log=lambda m: print(time.strftime("%H:%M:%S"), m, flush=True))
    h = hashlib.sha256("\n".join(sorted(p.name for p in files)).encode()).hexdigest()
    print(f"done: {len(files)} leaf files; name-list sha256 {h}")


if __name__ == "__main__":
    main(sys.argv[1:])
