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


def _free_gb(p):
    return shutil.disk_usage(p).free / 1e9


def run_job(adql, out):
    import pyvo
    svc = pyvo.dal.TAPService(TAP)
    job = svc.submit_job(adql, maxrec=MAXREC, uploads=None)
    job.execution_duration = JOB_TIMEOUT_S
    job.run()
    t0 = time.time()
    try:
        while True:
            ph = job.phase
            if ph in ("COMPLETED", "ERROR", "ABORTED"):
                break
            if time.time() - t0 > JOB_TIMEOUT_S:
                job.abort()
                raise Split(f"timeout after {JOB_TIMEOUT_S} s")
            time.sleep(15)
        if ph != "COMPLETED":
            msg = (job.job.errorsummary.message.content if job.job.errorsummary else "") or ph
            if "time" in str(msg).lower() or "abort" in str(msg).lower() or ph == "ABORTED":
                raise Split(f"{ph}: {str(msg)[:200]}")
            raise RuntimeError(f"job {ph}: {str(msg)[:300]}")
        tbl = job.fetch_result().to_table()
    finally:
        try:
            job.delete()
        except Exception:
            pass
    if len(tbl) >= MAXREC:
        raise Split(f"{len(tbl)} rows = row cap")
    tmp = out.with_suffix(".tmp")
    with gzip.open(tmp, "wt") as fh:
        tbl.to_pandas().to_csv(fh, index=False)
    tmp.rename(out)
    return len(tbl)


def fetch(name, lo, hi, root=ROOT_DIR, log=print):
    """Leaf files for [lo, hi): <root>/field_<name>.csv.gz (+ .adql). Splits on timeout / row cap."""
    root.mkdir(parents=True, exist_ok=True)
    out, qf = root / f"field_{name}.csv.gz", root / f"field_{name}.adql"
    if (root / f"field_{name}a.adql").exists():                       # split on an earlier run
        mid = (lo + hi) // 2
        return fetch(name + "a", lo, mid, root, log) + fetch(name + "b", mid, hi, root, log)
    q = query(lo, hi)
    if out.exists() and qf.exists() and qf.read_text() == q:
        return [out]
    if _free_gb(root) < MIN_FREE_GB:
        raise RuntimeError(f"less than {MIN_FREE_GB} GB free on {root}; stopping before range {name}")
    qf.write_text(q)
    t0 = time.time()
    try:
        n = run_job(q, out)
    except Split as e:
        qf.unlink()
        if hi - lo < 2:
            raise
        log(f"{name}: {e}; splitting")
        mid = (lo + hi) // 2
        (root / f"field_{name}a.adql").write_text(query(lo, mid))     # marks the split before fetching
        return fetch(name + "a", lo, mid, root, log) + fetch(name + "b", mid, hi, root, log)
    log(f"{name}: {n} rows in {time.time() - t0:.0f} s; free {_free_gb(root):.1f} GB")
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
