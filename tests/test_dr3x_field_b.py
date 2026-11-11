"""Method B field downloader robustness (DECISIONS E57): retries on dropped connections and DNS failures with
increasing waits, no retry on a bad query, and the completeness check. Offline: the ARI job is a fake."""
import gzip
import http.client
import pathlib
import socket
import tempfile

import requests
from astropy.table import Table
from pyvo.dal.exceptions import DALQueryError, DALServiceError

from src import dr3x_field_b as B


def _dns_error():
    try:
        try:
            raise socket.gaierror(8, "nodename nor servname provided, or not known")
        except socket.gaierror as g:
            raise requests.exceptions.ConnectionError(
                "Failed to resolve 'gaia.ari.uni-heidelberg.de'") from g
    except requests.exceptions.ConnectionError as e:
        return e


def _dropped():
    try:
        try:
            raise http.client.RemoteDisconnected("Remote end closed connection without response")
        except http.client.RemoteDisconnected as r:
            raise requests.exceptions.ConnectionError(("Connection aborted.", r)) from r
    except requests.exceptions.ConnectionError as e:
        return e


def _wrapped(e):                                   # pyvo wraps some network errors
    try:
        raise DALServiceError("could not submit") from e
    except DALServiceError as w:
        return w


class FakeJob:
    def __init__(self, errors, rows=5):
        self.errors, self.rows, self.calls = list(errors), rows, 0

    def __call__(self, adql):
        self.calls += 1
        if self.errors:
            raise self.errors.pop(0)
        return Table({"source_id": list(range(self.rows)), "ra": [1.5] * self.rows})


def _run(errors, rows=5):
    sleeps, logs = [], []
    job = FakeJob(errors, rows)
    with tempfile.TemporaryDirectory() as t:
        root = pathlib.Path(t)
        try:
            files = B.fetch("07", 0, 1000, root=root, log=logs.append, job=job, sleep=sleeps.append)
            n = B.verify(files[0])
            ok = files[0].with_suffix(".ok").exists()
        except Exception as e:                        # noqa: BLE001
            return job, sleeps, logs, e
    return job, sleeps, logs, (n, ok)


def test_classification():
    assert B.is_transient(_dns_error()) and B.is_transient(_dropped())
    assert B.is_transient(_wrapped(_dns_error())) and B.is_transient(_wrapped(_dropped()))
    assert B.is_transient(DALServiceError("503 Server Error: Service Unavailable"))
    assert not B.is_transient(DALQueryError("syntax error near FROM"))
    assert not B.is_transient(B.BadQuery("job ERROR: column not found"))
    assert not B.is_transient(B.Split("ABORTED: time limit"))
    assert not B.is_transient(ValueError("x"))


def test_dns_failure_and_dropped_connection_are_retried_with_increasing_waits():
    job, sleeps, logs, res = _run([_dns_error(), _dropped(), _wrapped(_dns_error())])
    assert res == (5, True), res
    assert job.calls == 4
    assert sleeps == [60, 120, 240], sleeps                         # 1, 2, 4 minutes
    retries = [m for m in logs if "retrying in" in m]
    assert len(retries) == 3 and "try 1/10" in retries[0] and "1 min" in retries[0], retries
    assert any("verified" in m for m in logs)


def test_gives_up_after_ten_tries():
    job, sleeps, logs, res = _run([_dropped() for _ in range(10)])
    assert isinstance(res, requests.exceptions.ConnectionError), res
    assert job.calls == 10 and len(sleeps) == 9
    assert sleeps == [60 * w for w in (1, 2, 4, 8, 16, 30, 30, 30, 30)]
    assert any("giving up" in m for m in logs)


def test_bad_query_is_not_retried():
    for err in (B.BadQuery("job ERROR: unknown column g.foo"), DALQueryError("syntax error")):
        job, sleeps, logs, res = _run([err])
        assert isinstance(res, type(err)) and job.calls == 1 and sleeps == [], (err, res)


def _write(path, text):
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def test_completeness_check():
    with tempfile.TemporaryDirectory() as t:
        d = pathlib.Path(t)
        good = d / "good.csv.gz"
        _write(good, "a,b\n1,2\n3,4\n")
        assert B.verify(good) == 2
        try:
            B.verify(good, max_rows=2)
        except B.Incomplete as e:
            assert "row limit" in str(e)
        else:
            raise AssertionError("row limit not enforced")
        trunc = d / "trunc.csv.gz"
        trunc.write_bytes(good.read_bytes()[:-12])                     # cut mid-stream
        nonl = d / "nonl.csv.gz"
        _write(nonl, "a,b\n1,2\n3,")
        wrong = d / "wrong.csv.gz"
        _write(wrong, "a,b\n1,2\n3\n")
        for bad in (trunc, nonl, wrong):
            try:
                B.verify(bad)
            except B.Incomplete:
                pass
            else:
                raise AssertionError(f"{bad.name} passed")


def test_incomplete_or_interrupted_range_is_rerun_from_scratch():
    with tempfile.TemporaryDirectory() as t:
        root = pathlib.Path(t)
        out, qf = root / "field_09.csv.gz", root / "field_09.adql"
        qf.write_text(B.query(0, 1000))
        out.write_bytes(b"\x1f\x8b\x08\x00garbage")                  # a broken leftover
        out.with_suffix(".tmp").write_bytes(b"partial")               # an interrupted write (field_09.csv.tmp)
        logs, job = [], FakeJob([])
        files = B.fetch("09", 0, 1000, root=root, log=logs.append, job=job, sleep=lambda s: None)
        assert job.calls == 1 and B.verify(files[0]) == 5
        assert not out.with_suffix(".tmp").exists()
        assert any("incomplete" in m for m in logs)
        # a complete range with a matching marker is not re-run
        job2 = FakeJob([])
        B.fetch("09", 0, 1000, root=root, log=logs.append, job=job2, sleep=lambda s: None)
        assert job2.calls == 0


def test_split_still_works():
    class SplitOnce(FakeJob):
        def __call__(self, adql):
            self.calls += 1
            if self.calls == 1:
                raise B.Split("ABORTED: time limit")
            return Table({"source_id": [1, 2], "ra": [1.0, 2.0]})
    with tempfile.TemporaryDirectory() as t:
        job = SplitOnce([])
        files = B.fetch("05", 0, 1000, root=pathlib.Path(t), log=lambda m: None, job=job, sleep=lambda s: None)
        assert [f.name for f in files] == ["field_05a.csv.gz", "field_05b.csv.gz"] and job.calls == 3
