"""Blocking privacy hooks (DECISIONS E49): a planted private value makes pre-commit, commit-msg and pre-push fail.

Runs the committed hook scripts in a throwaway repository against a FAKE private reference directory (invented
id, name and RV value; no real private data is used or written here).
"""
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

from src.config import ROOT

FAKE_ID, FAKE_NAME, FAKE_RV = "9876543210987654321", "FAKESTAR-00042", "-12.34567"
HOOKS = ROOT / "scripts" / "hooks"


def _env(priv):
    return {"PATH": "/usr/bin:/bin:/usr/local/bin:/opt/homebrew/bin", "HOME": str(priv.parent),
            "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t", "GIT_COMMITTER_NAME": "t",
            "GIT_COMMITTER_EMAIL": "t@t", "GIT_CONFIG_NOSYSTEM": "1",
            "AST192_PRIVATE_DIR": str(priv), "AST192_HOOK_PYTHON": sys.executable}


def _fake_private(d):
    priv = d / "private"
    (priv / "processed").mkdir(parents=True)
    (priv / "processed" / "extended_sample_all.csv").write_text(f"name,gaia_dr3_source_id\n{FAKE_NAME},{FAKE_ID}\n")
    (priv / "processed" / "extended_sample_pool_k.csv").write_text("our_dr3\n1111111111111111111\n")
    (priv / "rv" / "coralie07" / "dbase").mkdir(parents=True)
    (priv / "rv" / "coralie07" / "dbase" / "FAKE_coralie.rdb").write_text(
        f"rjd\tvrad\tsvrad\n---\t----\t-----\n55000.1\t{FAKE_RV}\t0.01\n")
    return priv


def _repo(d, env):
    remote, work = d / "remote.git", d / "work"
    run = lambda *a, cwd=d, check=True: subprocess.run(["git", *a], cwd=cwd, env=env, capture_output=True,
                                                        text=True, check=check)
    run("init", "-q", "--bare", str(remote))
    run("init", "-q", str(work))
    shutil.copytree(HOOKS, work / "scripts" / "hooks")
    r = subprocess.run([str(work / "scripts" / "hooks" / "install.sh")], cwd=work, env=env, capture_output=True,
                       text=True)
    assert r.returncode == 0, r.stderr
    run("remote", "add", "origin", str(remote), cwd=work)
    (work / "README").write_text("clean 0.5 1234\n")
    run("add", ".", cwd=work)
    run("commit", "-qm", "init", cwd=work)
    run("branch", "-M", "main", cwd=work)
    run("push", "-q", "origin", "main:victoria-barriola-08", cwd=work)
    run("fetch", "-q", "origin", cwd=work)
    return work, lambda *a, check=False: run(*a, cwd=work, check=check)


def _commit_fails(git, work, fname, text, msg="add"):
    (work / fname).write_text(text)
    git("add", "-A")
    r = git("commit", "-qm", msg)
    git("reset", "-q", "--hard", "HEAD")
    git("clean", "-qfd")
    return r.returncode != 0


def test_hooks_block_planted_private_values():
    with tempfile.TemporaryDirectory() as t:
        d = pathlib.Path(t)
        env = _env(_fake_private(d))
        work, git = _repo(d, env)
        head0 = git("rev-parse", "HEAD").stdout
        # pre-commit
        assert _commit_fails(git, work, "notes.md", f"| WASP-x | {FAKE_ID} |\n")             # id
        assert _commit_fails(git, work, "notes.md", f"see {FAKE_NAME} here\n")                # name
        assert _commit_fails(git, work, "notes.md", f"gamma = {FAKE_RV} km/s\n")              # rv as written
        assert _commit_fails(git, work, "notes.md", "gamma = -12.346 km/s\n")                 # rv rounded (3 dp)
        assert _commit_fails(git, work, "data_private_rv.rdb", "x\n")                         # private path
        assert git("rev-parse", "HEAD").stdout == head0, "a blocked commit was created"
        # commit-msg
        assert _commit_fails(git, work, "ok.md", "fine\n", msg=f"note {FAKE_ID}")
        # a clean commit passes
        assert not _commit_fails(git, work, "ok.md", "fine 0.5\n"), "clean commit was blocked"
        # pre-push: a planted id that skipped the commit hooks (--no-verify) is still stopped at push
        (work / "leak.md").write_text(f"{FAKE_ID}\n")
        git("add", "-A")
        git("commit", "-q", "--no-verify", "-m", "leak", check=True)
        (work / "later.md").write_text("clean\n")
        git("add", "-A")
        git("commit", "-qm", "later", check=True)                 # the leak is one commit below the tip
        before = git("ls-remote", "origin", "victoria-barriola-08").stdout
        r = git("push", "origin", "main:victoria-barriola-08")
        assert r.returncode != 0 and "extended-sample id" in r.stderr, r.stderr
        assert FAKE_ID not in r.stderr, "the hook printed the private value"
        assert git("ls-remote", "origin", "victoria-barriola-08").stdout == before, "push went through"
        # pre-push: a planted value in a commit message is stopped too
        git("reset", "-q", "--hard", "origin/victoria-barriola-08")
        git("commit", "-q", "--no-verify", "--allow-empty", "-m", f"msg {FAKE_RV}", check=True)
        assert git("push", "origin", "main:victoria-barriola-08").returncode != 0
        # clean history pushes
        git("reset", "-q", "--hard", "origin/victoria-barriola-08")
        (work / "ok2.md").write_text("clean\n")
        git("add", "-A")
        git("commit", "-qm", "ok", check=True)
        r = git("push", "origin", "main:victoria-barriola-08")
        assert r.returncode == 0, r.stderr


def test_hooks_fail_closed_without_reference_files():
    with tempfile.TemporaryDirectory() as t:
        d = pathlib.Path(t)
        env = _env(_fake_private(d))
        work, git = _repo(d, env)
        env["AST192_PRIVATE_DIR"] = str(d / "missing")
        (work / "ok.md").write_text("clean\n")
        subprocess.run(["git", "add", "-A"], cwd=work, env=env, check=True)
        r = subprocess.run(["git", "commit", "-qm", "x"], cwd=work, env=env, capture_output=True, text=True)
        assert r.returncode != 0 and "fail closed" in r.stderr, r.stderr
        r = subprocess.run(["bash", str(work / "scripts/hooks/pre-push")], cwd=work, env=env, input="",
                           capture_output=True, text=True)
        assert r.returncode != 0 and "fail closed" in r.stderr, r.stderr


def test_hook_wrappers_use_strict_mode():
    for h in ("pre-commit", "commit-msg", "pre-push", "install.sh"):
        txt = (HOOKS / h).read_text()
        assert "set -euo pipefail" in txt, h
        assert os.access(HOOKS / h, os.X_OK), h
