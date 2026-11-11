"""Seed labels (DECISIONS D19, E36): no two registered labels share a stream; amendment-1 labels share no
code with any label; wasp_a1_* streams refuse to start unless the frozen amendment matches its hash."""
import hashlib
import pathlib
import shutil
import subprocess
import tempfile
from itertools import combinations

import numpy as np

from src import seeds
from src.config import ROOT


def test_no_two_registered_labels_produce_the_same_seed():
    owner = {}
    for label, tails in seeds.REGISTRY.items():
        for t in tails:
            e = seeds.stream_entropy(label, t)
            assert e not in owner or owner[e] == label, (label, t, owner.get(e))
            owner[e] = label
    # the generated states, not just the entropy tuples, are distinct
    states = {tuple(np.random.SeedSequence(list(e)).generate_state(4)) for e in owner}
    assert len(states) == len(owner)


def test_trailing_zero_is_the_same_stream():
    a = np.random.default_rng(np.random.SeedSequence([seeds.MASTER, 1234])).random(4)
    b = np.random.default_rng(np.random.SeedSequence([seeds.MASTER, 1234, 0])).random(4)
    assert np.array_equal(a, b)
    assert seeds.stream_entropy("x", 0) == seeds.stream_entropy("x", ())


def test_amendment1_labels_share_no_code_with_any_label():
    codes = {}
    for label in seeds.REGISTRY:
        codes.setdefault(sum(map(ord, label)), set()).add(label)
    for label in seeds.AMENDMENT1:
        assert codes[sum(map(ord, label))] == {label}, (label, codes[sum(map(ord, label))])
    shares = {frozenset(v) for v in codes.values() if len(v) > 1}
    for a, b in (tuple(s) for s in shares):            # only the documented pre-existing pairs
        assert len({a, b}) == 2
    assert shares == seeds.KNOWN_CODE_SHARES, shares


def test_amendment1_label_list_matches_the_draft():
    txt = (ROOT / "prereg_drafts" / "DRAFT_AMENDMENT_wasp_membership_2026-11-09.txt").read_text()
    import re
    in_draft = set(re.findall(r"\bwasp_a1_[A-Za-z0-9_]+", txt))
    assert in_draft == set(seeds.AMENDMENT1), in_draft ^ set(seeds.AMENDMENT1)


def test_checkpointed_labels_are_registered():
    on_disk = {p.name.rsplit("_", 1)[0] for p in (ROOT / "data").rglob("mc_checkpoints/*") if p.is_dir()}
    assert on_disk <= set(seeds.REGISTRY), on_disk - set(seeds.REGISTRY)


def test_label_code_is_the_character_sum():
    for label in seeds.REGISTRY:
        if not label.startswith(seeds.A1_PREFIX):
            assert seeds.label_code(label) == sum(map(ord, label))


def _git(root, *a):
    subprocess.run(["git", "-C", str(root), *a], check=True, capture_output=True,
                   env={"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t", "GIT_COMMITTER_NAME": "t",
                        "GIT_COMMITTER_EMAIL": "t@t", "PATH": "/usr/bin:/bin:/usr/local/bin:/opt/homebrew/bin"})


def _refused(root):
    try:
        seeds.label_code("wasp_a1_primary", root=root)
    except RuntimeError:
        return True
    return False


def test_wasp_a1_streams_refuse_without_a_valid_frozen_amendment():
    assert _refused(ROOT) == (not list(ROOT.glob(seeds.A1_GLOB)))   # today: not frozen -> refused
    with tempfile.TemporaryDirectory() as d:
        remote = pathlib.Path(d) / "remote.git"
        root = pathlib.Path(d) / "work"
        _git(d, "init", "-q", "--bare", str(remote))
        _git(d, "init", "-q", str(root))
        _git(root, "remote", "add", seeds.REMOTE, str(remote))
        name = next(iter(seeds.FROZEN_PREREG_SHA256))
        shutil.copy(ROOT / name, root / name)
        _git(root, "add", name)
        _git(root, "commit", "-qm", "init")
        _git(root, "branch", "-M", "main")
        assert _refused(root)                                             # no amendment
        txt = root / "PREREG_wasp_kinematics_amendment1_2026-12-01.txt"
        txt.write_text("Amends ... sha256 " + seeds.FROZEN_PREREG_SHA256[name] + "\n")
        assert _refused(root)                                             # no recorded hash
        side = txt.with_suffix(".sha256")
        side.write_text("0" * 64 + "  " + txt.name + "\n")
        assert _refused(root)                                             # wrong hash
        side.write_text(hashlib.sha256(txt.read_bytes()).hexdigest() + "  " + txt.name + "\n")
        assert _refused(root)                                             # right hash, not committed
        _git(root, "add", ".")
        _git(root, "commit", "-qm", "freeze")
        assert _refused(root)                                             # committed, remote branch missing
        _git(root, "commit", "-q", "--allow-empty", "-m", "later")
        _git(root, "push", "-q", seeds.REMOTE, "main~1:refs/heads/other")
        assert _refused(root)                                             # pushed elsewhere; branch missing
        _git(root, "push", "-q", seeds.REMOTE, "main~2:refs/heads/" + seeds.REMOTE_BRANCH)
        assert _refused(root)                                             # branch exists, freeze not on it
        try:
            seeds.require_frozen_amendment(root)
        except RuntimeError as e:
            assert "not an ancestor" in str(e), e
        _git(root, "push", "-q", seeds.REMOTE, "main~1:refs/heads/" + seeds.REMOTE_BRANCH)
        assert not _refused(root)                                         # freeze commit on remote branch
        assert seeds.label_code("wasp_a1_primary", root=root) == sum(map(ord, "wasp_a1_primary"))
        _git(root, "remote", "set-url", seeds.REMOTE, str(pathlib.Path(d) / "missing.git"))
        assert _refused(root)                                             # fetch fails -> refused
        _git(root, "remote", "set-url", seeds.REMOTE, str(remote))
        assert not _refused(root)                                         # valid: allowed
        assert seeds.label_code("wasp_a1_primary", root=root) == sum(map(ord, "wasp_a1_primary"))
        txt.write_text(txt.read_text() + "edited after freeze\n")
        assert _refused(root)                                             # edited after freeze
        _git(root, "checkout", "--", txt.name)
        (root / name).write_text("tampered")
        assert _refused(root)                                             # frozen pre-registration changed
        _git(root, "checkout", "--", name)
        (root / "PREREG_wasp_kinematics_amendment1_2026-12-02.txt").write_text("second")
        assert _refused(root)                                             # two candidate files


def test_wasp_a1_refusal_reaches_the_monte_carlo_entry_points():
    if list(ROOT.glob(seeds.A1_GLOB)):
        return
    import pandas as pd
    from src import config, eblm_kinematics as ek, pipeline, wasp_kinematics as wk
    cfg = config.load()
    d = pd.DataFrame({"U": [1.0, 2.0], "V": [0.0, 1.0], "W": [3.0, 1.0]})
    for call in (lambda: wk.boot_draws(d, cfg, 10, "wasp_a1_C_W", wk.total),
                 lambda: ek.bootstrap(d, 0.0, cfg, 10, "wasp_a1_primary", "mean_abs"),
                 lambda: wk.mc_sigma_w(d, d, cfg, "wasp_a1_E_W", 10, ["z"]),
                 lambda: pipeline._run_chunked(None, d, d, cfg, "wasp_a1_D1", 10, "mean_abs")):
        try:
            call()
        except RuntimeError as e:
            assert "amendment 1" in str(e) or "refuses" in str(e)
        else:
            raise AssertionError("wasp_a1_* stream started without a frozen amendment")
