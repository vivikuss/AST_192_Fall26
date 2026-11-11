"""Blocking privacy scan for git hooks (DECISIONS E49). Standard library only; never prints a private value.

    python3 -I scripts/hooks/privacy_scan.py staged                 # pre-commit: lines added in the index
    python3 -I scripts/hooks/privacy_scan.py message <file>         # commit-msg: the commit message
    python3 -I scripts/hooks/privacy_scan.py push <base> <tip>...   # pre-push: every commit in <base>..<tip>

What is scanned: added lines of each diff (new files in full), the file paths touched, and (push mode) every
commit message. Removed lines are not scanned (they were scanned when they were added). Binary files are checked
by path only.

What counts as a match (same rules as the manual scan used since E36):
  - an extended-sample numeric id: a 12-20 digit token equal to any such value in
    <private>/processed/extended_sample_all.csv or extended_sample_pool_k.csv;
  - an extended-sample name-like string (letters + digits, 7-40 chars) occurring anywhere in the text;
  - a CORALIE RV token: a decimal equal to a vrad value of <private>/rv/coralie*/dbase/*_coralie.rdb, as written
    or rounded to 3 or 4 decimals;
  - a private path: data/private/, data/eblm/, *.rdb, or a path containing "coralie" or "extended_sample".

<private> is $AST192_PRIVATE_DIR, default ~/Documents/Astro2026/AST_192_Fall26/07_EBLM_vs_HJ_Tertiaries/data/private.
Fails closed: exit 2 if the reference files are missing or yield nothing, 1 on any match, 0 if clean.
"""
import csv
import glob
import os
import re
import subprocess
import sys
from pathlib import Path

DEFAULT = Path.home() / "Documents/Astro2026/AST_192_Fall26/07_EBLM_vs_HJ_Tertiaries/data/private"
PRIV_PATH = re.compile(r"(^|/)data/private/|(^|/)data/eblm/|\.rdb$|coralie|extended_sample", re.I)
ID_RE, NUM_RE = re.compile(r"\d{12,20}"), re.compile(r"-?\d+\.\d+")


def references():
    d = Path(os.environ.get("AST192_PRIVATE_DIR", DEFAULT)).expanduser()
    ids, names, rv = set(), set(), set()
    csvs = [d / "processed" / "extended_sample_all.csv", d / "processed" / "extended_sample_pool_k.csv"]
    missing = [str(p) for p in csvs if not p.is_file()]
    if missing:
        raise SystemExit(f"privacy scan: reference file(s) missing: {missing}; refusing (fail closed)")
    for p in csvs:
        with open(p, newline="", errors="ignore") as fh:
            rows = csv.reader(fh)
            next(rows, None)                                   # header: column names are not private values
            for row in rows:
                for v in row:
                    v = v.strip()
                    if re.fullmatch(r"\d{12,20}", v):
                        ids.add(v)
                    elif re.search(r"[A-Za-z]", v) and re.search(r"\d", v) and 7 <= len(v) <= 40:
                        names.add(v)
    rdb = glob.glob(str(d / "rv" / "coralie*" / "dbase" / "*_coralie.rdb"))
    for f in rdb:
        with open(f, errors="ignore") as fh:
            for line in fh:
                p = line.split()
                if len(p) >= 2 and re.fullmatch(r"-?\d+\.\d+", p[1]):
                    x = float(p[1])
                    rv.update((p[1], f"{x:.3f}", f"{x:.4f}"))
    if not ids or not rdb or not rv:
        raise SystemExit(f"privacy scan: empty reference set (ids {len(ids)}, rdb files {len(rdb)}, "
                         f"rv tokens {len(rv)}) under {d}; refusing (fail closed)")
    return ids, names, rv


def git(*a):
    r = subprocess.run(["git", *a], capture_output=True, text=True, errors="ignore")
    if r.returncode != 0:
        raise SystemExit(f"privacy scan: git {' '.join(a[:3])} failed: {r.stderr.strip()[:300]}")
    return r.stdout


def diff_items(diff, where):
    """(location, text) for every added line, plus ('path', name) for every file touched."""
    path, ln = None, 0
    for line in diff.splitlines():
        if line.startswith("diff --git "):
            path = line.split(" b/", 1)[-1]
            yield f"{where}{path}", "PATH", path
        elif line.startswith("+++ "):
            if line != "+++ /dev/null":
                path = line[6:] if line.startswith("+++ b/") else line[4:]
        elif line.startswith("@@"):
            m = re.search(r"\+(\d+)", line)
            ln = int(m.group(1)) if m else 0
        elif line.startswith("+"):
            yield f"{where}{path}:{ln}", "TEXT", line[1:]
            ln += 1
        elif line.startswith("Binary files"):
            yield f"{where}{path}", "PATH", path


def scan(items, refs):
    ids, names, rv = refs
    hits = []
    for loc, kind, text in items:
        if kind == "PATH":
            if PRIV_PATH.search(text):
                hits.append((loc, "private path"))
            continue
        if set(ID_RE.findall(text)) & ids:
            hits.append((loc, "extended-sample id"))
        if set(NUM_RE.findall(text)) & rv:
            hits.append((loc, "CORALIE RV token"))
        if any(n in text for n in names):
            hits.append((loc, "extended-sample name"))
    return hits


def items_staged():
    yield from diff_items(git("diff", "--cached", "-U0", "--no-color", "--no-ext-diff", "--no-renames"), "")


def items_message(f):
    for i, line in enumerate(Path(f).read_text(errors="ignore").splitlines(), 1):
        if not line.startswith("#"):
            yield f"commit message:{i}", "TEXT", line


def items_push(base, tips):
    have_base = subprocess.run(["git", "rev-parse", "-q", "--verify", base + "^{commit}"],
                               capture_output=True).returncode == 0
    commits = git("rev-list", *tips, *(["--not", base] if have_base else [])).split()
    if not have_base:
        print(f"privacy scan: {base} not found; scanning the full history of {' '.join(tips)}", file=sys.stderr)
    for c in commits:
        short = c[:7]
        for i, line in enumerate(git("log", "-1", "--format=%B", c).splitlines(), 1):
            yield f"{short} message:{i}", "TEXT", line
        d = git("show", "--format=", "-m", "--first-parent", "-U0", "--no-color", "--no-ext-diff", "--no-renames", c)
        yield from diff_items(d, f"{short} ")
    yield "summary", "COUNT", str(len(commits))


def main(argv):
    if not argv or argv[0] not in ("staged", "message", "push"):
        raise SystemExit(__doc__)
    refs = references()
    if argv[0] == "staged":
        items = list(items_staged())
    elif argv[0] == "message":
        items = list(items_message(argv[1]))
    else:
        if len(argv) < 3:
            raise SystemExit("usage: push <base> <tip>...")
        items = list(items_push(argv[1], argv[2:]))
    n_commits = [t for _, k, t in items if k == "COUNT"]
    items = [i for i in items if i[1] != "COUNT"]
    hits = scan(items, refs)
    scope = f"{argv[0]}" + (f", {n_commits[0]} commit(s)" if n_commits else "")
    if hits:
        print(f"privacy scan ({scope}): {len(hits)} match(es); aborting. Values are not shown.", file=sys.stderr)
        for loc, what in hits:
            print(f"  {loc}: {what}", file=sys.stderr)
        return 1
    print(f"privacy scan ({scope}): clean ({len(items)} items)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
