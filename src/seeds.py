"""Monte Carlo seed labels: derivation, registry, and the amendment-1 guard (DECISIONS D19, E36).

Derivation (D19): every random stream is np.random.SeedSequence([20191021, code, tail]), where
code = sum of the label's character codes (`label_code`) and tail is the stream tag of the call
site: the chunk index c for control samples (0..7 for 40,000 iterations in chunks of 5,000),
99 for bootstraps, 2 / 3 for eblm_kinematics C-2 / evolved fraction, 11 / 12 for eblm_width, and
no tail for hu2026_checks / eblm_compare.control_width. numpy's SeedSequence ignores trailing zero
words, so [seed, code] and [seed, code, 0] are the same stream.

REGISTRY lists every label the project uses, with its tails. tests/test_seeds.py checks that no
two labels produce the same stream, and that the amendment-1 labels share no code with any label.

Labels starting with "wasp_a1_" belong to the rerun of amendment 1. `label_code` refuses them unless
the frozen amendment exists at the repository root and matches its recorded SHA-256 (E36, item g).
"""
import hashlib
import subprocess

from .config import ROOT

MASTER = 20191021
CHUNKS = tuple(range(8))            # 40,000 / 5,000
A1_PREFIX = "wasp_a1_"
FROZEN_PREREG_SHA256 = {"PREREG_wasp_kinematics_2026-11-04.txt":
                        "fece9b6608624c8213d1fb5fa96a6223bc07f5328adf820bd133db4a3490ecc4"}
A1_GLOB = "PREREG_wasp_kinematics_amendment1_*.txt"

_POWER_SETS = ["WASP, all distances", "WASP, d < 600 pc", "WASP, d < 600 pc, 0.5-1.5",
               "WASP, d < 600 pc, Dec < 0 (proxy for a South-only sample)"]


def _reg():
    r = {}
    for k in ("hj_table1", "hj", "longp"):                                  # phase1 / pipeline
        r[k] = CHUNKS
    for k in ("eblm_precision", "hj_precision"):                            # eblm_compare.control_width
        r[k] = ((),)
    for k in ["eblm_primary"] + [f"eblm_G{i}" for i in range(1, 12)]:      # eblm_kinematics.primary_block
        r[k] = CHUNKS + (99,)
    for k in ("eblm_c2", "eblm_G9c2"):
        r[k] = (2,)
    for k in ("I_main", "I_bw005", "I_bw02", "I_d05", "I_now"):
        r[k] = (3,)
    for k in ("w_e_zd", "w_e_z", "w_h_z", "w_h_zd"):                        # eblm_width (exploratory)
        r[k] = (11,)
    for k in ("w_sub_z", "w_sub_zd"):
        r[k] = (12,)
    for n in _POWER_SETS:                                                   # hu2026_checks (Stage E)
        r[f"hu_pow_hs19_{n}"] = r[f"hu_pow_hu_{n}"] = ((),)
    for k in ("hu_contrast_eblm", "a1_pow_hs19_all", "a1_pow_hs19_d600", "a1_contrast_eblm"):
        r[k] = ((),)
    # frozen WASP-only run (wasp_kinematics, E31)
    r["wasp_primary"] = CHUNKS + (99,)
    for i in range(1, 6):
        r[f"wasp_D{i}"] = CHUNKS + (99,)
    for k in ("eblm_wD3", "eblm_wD4", "eblm_wD5", "eblm_wD6"):
        r[k] = CHUNKS
    for k in ["wasp_C"] + [f"wasp_D{i}C" for i in range(1, 7)] + ["wasp_Eii"]:
        r[k + "_W"] = r[k + "_E"] = (99,)
    r["wasp_E_W"] = r["wasp_E_E"] = CHUNKS
    r["wasp_Ei"] = r["wasp_Eiii"] = (99,)
    r.update(AMENDMENT1)
    return r


# amendment 1, section 5 (draft revision 3)
AMENDMENT1 = {"wasp_a1_primary": CHUNKS + (99,), "wasp_a1_E_W": CHUNKS, "wasp_a1_Ei": (99,)}
for _i in range(1, 6):
    AMENDMENT1[f"wasp_a1_D{_i}"] = CHUNKS + (99,)
for _k in ["wasp_a1_C"] + [f"wasp_a1_D{_i}C" for _i in range(1, 7)] + ["wasp_a1_Eii"]:
    AMENDMENT1[_k + "_W"] = AMENDMENT1[_k + "_E"] = (99,)

REGISTRY = _reg()
# Two pairs of existing labels share a code (found 2026-11-10, E36). Their streams differ by tail
# (chunk 0 / bootstrap 99 vs eblm_width 11), so no stream is shared. Kept as is: changing them would
# change results already reported.
KNOWN_CODE_SHARES = {frozenset({"eblm_G2", "w_e_zd"}), frozenset({"eblm_G5", "w_h_zd"})}


def stream_entropy(label, tail=()):
    """The SeedSequence entropy of a stream (trailing zeros dropped, as numpy does)."""
    tail = tail if isinstance(tail, tuple) else (tail,)
    e = [MASTER, sum(map(ord, label)), *tail]
    while e and e[-1] == 0:
        e.pop()
    return tuple(e)


def _sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _git_clean_and_tracked(root, paths):
    for p in paths:
        rel = str(p.relative_to(root))
        if subprocess.run(["git", "-C", str(root), "ls-files", "--error-unmatch", rel],
                          capture_output=True).returncode != 0:
            return f"{rel} is not committed"
        if subprocess.run(["git", "-C", str(root), "diff", "--quiet", "HEAD", "--", rel],
                          capture_output=True).returncode != 0:
            return f"{rel} differs from the committed version"
    return None


def require_frozen_amendment(root=ROOT):
    """Raise unless exactly one frozen amendment 1 exists, matches its committed .sha256, and cites the
    frozen WASP pre-registration's hash. Returns the amendment path."""
    found = sorted(root.glob(A1_GLOB))
    if len(found) != 1:
        raise RuntimeError(f"amendment 1 is not frozen: expected one {A1_GLOB} at the repository root, "
                           f"found {len(found)}. The wasp_a1_* rerun refuses to run (DECISIONS E36).")
    txt = found[0]
    side = txt.with_suffix(".sha256")
    if not side.exists():
        raise RuntimeError(f"{txt.name} has no recorded hash ({side.name}); the rerun refuses to run.")
    rec = side.read_text().split()[0].lower()
    if _sha256(txt) != rec:
        raise RuntimeError(f"{txt.name} does not match its recorded SHA-256; the rerun refuses to run.")
    text = txt.read_text()
    if not all(h in text for h in FROZEN_PREREG_SHA256.values()):
        raise RuntimeError(f"{txt.name} does not cite the frozen pre-registration's hash; refusing.")
    for name, h in FROZEN_PREREG_SHA256.items():
        if _sha256(root / name) != h:
            raise RuntimeError(f"{name} no longer matches its frozen hash; refusing.")
    err = _git_clean_and_tracked(root, [txt, side])
    if err:
        raise RuntimeError(f"{err}; the rerun refuses to run until the frozen amendment is committed.")
    return txt


def label_code(label, root=ROOT):
    """Seed code of a label (sum of character codes; D19). wasp_a1_* labels require the frozen amendment."""
    if label.startswith(A1_PREFIX):
        require_frozen_amendment(root)
    return sum(map(ord, label))
