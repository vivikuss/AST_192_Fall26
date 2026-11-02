"""Load config.yaml and resolve project paths."""
from pathlib import Path
import copy
import yaml

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"


def load(path=None, **overrides):
    """Return the config dict; `overrides` are dotted keys, e.g. {"cuts.uwe.form": "literal"}."""
    with open(path or ROOT / "config.yaml") as f:
        cfg = yaml.safe_load(f)
    for key, val in overrides.items():
        set_key(cfg, key, val)
    return cfg


def set_key(cfg, dotted, val):
    d = cfg
    *head, last = dotted.split(".")
    for k in head:
        d = d[k]
    d[last] = val
    return cfg


def variant(cfg, **overrides):
    """Deep copy of cfg with dotted-key overrides."""
    new = copy.deepcopy(cfg)
    for key, val in overrides.items():
        set_key(new, key.replace("__", "."), val)
    return new


def out_dir(cfg, kind="processed"):
    """data/processed/<run_name>/ or figures/<run_name>/."""
    base = ROOT / ("figures" if kind == "figures" else "data/processed")
    p = base / cfg["run_name"]
    p.mkdir(parents=True, exist_ok=True)
    return p
