"""Generic target lists (Phase 3 prep): a CSV of Gaia DR2 source_ids through the host pipeline.

CSV columns:
  source_id   Gaia DR2 source_id (required)
  name        label (required)
  is_binary   true/false (required). Binaries get an RV ONLY from `rv` with
              rv_kind == "gamma" (systemic velocity from an orbit fit). A Gaia mean RV
              or a single-epoch RV is never used for them; one without a gamma is
              flagged and has no UVW.
  rv, rv_err  optional per-target systemic-RV override [km/s]
  rv_kind     "gamma" (orbit fit) | "systemic" (other systemic value) | empty
  rv_ref      optional provenance
Non-binaries: the override if given, else the Gaia DR2 RV (as for hosts).
"""
import numpy as np
import pandas as pd

from . import hosts
from .config import ROOT

REQUIRED = ("source_id", "name", "is_binary")


def load(path):
    t = pd.read_csv(path, dtype={"source_id": "Int64"})
    missing = [c for c in REQUIRED if c not in t]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    for c in ("rv", "rv_err"):
        t[c] = pd.to_numeric(t.get(c), errors="coerce") if c in t else np.nan
    for c in ("rv_kind", "rv_ref"):
        t[c] = t[c].fillna("").astype(str) if c in t else ""
    t["is_binary"] = t["is_binary"].astype(str).str.lower().isin(("true", "1", "yes"))
    if t["source_id"].isna().any() or t["source_id"].duplicated().any():
        raise ValueError(f"{path}: source_id missing or duplicated")
    return t


def assign_rv(t):
    """rv / rv_source following the rules in the module docstring."""
    gamma = t["rv"].notna() & (t["rv_kind"] == "gamma")
    other = t["rv"].notna() & ~gamma
    t["rv_used_input"] = np.where(t["is_binary"], np.where(gamma, t["rv"], np.nan),
                                  np.where(t["rv"].notna(), t["rv"], t["radial_velocity"]))
    t["rv_source"] = np.where(t["is_binary"],
                              np.where(gamma, "gamma", np.where(other, "rejected_not_gamma", "none_binary_no_gamma")),
                              np.where(t["rv"].notna(), "override", np.where(t["radial_velocity"].notna(),
                                                                             "gaia_dr2", "none")))
    t["flag_binary_no_gamma"] = t["is_binary"] & ~gamma
    return t


def build(cfg, path, tag):
    """Target table with Gaia DR2 rows + BJ18 distances, same shape as hosts.build()."""
    t = load(path)
    g = hosts.gaia_rows(t["source_id"], cfg, f"targets_{tag}_{cfg['gaia']['release']}")
    g["source_id"] = g["source_id"].astype("Int64")
    t = t.merge(g, on="source_id", how="left", suffixes=("_input", ""))
    t.loc[t["is_binary"], ["radial_velocity", "radial_velocity_error"]] = np.nan   # never used for binaries
    t["in_gaia"] = t["parallax"].notna() | t["phot_g_mean_mag"].notna()
    t["hostname"] = t["name"]
    t = assign_rv(t)
    t["rv"] = t.pop("rv_used_input")
    return t


def default_path(tag):
    return ROOT / "data" / "eblm" / f"targets_{tag}.csv"
