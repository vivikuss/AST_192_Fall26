"""Stage D: systemic velocity (gamma) of the EBLM IV primaries from the CORALIE RVs.

    python -m src.gamma_fit             # fit all EBLM IV systems; per-system table to data/private only
    python -m src.gamma_fit --summary   # aggregates only (STOP D1)

Rules (frozen; project 07 PREREG_amendment_rv_gamma_2026-11-02.txt and this project's
PREREG_eblm_kinematics_2026-11-02.txt, A):
- Only the columns jdb, vrad, svrad of the curated <star>_coralie.rdb files are read.
- P, T0 (tab:obs), e, omega (tab:params; e = 0 for an upper limit) and any published drift or
  outer orbit are fixed. RV = gamma + K [cos(nu + omega) + e cos(omega)], weighted linear least squares.
- One global CORALIE07->14 offset (0 km/s); 0.020 km/s added in quadrature to sigma_gamma.
- fit_gamma() returns gamma, sigma_gamma, K, gamma - simple average and a flag, and nothing else.
  No residual, chi^2, per-epoch model value or per-era quantity is computed, saved or printed.
- Per-system results go only to the gitignored eblm.gamma.out_dir (C7).
"""
import re
import sys
import tarfile

import numpy as np
import pandas as pd

from . import config
from .config import ROOT
from .eblm import SRC, TEX, _num, systems

OUT_PUB = ROOT / "data" / "processed" / "eblm"
PERMITTED = ("gamma", "sigma_gamma", "K", "gamma_minus_mean", "flag")


# ---------- orbit geometry (no RVs) ----------
def tperi_from_t0(t0, P, e, omega_deg):
    """Periastron time from the primary-eclipse time: eclipse at nu = pi/2 - omega."""
    w = np.radians(omega_deg)
    nu = np.pi / 2 - w
    E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu / 2))
    M = E - e * np.sin(E)
    return t0 - M * P / (2 * np.pi)


def true_anomaly(t, P, tp, e):
    M = np.mod(2 * np.pi * (t - tp) / P, 2 * np.pi)
    E = M + e * np.sin(M) if e < 0.8 else np.full_like(M, np.pi)
    for _ in range(100):
        dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E = E - dE
        if np.max(np.abs(dE)) < 1e-12:
            break
    return 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))


def rv_shape(t, P, tp, e, omega_deg):
    """cos(nu + omega) + e cos(omega): the RV curve per unit K."""
    w = np.radians(omega_deg)
    return np.cos(true_anomaly(t, P, tp, e) + w) + e * np.cos(w)


def drift(t, lin, quad, cubic, t_ref):
    """Published drift polynomial [km/s]; coefficients in m/s per yr^j, tau in years from t_ref."""
    tau = (t - t_ref) / 365.25
    c = [0.0 if not np.isfinite(x) else x for x in (lin, quad, cubic)]
    return (c[0] * tau + c[1] * tau ** 2 + c[2] * tau ** 3) / 1e3


# ---------- the fit: returns only the permitted outputs ----------
def fit_gamma(t, rv, err, shape, fixed, sig_add, sys_kms, min_rv=3, max_cond=1e12):
    """Weighted linear LSQ of (rv - fixed) = gamma + K * shape.

    sigma_gamma comes from the parameter covariance with fixed weights 1/(err^2 + sig_add^2),
    plus sys_kms in quadrature. It is not rescaled by the scatter. Only PERMITTED keys return;
    intermediate arrays are local and discarded.
    """
    nan = dict(gamma=np.nan, sigma_gamma=np.nan, K=np.nan, gamma_minus_mean=np.nan)
    if len(rv) < min_rv:
        return dict(nan, flag="too_few_rv")
    y = rv - fixed
    w = 1.0 / (err ** 2 + sig_add ** 2)
    A = np.column_stack([np.ones_like(shape), shape])
    N = A.T @ (A * w[:, None])
    if not np.all(np.isfinite(N)) or np.linalg.cond(N) > max_cond:
        return dict(nan, flag="singular")
    cov = np.linalg.inv(N)
    gamma, K = cov @ (A.T @ (w * y))
    return dict(gamma=float(gamma), sigma_gamma=float(np.sqrt(cov[0, 0] + sys_kms ** 2)), K=float(K),
                gamma_minus_mean=float(gamma - y.mean()), flag="ok")


# ---------- inputs ----------
def _err(s):
    """'16.6460(60)' -> 0.0060 (error in the last digits)."""
    m = re.match(r"^\s*\$?([-+]?\d*\.?\d+)\((\d+(?:\.\d+)?)\)", str(s))
    if not m:
        return np.nan
    val, e = m.group(1), m.group(2)
    if "." in e:                                  # '4.080(53)' style vs '23.4(2.2)' style
        return float(e)
    dec = len(val.split(".")[1]) if "." in val else 0
    return float(e) * 10 ** -dec


def triples():
    """Outer orbits from tab:triple (P2, K2, e2, omega2, T_peri2 in BJD)."""
    with tarfile.open(SRC) as t:
        tex = t.extractfile(TEX).read().decode("utf-8", "replace")
    i = tex.index(r"\label{tab:triple}")
    j = tex.index(r"\end{tabular}", i)
    rows = []
    for line in tex[i:j].splitlines():
        if "tertiary orbit" in line:
            c = [x.strip() for x in line.rstrip("\\").split("&")]
            name = c[0].replace("EBLM ", "").replace("tertiary orbit", "").strip()
            rows.append(dict(name=name, P2=_num(c[1]), K2=_num(c[3]), e2=_num(c[4]), omega2=_num(c[5]),
                             tperi2=_num(c[6]) + 2455000.0))
    return pd.DataFrame(rows)


def k_errors():
    """Published K and its error from tab:params (parsed again, keeping the parenthetical error)."""
    with tarfile.open(SRC) as t:
        tex = t.extractfile(TEX).read().decode("utf-8", "replace")
    i = tex.index(r"\label{tab:params}")
    j = tex.index(r"\end{longtable}", i)
    out = {}
    for line in tex[i:j].splitlines():
        line = line.strip()
        if line.startswith("EBLM J"):
            c = [x.strip() for x in line.rstrip("\\").split("&")]
            out[c[0].replace("EBLM ", "")] = _err(c[3])
    return pd.Series(out, name="K_err_pub")


def rv_files(cfg, sysdf):
    """Map curated per-exposure files to EBLM IV systems by file-name coordinates. No RVs read."""
    g = cfg["eblm"]["gamma"]
    rows = []
    for d in g["rv_dirs"]:
        for f in sorted((ROOT / d).resolve().glob("*_coralie.rdb")):
            m = re.match(r"^(\d\d)(\d\d)(\d\d)([SN])(\d\d)(\d\d)(\d\d)_coralie\.rdb$", f.name)
            if not m:
                continue                          # DRS summary files
            h, mi, s, sg, dd, dm, ds = m.groups()
            ra = 15 * (int(h) + int(mi) / 60 + int(s) / 3600)
            dec = (-1 if sg == "S" else 1) * (int(dd) + int(dm) / 60 + int(ds) / 3600)
            rows.append(dict(file=f, ra=ra, dec=dec))
    files = pd.DataFrame(rows)
    out = {}
    for _, s in sysdf.iterrows():
        sep = 3600 * np.hypot((files["ra"] - s["ra_deg"]) * np.cos(np.radians(s["dec_deg"])),
                              files["dec"] - s["dec_deg"])
        hit = files[sep < g["match_radius_arcsec"]]
        out[s["name"]] = list(hit["file"])
    return out


def read_rv(paths):
    """jdb (-> BJD), vrad, svrad from the curated files, concatenated. No other column is read."""
    parts = [pd.read_csv(p, sep="\t", skiprows=[1], usecols=["jdb", "vrad", "svrad"]) for p in paths]
    d = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["jdb", "vrad", "svrad"])
    d = d.apply(pd.to_numeric, errors="coerce").dropna()
    return d["jdb"].to_numpy() + 2400000.0, d["vrad"].to_numpy(), d["svrad"].to_numpy()


# ---------- run ----------
VARIANTS = {                                       # name -> (period/epoch source, drift t_ref)
    "primary": ("obs", None),                      # drift t_ref from config (D2)
    "drift_none": ("obs", "none"),
    "drift_bjd2455000": ("obs", "bjd2455000"),
    "drift_mid": ("obs", "mid"),
    "params_orbit": ("params", None),              # tab:params P and T_peri instead of tab:obs P, T0
}


def run(cfg):
    g = cfg["eblm"]["gamma"]
    s = systems().set_index("name")
    tr = triples().set_index("name")
    files = rv_files(cfg, s.reset_index())
    out = []
    for name, r in s.iterrows():
        e = 0.0 if not np.isfinite(r["e"]) else float(r["e"])   # upper limit or circular model -> 0
        w = float(r["omega"]) if np.isfinite(r["omega"]) else (0.0 if e == 0 else np.nan)
        base = dict(name=name, n_files=len(files[name]), e_used=e)
        t, rv, err = read_rv(files[name])
        rv = rv + np.where(t >= 2456977.5, g["offset_kms"], 0.0)  # global offset (0); CORALIE14 from 2014-11-16
        sig_add = (r["sig_add"] if np.isfinite(r["sig_add"]) else 0.0) / 1e3
        for vname, (src, dref) in VARIANTS.items():
            dref = dref or g["drift_ref"]
            res = dict(base, variant=vname)
            if not np.isfinite(w):
                out.append(dict(res, gamma=np.nan, sigma_gamma=np.nan, K=np.nan, gamma_minus_mean=np.nan,
                                flag="no_omega"))
                continue
            if src == "obs":
                P, tp = r["P_obs"], tperi_from_t0(r["T0_pri_bjd"], r["P_obs"], e, w)
            else:
                P, tp = r["P"], r["T_peri"] + 2455000.0
            shape = rv_shape(t, P, tp, e, w)
            fixed = np.zeros_like(t)
            if len(t):
                t_ref = {"mid": 0.5 * (t.min() + t.max()), "bjd2455000": 2455000.0, "none": np.nan}[dref]
                if dref != "none":
                    fixed = fixed + drift(t, r["drift_lin"], r["drift_quad"], r["drift_cubic"], t_ref)
                if name in tr.index:
                    o = tr.loc[name]
                    fixed = fixed + o["K2"] * rv_shape(t, o["P2"], o["tperi2"], o["e2"], o["omega2"])
            fit = fit_gamma(t, rv, err, shape, fixed, sig_add, g["offset_sys_kms"], g["min_rv"], g["max_cond"])
            assert set(fit) == set(PERMITTED)
            out.append(dict(res, **fit))
            del shape, fixed
        del t, rv, err
    return pd.DataFrame(out)


def ephemeris_check():
    """Published-only check of the omega convention: T_peri from (T0_pri, e, omega) vs tab:params T_peri."""
    s = systems()
    ecc = s[(s["e"] > 0.05) & s["omega"].notna()].copy()
    tp = tperi_from_t0(ecc["T0_pri_bjd"], ecc["P_obs"], ecc["e"], ecc["omega"])
    d = np.mod(tp - (ecc["T_peri"] + 2455000.0) + ecc["P_obs"] / 2, ecc["P_obs"]) - ecc["P_obs"] / 2
    ecc["dphase"] = d / ecc["P_obs"]
    return ecc[["name", "P_obs", "e", "omega", "dphase"]]


def main():
    cfg = config.load()
    od = ROOT / cfg["eblm"]["gamma"]["out_dir"]
    if "--summary" in sys.argv:
        from .gamma_summary import summary
        summary(cfg, od)
        return
    od.mkdir(parents=True, exist_ok=True)
    ec = ephemeris_check()
    print(f"omega convention check (published values only; e > 0.05, n = {len(ec)}): "
          f"|dphase| median {ec['dphase'].abs().median():.4f}, max {ec['dphase'].abs().max():.4f}")
    res = run(cfg)
    res.to_csv(od / "gamma_fits.csv", index=False)
    print(f"fits written for {res['name'].nunique()} systems x {res['variant'].nunique()} variants "
          f"(per-system table in {cfg['eblm']['gamma']['out_dir']}; not printed)")


if __name__ == "__main__":
    main()
