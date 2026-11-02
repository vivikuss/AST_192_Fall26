"""Step 10: Monte Carlo z- and colour-matched control samples (HS19 Sec. 3).

Each iteration:
  (a) draw a subsample of the field whose z distribution matches the hosts'
      (histogram weights w_k = host fraction in bin k / field count in bin k);
  (b) for every host, pick one star from that subsample within color_tol in
      (BP-RP)_0 (without replacement within the iteration).
The control sample's velocity dispersion is computed about its own mean.

Matching variables beyond z (e.g. [Fe/H], Phase 3) enter step (a) through a
multi-dimensional histogram over cfg['monte_carlo']['match_vars'].
"""
import numpy as np

from .kinematics import dispersion


# ------------------------------------------------------------------ binning
def bin_edges(host_vals, mc_var_cfg):
    """Bin edges for one matching variable, spanning the host range."""
    x = np.asarray(host_vals, dtype=float)
    lo, hi = np.nanmin(x), np.nanmax(x)
    if mc_var_cfg["bins"] == "none":             # one bin over everything = no matching
        return np.array([-np.inf, np.inf])
    if mc_var_cfg["bins"] == "fixed":
        w = mc_var_cfg["bin_width"]
        edges = np.arange(np.floor(lo / w) * w, np.ceil(hi / w) * w + 0.5 * w, w)
        if edges[-1] <= hi:                      # hi exactly on an edge
            edges = np.append(edges, edges[-1] + w)
        return edges
    if mc_var_cfg["bins"] == "quantile":
        e = np.quantile(x, np.linspace(0, 1, mc_var_cfg["n_quantile_bins"] + 1))
        e[-1] = np.nextafter(e[-1], np.inf)      # include the max host
        return np.unique(e)
    raise ValueError(mc_var_cfg["bins"])


def var_cfg(mc, var):
    """Binning settings for matching variable `var` (z keeps its legacy keys)."""
    if var == "z":
        zc = mc["z_match"]
        return {"bins": zc["bins"], "bin_width": zc["bin_width_pc"],
                "n_quantile_bins": zc["n_quantile_bins"]}
    return mc["match_var_bins"][var]


def histogram_weights(host_X, field_X, edges_list):
    """Per-field-star weight w = p_host(bin) / n_field(bin) over the joint
    histogram of the matching variables. Field stars outside the host range,
    or in bins with no hosts, get w = 0.

    Returns (weights, info dict). Host bins with no field stars are reported
    (those hosts' z cannot be reproduced; the matched distribution is then
    renormalised over populated bins).
    """
    shape = tuple(len(e) - 1 for e in edges_list)

    def flat_index(X):
        idx = np.zeros(len(X), dtype=np.int64)
        ok = np.ones(len(X), dtype=bool)
        for j, e in enumerate(edges_list):
            k = np.searchsorted(e, X[:, j], side="right") - 1
            ok &= (k >= 0) & (k < len(e) - 1)
            idx = idx * (len(e) - 1) + np.clip(k, 0, len(e) - 2)
        return np.where(ok, idx, -1)

    hk = flat_index(np.asarray(host_X, dtype=float))
    fk = flat_index(np.asarray(field_X, dtype=float))
    nb = int(np.prod(shape))
    h_counts = np.bincount(hk[hk >= 0], minlength=nb).astype(float)
    f_counts = np.bincount(fk[fk >= 0], minlength=nb).astype(float)
    empty = (h_counts > 0) & (f_counts == 0)
    p_host = np.where(f_counts > 0, h_counts, 0.0)
    p_host /= p_host.sum()
    ratio = np.divide(p_host, f_counts, out=np.zeros(nb), where=f_counts > 0)
    w = np.where(fk >= 0, ratio[np.clip(fk, 0, None)], 0.0)
    info = {"n_bins": nb, "host_counts": h_counts.reshape(shape),
            "field_counts": f_counts.reshape(shape), "n_hosts_in_empty_bins": int(h_counts[empty].sum()),
            "n_hosts_outside": int((hk < 0).sum())}
    return w, info


# ------------------------------------------------------------------ Monte Carlo
class ControlSampler:
    """Precomputes everything that does not change between iterations."""

    def __init__(self, host, field, cfg, match_vars=None):
        mc = cfg["monte_carlo"]
        self.tol = mc["color_tol"]
        self.method = mc["z_match"]["method"]
        self.without_replacement = mc["control_without_replacement"]
        self.match_vars = list(match_vars or mc["match_vars"])
        hX = np.column_stack([np.asarray(host[v], dtype=float) for v in self.match_vars])
        fX = np.column_stack([np.asarray(field[v], dtype=float) for v in self.match_vars])
        self.edges = [bin_edges(hX[:, j], var_cfg(mc, v)) for j, v in enumerate(self.match_vars)]
        w, self.bin_info = histogram_weights(hX, fX, self.edges)

        hc = np.asarray(host["bp_rp0"], dtype=float)
        fc = np.asarray(field["bp_rp0"], dtype=float)
        # only field stars that can ever be used: w > 0 and within colour reach of a host
        usable = (w > 0) & (fc >= hc.min() - self.tol) & (fc <= hc.max() + self.tol)
        idx = np.flatnonzero(usable)
        order = np.argsort(fc[idx], kind="stable")
        self.field_idx = idx[order]                         # into the original field table
        self.fc = fc[self.field_idx]
        self.w = w[self.field_idx]
        self.p_keep = self.w / self.w.max()                 # thinning probabilities
        self.fU = np.asarray(field["U"], dtype=float)[self.field_idx]
        self.fV = np.asarray(field["V"], dtype=float)[self.field_idx]
        self.fW = np.asarray(field["W"], dtype=float)[self.field_idx]
        self.hc = hc
        self.hU, self.hV, self.hW = (np.asarray(host[c], dtype=float) for c in ("U", "V", "W"))
        self.n_host = len(hc)
        self.n_usable = len(self.field_idx)
        self.mean_subsample = float(self.p_keep.sum())

    def subsample(self, rng):
        """Step (a): positions (into the colour-sorted usable field) of the z-matched subsample."""
        if self.method == "thinning":
            return np.flatnonzero(rng.random(self.n_usable) < self.p_keep)
        if self.method == "weighted_resample":
            m = int(round(self.mean_subsample))
            pos = rng.choice(self.n_usable, size=m, replace=True, p=self.w / self.w.sum())
            return np.sort(pos)
        raise ValueError(self.method)

    def control(self, rng):
        """Steps (a)+(b): positions of one control sample (-1 = host unmatched)."""
        sub = self.subsample(rng)
        cs = self.fc[sub]
        lo = np.searchsorted(cs, self.hc - self.tol, side="left")
        hi = np.searchsorted(cs, self.hc + self.tol, side="right")
        n = hi - lo
        pick = np.full(self.n_host, -1, dtype=np.int64)
        has = n > 0
        pick[has] = sub[lo[has] + (rng.random(has.sum()) * n[has]).astype(np.int64)]
        if self.without_replacement:
            pick = self._dedupe(pick, sub, lo, hi, rng)
        return pick

    def _dedupe(self, pick, sub, lo, hi, rng):
        """Re-draw hosts whose pick duplicates an earlier host's (random host order)."""
        order = rng.permutation(self.n_host)
        p = pick[order]
        _, first = np.unique(p, return_index=True)
        dup = np.ones(len(p), dtype=bool)
        dup[first] = False
        dup &= p >= 0
        if not dup.any():
            return pick
        used = set(p[~dup & (p >= 0)].tolist())
        for j in np.flatnonzero(dup):
            h = order[j]
            cand = [c for c in sub[lo[h]:hi[h]].tolist() if c not in used]
            if cand:
                c = cand[int(rng.random() * len(cand))]
                pick[h] = c
                used.add(c)
            else:
                pick[h] = -1
        return pick

    def run(self, n_iter, rng, statistic="rms", chunk=500):
        """Per-iteration control dispersions, and the host dispersion over the same
        hosts: if some hosts get no match in an iteration, both the control and the
        host dispersion of that iteration use only the matched hosts (D13b)."""
        stats = ("rms", "mean_abs")
        disp = {k: np.empty(n_iter) for k in stats}
        hfull = {k: float(dispersion(self.hU, self.hV, self.hW, k)) for k in stats}
        hdisp = {k: np.full(n_iter, hfull[k]) for k in stats}
        n_unmatched = np.zeros(n_iter, dtype=np.int32)
        host_unmatched = np.zeros(self.n_host, dtype=np.int64)    # iterations each host went unmatched
        for s in range(0, n_iter, chunk):
            m = min(chunk, n_iter - s)
            picks = np.stack([self.control(rng) for _ in range(m)])
            miss = picks < 0
            n_unmatched[s:s + m] = miss.sum(axis=1)
            host_unmatched += miss.sum(axis=0)
            complete = n_unmatched[s:s + m] == 0
            pc = picks[complete]
            for stat in stats:
                vals = np.empty(m)
                vals[complete] = dispersion(self.fU[pc], self.fV[pc], self.fW[pc], stat)
                for i in np.flatnonzero(~complete):
                    ok = ~miss[i]
                    p = picks[i][ok]
                    vals[i] = dispersion(self.fU[p], self.fV[p], self.fW[p], stat)
                    hdisp[stat][s + i] = dispersion(self.hU[ok], self.hV[ok], self.hW[ok], stat)
                disp[stat][s:s + m] = vals
        return {"dispersion": disp[statistic], "dispersion_all": disp,
                "host_dispersion": hdisp[statistic], "host_dispersion_all": hdisp,
                "host_dispersion_full": hfull, "n_unmatched": n_unmatched,
                "host_unmatched_count": host_unmatched}


def candidates(host, field, tol):
    """Number of field stars within +-tol in colour of each host."""
    fc = np.sort(np.asarray(field["bp_rp0"], dtype=float))
    hc = np.asarray(host["bp_rp0"], dtype=float)
    return np.searchsorted(fc, hc + tol, side="right") - np.searchsorted(fc, hc - tol, side="left")


def p_value(host_disp, control_disp):
    """Fraction of controls with dispersion <= the hosts' (and the raw count).
    `host_disp` may be per-iteration (paired with each control)."""
    k = int(np.sum(np.asarray(control_disp) <= host_disp))
    return k / len(control_disp), k
