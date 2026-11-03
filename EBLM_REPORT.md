# EBLM primaries vs the field and vs hot-Jupiter hosts: Galactic velocity dispersion (Stage D)

Project 08, 2026-11-02. This report follows the frozen pre-registration `PREREG_eblm_kinematics_2026-11-02.txt`
(sha256 `2d7ba27a…9e672`; section letters below match it). The systemic velocities come from fits permitted by project
07's frozen amendment, `PREREG_amendment_rv_gamma_2026-11-02.txt` (07 commit dd73827).
- **Labels:** anything not in the pre-registration is marked **[not pre-registered]**; the one test added after STOP D1
  is marked **[post-D1]**.
- **Privacy:** only aggregates appear here. Per-system values derived from the CORALIE RVs stay in `data/private/`
  (C7).

## 1. Result in one paragraph

**The 79 EBLM IV primaries are kinematically indistinguishable from colour-, z- and distance-matched field stars:
Δσ = σ_EBLM − σ_field = +0.50 km/s, 95 % CI [−4.23, +4.71] km/s.**
- The CI includes zero. Its lower bound excludes a deficit as large as the HJ hosts' (−6.6 km/s in Phase 1).
- Compared directly with colour- and distance-matched HJ hosts, the EBLM primaries are **hotter** by +7.9 km/s,
  95 % interval [+2.1, +13.7].
- This is the first outcome listed in the pre-registration (H). The cold kinematics of HJ hosts are not shared by
  close M-dwarf companions to similar stars from the same WASP survey. That supports a planet-specific explanation,
  such as HS19's tidal destruction, over a survey-selection origin.
- The evidence is at the 95 % level, and power was about 65 % for an HJ-sized effect. Section 6 gives the caveats.

## 2. Inputs

### Sample (pre-registration A)
- **EBLM IV (Triaud et al. 2017):** 118 → cuts 1–4 → 107 → MS cut → **79** (F 54, G 24, K 1).
- **Gaia DR2 and Phase 1 settings, unchanged:**
  - Bailer-Jones+2018 distances;
  - STILISM v4.1 extinction (Lallement+2018 map) with the CV18 coefficients;
  - the published Pleiades MS polynomial;
  - astropy UVW relative to the LSR;
  - the mean residual speed σ = (1/N) Σ|v − v̄|.
- **Field:** the Phase 1 field (376,791 MS stars with UVW), with the 5 EBLM primaries it contained removed: **376,786**
  (DECISIONS D9).
- **HJ comparison sample:** the Phase 1 primary sample (HS19 Table 1 hosts): **297**, σ = 36.18 km/s.

### Systemic velocities (STOP D1; DECISIONS D1–D8)
- **Model:** RV = γ + K[cos(ν + ω) + e cos ω], fitted by weighted linear least squares to the CORALIE RVs (columns
  `jdb`, `vrad`, `svrad` only).
- **Fixed inputs:** P and T0 (`tab:obs`), e and ω (`tab:params`), and the published drifts and outer orbits.
- **Drift polynomials:** referenced to the midpoint of each system's RV span. EBLM IV does not state its reference
  epoch; this was the user's choice at D1, before any kinematics (D7).
- **Instrument offset:** one global CORALIE07→14 offset of 0 m/s, with 0.020 km/s added in quadrature to every σ_γ.
- **Fit quality:**
  - All 118 fits could be computed.
  - Median σ_γ 0.028 km/s for the 79 (max 2.1); fitted K matches EBLM IV's to a median 0.03 %.
  - One system's fitted K is half the published value. Its file header confirms it is the right star, so it is kept
    (D8) and tested by G11.
- **Cross-checks:**
  - **Gaia DR3 two-body orbits** (65 systems): γ − v_com median −0.08 km/s, robust SD 0.73 km/s. The spread is mostly
    Gaia's precision.
  - **Published γ** (15 systems): median +0.007 km/s, robust SD 0.14 km/s, none off by more than 0.7 km/s.

## 3. Primary result (pre-registration B)

| | Value |
|---|---|
| σ_EBLM (N = 79) | 42.27 km/s |
| Matched-field control median (40,000 samples; z 50 pc bins × 5 EBLM-distance quantiles; colour ±0.025 mag) | 41.77 km/s (SD 2.84; 2–98 % 36.54–48.34) |
| **Δσ̂ = σ_EBLM − median(σ_control)** | **+0.50 km/s** |
| **95 % CI** (percentile, 40,000 EBLM bootstrap resamples) | **[−4.23, +4.71] km/s** |
| Bootstrap bias (BCa threshold 0.2 × half-width = 0.89) | −0.27 km/s, so no BCa needed |
| Pre-registered reading | **consistent with the field** |

No EBLM was ever without a colour match, and none fell in an empty bin, so no star was dropped from either side.
The expected CI half-width was about ±5.6 km/s (D); the achieved half-width is ±4.5 km/s.

## 4. Secondary results (C)

- **C-1, as in HS19:**
  - 22,864 of 40,000 controls (0.572) are as cold as or colder than the EBLMs: a Gaussian-equivalent −0.18σ.
  - For comparison, the HJ hosts gave 0 of 40,000 controls as cold (4.5σ).
- **C-2, EBLM vs HJ hosts:**
  - Each EBLM resample is matched to HJ hosts within ±0.05 mag in colour and in the same EBLM-distance quintile.
  - Δσ_HJ = σ_EBLM − σ_HJ: **median +7.90 km/s, 95 % interval [+2.07, +13.74]**.
  - 2 EBLMs have no HJ host in their cell and are dropped in every iteration.

## 5. Sensitivity analyses (G; 5,000 iterations and 5,000 bootstrap resamples each)

| Test | N | σ_EBLM | Control median | Δσ̂ | 95 % CI | C-1 fraction |
|---|---|---|---|---|---|---|
| Primary (z + distance) | 79 | 42.27 | 41.77 | +0.50 | [−4.23, +4.71] | 0.572 |
| G1 z only (as Phase 1) | 79 | 42.27 | 41.71 | +0.56 | [−4.25, +4.91] | 0.582 |
| G2 z and distance, 5 × 5 quantile bins | 79 | 42.27 | 41.71 | +0.57 | [−4.18, +4.89] | 0.581 |
| G3 \|z\| (50 pc) and distance | 79 | 42.27 | 41.98 | +0.29 | [−4.34, +4.44] | 0.544 |
| G4 \|z\| only | 79 | 42.27 | 41.79 | +0.48 | [−4.27, +4.62] | 0.576 |
| G5 field restricted to Dec < 0 | 79 | 42.27 | 41.77 | +0.50 | [−4.33, +4.76] | 0.573 |
| G6 without the 11 drift or triple systems | 68 | 42.94 | 41.72 | +1.23 | [−3.94, +5.92] | 0.652 |
| G7 without the 2 "active" systems | 77 | 41.83 | 41.55 | +0.29 | [−4.46, +4.53] | 0.543 |
| G8 RMS statistic | 79 | 46.81 | 47.71 | −0.91 | [−5.62, +3.28] | 0.399 |
| G9 F stars only | 54 | 38.12 | 38.97 | −0.85 | [−6.52, +4.04] | 0.389 |
| G10 simple average instead of γ | 79 | 42.34 | 41.89 | +0.44 | [−4.22, +4.50] | 0.573 |
| G11 without J1500-33 **[post-D1]** | 78 | 41.65 | 41.84 | −0.19 | [−4.72, +3.86] | 0.476 |

- **G9, C-2 against F-type HJ hosts (102):** +6.27 km/s, 95 % interval [−0.08, +12.54].
- **Overall:** every variant gives a CI that includes 0, and Δσ̂ stays between −0.9 and +1.2 km/s. The result does not
  depend on how z is matched, on the southern sky, on the drift and triple systems, on the statistic, on the γ method
  or on the one system with a K discrepancy.
- **The weakest point is the F-star subset.** The lower bound (−6.52) barely excludes an HJ-sized deficit, and the C-2
  interval just touches 0.

## 6. Interpretation and caveats (H)

**Outcome.** The pre-registered first case holds:
- the CI includes 0;
- its lower bound (−4.23) is above −6.6;
- C-2 is positive.

The cold kinematics of HJ hosts (Phase 1: −6.6 km/s, 4.5σ) are not reproduced by close binaries with M-dwarf
companions around similar F/G stars from the same transit survey. This argues against a shared survey-selection
origin, for example WASP-South's magnitude and colour windows. It supports a planet-specific cause, consistent with
HS19's interpretation that HJ hosts are younger because HJs are destroyed by tides during the main sequence.

**Caveats (all stated in advance or recorded at D1):**
1. **Power.** About 65 % for an HJ-sized offset (D). The exclusion of −6.6 km/s rests on the 95 % CI, and only
   narrowly for the F-star subset (G9).
2. **Different matching for the two samples.** The EBLMs are matched in z and distance, while the Phase 1 HJ result is
   z-only. G1 (z-only) gives the same EBLM answer. C-2 compares the two directly with colour and distance matching.
3. **The samples differ in composition.** The EBLMs are bluer (68 % F stars, against 34 % (102 of 297) of the HJ hosts) and mostly
   southern (C8, C10). Colour matching and G5/G9 address this, but [Fe/H] could not be matched (E, below).
4. **γ for drift and triple systems** can differ from the centre-of-mass velocity of the whole system by up to a few
   km/s (D7). G6 removes these systems and gives Δσ̂ = +1.23.
5. **Outlier and eclipse handling.** EBLM IV's bisector outlier rejection and removal of in-eclipse points were not
   applied, because the amendment allows only `jdb`/`vrad`/`svrad`. The published-K, Gaia and published-γ
   cross-checks show no resulting bias at the km/s level.
6. **Gaia DR2 throughout,** like Phase 1. A DR3 update is Phase 2, which is on hold.

## 7. [Fe/H] (E): dropped
Spectroscopic [Fe/H] could be assembled for only **27 of the 79** primaries from the cached literature, below the
pre-registered 50 % (40). The analysis is therefore not run (DECISIONS D12). Five systems also have discrepant values
between sources.

## 8. Measurement errors (F)
The fitted median σ_γ (0.028 km/s, including the systematic) is far below the 1 km/s trigger. γ and DR2 astrometric
errors inflate a ~42 km/s dispersion by ≲ 0.01 km/s, so no correction is applied.

## 9. EXPLORATORY (pre-registered, I): fraction of evolved primaries at matched colour

**Method.** "Evolved" means M_G0 more than 1 mag above the Pleiades MS. The samples are the EBLMs after cuts 1–4 (107)
and the archive HJ hosts (discovered ≤ 2019) after cuts 1–4 that fall in the EBLM colour range (398 of 457). The HJ
hosts are reweighted to the EBLM colour distribution, and the CI comes from 10,000 joint bootstrap resamples.

| Variant | f_EBLM | f_HJ (colour-weighted) | Δf | 95 % CI |
|---|---|---|---|---|
| **Primary, 0.1 mag bins** | 0.262 | 0.126 | **+0.136** | **[+0.038, +0.227]** |
| 0.05 mag bins (1 EBLM in a bin with no HJ host dropped) | 0.255 | 0.128 | +0.127 | [+0.026, +0.213] |
| 0.2 mag bins | 0.262 | 0.130 | +0.132 | [+0.039, +0.225] |
| HJ hosts from the current archive (no year limit; 603 in range) | 0.262 | 0.159 | +0.103 | [+0.011, +0.195] |
| Threshold 0.5 mag above the MS | 0.617 | 0.407 | +0.210 | [+0.091, +0.322] |

- **Result:** colour reweighting barely changes the raw HJ fraction (0.126 vs 0.126), and the EBLM excess persists in
  every variant.
- **Possible explanations (pre-registered; none is tested here):**
  - (i) HJs being destroyed as their hosts evolve;
  - (ii) detection selection: shallower transits and harder RV confirmation for planets around evolved stars, while
    M-dwarf eclipses stay deep;
  - (iii) mass and metallicity differences at fixed colour;
  - (iv) the cut acting equally on both samples.
- **No conclusion about tidal survival is drawn from this result alone.**

## 10. Figures (`figures/eblm/`)
- **`fig3_uvw_eblm.png`** (as HS19 Fig. 3): 1, 2 and 3σ velocity ellipses of the EBLM primaries over the field, in
  UV/UW/VW, for F, G and K.
  - Per-star points are omitted because the γ values are private (D9).
  - The K panel has N = 1, so no EBLM ellipse is drawn there.
- **`fig4_mc_eblm.png`** (as HS19 Fig. 4): the control-dispersion histogram for the EBLM primaries (z and distance
  matched), with their value and percentiles. Below it, for reference, the Phase 1 HJ-host panel (z matched). Both
  panels share one horizontal range (30–58 km/s), so the two distributions can be compared directly (redrawn
  2026-11-02 at the collaborator's request).

## 11. EXPLORATORY [not pre-registered]: why the EBLM control distribution is wider
Requested by the collaborator after STOP D2 (`python -m src.eblm_width`; `data/processed/eblm/d2/width_exploratory.json`).
The EBLM control SD is 2.84 km/s against 1.48 km/s for the 297 HJ hosts. Pure 1/√N scaling predicts
1.48 × √(297/79) = **2.87 km/s**.

**(a) Spread of individual control-star velocities.** Residual speeds |v − v̄| of the control stars actually picked,
over 2,000 control samples per setup:

| Control pool | N hosts | Per-star mean | Per-star SD | SD/√N | Control SD measured | Pool median d, \|z\| |
|---|---|---|---|---|---|---|
| EBLM, z + distance (primary) | 79 | 42.1 | 25.0 | 2.81 | 2.84 | 486 pc, 211 pc |
| EBLM, z only | 79 | 42.1 | 24.7 | 2.78 | 2.81 | 549 pc, 176 pc |
| HJ, z only (Phase 1) | 297 | 42.9 | 25.1 | 1.46 | 1.45 | 412 pc, 130 pc |
| HJ, z + distance | 297 | 43.6 | 25.8 | 1.50 | 1.48 | 412 pc, 157 pc |

- The per-star spreads of the two pools are the same to within 0.1 km/s (25.0 vs 25.1).
- In every setup, SD/√N reproduces the measured control SD.

**(b) Random subsamples of 79 HJ hosts.** 200 subsamples, each with its own matched controls and 1,000 iterations:

| Matching | Control SD: median (16–84 %) | Control median |
|---|---|---|
| z only (as Phase 1) | **2.86** (2.80–2.93) km/s | 42.6 km/s |
| z + distance (as the EBLM primary) | **2.95** (2.87–3.05) km/s | 43.2 km/s |

**Conclusion.** The wider EBLM distribution comes from sample size. (b) gives 2.86 km/s for 79 HJ hosts, and (a)
shows the same per-star spread.

**How much comes from distance.** Very little, a few per cent of the width:
- 79 HJ hosts: +0.08 km/s (2.86 → 2.95);
- the EBLMs: about +0.03 to +0.10 km/s (2.81 → 2.84 here; 2.74 for G1 at 5,000 iterations, against 2.84 for the
  primary);
- per star: matching on distance raises the per-star SD by 0.3–0.7 km/s, because it pulls in farther, higher-|z|
  stars.

Of the 1.36 km/s difference between the two widths, essentially all is 1/√N.

**Aside.** For the 200 subsamples of 79 HJ hosts, σ_host − control median has median −6.6 km/s (16–84 %: −8.4 to
−5.0) with z matching, and −7.5 km/s with z + distance. An HJ sample the size of the EBLM sample would therefore
typically show the HJ deficit, whereas the EBLMs give +0.5 km/s. This is exploratory and consistent with C-2.

## 12. Reproduce
```
python -m src.gamma_fit             # γ fits (needs the private CORALIE copies; outputs to data/private/)
python -m src.gamma_fit --summary   # STOP D1 aggregates
python -m src.gamma_crosscheck nss  # Gaia DR3 two-body orbits (ARI)
python -m src.eblm_kinematics       # sections 1-10 (~4 min)
python -m src.eblm_width            # section 11, exploratory (~12 min)
python tests/run_tests.py           # 39 tests
```
Seeds are fixed (config `seed`; one stream per label and chunk). Monte Carlo chunks are checkpointed in
`data/private/eblm_d2/mc_checkpoints/`.
