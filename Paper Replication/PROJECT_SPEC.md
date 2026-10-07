# Project Spec: Morton & Winn (2014) Replication with Corrected Inclination Inference

## How to use this file

Save this at the root of a fresh git repository. Start Claude Code in that
directory and open with:

> Read PROJECT_SPEC.md. Work through the stages in order. Do not begin a stage
> until the previous stage's validation targets pass. Stop and report at the end
> of each stage rather than continuing automatically.

Attach the four PDFs (Morton & Winn 2014, Masuda & Winn 2020, Mazeh et al. 2015,
Morgan et al. 2024) to the session or place them in `papers/`.

---

## 1. Scientific goal

Morton & Winn (2014, ApJ 796, 47) found that Kepler stars with multiple
transiting planets have systematically lower obliquities than stars with a
single transiting planet. They inferred a Fisher concentration parameter
κ = 4.8 (+2.0/−1.6) for singles and κ = 19.1 (+73.4/−12.1) for multis, with a
null-hypothesis probability of 1.5% (total variation distance) or 3.6% (squared
Hellinger distance).

Masuda & Winn (2020, AJ 159, 81) later showed that the marginalization Morton &
Winn used to obtain p(cos i | D) is mathematically incorrect. It treats v and
v sin i as statistically independent when they are not. Masuda & Winn state in
their Section 6 that they *suspect* correcting such errors "would not change the
qualitative conclusions of those studies, although it would make quantitative
differences." This has not been tested for Morton & Winn specifically.

This project answers two questions:

1. Does the single-versus-multi obliquity difference survive the corrected
   inference on the original 70-star sample?
2. Does it survive on a sample several times larger, built from catalogs
   published since 2014?

Both outcomes are publishable as a result. A confirmation is a validation; a
reversal is a correction. Do not optimize toward either.

---

## 2. The mathematics (implement exactly as written)

### Notation

| Symbol | Meaning |
|---|---|
| `v` | stellar equatorial rotation velocity, 2πR/P_rot |
| `u` | projected rotation velocity, v sin i |
| `L_v(v)` | likelihood for v, derived from R and P_rot measurements |
| `L_u(u)` | likelihood for u, from spectroscopic line broadening |
| `P_v(v)` | prior on v (uniform) |
| `i` | stellar inclination, restricted to [0°, 90°] |

Both papers place a uniform prior on `cos i` over [0, 1]. This is identical to
an isotropic prior, and identical to a `sin i` prior on `i`. Do not apply it
twice.

### Correct marginal likelihood (Masuda & Winn Eq. 9)

```
p(D | cos i)  ∝  ∫  L_v(v) · L_u( v·√(1 − cos²i) ) · P_v(v)  dv
```

The posterior follows by Bayes' theorem with the uniform `cos i` prior
(their Eq. 10). The essential feature is that **the same `v` appears inside both
likelihood functions**. That is the correlation Morton & Winn dropped.

### Incorrect marginal likelihood (Morton & Winn Eq. 13 = Masuda & Winn Eq. 18)

```
p_incorrect(D | cos i)  ∝  ∫  L_v( u / √(1 − cos²i) ) · L_u(u)  du
```

Missing relative to the correct form: the `1/sin i` Jacobian and the `P_v(v)`
prior weighting.

Implement **both**. They are the central comparison of the project. Name them
unambiguously in code (`marginal_likelihood_correct`, `marginal_likelihood_mw14`)
and never let a default argument silently select one.

### Expected direction of the error

Masuda & Winn Figure 1 shows the incorrect marginalization overestimates
probability density at small `cos i` and underestimates it at large `cos i`,
biasing toward apparent alignment. It does this to **both** subsamples, so the
single-versus-multi *difference* may partially cancel. Whether it does is the
open question. Do not assume.

### Hierarchical layer (Morton & Winn Section 2)

Obliquity θ is modeled as a Fisher distribution:

```
f_θ(θ | κ) = [κ / (2 sinh κ)] · exp(κ cos θ) · sin θ
```

The induced distribution on `cos i` (their Eq. 11):

```
f_cosI(z | κ) = [2κ / (π sinh κ)] · ∫_z^1  cosh(κ√(1−y²))/√(1−y²) · 1/√(1−(z/y)²)  dy
```

Population likelihood via Hogg et al. (2010) importance sampling (their Eq. 4):

```
L_κ  ≈  ∏_n  (1/K) · Σ_k  f_κ(c_nk) / π_0(c_nk)
```

where `π_0` is the prior under which the per-star posterior samples were drawn
(uniform on [0,1], so `π_0 = 1`). Use K = 1000 samples per star.

Prior on κ (Fabrycky & Winn 2009, their Eq. 3):

```
π_κ(κ) ∝ (1 + κ²)^(−3/4)
```

### Numerical hazards (handle these explicitly)

- `f_cosI` has an integrable singularity at `y = z`. Use a substitution
  (e.g. `y = z / cos φ`) or a quadrature routine that tolerates endpoint
  singularities. Do not use naive fixed-grid trapezoid integration.
- `sinh κ` overflows in float64 for κ ≳ 710. Work with `κ/sinh κ` via
  `np.expm1` or in log space.
- `f_cosI(·|κ)` must integrate to 1 on [0,1]. Assert this numerically for a grid
  of κ values as a test.

---

## 3. Stages

### Stage 0 — Environment and unit tests

No science data yet. Verify the inference engine against published figures.

1. Set up the environment. Report Python version and package versions. Verify
   network access by fetching one small table and printing its shape.
2. Reproduce **Masuda & Winn Figure 1**, all four panels:
   - v = 10.0 ± 2.0 km/s, v sin i = 9.8 ± 2.0 km/s
   - v = 10.0 ± 2.0 km/s, v sin i = 4.0 ± 2.0 km/s
   - v = 10.0 ± 2.0 km/s, v sin i < 5.0 km/s (flat below threshold, zero above)
   - v = 10.0 ± 2.0 km/s, v sin i < 30.0 km/s
   Each panel must show the correct analytic curve, an independent MCMC
   sampling of the same posterior, and the incorrect-method curves. The analytic
   and MCMC results must agree.
3. Confirm the bottom-right panel recovers a flat posterior in `cos i`. A
   `v sin i` measurement with no constraining power must return the prior. This
   is the single best test that the implementation is correct.
4. Reproduce **Morton & Winn Figure 2**: V sin I = V = 10 ± 0.5 km/s, showing
   the correct posterior against the "divide-and-clip" Monte Carlo that produces
   a spurious delta function at cos i = 0.
5. Reproduce **Morton & Winn Figure 1**: analytic `f_cosI(z|κ)` against direct
   Monte Carlo, for several κ.

**Validation target:** all five reproductions visually match the published
figures. Do not proceed otherwise.

### Stage 1 — Reproduce the original result

Input: `data/morton_winn_table1.csv`, transcribed from Table 1 of the paper
(70 KOIs: `koi`, `R_star`, `R_star_err_hi`, `R_star_err_lo`, `P_rot`,
`P_rot_err`, `vsini`, `vsini_err`, `n_cands`). If this file is not present,
stop and ask for it rather than reconstructing it from catalogs.

1. Implement the differential-rotation treatment from their Section 4.2:
   spot latitudes ℓ = 20° ± 20°, shear α = 0.23, with
   `P_rot(ℓ) = P_eq / (1 − α sin²ℓ)`. Propagate by Monte Carlo over spot
   latitude and stellar radius to obtain `L_v(v)` per star.
   Note: the V sin I values in Table 1 are already corrected by the factor
   `(1 − α/2)`. Do not apply that correction a second time.
2. Compute per-star `cos i` posteriors using the **incorrect** formula.
3. Run the hierarchical inference for the full sample, the 45 singles, and the
   25 multis.

**Validation targets:**
- κ_single maximum-posterior ≈ 4.8, κ_multi ≈ 19.1
- Reproduce their Figure 5 (stacked `cos i` posteriors) and Figure 6 (κ posteriors)
- Reproduce the 95% confidence upper limits on `i` in the last column of Table 1

If these do not match, the discrepancy is itself a finding, but diagnose it
before moving on. Report which stars drive any difference.

### Stage 2 — The corrected inference, same data

Re-run Stage 1 unchanged except for swapping in the correct marginal likelihood.

Produce:
- Per-star `cos i` posteriors under both methods, overlaid, for all 70 stars
- κ posteriors for singles and multis under both methods, on one figure
- A table of maximum-posterior κ and credible intervals, both methods, both
  subsamples
- The shift in the 95% upper limit on `i` per star, so it is clear which stars
  move and by how much

This is the core result. Treat the figure comparing the four κ posteriors as the
paper's money plot.

### Stage 3 — Extended sample

Build a new sample from post-2014 catalogs. **Do not hard-code VizieR catalog
identifiers.** Search by author and year, print the returned catalog
description and row count, and confirm the column names before use. Report any
catalog you cannot locate rather than substituting one.

Target catalogs:

| Quantity | Sources |
|---|---|
| Planet multiplicity, dispositions | NASA Exoplanet Archive, Kepler DR25 KOI table (`q1_q17_dr25_koi`) and cumulative table |
| Stellar radii | Berger et al. 2020 (AJ 159, 280), Gaia DR2 based |
| Rotation periods | McQuillan et al. 2014 (ApJS 211, 24); Santos et al. 2019 (ApJS 244, 21); Santos et al. 2021 (ApJS 255, 17); Reinhold et al. 2023 |
| v sin i | Petigura et al. 2017 (AJ 154, 107, California-Kepler Survey); Petigura et al. 2022 (AJ 163, 179); Brewer & Fischer 2018 (ApJS 237, 38); Buchhave et al. 2012; APOGEE (Jönsson et al. 2020); LAMOST-Kepler |

Requirements:

1. **Re-derive multiplicity from DR25, do not inherit it from the 2014 paper.**
   Dispositions and planet counts changed. Some of Morton & Winn's singles are
   now multis; some candidates are now false positives. Report explicitly how
   many of the original 70 change category or drop out. This alone is a result
   worth a figure.
2. Apply the same main-sequence cuts the original sample inherited from
   McQuillan et al. (Teff and log g). *Decision point: confirm with the user
   whether to keep these cuts or widen them.*
3. When multiple v s i n i measurements exist for a star, adopt a weighted mean
   following the procedure in Bowler et al. (2023) as described in Morgan et al.
   (2024) Section 2. Record the provenance of every value.
4. **Handle v sin i upper limits properly.** Many cool KOIs have only
   `v sin i < 2 km/s`. These are not missing data and must not be dropped or
   treated as detections. Use the error-function likelihood from Bowler et al.
   (2023) Eq. A17, as applied in Morgan et al. (2024) Eq. 3. The bottom two
   panels of Masuda & Winn Figure 1 are the Stage 0 test for this case.
   Mishandling upper limits will spuriously inflate the misaligned fraction
   among slow rotators, which correlates with host star temperature, which
   correlates with multiplicity. Getting this wrong invalidates the result.
5. Flag any star where `v sin i` exceeds `v` by more than 2σ. This is
   unphysical and usually signals a bad rotation period or a blend.
   Morton & Winn dropped KOI-244 for exactly this reason.

Then run Stages 1 and 2 again on the extended sample, both methods.

**Report:** final sample size, number of singles and multis, and the κ
posteriors. Morton & Winn predicted the null-hypothesis probability would fall
below 1% with a larger sample. State whether it does.

### Stage 4 — Significance and controls

1. Implement the permutation test from their Section 4.3: 1000 random splits of
   the sample into groups matching the observed subsample sizes, computing both
   metrics between the resulting κ posteriors.
   ```
   δ_TVD(p_A, p_B) = max |p_A(κ) − p_B(κ)|
   δ_H²(p_A, p_B)  = 1 − ∫ √(p_A(κ) p_B(κ)) dκ
   ```
   Reproduce their Figure 7 for the original sample, then produce the equivalent
   for the extended sample under the corrected method.
2. Repeat their Section 4.4 planet-radius control (1.3 R⊕ < R_p < 10 R⊕). They
   could not draw conclusions at N=70. With a larger sample, determine whether
   multiplicity or planet radius is the operative variable.
3. Write a caveats section covering, at minimum:
   - The Masuda & Winn Assumption II selection effect: a sample with measured
     rotation periods is biased against low-inclination systems, because
     photometric modulation amplitude goes to zero as the star approaches
     pole-on. This is the same geometric effect Mazeh et al. (2015) exploit, and
     it works against detecting exactly the misaligned systems of interest.
   - Heterogeneity of v sin i measurements across instruments and pipelines.
   - Sensitivity of results to the differential-rotation assumptions
     (α = 0.23, ℓ = 20° ± 20°). Run at least one alternative.

---

## 4. Engineering requirements

- **Structure:** `src/` for importable modules, `notebooks/` for one notebook
  per stage, `data/raw/` for untouched downloads, `data/processed/` for derived
  tables, `figures/` for output. Notebooks orchestrate and plot; they do not
  define the physics.
- **Caching:** every catalog download writes to `data/raw/` and is loaded from
  disk on subsequent runs. Never re-query in a loop.
- **Reproducibility:** fix random seeds, pin package versions in
  `requirements.txt`, log sample sizes at every filtering step so attrition is
  auditable.
- **Testing:** `pytest` suite covering the Stage 0 validations, normalization of
  `f_cosI`, and recovery of the flat prior under an uninformative `v sin i`.
- **Provenance:** every derived quantity carries a column recording its source
  catalog.
- **Honesty:** if a validation target is missed, say so in the notebook output
  and stop. Do not tune parameters to hit a published number.

---

## 5. Things that will go wrong

- Confusing `i` (stellar inclination to line of sight) with `θ` (obliquity, the
  star-planet angle). They are related through Morton & Winn Eq. 5,
  `cos i = sin θ cos φ`, not equal.
- Double-applying the `(1 − α/2)` v sin i correction in Stage 1.
- Applying the isotropic prior twice: once as `sin i` on `i`, again as uniform
  on `cos i`.
- Using importance-sampling weights whose `π_0` does not match the prior the
  per-star posteriors were actually drawn under.
- Dropping upper-limit `v sin i` stars, which removes slow rotators
  non-randomly.
- Trusting a VizieR catalog ID from memory.
