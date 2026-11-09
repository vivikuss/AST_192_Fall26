# Cuts and analysis choices, by analysis

Written 2026-11-09 to compare with a teammate's implementation line by line. **Read-only:** no
code was changed, and no velocity or dispersion was computed for this note. File:line references
are to commit `55aa85d`, which had no uncommitted changes in `src/` or `config.yaml`. Counts come
only from what is already logged: `logs/counts.log` /
`data/processed/phase1_dr2/count_table.csv` (= REPRODUCTION_REPORT.md §2), `logs/eblm_d2.log`,
`logs/wasp_kinematics.log`, EBLM_REPORT.md, DECISIONS.md and the two frozen pre-registrations.
"–" means no logged count.

Codes: **HS19 n** = cut n in the HS19 appendix. **D/A/C/E + number** = DECISIONS.md entry.
**Hu+** = Hu, Zhu, Huang & Zhang 2026 (arXiv:2608.12790). "Field" = Gaia DR2 comparison stars.
All thresholds live in `config.yaml`. The code reads them from there, so both locations are given.

Contents: [0. Shared choices](#0-shared-easy-to-miss-choices-all-three-analyses) ·
[1. Phase 1 HS19 reproduction](#1-phase-1-hs19-reproduction-gaia-dr2) ·
[2. EBLM](#2-eblm-analysis-frozen-prereg_eblm_kinematics_2026-11-02) ·
[3. WASP-only](#3-wasp-only-analysis-frozen-prereg_wasp_kinematics_2026-11-04) ·
[4. Most likely to differ](#4-choices-most-likely-to-differ-from-another-implementation)

---

## 0. Shared, easy-to-miss choices (all three analyses)

The EBLM and WASP analyses reuse the Phase 1 pipeline. These rows apply to all three unless a
later table says otherwise.

| # | Name | Exact definition as implemented (file:line) | Applied to | Source | Differs from paper? |
|---|---|---|---|---|---|
| S1 | Gaia release / table | Gaia **DR2** `gaiadr2.gaia_source`, queried on **ARI Heidelberg** TAP (`config.yaml:10,15-19`; field `src/field.py:15-22`; hosts `src/hosts.py:154-159`). Pleiades and the DR3-RV variant query ESA instead (`src/pleiades.py:48-57`, `src/hosts.py:200-215`). | both | HS19 (DR2); ARI is ours (A6, A8) | No. ARI was checked bit-identical to ESA on 4 field regions and on BJ18 distances for 10 hosts (A6, `logs/ari_validation.json`). |
| S2 | Columns | `source_id, ra, dec, l, b, parallax(_error), pmra(_error), pmdec(_error), radial_velocity(_error), phot_g/bp/rp_mean_mag` + the cut inputs + BJ18 `r_est, r_lo, r_hi, result_flag` (`src/gaia_columns.py:4-16`). The Gaia DR2 `l, b` columns are used for extinction. | both | ours | – |
| S3 | Distance | **Bailer-Jones et al. (2018) `r_est`** (posterior mode), **not 1/parallax**. ARI `gaiadr2_complements.geometric_distance` (= ESA `external.gaiadr2_geometric_distance`), LEFT JOIN on source_id (`config.yaml:13-14,20`; `src/field.py:19-22`; `src/hosts.py:156-158`). Copied to `dist_pc` at `src/pipeline.py:42`, then used for M_G, extinction, z and UVW. A star with no `r_est` gets no M_G/UVW and drops out. **Exception:** the Pleiades check uses 1/plx with the cluster's optimised parallaxes, `M_G0 = G + 5 log10(plx/100) − R_G·0.045` (`src/pleiades.py:71`). | both | D16 (ours). DECISIONS does not record whether HS19 name BJ18 `r_est`. | Not established. D8c: 1/parallax would give 319 Table 1 hosts on the MS instead of 322. |
| S4 | Extinction map | **STILISM cube v4.1 (2018-03-19)**, `data/raw/stilism/stilism_cube.h5` (sha256 `a768d80a…c559d906`), from the Wayback Machine. It is the **Lallement et al. (2018)** map, which builds on Capitanio+2017. The config label `stilism_capitanio2017` is only a name (`config.yaml:93-94`). | both | HS19 cite Capitanio+2017; A1 | **Yes, in version:** HS19 cite Capitanio+2017. The cube STILISM served while HS19 worked is the Lallement+2018 v4.1 map (A1, checked 2026-11-02). |
| S5 | Line-of-sight integration | E(B−V) = Σ dE/ds·h, midpoint rule with n = ceil(d / 1 pc) equal steps (h = d/n), trilinear interpolation (`map_coordinates order=1`). **Voxels outside the cube (±2 kpc in x,y; ±300 pc in z) contribute 0** (`src/extinction.py:27-55`; `config.yaml:95-96`). | both | D9 (ours; paper: "integrate the interpolated grid") | Paper gives no scheme |
| S6 | Extinction coefficients | A_G = 2.740·E(B−V); E(BP−RP) = (3.374 − 2.035)·E(B−V) = 1.339·E(B−V) (`src/extinction.py:58-62`; `config.yaml:98-100`). Casagrande & VandenBerg 2018 Table 2 **mean** coefficients. No Teff dependence. | both | HS19 ("mean extinction coefficients"); D9b | No |
| S7 | Dereddening | `bp_rp0 = bp_rp − E_BPRP`; `M_G0 = phot_g_mean_mag − 5 log10(dist_pc/10) − A_G` (`src/extinction.py:71-72`) | both | ours | – |
| S8 | Solar motion / LSR | astropy `GalacticLSR(v_bary = (11.1, 12.24, 7.25) km/s)` (Schönrich+2010) (`src/kinematics.py:16-24`; `config.yaml:119`) | both | D15 | Paper does not specify. The value does not affect a dispersion about the sample mean. It **does** affect the D5 halo cut (v_tot) in the WASP analysis. |
| S9 | UVW conversion | ICRS SkyCoord(ra, dec, `dist_pc`, pmra, pmdec, RV) → GalacticLSR. U is positive **toward the Galactic centre**, V toward rotation, W toward the NGP (`src/kinematics.py:8-24,42-47`). **No uncertainty propagation.** | both | D15 | Paper propagates none either |
| S10 | Height z | astropy `Galactocentric(galcen_distance = 8.3 kpc, z_sun = 27 pc)` z, with the other frame parameters at astropy defaults (pre-v4.0 values) (`src/kinematics.py:32-39`; `config.yaml:120-124`). Stars with a distance but no RV still get z (`src/pipeline.py:53-59`). | both | D11 (ours; paper silent) | Paper silent. Heliocentric d·sin b is a variant. |
| S11 | Main-sequence cut | Keep `M_G0 > P8(bp_rp0) − 1.0`, with P8 = HS19's **published polynomial** in increasing order: `[1.07857, 6.23258, −10.85944, 21.65561, −20.27879, 9.78665, −2.52543, 0.33004, −0.0170556]`. Strict `>`. **No lower limit** (stars below the MS are kept) and **no colour limit** (`src/main_sequence.py:11-24`; `config.yaml:105-108`) | both | HS19 polynomial; D8, D8b | **Possibly:** HS19's Fig. 1 caption says the cut used the *spline*. Our refit polynomial matches the paper's field count better (382,499 vs 376,791; D10, REPRODUCTION_REPORT §4.2). |
| S12 | Velocity measure (eq. 1) | **Mean residual speed** `σ = (1/N) Σ sqrt[(U−Ū)² + (V−V̄)² + (W−W̄)²]`, each sample about its **own** mean (`src/kinematics.py:50-66`, `mean_abs`; `config.yaml:138`). RMS `sqrt((1/N) Σ[…])` is computed alongside. | both | HS19 eq. 1 as typeset (root inside the sum); D14; Hu+ fn. 1 confirm it is not a typo | Matches the paper (36.18 vs 36.11 km/s). It differs from the project spec's RMS. |
| S13 | Missing values | Any NaN fails every cut (`src/cuts.py:103-113`) | both | ours | – |
| S14 | RNG | numpy `default_rng` (PCG64). Master `seed: 20191021` (`config.yaml:6`). Monte Carlo chunk c of 5,000 uses `SeedSequence([seed, sum(ord(c) for c in label), c])` (`src/pipeline.py:145,161`; `config.yaml:137`). Bootstraps use `SeedSequence([seed, label code, 99])`, C-2 uses stream `2`, I uses stream `3`. | – | D19 (ours) | Not in paper. Note: the label code is a character sum, so two labels that are anagrams would share a stream (none in use do). |
| S15 | Percentiles | `np.percentile` with default linear interpolation, everywhere | – | ours | – |

---

## 1. Phase 1 HS19 reproduction (Gaia DR2)

Driver: `python -m src.phase1` → `src/pipeline.py:run` (`src/pipeline.py:226-259`). Run `phase1_dr2`.
Counts, in the column order **field / Table 1 HJ / archive HJ (≤ 2019) / long-period**, are cumulative.

### 1a. Samples and cuts

| # | Name | Exact definition as implemented (file:line) | Applied to | Source | Differs from paper? (how / why) | Remaining (field / T1 / arch HJ / LP) |
|---|---|---|---|---|---|---|
| P1 | HJ host list (primary) | HS19 **machine-readable Table 1** (338 rows, ArchiveName + SimbadName), `data/hs19_table1.txt` (`config.yaml:61`; `src/hosts.py:182-197`). Archive host columns are joined on ArchiveName; `hostname` is set to **SimbadName** (`src/hosts.py:193`). | hosts | HS19 Table 1; A3 | No | – / 338 / – / – |
| P2 | Archive HJ list (run alongside) | NASA Exoplanet Archive `pscomppars`, snapshot **2026-10-28** (`src/hosts.py:26-38`). `pl_orbper < 10` and `pl_bmassj > 0.1` and `pl_bmassprov ∈ {Mass, Msini, Msin(i)/sin(i)}` (drops "M-R relationship"), with **upper-limit masses kept**, and `disc_year ≤ 2019` (`src/hosts.py:22,41-60`; `config.yaml:43-49`). One row per hostname. | hosts | HS19 (P < 10 d, Mp sin i > 0.1 MJup); D4, D5 | Today's archive, not 2019's. The `pl_bmassprov` filter undoes today's M–R-filled masses (D4). | – / – / 505 / – |
| P3 | Long-period host list | Same as P2 with `pl_orbper > 10`. **Hosts that also have an HJ are kept** (`exclude_hj_hosts_from_longp: false`) (`src/hosts.py:54-55,66-68`; `config.yaml:50-53`) | hosts | HS19; D7 | Paper silent on overlap. 264 vs 367 (D7). | – / – / – / 575 |
| P4 | Gaia DR2 ID | SIMBAD (Harvard mirror, sim-script) by name (SimbadName for Table 1). A record is accepted only if the queried name is in its ID list; fallback is the archive `gaia_dr2_id` (`src/hosts.py:85-146`; `config.yaml:54-55`). **HATS-67:** the paper's SimbadName is the wrong star; that star is kept in the Table 1 run (D6b). | hosts | HS19 (SIMBAD); D6b | No (the wrong HATS-67 star is reproduced) | – / 338 / 500 / 562 |
| P5 | Field download (server prefilter) | One ARI job per HEALPix level-2 source_id range (192). WHERE: `radial_velocity IS NOT NULL AND parallax IS NOT NULL AND` cuts 1–3, 5–7, 9, 10 exactly `AND` cut 8 with the looser `0.0015` (`src/cuts.py:126-149`; `src/field.py:15-53`; `config.yaml:23,85`). Count checked by COUNT(*) (A7). | field | HS19 (field = DR2 with RV + 5-par.); A6, A7 | No | 1,372,060 / – / – / – |
| P6 | 5-parameter solution | `parallax IS NOT NULL` (`src/cuts.py:138`) | field | HS19; D3b | Implementation of "5-parameter" | (in P5) |
| P7 | Field RV | Gaia DR2 `radial_velocity IS NOT NULL` (`src/cuts.py:137`), used directly as RV (`src/pipeline.py:87`) | field | HS19 | No | (in P5) |
| P8 | Cut 1: parallax quality | `parallax_over_error > 10` (`src/cuts.py:48-49`; `config.yaml:65`) | both | HS19 1 | No | 1,372,060 / 338 / 487 / 556 |
| P9 | Cut 2: parallax factor | `mean_varpi_factor_al > −0.23 AND < 0.36` (`src/cuts.py:52-55`; `config.yaml:66`) | both | HS19 2; D3 | Paper names the column `mean_varpi_factor`; DR2 column is `…_al` | 1,372,060 / 338 / 485 / 551 |
| P10 | Cut 3: visibility periods | `visibility_periods_used > 8` (`src/cuts.py:58-59`; `config.yaml:67`) | both | HS19 3 | No | 1,372,060 / 338 / 458 / 512 |
| P11 | Cut 4: unit weight error | `u = sqrt(astrometric_chi2_al / (astrometric_n_good_obs_al − 5))`; keep `u < 1.2 · max(1, exp(−0.2 · (phot_g_mean_mag − 19.5)))` (`src/cuts.py:27-45,62-65`; `config.yaml:68-72`, `form: lindegren`). Applied **locally** (field) | both | HS19 4 = Lindegren+2018 C.1; D1 | **Yes, typo fix:** paper prints `exp(−0.2*G − 19.5)` (no parentheses). Literally that is ≈ `u < 1.2` and would remove 173/338 Table 1 hosts (D1). Hu+ fn. 7 agree. | 1,372,060 / 338 / 457 / 510 (no local field loss) |
| P12 | Cut 5: astrometric GOF | `astrometric_gof_al < 3` (`src/cuts.py:68-69`; `config.yaml:73`) | field | HS19 5 | No | 1,372,060 |
| P13 | Cut 6: excess noise | `astrometric_excess_noise_sig < 2` (`src/cuts.py:72-73`; `config.yaml:74`) | field | HS19 6 | No (paper skips it for hosts, as we do) | 1,372,060 |
| P14 | Cut 7: RV transits | `rv_nb_transits > 5` (`src/cuts.py:76-77`; `config.yaml:75`) | field | HS19 7 | No (skipped for hosts, as paper) | 1,372,060 |
| P15 | Cut 8: BP/RP excess | `1.0 + 0.015·bp_rp² < phot_bp_rp_excess_factor < 1.3 + 0.06·bp_rp²` (`src/cuts.py:80-89`; `config.yaml:76-80`). The server uses `0.0015`, then 0.015 is applied locally. | field | HS19 8 = Lindegren C.2; D2 | **Yes, typo fix:** paper prints `0.0015`. No effect (min excess 1.15; REPRODUCTION_REPORT §2 fn. 1) | 1,372,060 |
| P16 | Cuts 9, 10: BP, RP flux S/N | `phot_bp_mean_flux_over_error > 10`; `phot_rp_mean_flux_over_error > 10` (`src/cuts.py:92-97`; `config.yaml:81-82`) | field | HS19 9, 10 | No | 1,372,060 |
| P17 | Host cut set | Hosts get **cuts 1–4 only** (`config.yaml:88`; applied `src/pipeline.py:118-121`) | hosts | HS19 ("cuts 1-4"); D3c | No (5, 8–10 are not mentioned for hosts) | (= P11) |
| P18 | Exclude hosts from field | Remove every source_id in the three host tables after cuts 1–4 (`src/pipeline.py:233-236`, applied at `:95`). Equivalent to "all listed hosts", because a host that fails 1–4 also fails the field cuts. | field | HS19 ("field stars without hot Jupiters"); D17 | Paper does not say whether long-period hosts are also removed; we remove them | 1,371,976 |
| P19 | Distance available | BJ18 `r_est` not null (S3) | both | D16 | see S3 | 1,371,976 / 338 / 457 / 510 |
| P20 | Extinction / dereddening | S4–S7 | both | see S4–S7 | see S4 | 1,371,976 / 338 / 457 / 510 |
| P21 | Main-sequence cut | S11 (`src/pipeline.py:45,98,127`) | both | HS19; D8 | Polynomial vs spline (S11). 16 Table 1 hosts fail (13 marginal; CoRoT-17/-23/-28 have high E(B−V)). No rescue (D8c). | **376,791 / 322 / 402 / 322** |
| P22 | Host RV | `rv = st_radv` if `st_radv` not null **and** `st_radvlim == 0` **and** the reference does not match "Gaia DR3"/"Gaia EDR3"; else Gaia DR2 `radial_velocity`; else none (`src/hosts.py:170-179`; `config.yaml:56-60`). **No DR3 fill** in the primary (`rv_fill_dr3: false`) | hosts | HS19 (archive systemic RV, else DR2); D6, D6c | **Yes, in effect:** today's archive (some values from post-2019 papers) replaces the unrecoverable 2019 archive. 25 Table 1 hosts have no RV today (D6c). | – / **297** / 360 / 264 |
| P23 | Final = UVW computable | `on_ms & U.notna()` (`src/pipeline.py:98,127`) | both | – | – | **376,791 / 297 / 360 / 264** (paper 385,036 / 338 / – / 367) |
| P24 | Spectral bins | F `0.438 < c < 0.742`, G `0.742 < c < 1.00`, K `1.00 < c < 2.02`, c = `bp_rp0`, **strict** (`src/pipeline.py:64-69`; `config.yaml:126-129`). Reporting only, not a cut. | both | HS19 Sec. 3; D20 | No | F/G/K field 89,496 / 179,780 / 103,636; T1 102 / 133 / 56 (paper 119 / 146 / 66) |
| P25 | Pleiades MS (check only) | Members from J/A+A/616/A10 tablea1a (1,326). **Cuts 1, 8, 9, 10** (`config.yaml:89`), EGGR 25 removed, E(B−V) = 0.045. Spline `s = N·σ²` (σ = MAD scatter) → 8th-order polynomial (`src/pleiades.py:60-90`). **The selection uses the published polynomial, not this refit.** | Pleiades | HS19 Sec. 2; D10, D10b | **Yes:** paper says cuts 1, 5, 8–10. Cut 5 removes 88 % of members, including EGGR 25, which HS19 removed by hand (D10b). | 1,178 → 1,177 members |

### 1b. Monte Carlo control matching

| # | Name | Exact definition as implemented (file:line) | Source | Differs from paper? |
|---|---|---|---|---|
| P26 | Matching variables | `match_vars: [z]` (`config.yaml:139`) | HS19 | No |
| P27 | z-matching (step a) | **Fixed 50 pc bins**, edges from `floor(min z_host/50)·50` to `ceil(max z_host/50)·50` (`src/matching.py:19-35`; `config.yaml:143-147`). Field weight `w = (host fraction in bin) / (field count in bin)`; 0 outside the host z range or in host-empty bins (`src/matching.py:47-80`). **Bernoulli thinning** keeps each field star with `p = w / max(w)`, a new subsample every iteration (`src/matching.py:107,117-120`) | D12 (ours; paper: "matched in its z distribution", method not given) | Paper gives no method. Variants: 25/100 pc bins, 5/10/20 quantile bins, weighted resampling, no matching. |
| P28 | Usable field pre-restriction | Only field stars with `w > 0` and `bp_rp0` within [min host − 0.025, max host + 0.025] (`src/matching.py:101`). Optimisation only. | ours | – |
| P29 | Colour matching (step b) | For each host, one star drawn **uniformly** from the thinned subsample with `abs(bp_rp0_host − bp_rp0_field) ≤ 0.025` (**inclusive**) (`src/matching.py:127-139`; `config.yaml:134`) | tolerance 0.025 (D13); draw rule ours | The tolerance is not flagged as a deviation in DECISIONS |
| P30 | Replacement within an iteration | **Without replacement:** a field star is used at most once per control sample. Collisions are resolved in random host order by re-drawing from the remaining candidates (`src/matching.py:141-161`; `config.yaml:148`) | D13 (ours) | Paper silent. Variant: with replacement. |
| P31 | Hosts with no colour match | Dropped **from both** the control and the host dispersion of that iteration; p-value paired per iteration (`src/matching.py:163-193`). None excluded permanently (`config.yaml:135`). | D13b (user, 2026-10-30) | Paper silent. KELT-9 is unmatched in 22 % of iterations, KELT-20 in 6 %. |
| P32 | Number of draws | **40,000** control samples (`config.yaml:133`). Sensitivity variants: 5,000 (`config.yaml:153`). | HS19 (40,000) | No |
| P33 | Seeds | S14. Labels `hj_table1`, `hj`, `longp`. | D19 | – |
| P34 | Test statistic and p | p = fraction of controls with `σ_control ≤ σ_host` (paired with that iteration's host value, P31) (`src/matching.py:203-207`; `src/pipeline.py:194`). Gaussian-equivalent = `(median_control − σ_host) / SD_control(ddof=1)` (`src/pipeline.py:210`). Percentiles 2/16/50/84/98 reported (`:211`). | D18; spec | HS19 quote "< 1 in 40,000". Their Fig. 4 shows the same percentiles (A9). |
| P35 | Interval type | **None.** Phase 1 reports the control distribution only. There is no bootstrap of the hosts. | HS19 | No |

Phase 1 primary result (logged): σ_host 36.18 km/s, control median 42.82, 0/40,000 (REPRODUCTION_REPORT §3).

---

## 2. EBLM analysis (frozen `PREREG_eblm_kinematics_2026-11-02`)

Driver: `python -m src.eblm_kinematics` (`src/eblm_kinematics.py:204-276`). The pipeline is the
Phase 1 pipeline (S1–S15, P8–P11, P17, P19–P21) except as listed. There is no published paper
for this analysis. "Differs" compares with the frozen pre-registration, or with EBLM IV where
EBLM IV data are used.

### 2a. Sample, RVs and cuts

| # | Name | Exact definition as implemented (file:line) | Applied to | Source | Differs? | Remaining (EBLM / field) |
|---|---|---|---|---|---|---|
| E1 | Sample source | EBLM IV (Triaud+2017), `tab:obs` / `tab:primary` / `tab:params` parsed from the arXiv source 1707.07521 (`src/eblm.py:19-20,54-83`; `config.yaml:157`). **EBLM XIX is excluded.** | EBLMs | C1, C8; Stage D user decision | – | 118 / – |
| E2 | Gaia DR2 match | One ARI cone search, **10″** around each EBLM IV J2000 position. Gaia positions are propagated **to J2000 with the DR2 proper motion** (`dt = −15.5 yr`); the **nearest** source is the match (`src/eblm.py:86-119`; `config.yaml:158-160`). Flags `sep < 3″`, `G − V ∈ [−1.0, 0.3]` and `dG_second < 2` are **logged only**: all 118 are kept (`src/eblm.py:122-131`; C11) | EBLMs | C11 (ours) | – | 118/118 matched (114 "ok", 4 inspected and kept) |
| E3 | Gaia rows | Same query as the hosts (DR2 + BJ18) via `targets.build` → `hosts.gaia_rows` (`src/targets.py:53-64`) | EBLMs | C12 A | – | 118 |
| E4 | Gaia RV never used | Gaia `radial_velocity` is set to NaN for binaries (`src/targets.py:59`). RV is taken only from `rv` with `rv_kind == "gamma"` (`src/targets.py:39-50`). | EBLMs | C11, C12 A | – | – |
| E5 | Cuts 1–4 | As P8–P11, host cut set (`config.yaml:88`; `src/pipeline.py:113-121`) | EBLMs | C12 A | No | 118 → 117 → 115 → 107 → **107** (C12 A) |
| E6 | Distance, extinction, MS | S3–S7, S11 | EBLMs | C12 A | No | **79** (F 54, G 24, K 1) |
| E7 | γ (systemic RV) model | `RV = γ + K[cos(ν+ω) + e cos ω] + fixed terms`, weighted linear least squares in (γ, K), weights `1/(svrad² + σ_add²)` with σ_add from `tab:obs` (m/s → km/s) (`src/gamma_fit.py:9-11,66-85,185`). Only the `jdb, vrad, svrad` columns of the CORALIE `.rdb` files are read (D1). | EBLMs | 07 amendment; C6, D1 | **vs EBLM IV:** no bisector-based outlier rejection and no removal of in-eclipse points (the amendment permits no other columns; D1) | 118/118 computable |
| E8 | Fixed orbit inputs | P and T0,pri from `tab:obs`; e and ω from `tab:params`. **e = 0 where only an upper limit is given** (65 systems); ω = 0 if e = 0 and ω is missing (`src/gamma_fit.py:181`). Eclipse at ν = π/2 − ω (`:32-36`) | EBLMs | C12 A; D1 | – | – |
| E9 | Drift / outer orbit | Published drift polynomial subtracted with **t_ref = midpoint of each system's RV span** (`src/gamma_fit.py:58,200-202`; `config.yaml:170`). The 4 `tab:triple` outer orbits are subtracted as fixed Keplerians (`:205`). | EBLMs | D2, D7 (user at STOP D1) | EBLM IV does not publish t_ref (D2) | – |
| E10 | CORALIE07→14 offset | **0 km/s** for all RVs at BJD ≥ 2456977.5 (`src/gamma_fit.py:184`). **0.020 km/s** added in quadrature to σ_γ (`config.yaml:168-169`). No per-system offset. | EBLMs | C6 (user 2026-11-01) | – | – |
| E11 | γ flags | `< 3 RVs` or `cond(normal matrix) > 1e12` → uncomputable, excluded (`config.yaml:171-172`). Only `flag == "ok"` gets an RV (`src/eblm_kinematics.py:36-39`). The "K disagrees" and "σ_γ > 0.5" flags are descriptive only (D3). J1500-33 is kept (D8). | EBLMs | C12 A; D3 | – | 79 (none excluded) |
| E12 | γ errors | σ_γ is **not** propagated into UVW. A deconvolved value would be reported only if the median σ_γ > 1 km/s (not triggered: 0.028) | EBLMs | C12 F; D13 | – | – |
| E13 | Field | Phase 1 field (P5–P23, same host exclusions) **minus all 118 EBLM IV source_ids** (`src/eblm_kinematics.py:59-67`) | field | D9 (pre-registration silent) | – | – / **376,786** (376,791 − 5) |

### 2b. Matching, statistic and intervals

| # | Name | Exact definition as implemented (file:line) | Source | Differs? |
|---|---|---|---|---|
| E14 | Matching variables (primary) | `match_vars: [z, dist_pc]`, a **joint 2-D histogram**: z in fixed 50 pc bins (P27), `dist_pc` in **5 quantile bins of the 79 EBLM distances** (the last edge nudged up to include the max) (`src/eblm_kinematics.py:223`; `config.yaml:140-141`; `src/matching.py:31-34`). Weights and thinning as P27. | C12 B | Phase 1 uses z only |
| E15 | Colour match / replacement / unmatched | As P29–P31 (0.025, without replacement, dropped from both sides) | C12 B | No (no EBLM was ever unmatched; D11) |
| E16 | Statistic | S12, mean residual speed | C12 A | – |
| E17 | MC draws, seed | 40,000; label `eblm_primary` (S14) | C12 B | – |
| E18 | Point estimate | `Δσ̂ = σ_EBLM − median(σ_control)` (`src/eblm_kinematics.py:82,116-123`) | C12 B | – |
| E19 | Bootstrap CI | **40,000 resamples of the 79, with replacement** (`rng.integers`). The control median is **held fixed**. 95 % **percentile** interval `np.percentile(d, [2.5, 97.5])`. Seed `SeedSequence([seed, label code, 99])` (`src/eblm_kinematics.py:75-86`) | C12 B; D10 | – |
| E20 | Bias / BCa rule | `bias = mean(d) − Δσ̂`. A BCa interval (z0 from `mean(d < est)`, jackknife acceleration) is computed **only if** `abs(bias) > 0.2 × half-width`. The percentile interval stays primary (`src/eblm_kinematics.py:84-97`) | C12 B; D10 | Not triggered (bias −0.27, threshold 0.89) |
| E21 | C-1 (HS19-style) | As P34 | C12 C | – |
| E22 | C-2 EBLM vs HJ | 40,000 iterations: bootstrap the 79. For each drawn EBLM, draw one Phase 1 HJ host (297) with `abs(Δbp_rp0) ≤ 0.05` **and** in the same EBLM-distance quintile, **with replacement**. HJ hosts outside the EBLM distance range are ineligible; EBLMs with no candidate are dropped from both sides. Seed stream 2 (`src/eblm_kinematics.py:126-153`) | C12 C; D10 | 2 EBLMs have no candidate (D11) |
| E23 | G sensitivity tests | **5,000** MC + **5,000** bootstrap each (`src/eblm_kinematics.py:231-252`): G1 z only; G2 5×5 quantile z × distance; G3 abs(z) (50 pc) + distance; G4 abs(z) only; G5 field `dec < 0`; G6 without `drift`/`triple` (n 68); G7 without `active` (77); G8 RMS; G9 F only (54) + C-2 vs F-type HJ; G10 simple average instead of γ; **G11 without J1500-33 (78), added after STOP D1** | C12 G; D8, D10 | G11 is not in the frozen file |
| E24 | I (exploratory) | Evolved = fails the MS cut (S11). Samples: EBLMs after cuts 1–4 (107) and **archive HJ list ≤ 2019 after cuts 1–4 (457)**, no RV requirement. 0.1 mag colour bins from the bluest EBLM; HJ hosts reweighted to the EBLM colour distribution; **10,000 joint bootstraps**, stream 3 (`src/eblm_kinematics.py:157-200,259-273`). Variants: 0.05/0.2 mag bins, 0.5 mag threshold, current archive. | C12 I; D10 | – |
| E25 | [Fe/H] (E) | **Dropped:** spectroscopic [Fe/H] for 27/79 < 40 required | C12 E; D12 | – |

EBLM primary result (logged): σ_EBLM 42.27, control median 41.77, Δσ +0.50 [−4.23, +4.71]; C-2 +7.90 [+2.07, +13.74] (`logs/eblm_d2.log`).

---

## 3. WASP-only analysis (frozen `PREREG_wasp_kinematics_2026-11-04`)

Driver: `python -m src.wasp_kinematics` (`src/wasp_kinematics.py:130-227`). Design principle
(pre-registration §0, E7): the **only** change from Phase 1 is keeping WASP-discovered hosts.
Everything in §0 and §1 applies. "Differs from paper" compares with Hu+, the paper this test
addresses.

| # | Name | Exact definition as implemented (file:line) | Applied to | Source | Differs from Hu+? | Remaining |
|---|---|---|---|---|---|---|
| W1 | Base sample | Phase 1 final Table 1 hosts, `hosts_hj_table1_final.csv` (297) (`src/wasp_kinematics.py:141`) | hosts | prereg A | Hu+ use their updated HS19 list, DR3 | 297 |
| W2 | WASP membership (primary) | `facility` = archive `pscomppars.disc_facility` (snapshot 2026-11-04, `data/raw/exoarchive/pscomppars_disc_2026-11-04.csv`) of the host's planets that pass the **P2 archive HJ selection (incl. disc_year ≤ 2019)**. Labelled WASP if the string contains "SuperWASP-North", "SuperWASP-South"/"WASP-South", or "WASP" (`src/hu2026_checks.py:110-147`). Keep `facility.notna()` (`src/wasp_kinematics.py:142-144`). **Join key: see the warning below.** | hosts | prereg A; E3, E9 | Hu+ separate surveys by their own labels | **98** |
| W3 | Membership by name (D1) | W2 **or** `hostname.startswith("WASP-")` (`src/hu2026_checks.py:146`; `src/wasp_kinematics.py:145`). Adds WASP-53, WASP-81. | hosts | prereg D1; E9 | – | 100 |
| W4 | No distance / colour cut (primary) | none | hosts | prereg A; E11, E12 | **Yes:** Hu+ cut `0.5 ≤ (BP−RP)0 ≤ 1.5` and use `d < 600 pc` subsamples | 98 |
| W5 | Field (WASP half) | Phase 1 `field_ms.csv.gz` (P23) (`src/wasp_kinematics.py:146`) | field | prereg A | Hu+: DR3, halo cut v_tot > 180 | 376,791 |
| W6 | Field (EBLM half) | E13: Phase 1 field minus the EBLM IV IDs (`:148-149`) | field | prereg A | – | 376,786 |
| W7 | Matching (primary) | z only, 50 pc fixed bins + colour 0.025, without replacement, D13b (P26–P31) (`src/wasp_kinematics.py:157`) | – | prereg A | **Yes:** Hu+ match per host on angular separation ≤ 10°, `d_HJ − 100 ≤ d_FS ≤ d_HJ + 50 pc`, `G_HJ − 1 ≤ G_FS ≤ G_HJ + 0.1` and `abs(Δ(BP−RP)0) ≤ 0.2`, **with** replacement, 10,000 iterations (E1). Their matching is dropped for DR2 (E8). | – |
| W8 | Statistic (primary) | S12, mean residual speed (`src/wasp_kinematics.py:31-32`) | – | prereg A | **Yes:** Hu+ use σ_W | – |
| W9 | MC draws, seed | 40,000; label `wasp_primary` (S14) | – | prereg A; E14 | Hu+ 10,000 | – |
| W10 | Δ_W and CI | `Δ_W = σ_WASP − median(control)`. **40,000** bootstrap resamples with replacement, control median fixed, **percentile** 95 % CI, BCa only if `abs(bias) > 0.2 × half-width` (the E19/E20 code, `src/wasp_kinematics.py:103-112`) | – | prereg B; E27 | – | – |
| W11 | δ_HS19 and decision rule | `δ = 36.1825 − 42.8234 = −6.6409 km/s` (`src/wasp_kinematics.py:25`). Six cases on [L, U] with closed ends (`:49-62`) | – | prereg §2, B; E18 | – | case 1 |
| W12 | Contrast C | `(σ_EBLM − ref_E) − (σ_WASP − ref_W)`, both control medians fixed. **Independent** 40,000-draw bootstraps of each sample (labels `wasp_C_W` / `wasp_C_E`, stream 99), percentile CI (`src/wasp_kinematics.py:115-126`). EBLM half = cached `eblm_primary` MC (z + distance), asserted median 41.77 (`:163-164`). | – | prereg C; E30 | – | 98 / 79 |
| W13 | D sensitivity | 5,000 MC + 5,000 bootstrap; a half the test does not change keeps its 40,000 reference (`src/wasp_kinematics.py:171-198`). **D2:** WASP on `[z, dist_pc]`, 5 quantiles of the WASP distances. **D3:** `dist_pc < 600` (strict, BJ18), both samples. **D4:** `0.5 ≤ bp_rp0 ≤ 1.5` (inclusive), both. **D5:** keep field stars with `sqrt(U²+V²+W²) ≤ 180` (LSR frame, S8), both fields. **D6:** EBLMs without J1500-33 (C only). | – | prereg D; E26, E30 | D3–D5 adopt Hu+'s cuts one at a time | D1 100; D3 97 / 71; D4 95 / 76; D5 788 field stars removed |
| W14 | σ_W (exploratory) | `σ_W = sqrt((1/N) Σ (W − W̄)²)` (`src/wasp_kinematics.py:35-37`), with the **matching of A**, 40,000 controls (`:82-99`) and 40,000 resamples. The HS19-sized deficit is `−(1 − 36.1825/42.8234) × median(σ_W control)` (`:26,212`) | – | prereg E; E16, E25 | Same formula as Hu+ eq. 2; matching not Hu+'s | – |

WASP primary result (logged): σ_WASP 37.23, control median 43.04, Δ_W −5.81 [−9.69, −2.09], case 1;
contrast +6.31 [+0.58, +12.02] (`logs/wasp_kinematics.log`).

> **⚠ Membership join (found while writing this note). Update 2026-11-09: fixed in code (source_id join, 107 hosts; DECISIONS E32); the frozen run stays pinned to the old join; amendment drafted (E33), not run.**
> In W2 the facility is looked up by `hostname` (`src/hu2026_checks.py:145`). For Table 1 rows,
> `hostname` is **HS19's SimbadName** (`src/hosts.py:193`), but the facility table is keyed by the
> **archive** hostname (`src/hu2026_checks.py:127`). Of the 297 final hosts, 36 have a SimbadName
> different from their ArchiveName. **Nine of them are WASP hosts** and therefore got no label:
> WASP-70 A (TYC 5782-970-1), WASP-147, WASP-151 (EPIC 246441449), WASP-156, WASP-158, WASP-161,
> WASP-168, WASP-172 (2MASS J13174412-4714152) and WASP-173 A. Their `disc_facility` in the cached
> 2026-11-04 file is SuperWASP or SuperWASP-South. The name variant (W3) misses them for the same
> reason. One caveat: W2 also requires the planet to pass the `disc_year ≤ 2019` archive
> selection, so not all 9 would necessarily be labelled even with the right key. Whether the
> intended sample is 98 or up to 107 has not been computed. The frozen pre-registration states
> 98, so any change needs a dated amendment (your call). The same join feeds
> `notes/hu2026_wasp_counts.md` and the power estimates (E3, E4).

---

## 4. Choices most likely to differ from another implementation

1. **Statistic.** Eq. 1 as typeset (mean of the 3-D residual speeds), not the RMS (S12). The RMS
   gives 40.5 instead of 36.2 km/s for the same hosts. Hu+ use σ_W.
2. **Cut 4 parenthesis and cut 8 coefficient.** Lindegren form `exp(−0.2(G − 19.5))` and 0.015
   (P11, P15). The printed cut 4 would drop about half the hosts.
3. **Distance.** BJ18 `r_est` throughout, not 1/parallax or BJ21/DR3 (S3).
4. **Extinction.** STILISM v4.1 cube (the Lallement+2018 map, not Capitanio+2017 as cited),
   1 pc midpoint integration, zero outside the cube, CV18 *mean* coefficients (S4–S6).
5. **Main sequence.** Published polynomial (not the spline), strict 1 mag, no lower and no colour
   limit (S11). This step alone sets the field count (376,791 vs 385,036).
6. **Host RVs and list.** HS19 Table 1 with today's archive `st_radv`, else DR2, with no DR3 fill.
   That gives 297 hosts, not 338 (P22). The archive snapshot date also matters (P2).
7. **Host cut set 1–4 only; Pleiades cuts 1, 8, 9, 10** (P17, P25).
8. **z definition and matching.** Galactocentric z with z_sun = 27 pc, R0 = 8.3 kpc; 50 pc
   histogram thinning (S10, P27). Many implementations use d·sin b, or nearest-neighbour matching.
9. **Colour matching details.** Tolerance inclusive at 0.025; without replacement within a
   control sample; hosts with no match dropped from both sides; p paired per iteration
   (P29–P31, P34).
10. **Intervals (EBLM/WASP).** Percentile bootstrap of the target sample only, with the control
    median held fixed (no Monte Carlo uncertainty added), 40,000 resamples (E19, W10). BCa only
    by rule.
11. **EBLM specifics.** Distance-quintile matching on top of z; field minus the EBLM IDs; γ from
    the fixed-ephemeris fit with e = 0 for upper limits and drift t_ref = mid-span (E7–E14).
12. **WASP membership.** `disc_facility` label, the disc_year ≤ 2019 archive selection, and the
    SimbadName join key flagged above (W2).
13. **Seeds.** Chunked `SeedSequence([seed, label-char-sum, chunk])`. Results will not match
    draw-for-draw anyway, but control medians should agree to about 0.01–0.02 km/s (S14).
