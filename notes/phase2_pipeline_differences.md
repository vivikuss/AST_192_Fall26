# Phase 2 groundwork: Kayla's two pipelines vs our frozen DR2 pipeline

**From reading the code only (DECISIONS E50).** Neither script was run, imported or executed. No number in this
note was computed from either script or from Kayla's data. Where an effect is not stated, it is not known.

- **Ours (frozen DR2):** the Phase 1 pipeline (`config.yaml`, DECISIONS D1–D20) as frozen in
  `PREREG_wasp_kinematics_2026-11-04.txt`, with the membership of amendment 1
  (`PREREG_wasp_kinematics_amendment1_2026-11-10.txt`).
- **K-DR2:** `velocity_dispersion_replication.py` (sha256 `7829cee3…`), her HS19 replication.
  - Its docstring says its choices were aligned with our `notes/cuts_and_choices.md` at commit 55aa85d.
  - It lists the settings where it was "not aligned".
- **K-DR3:** `hu2026_replication.py` (sha256 `6f7bd10a…`), her Hu et al. (2026) replication. It uses Gaia DR3.
  - With `USE_PUBLISHED_SAMPLES = True` (as received), the HS19 and TESS analyses use the hosts of Hu et al.
    Tables 2 and 3. Its own reconstruction is only compared against them.

Both scripts live at `~/data/hu2026_dr3_Kayla/code/` (outside every repository).

**Legend:** "= ours" means the same setting as ours; **bold** marks a difference from ours.

## Quality cuts

| # | Setting | Ours (frozen DR2) | K-DR2 | K-DR3 |
|---|---|---|---|---|
| 1 | BP/RP excess-factor lower bound (HS19 cut 8) | 1.0 + **0.015**·(BP−RP)², Lindegren C.2 (D2). The server prefilter uses 0.0015; 0.015 is applied locally | **0.0015** (as printed in HS19; listed as "not aligned") | **0.0015** (Hu App. A cut 7, hard-coded) |
| 2 | mean_varpi_factor_al in [−0.23, 0.36] (cut 2) | field and hosts | = ours | **not applied** (no such column in DR3) |
| 3 | astrometric_excess_noise_sig < 2 (cut 6) | field | = ours | **not applied** |
| 4 | RUWE | no cut | = ours | **ruwe < 1.4** (field) |
| 5 | Field parallax prefilter | none | = ours | **parallax > 0.4 mas** (d < 2500 pc) |
| 6 | Host cuts | HS19 cuts 1–4 (parallax_over_error, mean_varpi_factor, visibility periods, UWE) (D3c) | = ours | **App. A cuts 1–3** (parallax_over_error, visibility periods, UWE; no mean_varpi_factor) |
| 7 | Cuts on the hosts actually analysed | every host passes cuts 1–4, the MS cut and has an RV | = ours | **the published hosts of Hu Table 2/3 are used as listed**: only a W value is required, and the cuts are applied only to her reconstruction |

The other quality cuts are written the same way in all three:
- parallax_over_error > 10, visibility_periods_used > 8, astrometric_gof_al < 3, rv_nb_transits > 5, BP/RP
  flux_over_error > 10;
- the UWE cut in the Lindegren C.1 form;
- an RV is required for field stars.

## Distances, extinction, main sequence

| # | Setting | Ours (frozen DR2) | K-DR2 | K-DR3 |
|---|---|---|---|---|
| 8 | Distance | Bailer-Jones (2018) `r_est` (D16) | = ours (ARI `gaiadr2_complements.geometric_distance`) | **1000 / parallax** |
| 9 | Extinction source | STILISM cube v4.1 (2018-03-19 Wayback copy; D9, A1) | STILISM `stilism_cube.h5`. Her methods string says "v4.1", but **the file was not provided, so the version is unverified**. **If the file is missing, E(B−V) = 0 with only a warning** | **Gaia GSP-Phot** `ag_gspphot` / `ebpminrp_gspphot`, else **Bayestar19** (median of samples) × 0.884, else **zero** (no Bayestar19 coverage: dec < −30°, d < 63 pc). **If `bayestar2019.h5` is missing, zero with only a warning** |
| 10 | Line-of-sight integration | ceil(d / 1 pc) equal steps, midpoint, trilinear, zero outside the cube | = ours | not used (Bayestar19 is read per distance). The STILISM option, off as received, would use **5 pc steps** |
| 11 | Colour limit on the main sequence | **none** (D8b) | = ours | **0.5 ≤ (BP−RP)₀ ≤ 1.5**, hosts and field (UHJ: 0.8–1.5) |
| 12 | Halo cut | none in the primary (D5 is a sensitivity) | = ours | **field stars with v_tot > 180 km/s removed** |

The same in all three: the Pleiades-ZAMS 8th-order polynomial with stars < 1 mag above it kept, and the
Casagrande & VandenBerg (2018) coefficients R_G, R_BP, R_RP = 2.740, 3.374, 2.035.

## Field sample

| # | Setting | Ours (frozen DR2) | K-DR2 | K-DR3 |
|---|---|---|---|---|
| 13 | Planet hosts removed from the field | every source_id in the HS19 Table 1, archive HJ and archive long-period host lists (disc_year ≤ 2019) (D17) | the same lists' hosts that pass cuts 1–4 (equivalent, since the field cuts include 1–4), **plus the warm-Jupiter hosts and Hu Table 3 TESS hosts** (DR2 IDs via the closest `dr2_neighbourhood` match) | **every NEA host with a DR3 ID (any planet, current NEA) plus the WJ/CJ hosts** |

## Kinematics

| # | Setting | Ours (frozen DR2) | K-DR2 | K-DR3 |
|---|---|---|---|---|
| 14 | Height z | Galactocentric, R0 = 8.3 kpc, z_sun = 27 pc (D11) | = ours | **not computed** (no z matching) |

The same in substance in all three: UVW relative to the LSR with solar motion (11.1, 12.24, 7.25) km/s. Ours and
K-DR2 use astropy `GalacticLSR`; K-DR3 adds the solar motion to the heliocentric Galactic velocity, which is
the same transformation.

## Matching

| # | Setting | Ours (frozen DR2) | K-DR2 | K-DR3 |
|---|---|---|---|---|
| 15 | Matching variables | z histogram (fixed 50 pc bins, Bernoulli thinning, D12), then one field star with \|Δ(BP−RP)₀\| ≤ 0.025 per host (D13) | = ours | **per-host pool** (Hu Table 1): separation ≤ 10°, −100 ≤ d_FS − d_HJ ≤ +50 pc, −1 ≤ G_FS − G_HJ ≤ +0.1, \|Δ(BP−RP)₀\| ≤ **0.2**. WJ/CJ: 20°, ±100 pc, −1 to +0.5 mag |
| 16 | Reuse of a field star within a control sample | **without** replacement; collisions resolved in random host order (D13) | = ours (same rule; the code was not compared line by line with ours) | **with replacement** |
| 17 | Host with no match | dropped from both sides **in that iteration**; P paired per iteration (D13b) | = ours | **dropped from the sample entirely** if its pool is empty |
| 18 | z-bin top edge | if the highest host z lies exactly on a bin edge, one bin is appended | that host falls in the last bin (`np.histogram`); a boundary detail | (no z bins) |

## Statistic and inference

| # | Setting | Ours (frozen DR2) | K-DR2 | K-DR3 |
|---|---|---|---|---|
| 19 | Statistic | HS19 eq. 1, (1/N) Σ\|v_i − v̄\| (D14); σ_W is exploratory only (prereg E) | = ours | **σ_W = std(W), ddof = 0** (Hu eq. 2) is primary. σ_tot and eq. 1 are also run, for "HS19 All" only |
| 20 | Control samples | 40,000 (sensitivity tests: 5,000) | 40,000 (warm Jupiters: 10,000) | **10,000** |
| 21 | Result and decision | Δ = σ_hosts − median(σ_control). 95 % **percentile bootstrap CI** over 40,000 host resamples, control median fixed. Six-case rule against δ_HS19 = −6.64 km/s (prereg B). P reported alongside | **P(control ≤ host)** paired per iteration, and (median − host) / SD. **No CI, no decision rule** | **Tail fraction on the hosts' side of the median, converted to a one-sided Gaussian σ**, plus the "age ratio" (host / median)². **No CI, no decision rule** |
| 22 | Random numbers | `SeedSequence([20191021, label code, tail])` (D19) | **`default_rng([42, crc32(sample key)])`** | **`default_rng([2026, crc32(key/statistic)])`** |

## Host sample, RVs and IDs

| # | Setting | Ours (frozen DR2) | K-DR2 | K-DR3 |
|---|---|---|---|---|
| 23 | Starting HJ list | HS19 Table 1 as published (338 → 297 after cuts, MS and RV); not re-selected | = ours (Table 1 as is; archive values joined on ArchiveName) | Table 1 **mapped to the current NEA and re-selected with current NEA values**: P < 10 d, Mp sin i > 0.1 MJ, **no mass-provenance filter**, discovered by transit or RV. **Replaced by Hu Table 2 (282 hosts)** as received |
| 24 | WASP membership | archive `disc_facility` of the host's HJ(s) contains a WASP facility, joined on the Gaia DR2 source_id (amendment 1): **107 hosts**, South and North, no distance cut | **no WASP sample** | `disc_facility` of the host's **earliest-discovered HJ** contains "WASP", **ground-transit label**, **d < 600 pc**, within Hu Table 2 (Hu et al.: 102) |
| 25 | Distance cut | none in the primary (D3 sensitivity: d < 600 pc) | = ours | **d < 600 pc** for the WASP / HATNet / HATSouth / TG / TESS runs. No cut for "HS19 All" and RV; UHJ uses d < 600 pc |
| 26 | Host exclusions | none beyond the cuts, MS and RV | = ours | **close stellar companion** (her "best guess" list: KELT-2 A, KELT-4 A, WASP-70 A, WASP-94 A, WASP-94 B, NGTS-3 A, WASP-173 A); **cluster members** (NGC 2682 YBP 1194/1514/401, Pr0201, Pr0211); **extreme G** (ups And, POTS-1); **not discovered by transit or RV** |
| 27 | Host RV | archive `st_radv` if `st_radvlim` is 0 **or missing**, and its reference does not contain "Gaia DR3", "Gaia EDR3", "gaia_dr3" or "GAIA DR3"; else Gaia DR2 RV. No fill (D6c): hosts without an RV are dropped | archive `st_radv` only if `st_radvlim` **equals 0 (missing → rejected)**; **broader reference filter** (also "Gaia Collaboration et al. 2021–23", "Katz et al. 2022/23", `Gaia\s*E?DR3`); else Gaia DR2 RV. TESS sample only: **Gaia DR3 RV fill** | **NEA `st_radv` with no limit or reference filter** (those columns are not downloaded), **else Gaia DR3 RV** |
| 28 | Host Gaia IDs | DR2 from SIMBAD (Harvard mirror) by SimbadName, else the archive `gaia_dr2_id` (D6b) | DR2 from SIMBAD (CDS) by SimbadName, else the archive `gaia_dr2_id`: the same rule | **DR3**: the Hu Table 2/3 IDs (**TOI-5261 differs from the NEA**), else the NEA `gaia_dr3_id`, else the closest `gaiadr3.dr2_neighbourhood` match (per E47 used only for nu Oct A) |
| 29 | NEA snapshot | `pscomppars` cached 2026-10-28 (facility labels: 2026-11-04 cache) | her download; **date not recorded** | her `nea_pscomppars.csv`; **date not recorded**. TESS, UHJ and WJ/CJ selections keep disc_year ≤ 2025 |

## Consistency with the files and answers already logged

The following were checked by reading only:
- **E45 (field ADQL):** K-DR3's field query (rv not null, parallax > 0.4, App. A cuts 1–9 with 0.0015) matches the
  WHERE clause recorded from the `.adql` files.
- **E47 (Kayla's answers):** the `.key` files are written by `cache_key` (SHA-1 of the build settings). The
  DR2→DR3 fallback runs only for hosts without an NEA DR3 ID. The colour tolerances are 0.025 (K-DR2) and 0.2
  (K-DR3).
- **Not checked:** whether the files in `~/data/hu2026_dr3_Kayla/` were produced by exactly these versions of
  the scripts.

## Open for the Phase 2 pre-registration (no decision here)

K-DR3 differs from ours in almost every block at once: cuts, distances, extinction, colour window, halo cut,
matching, statistic, inference, host list, exclusions and RVs. A K-DR3-vs-ours difference in a result could not
be attributed to any single setting. K-DR2 differs from ours mainly in:
- the cut-8 coefficient (row 1);
- the host-RV filters (row 27);
- the extra field exclusions (row 13);
- the inference reported (row 21);
- the RNG (row 22).

This adds to the open questions of E46 and E47.
