# DR3 WASP analysis: EXPLORATORY plan (Stage 1 of 2)

**EXPLORATORY — not pre-registered.** This was decided on 2026-11-11, before any WASP DR3 kinematics were computed
(DECISIONS E52). Results will be preliminary, and this analysis cannot later be pre-registered as confirmatory.
The Stage 1 record is DECISIONS E53.

**Stage 1 computed no WASP host velocity, dispersion or result.** Nothing was queried for the WASP hosts in
Stage 1.

## 0. Status at the Stage 1 stop

| Item | Status |
|---|---|
| E52 exploratory declaration | committed (7f40c4b) before any other step |
| Seed labels `dr3x_*` | registered. The collision tests pass, and no new label shares a code with any label |
| Method A code (`src/dr3x_method_a.py`) | written; unit tests on synthetic data only |
| **Method A run and validation** | **Validation done (E60):** HS19-all N 279, σ_W 15.77 vs field 18.83 (paper 15.8 vs 18.7); TESS d<600 N 100, 18.07 vs 17.65 (paper 18.1 vs 17.5). Bayestar19 from the local file in Kayla's order (E59, E60) |
| Hu et al. published host tables (validation input) | not available locally. Validation uses only the name and DR3-ID columns of Kayla's `hosts_hs19.csv` / `hosts_tess.csv` (E54) |
| Method B field code (`src/dr3x_field_b.py`) | written |
| **Method B field download** | **started 2026-11-11 at ARI** (detached, pid 47363; log `~/data/dr3x_field_b/download.log`). The server query fixes none of the open choices (§3.1) |
| Host step counts, 107 WASP hosts (E55, E60) | Method A main 97, without exclusions 101 (6 hosts resolved with Bayestar19, E60); Method B 101 (`astrometric_params_solved = 31` removes WASP-29, -34, -120, -94 A, -100) |
| Host kinematics | none computed |

## 1. Method A: Hu et al. (2026) as implemented by Kayla

- **Source:** our own code, written from reading `hu2026_replication.py` as text (sha256 `6f7bd10a…`, E50). Her
  script was not imported, run or copied into the repository; the docstring cites it.
- **Field:** her `field_main_sequence_uvw.csv`. Its build settings were reproduced from her source text (cache
  key match, E51). Every NEA host with a DR3 ID (her `nea_pscomppars.csv`) is removed at load, and the count is
  reported.
- **Hosts:** the 107 WASP hosts of frozen amendment 1, with DR3 IDs from `data/processed/dr3_crossmatch/`
  (both routes agree for all 107, E48).
  - Gaia DR3 query: her host columns plus `ruwe`.
  - Steps, in her order, each reported by name:
    1. no Gaia DR3 row;
    2. App. A cuts 1–3 (parallax_over_error > 10, visibility_periods_used > 8, UWE in the Lindegren C.1 form);
    3. RV from NEA `st_radv`, else Gaia DR3 RV (only the presence of an RV is used);
    4. extinction;
    5. MS (< 1 mag above the HS19 ZAMS) plus 0.5 ≤ (BP−RP)₀ ≤ 1.5;
    6. d = 1000/parallax < 600 pc.
- **Extinction:**
  - GSP-Phot `ag_gspphot` / `ebpminrp_gspphot` if both exist;
  - else Bayestar19: the per-bin median of the samples × 0.884, with CV18 coefficients, interpolated linearly
    in distance modulus and held at the last bin;
  - **else zero only where the map has no value by design**: south of dec −30° (outside the map) or closer than
    its first distance bin (d < 63 pc).
  - Each host's `ext_source` records which case applied, and the zero cases are listed by name.
  - **Bayestar19 source (E54):** the full map is not downloaded. Hosts without GSP-Phot only are queried through
    the dustmaps web query (bayestar2019, median mode, × 0.884). Positions without coverage get zero, as in
    Kayla's method. **If the web service fails, STOP** (it did on 2026-11-11). Her script, by contrast, warns
    and uses zero when the map file is missing (E50).
  - Differences from her file-based reader: the web service interpolates in distance itself. It may also
    differ in whether it takes the median of the samples before or after interpolating, and in its handling
    of distances outside the map's range. To match her rules, we set zero below DM 4 (d < 63 pc) and query at
    the last bin beyond it. This is to be checked against the service when it is available.
  - Note: many WASP-South hosts lie south of dec −30°. Those without GSP-Phot get zero extinction under her
    rules. The paper does this too; it is reported, not hidden.
- **Matching and statistics:**
  - Per-host pool, all bounds inclusive: separation ≤ 10°, −100 ≤ d_FS − d_HJ ≤ +50 pc,
    −1 ≤ G_FS − G_HJ ≤ +0.1, |Δ(BP−RP)₀| ≤ 0.2.
  - 10,000 draws of one pool star per host, with replacement. Hosts with an empty pool are dropped and named.
  - Statistics: σ_W = std(W) (ddof 0; primary), σ_tot, HS19 eq. 1.
  - Reported: host value, control median, 1σ = [16, 84] and 2σ = [2.5, 97.5] percentiles, her tail fraction
    and one-sided Gaussian σ (floor 1/N).
- **Adaptations and choices** (her script vs ours):
  - A1. *Membership:* our 107 (frozen amendment 1, DR2 source_id facility join; South and North), as the user
    specified. Hers is Hu Table 2 hosts whose earliest HJ has a WASP facility, ground-transit label,
    d < 600 pc.
  - A2. *Exclusions (E54):* her close-companion list **is applied for the main Method A result**, and Method A
    without it is a labelled extra line. Four of our 107 are on it: **WASP-70 A, WASP-94 A, WASP-94 B,
    WASP-173 A**. None is on her cluster or extreme-G lists. Method B: no exclusions.
  - A3. *NEA st_radv:* taken from her `nea_pscomppars.csv` (her snapshot, date not recorded), matched on the DR3
    ID. It is the first non-null value over the host's HJ rows sorted by (disc_year, pl_name), as her pandas
    `first` does. No limit or reference filter, as in her code.
  - A4. *Draws:* σ_W, σ_tot and eq. 1 are computed on the same 10,000 draws (her script draws separately per
    statistic). Only the RNG stream differs; the distribution is the same.
  - A5. *HEALPix:* no healpy or astropy-healpix in the environment. We use our own NESTED `ang2pix`
    (Górski et al. 2005 / chealpix algorithm), tested for equal-area pixels and the parent = child // 4
    hierarchy.
  - A6. *Bayestar19:* web query only (E54), not a local file. `bayestar_ebv` (file reader) stays in the code
    for tests only.
- **Validation, before any WASP run:**
  - Reproduce Hu et al. **TESS d < 600 pc** (100 hosts, σ_W 18.1 vs field 17.5) and **HS19-all** (282 hosts,
    15.8 vs 18.7) from their published host tables. Not run (Bayestar19 missing).
  - **Input (E54: option b chosen):** Hu et al. Tables 2 and 3 (name + DR3 ID) were not available locally. Options were:
    - (a) take them from the paper (arXiv:2608.12790 source or journal version);
    - (b) read only the name and DR3-ID columns of Kayla's `hosts_hs19.csv` / `hosts_tess.csv`, which her
      script pins to the published tables (E50).
  - Seed labels: `dr3x_A_val_hs19`, `dr3x_A_val_tess`.

## 2. Seeds

| Label | Tails | Use |
|---|---|---|
| `dr3x_A_val_hs19` | () | Method A validation, HS19-all |
| `dr3x_A_val_tess` | () | Method A validation, TESS d < 600 |
| `dr3x_A_wasp` | () | Method A, WASP (Stage 2) |
| `dr3x_B_wasp` | chunks 0–7, 99 | Method B, WASP control samples and bootstrap (Stage 2) |
| `dr3x_B_dr2match` | chunks 0–7, 99 | Stage 2 comparison line: frozen test B on DR2 data for the 101 Method B survivors (§3.6) |

The derivation is unchanged (D19). The registry-wide tests pass: no shared stream, and no new code shared with
any label.

## 3. Method B: our frozen DR2 primary, ported to DR3

### 3.1 Field cuts: frozen DR2 → DR3

| HS19 cut | Frozen DR2 (config `cuts`) | DR3 form | Exact? |
|---|---|---|---|
| RV required | `radial_velocity` not null (DR2) | DR3 `radial_velocity` not null | exact in form. **Open:** Katz et al. (2023) correction for faint stars (`grvs_mag`) and hot templates (`rv_template_teff`); see B4 |
| 5-parameter solution | implicit in DR2 (a parallax implies 5p) | DR3 also has **6-parameter** solutions (`astrometric_params_solved` = 95) | **not exact.** Closest: require `astrometric_params_solved = 31`; see B1 |
| 1 | parallax_over_error > 10 | same | exact |
| 2 | −0.23 < mean_varpi_factor_al < 0.36 | **no DR3 column** (dropped after DR2) | **no equivalent.** Closest: drop it (as Hu / Kayla). Scan-direction columns exist but have no published equivalent threshold; see B2 |
| 3 | visibility_periods_used > 8 | same | exact |
| 4 | UWE = √(χ²/(n−5)) < 1.2·max(1, e^{−0.2(G−19.5)}) | same columns and formula | exact for 5p solutions. For 6p solutions n−5 is not the degrees of freedom (B1 resolves this). RUWE is not used, as frozen |
| 5 | astrometric_gof_al < 3 | same | exact in form (the DR3 calibration differs) |
| 6 | astrometric_excess_noise_sig < 2 | same | exact |
| 7 | rv_nb_transits > 5 | same column name | exact in form |
| 8 | 1.0 + **0.015**·c² < E < 1.3 + 0.06·c² | `phot_bp_rp_excess_factor`, same formula with **0.015** | exact on the raw factor. Not the DR3-recommended corrected C* (Riello et al. 2021); see B3 |
| 9, 10 | BP / RP flux_over_error > 10 | same | exact |

- **Server-side WHERE:** RV not null; cuts 1, 3, 5, 6, 7, 9, 10; cut 8 with the **looser** lower coefficient
  0.0015. This mirrors the frozen DR2 practice (`cuts.server_lo_b`).
- **Applied locally:** 0.015 for cut 8, cut 4, and every open choice. All their input columns are downloaded.

### 3.2 Distances (proposal)

- **Bailer-Jones et al. (2021) EDR3 distances, `r_med_geo`** (median of the geometric posterior).
  - Tables: ARI `gaiaedr3_complements.distances` = ESA `external.gaiaedr3_distance`. EDR3 and DR3 share
    source_ids and astrometry.
  - It is the closest DR3-era analogue of the frozen BJ18 `r_est` (D16): same author, geometric (parallax-only)
    prior.
  - Differences: BJ21 gives the posterior **median**, not the mode, and its prior is a GGD fitted per HEALPix,
    not BJ18's EDSD.
- **Alternative:** `r_med_photogeo` (uses the photometry, so the distance depends on colour and magnitude).
- **Stars without a BJ21 distance are dropped and counted, as in D16.** Both columns and `flag` are downloaded,
  so the choice stays open.

### 3.3 The rest of the frozen primary (unchanged unless stated)

- **Extinction:** STILISM cube v4.1 (the frozen file), ceil(d / 1 pc) steps, zero outside the cube, CV18
  coefficients; d = BJ21 distance.
- **Main sequence:** the published ZAMS polynomial, 1 mag, **no colour limit** (D8b).
- **Kinematics:** UVW relative to the LSR (11.1, 12.24, 7.25); z Galactocentric with R0 = 8.3 kpc,
  z_sun = 27 pc (D11).
- **Matching:**
  - z histogram, 50 pc fixed bins, Bernoulli thinning;
  - one field star within 0.025 mag in (BP−RP)₀ per host, without replacement;
  - unmatched hosts dropped from both sides per iteration (D12, D13, D13b).
- **Statistic and inference:**
  - eq. 1;
  - 40,000 control samples;
  - Δ = σ_hosts − median(control), with a 95 % percentile bootstrap over 40,000 host resamples (test B of
    the frozen pre-registration);
  - seed label `dr3x_B_wasp`.
- **B5. Reading of the six cases:** exploratory, so the case is reported descriptively only, not as a test
  outcome. **Proposal; user to confirm.**
- **Field planet-host removal (D17):** the DR3 IDs of the HS19 Table 1, archive HJ and archive long-period
  host lists. **Adaptation:** the NEA `gaia_dr3_id` of each host, else the closest `dr2_neighbourhood` match of
  its DR2 ID (as E48).

### 3.4 Hosts: the same 107, frozen host rules adapted to DR3

| Frozen DR2 rule | DR3 adaptation |
|---|---|
| Gaia DR2 source_id from SIMBAD | DR3 ID from E48 (NEA route = dr2_neighbourhood route for all 107) |
| host cuts 1–4 (D3c) | cuts **1, 3, 4** (cut 2 has no DR3 column); 5-parameter restriction per B1 |
| BJ18 `r_est` | BJ21 `r_med_geo` (§3.2) |
| RV: archive `st_radv` (limits and Gaia-DR3-referenced values rejected), else Gaia DR2 RV (D6) | archive `st_radv` from our 2026-10-28 snapshot under the same D6 rule, **else Gaia DR3 RV**. **Open (B6):** whether to keep rejecting archive RVs that cite Gaia DR3, now that DR3 is the fallback anyway |
| MS cut, STILISM, UVW, z | as §3.3, with DR3 photometry and BJ21 distances |

### 3.5 Field download (Kayla's chunking approach)

- **Chunking:** 48 equal source_id ranges up to 12·4¹²·2³⁵. A range is split in half when its job times out or
  fills the row cap. Leaves are written as gzipped CSV with their ADQL; the download is resumable, never keeps a
  partial file, and has a free-disk guard (stops below 4 GB).
- **Location:** `~/data/dr3x_field_b/`, outside every repository.
- **ESA does not work anonymously.** The smoke test on 1/196,608 of the source_id range: ESA aborts anonymous
  async jobs after **60 s**. After three splits it answered **HTTP 429 (rate limit)**. Kayla's 4 h limit
  suggests she was logged in.
- **ARI Heidelberg:** has `gaiadr3.gaia_source` and the BJ21 table (4 h jobs, 10M rows); the frozen DR2 field
  also ran there (A6). The downloader now points at ARI.
- **ARI smoke test, same slice:**
  - 50 rows in 17 s, 32 columns; every star has a BJ21 `r_med_geo`;
  - 3 of the 50 have 6-parameter solutions;
  - 180 bytes per row gzipped.
  - Extrapolation from this one slice: about 10M rows, about 1.8 GB. Densities vary strongly with sky
    position, so this is an order of magnitude only.
- **Download started 2026-11-11, detached,** sequential (one job at a time). It is resumable: rerun
  `python -m src.dr3x_field_b download` and finished leaves are kept.

### 3.6 Exploratory comparison line: the same stars in DR2 (Stage 2; user request, E56)

- **What:** the frozen DR2 primary method (test B) run **on DR2 data** for exactly the **101 hosts that survive
  Method B** (E55). A DR2-vs-DR3 difference then compares the same stars, separating the change of data release
  from the change of sample.
- **Unchanged from frozen test B:** Gaia DR2, BJ18, STILISM, MS, z + colour matching, eq. 1, 40,000 control
  samples, 40,000-resample percentile bootstrap, Phase 1 field.
- **Hosts:** the 101 Method B survivors, identified by their DR2 source_ids (E48), with their frozen Phase 1
  DR2 rows.
- **Seed label:** `dr3x_B_dr2match`. **The output is labelled EXPLORATORY.**
- **Not run in Stage 1.** It belongs to Stage 2.

## 4. Outputs

- Results go under `data/processed/dr3_exploratory/` (gitignored).
- Every results file is headed **"EXPLORATORY — not pre-registered"**: JSON has a `header` key, CSV a first
  comment line.

## 5. User decisions (2026-11-11, DECISIONS E54)

1. Bayestar19: web query (dustmaps, bayestar2019, median, × 0.884) for hosts without GSP-Phot only; no
   coverage → zero; **service failure → STOP**.
2. Validation: name and DR3-ID columns only of Kayla's `hosts_hs19.csv` / `hosts_tess.csv`. Hu Tables 2/3 were
   not available locally.
3. Method A main result with her close-companion exclusions; without them as a labelled extra line. Method B:
   no exclusions.
4. Method B only: `astrometric_params_solved = 31` (field and hosts); WASP hosts removed by it are reported.
5. Cut 2 dropped (no DR3 equivalent).
6. Cut 8: raw excess factor with 0.015. The Riello et al. (2021) corrected excess is a possible later variant,
   not used.
7. No Katz et al. (2023) RV corrections (known limitation).
8. Distances: BJ21 `r_med_geo`.
9. Six cases: descriptive only, labelled "for comparison only".
10. Archive RVs citing Gaia DR3 still rejected; the Gaia DR3 RV is then used directly.

**Still open:** the Bayestar19 web service is down (E54). Method A extinction for hosts without GSP-Phot, and
the validation, wait for it, or for a new instruction from the user.

## 6. Follow-up instructions (2026-11-11, DECISIONS E56)

- **Bayestar19:** wait.
  - When the Method B field download finishes, retry the web query **once**.
  - If it is still down: check the size of the full bayestar2019 file on Harvard Dataverse and the free disk
    space, and **report before downloading anything**.
  - **The validation and Stage 2 are not run while the 6 hosts are pending.**
- **Comparison line:** §3.6 (`dr3x_B_dr2match`; Stage 2).
- **Software environment:** future package installs go in a **separate virtual environment**, not the main
  `ast192-hj` environment (after the numpy change of E54).
