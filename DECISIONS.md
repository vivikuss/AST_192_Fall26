# DECISIONS — HS19 reproduction

Every judgment call, with what the paper says, the options, the choice, and why.
Paper = Hamer & Schlaufman 2019, AJ 158, 190 (arXiv:1908.06998v2), copy in `Astro2026/Papers/Hamer_2019_AJ_158_190.pdf`.
Settings referenced here live in `config.yaml`. Phase 1 follows the paper; where
the paper is silent or self-inconsistent the choice below is the default and the
alternatives are run as sensitivity checks where they could matter.

---

## Checkpoint 1 decisions (user, 2026-10-30)

- **Phase 1 primary unchanged:** HS19 Table 1 hosts, DR2 only (archive or DR2
  RV), the published MS polynomial, the mean-residual-speed statistic.
- **The refit MS polynomial stays a sensitivity test.**
- **Pleiades cut set 1, 8, 9, 10 accepted**, with the reasons and the
  printed-set refit recorded (D10b).
- **Phase 2 follows the original one-change-at-a-time plan:** 2a DR3
  astrometry + RVs (which also brings DR3 RVs to the hosts), 2b BJ21
  distances, 2c current archive list, 2d dust map (needs approval), 2e z_max.
  No new primary sample is chosen now.
- **Next is Stage C:** the EBLM sample and pre-registration only, no analysis.

---

## Data access (2026-10-28)

**A1. STILISM cube source.** stilism.obspm.fr no longer resolves (NXDOMAIN). The
site's "About" page (Wayback, 2024-12-11) lists the downloadable cube as
`static/cube/stilism_cube.h5`, version 4.1, created 2018-03-19, 5 pc grid,
x,y ∈ [−1997.5, 1997.5] pc, z ∈ [−297.5, 297.5] pc, values in mag/pc. The Wayback
Machine holds that file (snapshot 2023-11-06, server Last-Modified 2018-03-19,
360,599,622 bytes, sha256 a768d80a…c559d906). Its internal attributes say
`cube_version 4.1`, `creation_date 2018-03-19`, from `CUBE_ebvparparsec.txt`
(2018-01-27). This is the cube STILISM distributed for Capitanio et al. (2017)
when HS19 did the work, i.e. **the same map, not a substitute**. Choice: use it.
Axis orientation verified (tests/test_extinction.py): dust at Ophiuchus
(l=353°, b=+17°), not at b=−17°; Taurus (b=−15°) not at b=+15° or with x↔y
swapped; Pleiades line of sight E(B−V)=0.042 (Gaia Collab. 2018a adopt 0.045).
**Provenance and map version (checked 2026-11-02, user request).**
- **File:** `data/raw/stilism/stilism_cube.h5`, 360,599,622 bytes, sha256
  `a768d80ad7dea541eb57fb9b240c4f112862eb4355580d1cd3be1068c559d906` (re-hashed 2026-11-02, unchanged).
- **Downloaded from:**
  `http://web.archive.org/web/20231106201424id_/https://stilism.obspm.fr/static/cube/stilism_cube.h5`
  (Wayback snapshot 2023-11-06 20:14:24 UTC; the only archived copy of the file, per the Wayback CDX index).
- **About page:** `http://web.archive.org/web/20241211053929/https://stilism.obspm.fr/about`.
- **Which map (correction).** The About page names no paper. Archived About pages from
  2018-05-05, 2018-09-04, 2019-11-19 and 2024-12-11 all describe the downloadable cube as
  **Version 4.1, created 2018-03-19**, with this grid. So v4.1 was served without change from
  at least May 2018 to Dec 2024, which covers HS19's work (received 2019-06-28).
  - The 2019-08-19 homepage says: "See Lallement et al. 2014, Capitanio et al. 2017, Lallement et al.,
    submitted (2018). The current map is based on the inversion of reddening estimates towards 71,000
    target stars."
  - Lallement et al. 2018 (A&A 616, A132; arXiv:1804.06060) describe exactly this grid ("a cube
    [−2000,2000]pc×[−2000,2000]pc×[−300,300]pc … voxels of 5pc side"). They merged new data
    "with those corresponding to the ∼30,000 targets used by Capitanio et al.", and the STILISM site
    was updated with their map.
  - **Conclusion:** the cube is the **Lallement et al. (2018) map**, i.e. Capitanio et al. (2017) data
    plus new APOGEE-DR14/Gaia data, released on STILISM as v4.1. It is not the Capitanio et al. 2017 map
    as published. HS19 cite Capitanio 2017, but the cube STILISM served while they worked was this one,
    so it is still the map HS19 queried, as far as the archive can show.
  - The earlier wording above ("the cube STILISM distributed for Capitanio et al. (2017)") is therefore
    imprecise. The config/cache label `stilism_capitanio2017` is kept as a name only (changing it would
    only invalidate caches).
  - Lallement et al. 2019 (6 × 6 × 0.8 kpc, hierarchical grid) does not match the cube. It is a
    candidate for Phase 2d (newer map; needs user approval).

**A2. CDS outage.** VizieR (CDS, CfA mirror and others), SIMBAD-CDS and
TAPVizieR timed out all day from this network. SIMBAD is reached via the Harvard
mirror (simbad.harvard.edu). The Gaia Collab. (2018a) ReadMe was recovered from
the Wayback Machine; the member table `tablea1a.dat` is **not** archived, so the
Pleiades refit (spec step 7 check, paper Fig. 1) is pending CDS. This does not
block the main pipeline, which uses the published polynomial.

**A3. HS19 Table 1.** The machine-readable Table 1 (338 hosts, ArchiveName +
SimbadName) was downloaded from the AAS/IOP content server
(`content.cld.iop.org/journals/1538-3881/158/5/190/revision1/ajab3c56t1_mrt.txt`)
into `data/hs19_table1.txt`. Per the spec it is the primary HJ list; the
archive-derived list (D5) is run alongside and the overlap reported.

**A4. Gaia archive limits.** Anonymous async jobs are aborted after 60 s, and
the archive's ADQL has no `CASE`/`GREATEST`. The field query is therefore split
into HEALPix pixels via source_id ranges (level 2, a pixel that fails is split
into its 4 children; cached per pixel in `data/raw/gaia/field_dr2_prefilter/`),
and cut 4 is applied locally. Host queries are batched 50 source_ids at a time.

**A5. Network robustness (2026-10-30).** A download stalled overnight on a dead
connection (astroquery's http.client has no timeout). Now every Gaia query has a
socket timeout (120 s) and a 300 s wall-clock limit per job. Timeouts and
dropped connections are retried 4×, and HTTP 500 / failed jobs once before the
pixel is split (a 500 is usually the 60 s cap, so retrying a too-big pixel is
wasted time) or 4× at the deepest level (7), with waits of 15/30/60/120 s.
A region that still fails is logged, skipped and listed in
`logs/field_dr2_prefilter_failed.json`. The pipeline refuses to run on an
incomplete field. After download, `python -u -m src.download_field --verify`
compares a server COUNT(*) per cached chunk with the cached rows and checks for
duplicate source_ids. None of this changes the query results.
**Truncated chunks found.** The old code fetched results without checking the
job phase. Jobs the server aborted at the 60 s cap came back with *partial*
results, and those were cached. Example COUNT(*) vs cached rows: L2 pixel 5
10,642 vs 8,193; pixel 12 14,430 vs 4,931; pixel 63 9,789 vs 9,193.
(Correction: the earlier 1,372,060 figure was *not* derived from the cache. It
was a whole-sky ESA COUNT(*) that completed in 43 s on 2026-10-28, so it is a
correct server count; an earlier version of this note wrongly called it an
undercount.) Fix: a job
counts only if its phase is COMPLETED and it returned fewer rows than the
anonymous 50,000-row cap (`max_rows`). Every cached chunk is checked by
`--verify`; mismatched chunks are deleted and downloaded again. Host and
Pleiades batches are unaffected (each returned every requested source_id).

**A6. Field query moved to ARI Heidelberg (2026-10-30, approved by the user).**
ESA anonymous limits (1 min per async job, 50,000 rows; ESA quota page) forced
~300 chunks, many of which came back truncated (A5). ARI's Gaia TAP
(`gaia.ari.uni-heidelberg.de/tap`) serves the same DR2 `gaiadr2.gaia_source`
and BJ18 distances as `gaiadr2_complements.geometric_distance`. It accepts
anonymous jobs of up to 4 h and 10M rows, so the field is one job with
identical ADQL (only the distance table name differs). Settings:
`gaia.field_service: ari`, `gaia.ari`. Hosts and Pleiades stay on ESA (small
queries, complete). DR2 is a static release, so the data should be identical;
this was checked before switching (`python -u -m src.validate_ari`,
`logs/ari_validation.json`):
- Field: L2 regions 0, 109, 110, 111 (1,568 / 1,372 / 412 / 831 rows; each
  ESA download equals its ESA COUNT(*)). ARI returns the same source_ids, and
  all 31 other selected columns (28 Gaia columns plus r_est, r_lo, r_hi,
  result_flag) are bit-identical (0 differing values).
- BJ18 distances for 10 Table 1 HJ hosts (random, seed 20191021): r_est, r_lo,
  r_hi and result_flag are identical between ESA `external.gaiadr2_geometric_distance`
  and ARI.
Completeness: one ARI COUNT(*) with the same WHERE vs the downloaded rows, and a
duplicate source_id check (`--verify`, `logs/field_dr2_prefilter_verify.json`).
A downloaded result with as many rows as the row limit (maxrec) is treated as truncated.
The ESA chunk cache is kept but no longer used (`field_service: esa` reverts).

**A7. ARI single job aborted; field split into source_id ranges (2026-10-30).**
The whole-sky ARI job (1793359593302, started 07:26:34 EDT) ran the full
4 h execution limit and was ABORTED at 11:26:34 EDT with no result. Our script
did not notice. pyvo's `job.wait(timeout=…)` long-polls with `WAIT=-1`; its
`timeout` caps each poll request, not the total wait, so the loop has no
deadline. Fixes in `src/ari_tap.py`:
- The script polls the job phase itself (plain GET every 10 s) and logs the
  status with a timestamp every 3 min.
- ABORTED or ERROR is a failure; it is not retried.
- It gives up and aborts the job after `max_wait_s` (4 h + 5 min).
- A dropped connection while polling retries the poll on the same job and
  never resubmits. A first timing attempt had resubmitted and left an orphan
  job (aborted by hand).
- Any exit, including SIGTERM or Ctrl-C, aborts the running ARI jobs.
- The result fetch is retried on the same job.
The single whole-sky query is not retried. The field query now runs one job per
HEALPix nested source_id range (`gaia.ari.range_level`; ranges contiguous in
source_id, so the server can use its index). Up to 3 run concurrently. Each
range is saved only if its job reached COMPLETED, so a restart resumes from the
completed ones. Range size is chosen by timing the densest range first (below).
Verification: COUNT(*) per range on ARI vs the saved rows, and a duplicate
source_id check (`python -u -m src.download_field --verify`).
(2026-11-02) `ari_tap.query(..., max_wait_s=...)` overrides the 4 h 5 min deadline per call; used
with 1,800 s for project 07's truncation audit. Default behaviour unchanged.
**Result:**
- **Timing:** the densest level-2 range (pixel 59, Galactic plane, 13,659 rows)
  took 31 s. A level-1 range had run for more than 7 min, so larger ranges
  scale badly. Choice: `range_level: 2` (192 ranges).
- **Download (2026-10-30):** 12:05–13:29 EDT, 77 min, 3 jobs at a time, 0
  failed ranges. Plane ranges took 3–5 min each.
- **Verify:** for all 192 ranges the ARI COUNT(*) equals the saved rows
  (1,372,060 = 1,372,060, 0 mismatched ranges); 0 duplicate source_ids
  (`logs/field_dr2_prefilter_verify.json`, per-range counts included).
- **ESA cross-check:** the independent whole-sky ESA COUNT(*) from 2026-10-28
  (A5) is also 1,372,060.

**A9. Paper values from Fig. 4 (2026-10-30).** HS19 give no dispersion values
in the text (only "less than one in 40,000 (or more than 4σ)"). The arXiv
source (`data/raw/arxiv/1908.06998_source.tar.gz`) contains Fig. 4 as a vector
PDF (`combined_hist.pdf`, matplotlib 3.0.2). `python -m src.paper_fig4` reads
the orange host line and the black 2/16/50/84/98 % lines, converting page
coordinates to km/s with the 36 and 38 km/s tick marks (27.14 pt per km/s, so
the reading error is ≲ 0.01 km/s). Results (`data/processed/paper_fig4_values.json`):
- HJ hosts: 36.11; controls 40.55 / 41.97 / 43.33 / 44.77 / 46.43 km/s.
- Long-period hosts: 43.10; controls 40.67 / 42.08 / 43.41 / 44.76 / 46.28 km/s.
An earlier session had quoted "≈ 36.1 km/s" without recording its source; this
extraction is now the reference.

**A8. Host Gaia queries on ARI (2026-10-30).** ESA host batches (50 IDs) began
to come back ABORTED at the 60 s cap (e.g. a long-period variant list, 14:12–14:36).
A changed host list also reshuffled every batch, so the cache was useless for
variants. Choice: `gaia.host_service: ari`, one ARI job per host list. Same
validated DR2 tables (A6); the rows are compared with the cached ESA batches
for the three primary lists (result below). The DR3 RV fill variant
(`gaiadr3.dr2_neighbourhood`) and the Pleiades members stay on ESA (their
cached results are complete).

---

## Samples

**D4. Mass column and what counts as "Mp sin i > 0.1 MJup".**
Paper: "minimum mass Mp sin i > 0.1 MJup" (Wright et al. 2012); column not named.
Options: `pl_bmassj` (best mass: true mass if known, else Msini), `pl_msinij`
(only planets with an RV Msini), or `pl_massj`.
Choice: `pl_bmassj` from `pscomppars`, the modern equivalent of the 2019
`exoplanets` table's `pl_bmassj` that a 2019 query would have used.
**Two filters are required to make it mean what it meant in 2019:**
(i) `pl_bmassprov ∈ {Mass, Msini, Msin(i)/sin(i)}`. Today's `pscomppars` fills
`pl_bmassj` from a **mass–radius relation** for ~3,000 planets without a measured
mass (`pl_bmassprov = "M-R relationship"`). The 2019 table did not; keeping them
would add transiting planets that never had a mass measurement.
(ii) Upper limits (`pl_bmassjlim = 1`). First choice was to drop them
(→ 468 HJ hosts). **Revised after comparing with Table 1:** 14 Table 1 hosts are
not in that list, and 10 of them qualify *only* through an upper-limit mass
(KELT-19 A, KELT-20, KELT-21, XO-6, K2-77, Kepler-24, -32, -49, -52, -63).
The other 4 (Kepler-23, -28, -54, -57) now have measured masses < 0.1 MJup; in
2019 their entries were presumably TTV upper limits. So HS19 treated
`pl_bmassj` at face value, limit or not. Faithful choice:
`drop_mass_upper_limits: false` (keep limits); `true` is the sensitivity variant.
Effect on disc_year ≤ 2019 HJ hosts, before Gaia: 527 (no filters) → 505 (i).

**D5. Approximating the 2019 archive.** Spec: `disc_year ≤ 2019`. HS19 was
received 2019-06-28, so planets announced in the second half of 2019 could not
be in it; `≤ 2018` is the stricter alternative (HJ hosts 429 vs 468). Choice:
`≤ 2019` (spec); since Table 1 exists (A3) it is the primary HJ list anyway, and
the archive list is used for the long-period sample and the overlap check.

**D6. Host RVs.** Paper: archive systemic RV if available, else Gaia DR2 RV.
`pscomppars.st_radv` today cites post-2019 papers for some hosts (e.g. Kokori et
al. 2023, Standing et al. 2026), which a 2019 query could not have returned; we
cannot reconstruct 2019 archive values, so the archive value is used as is.
A reference to a post-2019 *Gaia* release (DR3/EDR3) would mix data releases, so
such values would be replaced by the DR2 RV (`archive_rv_reject_ref_patterns`);
in today's table no `st_radv` cites DR3/EDR3 (952 cite Gaia DR2), so this rule
currently changes nothing. RV upper/lower limits (`st_radvlim ≠ 0`) are not used.

**D6b. Gaia DR2 identifiers.** Paper: SIMBAD. Choice: SIMBAD (Harvard mirror,
sim-script), query name = paper's SimbadName for Table 1 hosts, archive
hostname otherwise (with simple spacing variants). A SIMBAD record is accepted
only if the queried name is in its identifier list. If SIMBAD has no DR2 ID,
fall back to the archive's `gaia_dr2_id` (flagged `id_source = fallback`).
Results: Table 1 336 SIMBAD + 2 fallback (Pr0201, WASP-70 A); archive HJ list
448 + 16 fallback + 5 none (pulsars, SWEEPS); long-period 502 + 21 + 13 none.
SIMBAD and archive DR2 IDs agree for every host where both exist, **except
HATS-67 in HS19 Table 1**: the paper's SimbadName "2MASS J01000141-5854172"
(RA 01h00m, Dec −58.9°, DR2 4904279261014267648) is not HATS-67 (RA 12h00m,
Dec −46.1°, DR2 6144060260072337024). HS19 evidently analysed the wrong star.
Faithful choice: keep the paper's star in the Table 1 run (that is what they
analysed); the archive-list run uses the correct HATS-67.

**D6c. Table 1 hosts with no RV.** 26 of the 338 Table 1 hosts have no
`st_radv` in today's archive (neither in `pscomppars` nor in any `ps` row) and
no Gaia DR2 RV, so the paper's rule gives them no UVW. HS19 must have had an RV
for them (all 338 are in their analysis), presumably from the 2019 archive
tables, which are no longer available. Options: (a) drop them (strict paper
rule with today's archive); (b) take the SIMBAD RV; (c) take the Gaia DR3 RV.
SIMBAD's values for these stars are post-2019 (mostly Gaia DR3 itself or
2022 compilations), so (b) is (c) in disguise. Choice: (a) for the primary run,
(c) as the labelled variant `host_rv_fill_dr3` (DR3 via `gaiadr3.dr2_neighbourhood`,
used only where the paper's rule gives nothing). The primary run mixes no data
releases; the variant shows whether the missing stars matter.
**Recovery attempt (2026-10-30, time-boxed, user request).** After cuts and
the MS cut, 25 Table 1 hosts lack an RV. The Exoplanet Archive has no systemic
RV for any of them: all 439 `stellarhosts` parameter sets (every reference,
not just the composite) and all 500 `ps` rows for these hosts have an empty
`st_radv` (queries saved in `data/raw/exoarchive/*missing_rv_2026-10-30.csv`).
Nothing was recovered, so no archive-based rerun; the DR3 fill variant remains
the check. **Hosts recovered with DR3 RVs (Checkpoint 1, user request to note):** the
`host_rv_fill_dr3` variant recovers 19 of the 25. These 19 are fainter (median
G 13.5, range 8.1–15.2, vs 11.9) and farther away (median 438 pc vs 294 pc),
with similar colour ((BP−RP)₀ 0.86 vs 0.81). A sample that includes them
therefore reaches deeper; Phase 2a (which brings DR3 RVs to every star) has
to keep that in mind. The missing-RV hosts are fainter and farther away than the rest:
median G 14.0 vs 11.9, distance 530 vs 294 pc, similar colour ((BP−RP)₀
0.86 vs 0.81). They are mostly Kepler, CoRoT and HATS hosts, plus WASP-33, a
pulsating A star (`data/processed/phase1_dr2/hj_table1_missing_rv.csv`).

**D7. Hosts with both an HJ and a long-period giant** (32 hosts with disc_year ≤ 2019,
e.g. HAT-P-13, HD 217107, HIP 14810, KELT-6). Paper: silent. HS19 Table 1
includes HD 217107 and HIP 14810, so they are in the HJ sample. Choice: remove
HJ hosts from the long-period sample (`exclude_hj_hosts_from_longp: true`) so
that the two tests are independent and the "control" test is not diluted by
HJ hosts. Sensitivity: rerun with overlap allowed.
**Revised (2026-10-30, user asked for the variant closest to 367):** final
long-period counts were 242 (HJ hosts excluded), **264 (included)**, 243 (mass
upper limits dropped), 232 (disc. ≤ 2018). Choice: include them
(`exclude_hj_hosts_from_longp: false`). It is the closest to 367, and HS19
describe no such removal, so it is also the more literal reading. Removing
them is now the sensitivity variant.
Remaining gap 264 vs 367 (−103, −28 %), by step: 575 listed → 510 after cuts 1–4
→ 322 on the MS → 264 with an RV. Likely causes:
(i) the MS cut removes 188 of 510 (37 %; for HJ hosts it is 5 %). Many
RV-survey long-period hosts are subgiants and giants, so the count is very
sensitive to the MS boundary and to extinction, distance and the 2019 host
list;
(ii) 58 MS hosts have no RV today (no archive `st_radv`, no DR2 RV; median
G ≈ 14, transit/microlensing hosts). The 2019 archive may have listed RVs for
some of them;
(iii) the 2019 archive's composition (planets later retracted, re-classified
masses) cannot be reconstructed.
Not pursued further (user: "don't spend much more time on this").

**D17. Field sample excludes planet hosts.** Paper: "field stars without hot
Jupiters". Choice: remove from the field every source_id in any of the host
lists (Table 1, archive HJ, archive long-period) before matching.

---

## Quality cuts

**D1. Cut 4 (C.1) form.** Paper prints `u < 1.2*MAX(1, exp(−0.2*phot_g_mean_mag − 19.5))`.
Lindegren et al. (2018) eq. C.1 (checked in arXiv:1804.09366, Appendix C):
`u < 1.2 × max(1, exp(−0.2 (G − 19.5)))`. As printed, exp(−0.2G − 19.5) ≈ 10⁻⁹
so the cut reduces to `u < 1.2` at all G, much stricter than C.1 for bright
stars (C.1 threshold is 1.2·e^1.5 = 5.4 at G = 12). The paper says "Cuts 4 and 8
are cuts C.1 and C.2 of Lindegren et al. (2018)", so the missing parentheses are
a typesetting error. Choice: Lindegren form (`cuts.uwe.form: lindegren`).
Sensitivity: `literal` form, to see which one reproduces the paper's counts.

**D2. Cut 8 (C.2) coefficient.** Paper prints `1.0 + 0.0015*bp_rp²`; Lindegren
C.2 has `1.0 + 0.015 (G_BP − G_RP)²` (checked). Same argument → 0.015. The lower
bound only matters for red stars (at BP−RP = 2: 1.06 vs 1.006), so this has
little effect; the server prefilter uses the looser 0.0015 so either can be
applied locally.

**D3. Cut 2 column name.** Paper: `mean_varpi_factor`; the DR2 column is
`mean_varpi_factor_al`. These bounds (−0.23, 0.36) and cut 3 (> 8 visibility
periods) are not in Lindegren et al. (2018) Appendix C (checked); presumably from
Marchetti et al. (2018). Used as printed.

**D3b. "5-parameter astrometric solution".** Implemented as `parallax IS NOT
NULL` (DR2 2-parameter sources have null parallax).

**D3c. Host cuts.** Paper: "astrometric quality cuts 1-4" for hosts; cut 6 and 7
explicitly skipped; 5, 8–10 not mentioned for hosts. Choice: exactly 1–4.

---

## Distances, extinction, main sequence

**D16. Distance estimate.** BJ18 `r_est` (mode of the posterior, the catalogue's
point estimate), from `external.gaiadr2_geometric_distance`, for absolute
magnitudes, extinction, z and UVW. Stars without a BJ18 distance are dropped
(counted in the logs).

**D9. Line-of-sight integration.** Paper: "integrate the interpolated grid along
the line of sight". Choice: trilinear interpolation of the dE(B−V)/ds cube,
midpoint rule with n = ceil(d / 1 pc) equal steps (`step_pc: 1.0`; 0.5 pc
changes E(B−V) by < 0.005 mag, test_step_convergence). Points outside the cube
contribute zero (`beyond_cube: zero`): the cube is ±2 kpc in x, y and ±300 pc in
z; beyond |z| = 300 pc the dust density is negligible, and stars beyond 2 kpc in
the plane are rare in this sample (parallax_over_error > 10).

**D9b. Extinction coefficients.** Casagrande & VandenBerg (2018, MNRAS 479,
L102) Table 2 mean values (checked in arXiv:1806.01953): ⟨R_G⟩ = 2.740,
⟨R_BP⟩ = 3.374, ⟨R_RP⟩ = 2.035, with A_ζ = R_ζ E(B−V). Hence A_G = 2.740 E(B−V)
and E(BP−RP) = 1.339 E(B−V). (Their Teff-dependent fit is not used: the paper
says "mean extinction coefficients".)

**D8. Main-sequence relation and cut.** Text: polynomial fitted to the densely
sampled spline, coefficients published. Fig. 1 caption: "We use the spline fit
… to exclude all stars more than one G magnitude above the main sequence".
Choice: the published polynomial (reproducible; the spline is not published).
Cut: keep M_G0 > MS((BP−RP)_0) − 1 ("less than one magnitude above"); no lower
limit (stars below the line are kept, the paper describes none). Sign checked in
test_sign_convention and on the CMD (figures/).
**D8b. Colour range.** The polynomial is fitted over the Pleiades colour range
(about −0.1 to 3.7) and nothing in the paper limits colour. Choice: no colour
limit (`color_range: null`); the F/G/K analysis only uses 0.438–2.02 anyway, and
the Monte Carlo is colour-matched.

**D8c. Table 1 hosts that fail the MS cut (not resolved; reported).** HS19 say
all 338 Table 1 hosts are main-sequence stars; with our pipeline 16 lie more
than 1 mag above the published polynomial. 13 are marginal (0.02–0.36 mag past
the edge); CoRoT-17, -23 and -28 have STILISM E(B−V) = 0.44–0.75. Tried:
no extinction (303 pass), half extinction (fewer pass), 1/parallax distances
(319 pass), our full pipeline (322 pass). Nothing tried reproduces 338.
Choice: apply the cut as written (the paper's rule) and keep the 322; no
ad hoc rescue. The likely causes are a different distance/extinction
implementation or the spline (not the polynomial) being used for the cut.
**Sensitivity (user request):** the threshold is shifted to MS − (1 ± 0.1,
0.2, 0.3) mag (`ms_shift_*` variants). Field and HJ counts and the HJ result
for each are in REPRODUCTION_REPORT.md.

**D10. Pleiades refit (check only; the selection uses the published polynomial).**
Members: J/A+A/616/A10 tablea1a, Cluster = Pleiades, 1,326 stars (fetched from the
CfA VizieR mirror once it recovered, 14:54 EDT), optimised parallaxes `plx`,
Gaia DR2 photometry/cut columns joined by source_id, E(B−V) = 0.045 with the
CV18 coefficients. Smoothing spline = `scipy.interpolate.UnivariateSpline`
(cubic). The paper does not give the smoothing; scipy's default (s = N, i.e.
unit errors) over-smooths to a single cubic. Choice: s = N·σ², σ = robust
(MAD) scatter of the members about the main sequence (0.244 mag). Then an
8th-order polynomial is fitted to 1000 points along the spline. Result: max
|refit − published| = 0.074 mag, RMS 0.046 mag over 0.44 < (BP−RP)₀ < 2.02
(s = N·σ²/2 or 2N·σ² give 0.12 / 0.07 max). Reproduced. The refit polynomial is
also run as a main-sequence variant (`ms_refit_poly`).

**D10b. Pleiades cut set.** Paper: "We apply cuts 1, 5, and 8-10 to the Pleiades
members". Applied literally, cut 5 (`astrometric_gof_al < 3`) keeps only 159 of
1,326 members (median gof = 6.3 for these bright, nearby stars), with
(BP−RP)₀ ≥ 0.84 and no A/F stars, unlike the paper's Fig. 1, which spans
−0.1 to 3.7. Cut 5 also removes the white dwarf EGGR 25 (gof = 5.86), yet the
paper says EGGR 25 had to be removed by hand *after* the quality cuts. Both
facts mean HS19 did not effectively apply cut 5 to the Pleiades. Choice: cuts
1, 8, 9, 10 → 1,178 members, 1,177 after removing EGGR 25, colour range
−0.15 to 3.69 (matches Fig. 1).
**Why cut 5 was dropped (user accepted this at Checkpoint 1, provided it is
recorded):**
1. It removes 88 % of the members (1,326 → 159). It removes every member
   bluer than (BP−RP)₀ = 0.84, so there would be no A/F stars to define the
   F-star main sequence. HS19's Fig. 1 shows members from −0.1 to 3.7.
2. It would remove EGGR 25 automatically, but HS19 say they removed EGGR 25 by
   hand after the cuts. That only makes sense if cut 5 was not effectively
   applied.
3. Cut 5 is an astrometric goodness-of-fit cut, and bright, nearby cluster
   stars have high gof. The paper's polynomial is reproduced to 0.07 mag
   without it.
**Refit with the printed set (1, 5, 8–10), for comparison**
(`python -m src.pleiades`; `pleiades_refit_printed_cuts.json`,
`figures/phase1_dr2/fig1_pleiades_printed_cuts.png`): 159 members, (BP−RP)₀
0.84–3.64, no members at 2.1–2.8. The spline + polynomial oscillates. Over
the colours the members cover within F/G/K (0.84–2.02) it deviates from the
published polynomial by up to **1.62 mag** (RMS 0.58). Bluer than 0.84 it is
unconstrained extrapolation (up to 97 mag). With the cuts used (1, 8, 9, 10)
the deviation is 0.074 mag max (RMS 0.046) over 0.44–2.02.
(Correction: an earlier version of this entry quoted "1.19 mag" for the printed
set. That number came from a script that was not kept and cannot be
reproduced; the values above replace it.)

---

## Kinematics

**D11. Height z.** Paper: "height above or below the Milky Way midplane z",
method not given. HS19 used pyia + astropy (2019, astropy ≤ 3.2), whose
Galactocentric defaults were R0 = 8.3 kpc, z_sun = 27 pc. Choice: astropy
Galactocentric z with those pre-v4.0 defaults. A constant offset in z does not
change the matching except through bin placement; heliocentric d·sin b is the
sensitivity alternative.

**D15. UVW and the LSR.** astropy `GalacticLSR` with v_bary = (11.1, 12.24,
7.25) km/s (Schönrich, Binney & Dehnen 2010, astropy's default, also what a
2018–19 pyia/astropy user would get). U positive toward the Galactic centre, V
toward rotation, W toward the NGP. The dispersion about the sample's own mean
does not depend on the solar motion (a constant offset), so this only affects
the figures. Verified against the IDL Astronomy Library gal_uvw.pro worked
example (HD 6755) and an independent Johnson & Soderblom implementation
(tests/test_kinematics.py). RV and astrometric uncertainties are not propagated
(the paper does not).

**D14. Dispersion statistic.** Spec: σ = sqrt((1/N) Σ[(U−Ū)² + (V−V̄)² + (W−W̄)²]).
Paper eq. 1 as typeset: (1/N) Σ [(U−Ū)² + (V−V̄)² + (W−W̄)²]^{1/2}, i.e. the mean
residual speed, which is ≈ 0.92× the RMS for a Gaussian. First choice: the spec's RMS
as primary, both computed every run, and the host value compared with the
paper's Fig. 4 (≈ 36.1 km/s) to see which one HS19 used.
**Resolved (2026-10-30, test run, 297 Table 1 hosts):** the mean residual
speed (average distance from the sample mean in UVW, eq. 1 as typeset) gives
**36.18 km/s**; the RMS gives 40.50 km/s. HS19's Fig. 4 host line is at
**36.11 km/s** (read from the vector figure, A9), which matches the average
distance, not the RMS. The arXiv TeX source confirms eq. 1 is
(1/N) Σ [(U−Ū)² + (V−V̄)² + (W−W̄)²]^{1/2}, i.e. the root is inside the sum. Choice (confirmed by the user): primary
statistic = mean residual speed (`statistic: mean_abs`); the RMS is still
computed and reported alongside.

---

## Monte Carlo

**D12. z-matching (step a).** Paper: "select a subsample of the field star
sample that is matched in its z distribution to the hot Jupiter host sample";
method unspecified. Choice: histogram weighting. Bins of fixed width (50 pc)
spanning the host z range; each field star in bin k gets
w = (host fraction in bin k) / (field count in bin k); field stars outside the
host range get 0. Subsample = Bernoulli thinning with p_i = w_i / max(w)
(`method: thinning`): the largest subsample whose expected z histogram equals
the hosts'. Sensitivity: bin widths 25 / 50 / 100 pc, 5 / 10 / 20 host-quantile
bins, and `weighted_resample` (draw, with replacement and ∝ w, as many stars as
the expected thinning size).

**D13. Colour matching (step b).** For every host, one star from the z-matched
subsample with |Δ(BP−RP)_0| ≤ 0.025, drawn uniformly among candidates; a field
star is used at most once per control sample (random host order resolves
collisions). If a host has no candidate in an iteration it is left unmatched
for that iteration and the dispersion is computed on the matched stars; the
number of such iterations is reported.

**D13b. Hosts without a colour match (2026-10-30, user's instruction).** Some
hosts have few or no field stars of their colour:
- A-type HJ hosts at (BP−RP)₀ ≈ 0–0.2 have 34–210 field stars within
  0.025 mag, so the ~19,000-star z-matched subsample often has none of them.
- M-dwarf hosts at (BP−RP)₀ > 2.2 have none: no field MS star with a DR2 RV
  is that red.
Previously the control dispersion of such an iteration used the N−k matched
stars while the host dispersion used all N hosts. That is an asymmetric
comparison. Fix: in each iteration both use only the hosts that got a match.
The p-value compares each control with its own iteration's host value
(`matching.ControlSampler.run`, test `test_unmatched_host_dropped_from_both_sides`).
How often each host goes unmatched is saved (`unmatched_hosts` in summary.json).
Sensitivity: `exclude_unmatchable_hosts: true` permanently drops hosts with
fewer than `min_field_candidates` (1) field MS stars within 0.025 mag.

**D18. p-value.** Fraction of control samples with dispersion ≤ the host
sample's (spec), with the raw count. With 0 of 40,000 we report p < 2.5×10⁻⁵.

**D19. Seeds.** `seed: 20191021`; each Monte Carlo stream is
`SeedSequence([seed, sum of the label's character codes])`, so each sample's result is reproducible
independently of run order.
**Update (2026-10-30, Stage B):** the Monte Carlo now runs in chunks of
`checkpoint_every` = 5,000 iterations. Chunk c uses its own stream
`SeedSequence([seed, label code, c])`, so a resumed run gives exactly the same
numbers as an uninterrupted one. Finished chunks are saved under
`data/processed/<run>/mc_checkpoints/<label>_<hash>/`. The hash covers the
seed, the Monte Carlo and main-sequence settings, the host source_ids and the
field, so changed inputs never reuse old chunks. Results therefore differ
slightly from the Stage A test run, which used a single stream.

**D20. Spectral bins.** Strict inequalities as printed (0.438 < c < 0.742, etc.).

---

## Stage C: EBLM sample (2026-10-30). Source inventory only; nothing compiled yet

**C1. Sample counts.** EBLM IV (Triaud+2017, arXiv:1707.07521) has **118**
binaries (the paper's text and all three tables; the group's note said 188).
Davis+2026 (arXiv:2610.18338) is EBLM XIX, the EBLM-North sample: 66 binaries,
SOPHIE RVs + TESS, accepted in MNRAS.

**C2. Published γ.** A survey of 27 EBLM/BEBOP arXiv sources (tables grepped
for γ, V0, v_sys, "systemic") found no γ table in EBLM IV, V, VIII–XI, XIX
or BEBOP I. About 25 single systems have published γ (EBLM I, II, VI, XII,
XVI, XVIII, XX; BEBOP), on mixed instrument zero points. These can serve
only as checks. Gaia DR3 `nss_two_body_orbit.center_of_mass_velocity` is an
orbit-fit γ with errors for every SB1 solution. Its overlap with the EBLMs
was not measured, EBLM XIV reports extra scatter in Gaia NSS K1, and using it
would bring DR3 into a DR2 analysis.

**C3. Where the raw RVs are.**
- **EBLM IV:** not public. VizieR J/A+A/608/A129 does not exist (07
  DECISIONS, 2026-10-14). CDS was unreachable again today, so this could not
  be re-checked.
- **The user's CORALIE files** (`AST_192_Fall26/coralie07`, `coralie14`;
  private copies in 07 `data/private/rv/`) contain **all 118 EBLM IV
  systems**: 118 in CORALIE07, 69 of them also in CORALIE14. The match was by
  file-name coordinates within 30″; no RV values were read. CORALIE14 ends on
  EBLM IV's cutoff date, so these are EBLM IV's RVs. Some data may be
  unpublished, and the 07 rules (private, never committed) apply.
- **EBLM XIX:** per-epoch RVs only for the worked example J0040+01 (80 RVs in
  the arXiv source). The data statement points to the SOPHIE archive. Its
  tables give P, K and e but no T0, ω, γ or coordinates (07 had already
  noted the missing coordinates). The fixed-ephemeris γ fit therefore cannot be
  done from the paper alone. Whether the journal's supplementary material has
  RVs was not verified.
- **Inputs for EBLM IV:** `tab:obs` gives P and the primary-eclipse epoch
  T0,pri; `tab:params` gives K, e, ω, T_peri and drift terms. 21/118 have
  a drift or a second Keplerian; 4 have published outer orbits (`tab:triple`).

**C4. Zero points (for the note the plan asks for).** EBLM XIX's J0040+01 fit
gives per-instrument γ of SOPHIE 17.0357, HARPS 17.0498, CORALIE07 17.0269,
CORALIE14 17.0353 km/s: differences of ≤ 23 m/s. 07 recorded a CORALIE07→14
offset of order 10–20 m/s. Against a ~40 km/s dispersion these are negligible
(< 0.1 %).

**C5. Conflict with project 07's blinding (needs a user decision).** 07's RV
arm forbids orbit fits, residuals and trends on these CORALIE RVs before its RV
pre-registration. The γ fit here is an orbit fit on the same data.

### Stage C user decisions (2026-10-30)
**C6. Blinding (option b).** γ fits on the CORALIE RVs run only after the user
adds a dated amendment to project 07's blinding rules. The amendment permits
outputs of γ, σ_γ, K and γ − simple average, with nothing else saved or
viewed. Draft: `docs/DRAFT_07_PREREG_amendment_rv_gamma_2026-10-30.txt`.
Nothing has been run. **No per-system CORALIE07/14 offset is fitted** (it could
absorb the long-term trends that 07 is blinded to). One fixed global offset is
used instead, and nothing is output per instrument era. The literature gives
18 ± 11 m/s (Neveu-VanMalle+2016, WASP-47; sign convention not stated) and
"of order 10 m/s" (EBLM IV). **User decision (2026-11-01): 0 m/s, with ±20 m/s
as a global systematic** added in quadrature to every σ_γ (amendment draft
updated). The effect is < 0.1 % of the dispersion.
**C7. CORALIE data.** Per-system tables stay out of git (`data/private/`,
gitignored); only aggregates go into reports. The user is confirming
permission with the data owners.
**C8. EBLM XIX excluded from the primary sample.** It is the only large
northern sample (SOPHIE). Without it the EBLM sample (EBLM IV,
WASP-South/CORALIE) is **mostly southern**, while the HJ hosts cover both
hemispheres. Sky position (RA/Dec and l/b) is therefore part of the
comparison plots (C10).
**C9. Gaia DR3 two-body orbits and the ~25 published γ values: cross-checks
only**, never a source of γ.

### Stage C results without RVs (2026-10-30)
**C11. Gaia DR2 IDs of the EBLM IV primaries** (`python -m src.eblm`).
- **Method:** a positional match to EBLM IV's J2000 positions (`tab:obs`), as
  one ARI job with the target list uploaded. Every DR2 source within 10″ was
  taken, and Gaia positions were propagated back to J2000 with the DR2 proper
  motions. The nearest source is the match.
- **Result:** 118/118 matched. 114 pass all checks (separation < 3″, G − V in
  [−1.0, 0.3], no second source within 2 mag). Median separation 0.13″,
  max 1.84″.
- **The 4 flagged systems were inspected and kept:**
  - J0040+01: a second source 8″ away and 1.5 mag fainter; the match is at
    0.08″.
  - J0621-46, J1509-10, J2353-10: G − V = +0.5 to +0.9, separation ≤ 0.35″
    and a single bright source. EBLM IV computed V from models for three
    systems.
- **Independent check:** DR2 → DR3 via `gaiadr3.dr2_neighbourhood` agrees with
  project 07's positional DR3 matches for **118/118**.
- **Outputs:** match details in `data/processed/eblm/eblm4_dr2_match.csv`;
  target list in `data/eblm/targets_eblm4.csv` (`src/targets.py` format).
- **No RVs are read or used.** The Gaia RV columns are neither queried nor
  kept for binaries.

**C10. EBLM primaries vs HJ hosts** (`python -m src.eblm_compare`;
`figures/eblm/fig_eblm_vs_hj_{distributions,sky}.png`,
`data/processed/eblm/compare_stats.csv`).
- **Samples:** the EBLMs after cuts 1–4 and the MS cut (118 → 117 → 115 → 107
  → 107 → **79**; F/G/K = 54/24/1), vs the Phase 1 HJ hosts (297; F/G/K
  102/133/56).

| | EBLM median [16–84 %] | HJ median [16–84 %] | KS p |
|---|---|---|---|
| Distance (pc) | 352 [206–545] | 294 [137–634] | 0.003 |
| G (mag) | 11.5 [10.1–12.2] | 11.9 [10.0–13.5] | < 0.001 |
| (BP−RP)₀ (mag) | 0.69 [0.58–0.79] | 0.81 [0.64–1.05] | < 0.001 |
| z (pc) | −124 [−279, +157] | +26 [−163, +194] | < 0.001 |
| \|z\| (pc) | 168 [100–302] | 127 [40–266] | < 0.001 |
| Dec < 0 | 87 % | 51 % | |

- **Distance:** the EBLMs are slightly farther away and span a *narrower*
  range: no very nearby bright stars, and no faint distant ones (Kepler/CoRoT
  hosts).
- **Colour:** they are bluer and dominated by F stars.
- **Sky:** they are mostly southern, at negative Galactic latitude (b median
  −30°, so z mostly negative).
- **Magnitude:** they lie in a narrow G window (10–12.5), the WASP-South +
  CORALIE regime. The HJ list mixes ground-based and Kepler/CoRoT hosts.

**C13. [Fe/H] availability.**
- **Field:** Gaia DR2 has no metallicities.
- **EBLM IV:** tabulates no [Fe/H].
- **Gaia DR3 GSP-Spec [M/H]:** exists for 38/79 EBLM primaries, but only 1 has
  all quality flags clean. Orbital smearing of the combined RVS spectra is the
  likely cause.
- **GSP-Phot [M/H]:** 74/79, but this photometric metallicity is known to be
  unreliable.
- A homogeneous [Fe/H] match is therefore not feasible. Proposal below
  (C12, E).

**C12. PRE-REGISTRATION (revision 2; APPROVED by the user 2026-11-02; FROZEN as
`PREREG_eblm_kinematics_2026-11-02.txt`, sha256 2d7ba27a4a1193bcf65d9da5ce8b76163e2ee9b1e7b748b95bac7b0ac869e672).
The frozen file is authoritative; this entry is not edited further.**
(Text as approved:)
Written before any EBLM γ, UVW, dispersion, or colour-matched evolved fraction has
been computed or seen. After approval it is frozen as
`PREREG_eblm_kinematics_<date>.txt` with a sha256.
Revision 2 (user, 2026-11-01): one explicit primary result (B); |z| matching added
(G); the evolved-fraction exploratory analysis added (I); the CORALIE offset fixed
at 0 ± 20 m/s (A).

A. **Pipeline.** Gaia DR2, identical to the Phase 1 HJ analysis: same config,
   cuts, BJ18 distances, STILISM extinction with CV18 coefficients, published MS
   polynomial, astropy UVW, and the mean-residual-speed statistic
   σ = (1/N) Σ |v_i − v̄|.
   - **EBLM sample:** EBLM IV only (EBLM XIX excluded; C8). Cuts 1–4 as for the
     HJ hosts (per-cut counts: 118 → 117 → 115 → 107 → 107), then the MS cut
     → **79**.
   - **RV:** γ only, from the amended fit (C6): fixed P and T0 (`tab:obs`), e
     and ω (`tab:params`; e = 0 where only an upper limit is given), published
     drift and outer orbit subtracted as fixed terms.
   - **CORALIE07→14 offset:** 0 m/s, with 0.020 km/s added in quadrature to every
     σ_γ as a systematic. No Gaia RV is used.
   - **Exclusions:** systems whose fit is flagged as uncomputable are excluded
     and listed. No other exclusion is made after the fit.

B. **THE PRIMARY RESULT.** Δσ = σ_EBLM − σ_field, the difference between the
   EBLM primaries' mean-residual-speed dispersion and that of the matched field,
   with a 95 % confidence interval. Everything else in this entry is secondary,
   sensitivity or exploratory.
   - **Matched field:** HS19's Monte Carlo, 40,000 control samples. One field
     star within 0.025 mag in (BP−RP)₀ per EBLM, none reused within a sample,
     from a field subsample matched in **z (50 pc fixed bins, as in Phase 1) and
     distance (5 quantile bins of the EBLM distances)** by joint histogram
     weighting. Hosts without a match in an iteration are dropped from both sides
     (D13b).
   - **Point estimate:** Δσ̂ = σ_EBLM − median(σ_control).
   - **95 % CI:** percentile interval over 40,000 bootstrap resamples of the
     EBLM sample, Δσ_b = σ(EBLM resample b) − median(σ_control). The control
     median is a fixed reference (its Monte Carlo uncertainty is ≈ SD/200,
     negligible). The bootstrap bias, mean(Δσ_b) − Δσ̂, is reported. If
     |bias| > 0.2 × the CI half-width, a bias-corrected (BCa) interval is
     reported alongside; the percentile interval stays primary.
   - **Interpretation:** "EBLM primaries are colder than the field" if the
     whole CI is < 0; "hotter" if the whole CI is > 0; otherwise "consistent".
   - **Expected precision:** under a field-like EBLM population, the SD of σ
     for N = 79 is ≈ 2.84 km/s (D). The CI half-width is therefore
     ≈ ±5.6 km/s. The HJ hosts' offset was −6.6 km/s.

C. **Secondary results.**
   - **C-1, HS19-style:** the fraction of the 40,000 controls with σ_control ≤
     σ_EBLM (paired as in D13b), and the Gaussian-equivalent significance.
   - **C-2, EBLM vs HJ directly:** colour- and distance-matched resampling.
     In each of 40,000 iterations:
     1. Bootstrap the 79 EBLMs.
     2. For each drawn EBLM, draw one HJ host (Phase 1 primary, 297) within
        0.05 mag in (BP−RP)₀ and in the same EBLM-distance quintile, with
        replacement.
     3. Compute Δσ_HJ = σ(EBLM resample) − σ(matched HJ resample).

     Report the median and the 95 % interval. EBLMs without an HJ match in an
     iteration are dropped from both sides.

D. **Expected precision** (`data/processed/eblm/precision.json`; field
   velocities only).
   - **Control SD:** for N = 79 with z + distance matching, 2.84 km/s (z only:
     2.72), against 1.46 for the 297 HJ hosts.
   - **Primary (B):** the smallest detectable Δσ (two-sided 5 %, 80 % power) is
     **≈ 7.9 km/s**. Power for an HJ-sized offset (6.6 km/s) is **≈ 65 %**.
   - **C-2:** SD ≈ 3.0 km/s, smallest detectable difference ≈ 8.5 km/s, power
     ≈ 59 % if the EBLMs are field-like.
   - **So:** a CI that includes 0 is weak evidence of no effect. It is reported
     as a constraint, not as a null result.

E. **[Fe/H]:** not feasible homogeneously (C13). Exploratory only: done if
   spectroscopic [Fe/H] can be assembled for ≥ 50 % of the 79 from the
   literature, with the field from GALAH DR4 and the HJ hosts from `st_met`.
   Otherwise it is dropped, with the reason recorded.

F. **Measurement errors** (estimated without RVs).
   - **γ errors:** σ_γ estimated from EBLM IV's *published* summaries has median
     0.025 km/s (max 2.5); with the 0.020 km/s offset systematic, about
     0.032 km/s. The inflation of a ~42 km/s dispersion is ≈ 0.001 km/s.
   - **DR2 astrometric errors:** the 3D velocity error has median 0.40 km/s
     (EBLM) and 0.31 km/s (HJ); the inflation is 0.011 vs 0.005 km/s.
   - **Correction:** none, because these are negligible. Fixed rule: if the
     fitted median σ_γ (with the systematic) exceeds 1 km/s, also report the
     deconvolved (σ² − ⟨σ_err²⟩)^½. The raw value stays primary.

G. **Sensitivity analyses** (5,000 iterations each; report Δσ̂, its CI and the
   C-1 fraction):
   1. z-only matching (as Phase 1).
   2. z and distance with quantile bins (5 × 5).
   3. **|z| instead of signed z** (50 pc fixed bins on |z|), with distance.
      The EBLMs lie mostly below the plane (median z = −124 pc), and the
      dispersion depends on |z|, not on its sign.
   4. **|z| only** (no distance).
   5. Field restricted to Dec < 0 (the EBLMs are 87 % southern).
   6. Excluding the 11 drift or triple systems (9 + 2 among the 79).
   7. Excluding the 2 "active" systems.
   8. RMS statistic instead of mean residual speed.
   9. F stars only (54), against the field and against F-type HJ hosts.
   10. The simple average used in place of γ (method check of the γ fit).

H. **What each outcome would mean** (primary result B, read with C-2).
   - **Δσ CI includes 0, its lower bound is above −6.6 km/s (so an HJ-sized
     deficit is excluded), and C-2 shows the EBLMs hotter than the HJ hosts.** The cold kinematics are specific to HJ
     hosts and not shared by close binaries from the same surveys. This
     supports the planet-specific (tidal) interpretation and argues against a
     survey-selection origin.
   - **Δσ CI entirely < 0, and C-2 consistent with 0.** EBLM primaries are as
     young as the HJ hosts. A shared bias in transit-survey targets is more
     likely (e.g. WASP-South magnitude and colour windows), or a common age
     dependence of close companions. That would weaken the tidal-destruction
     interpretation of the HJ signal.
   - **In between or inconclusive.** Likely: power is ≈ 60 % even for an
     HJ-sized effect (D). Report the intervals and make no claim. More power
     needs Phase 2 (DR3) and more EBLMs with γ.
   - **Scope:** conclusions stay within the stated precision. Anything not
     listed here is labelled exploratory.

I. **EXPLORATORY (pre-registered now, before any colour-matched value has been
   computed): fraction of evolved primaries, EBLM vs HJ hosts at matched colour.**
   - **Question:** at the same (BP−RP)₀, is the fraction of primaries that fail
     the MS cut (M_G0 more than 1 mag above the Pleiades MS) different for
     EBLMs and HJ hosts?
   - **Disclosure:** the *raw* counts are already known from the count tables:
     EBLMs 28 of 107 fail the MS cut (26 %); archive HJ hosts 55 of 457 (12 %).
     The colour-matched comparison has not been computed.
   - **Samples:** stars passing cuts 1–4 with a distance and finite
     extinction-corrected photometry. No RV requirement for either sample,
     because the MS cut does not need one.
     - **EBLMs:** EBLM IV, 107.
     - **HJ hosts:** the **archive HJ list (disc. ≤ 2019), 457**, not HS19
       Table 1. Table 1 was selected as MS hosts by HS19, so its evolved
       fraction is near zero by construction.
     - Same colour range for both: EBLM min–max (BP−RP)₀.
   - **Method:**
     1. Bin (BP−RP)₀ in 0.1 mag bins over the EBLM colour range.
     2. Reweight the HJ hosts to the EBLM colour distribution
        (w = n_EBLM,bin / n_HJ,bin, normalised; bins with no HJ host are
        dropped from both sides and counted).
     3. Compute f_evolved for each sample and Δf = f_EBLM − f_HJ,weighted.
     4. 95 % CI from 10,000 joint bootstrap resamples of both samples.
   - **Sensitivity:** 0.05 and 0.2 mag bins; HJ hosts from the current archive
     list (no year limit); a threshold of 0.5 mag above the MS instead of 1.
   - **Interpretation (several explanations, none tested here):**
     - (i) **Tidal survival.** If HJs are destroyed as their hosts evolve and
       expand (HS19's argument), evolved hosts should be rarer among HJ hosts.
       Close M-dwarf companions are not destroyed the same way.
     - (ii) **Detection selection.** Transit depth scales as (R_p/R_*)², so a
       Jupiter around a subgiant gives a shallower, harder-to-detect transit.
       An M-dwarf companion (R ≈ 0.1–0.5 R☉) still gives a deep eclipse.
       Planet RV confirmation is also harder for evolved, active or broadened
       stars. Both effects favour evolved EBLM primaries.
     - (iii) **Colour and mass.** The EBLMs are bluer (68 % F stars). At fixed
       colour, more massive primaries live shorter on the MS, and the MS turnoff
       region crosses the F colours. Colour matching removes the first-order
       effect but not mass differences at a given colour (metallicity).
     - (iv) **The cut itself.** The MS cut uses the Pleiades polynomial. Binary
       light from the M dwarf is negligible in G, but extinction or distance
       errors move stars across the cut equally in both samples.
   - **Reported as exploratory only.** No conclusion about tidal survival is
     drawn from it alone.

### Status (2026-11-02)
- **C12 approved and frozen** (above).
- **Stage D waits on** (1) the user's sign-off of the project 07 blinding amendment
  (`docs/DRAFT_07_PREREG_amendment_rv_gamma_2026-10-30.txt`) and (2) the CORALIE data-owner permission.
  Both are being raised at the group meeting.
- **Phase 2 (DR3 update) is on hold** (user): a teammate is doing a DR2 → DR3 reproduction, and the split
  of work will be decided at the meeting.

### Storage clean-up (user decision, 2026-11-02)
Deleted superseded data, about 214 MB (project 08: 1.3 GB → 1.1 GB). Nothing used by current results was
removed.
- `data/raw/gaia/field_dr2_prefilter/` (146 MB): the old ESA field chunks, partly truncated (A5) and
  superseded by the COUNT-verified ARI field (A6, A7).
  - `src/validate_ari.py` compared ARI against 4 of these chunks (L2 regions 0, 109, 110, 111). Its result
    stays in `logs/ari_validation.json`; a rerun would re-fetch those regions from ESA.
  - `field_service: esa` would re-download everything.
- `data/raw/gaia/hosts_{hj_table1,hj_dy2019,longp_dy2019}_dr2/` (1.2 MB): ESA 50-ID host batches,
  superseded by the ARI host queries (A8; identical rows).
- `data/processed/phase1_dr2_test1000/`, `phase1_dr2_test1000_excl/` (64 MB): Stage A test runs,
  superseded by the 40,000-iteration run. Their numbers remain in the conversation summary at STOP A;
  `logs/test_mc_1000*.log` are kept.
- `data/raw/ari/audit07/` (0.8 MB): the cache of project 07's truncation audit, copied into 07's
  `audit/truncation_2026-11/` (public) and `data/private/processed/audit_truncation/` (extended sample).

**Kept on purpose:**
- `data/raw/gaia/field_dr2_prefilter_count/`: the ESA COUNT(*) results that documented the truncation (A5).
- The `*_dr3rv` batches (DR3 RV-fill sensitivity variant), `pleiades_dr2`, and `validate` (BJ18 ESA vs
  ARI check).
- The sensitivity-variant folders (≈ 440 MB of regenerable `field_ms.csv.gz` copies), at the user's
  request.

## Stage D: EBLM γ fits and kinematics (2026-11-02)

### User decisions (2026-11-02)
- **CORALIE data permission granted** by the data owners (via the user). C7 still applies: per-system tables stay in
  `data/private/` and only aggregates go into reports or git.
- **EBLM XIX (C8):** the northern EBLM XIX systems resemble the BEBOP sample, which was selected for circumbinary-planet
  searches, so they are not merged with EBLM IV. If they are added later, they form a separate sample with its own
  selection caveat.
- **07 amendment adopted** and frozen in project 07 as `PREREG_amendment_rv_gamma_2026-11-02.txt` (sha256
  `36c49d2f…2746d`), committed there as dd73827 **before any γ fit**. The text is the approved draft
  (`docs/DRAFT_07_PREREG_amendment_rv_gamma_2026-10-30.txt`) with only the status, permission and adoption lines
  changed.
- Stage D runs under the frozen `PREREG_eblm_kinematics_2026-11-02.txt`. STOP D1 after the γ fits, with no
  kinematics; STOP D2 after the pre-registered comparisons.

### Implementation (`src/gamma_fit.py`; settings in config `eblm.gamma`; tests `tests/test_gamma_fit.py`)
**D1. Inputs as the amendment states.**
- **Data:** only the columns `jdb`, `vrad` and `svrad` of the curated `<star>_coralie.rdb` files, read from 07's
  private copies (identical to the originals, which are not touched).
- **Matching:** files are matched to systems by file-name coordinates within 30″ (C3). All 118 systems have files:
  69 have CORALIE07 and CORALIE14 files, 49 have CORALIE07 only.
- **Fixed orbit inputs:**
  - P and T0,pri from `tab:obs`; e and ω from `tab:params`. For the 65 systems with only an upper limit on e, e = 0.
  - Weights are 1/(svrad² + σ_add²), with σ_add from `tab:obs` in m/s.
  - The 4 outer orbits in `tab:triple` are subtracted as fixed Keplerians.
- **What the fit returns:** `fit_gamma()` returns exactly γ, σ_γ, K, γ − mean and a flag. A test enforces this.
  The per-system table is written only to `data/private/eblm_gamma/`.
- **ω convention, checked before any fit using published values only:**
  - **Eccentric orbits:** T_peri computed from (T0,pri, e, ω), with the eclipse at ν = π/2 − ω, matches the
    `tab:params` T_peri to |Δphase| ≤ 0.0002 for all 36 orbits with e > 0.05.
  - **Circular orbits:** for all 65 upper-limit systems, (T_peri − T0)/P = 0.75 exactly, as expected for ω = 0.
- **Not done, because the amendment does not list it:** EBLM IV's bisector-based outlier rejection and removal of
  in-eclipse points. The amendment allows no other column, and the rejection rule needs the bisectors. The effect on
  γ is expected to be small: a Rossiter–McLaughlin anomaly is ≲ 0.1 km/s on a few points, against K ~ 10–40 km/s.
  It is part of what the published-K and Gaia checks (D3) test.

**D2. Drift polynomials: the reference epoch is not published (judgment call; flagged for the user at STOP D1).**
- **The problem:** EBLM IV's `tab:params` gives lin/quad/cubic coefficients (m/s/yr^j) but no reference epoch t_ref.
  For a quadratic or cubic term, t_ref changes the shape of the polynomial, not only its zero point.
- **Options:**
  - (a) t_ref = the midpoint of the system's RV span;
  - (b) t_ref = BJD 2,455,000, the zero point of the table's time columns;
  - (c) no drift subtracted, so γ averages over the drift.
- **Evidence, from published values only:**
  - Under (a), the polynomials span 0.06–10 km/s over each system's RV baseline for the 9 drift systems among the
    79. That is plausible for tertiaries.
  - Under (b), the cubic terms would reach thousands of km/s (for example J0629-67: cubic −1.1 × 10⁵ m/s/yr³, with
    data years from 2009.5).
- **Choice: (a) is primary** (config `drift_ref: mid`); (b) and (c) are computed as variants. Pre-registration G6
  (excluding the drift and triple systems) covers the residual risk.
- **Affected systems:** 9 drift and 2 triple systems among the 79, plus J0021-16 ("active", with a linear term).

**D3. Cross-checks and descriptive "poor fit" criteria (fixed before the first summary was run).**
- **Uncomputable** (fewer than 3 RVs, or a singular design) → excluded and listed (pre-registration A).
- **Descriptive flags only; no exclusion (pre-registration A):**
  - **K disagrees:** |K_fit − K_pub|/K_pub > 5 % **and** more than 5 published σ_K away.
  - **Poorly constrained:** σ_γ > 0.5 km/s.
- **Variants** (only the permitted outputs): drift (b) and (c); and the `tab:params` P and T_peri in place of the
  `tab:obs` P and T0, a check on the 4-decimal `tab:obs` period.
- **Gaia DR3** `nss_two_body_orbit`: SB1 solutions with a centre-of-mass velocity and a period within 1 % of EBLM IV's
  (the inner orbit). Reported: γ_fit − v_com, its robust SD and pulls. A zero-point difference between the CORALIE and
  Gaia RV scales is expected, so the median is reported and not removed.
- **Published γ (C2):** compiled from the cached arXiv sources (`data/processed/eblm/published_gamma.csv`, public
  values) and compared in the same way.

### STOP D1 results (2026-11-02; aggregates only, from `data/processed/eblm/gamma_d1_summary.json`)
**D4. γ fits.**
- **Computable:** 118/118, so pre-registration A excludes no system.
- **σ_γ (including the 0.020 km/s systematic), kinematics sample of 79:**

  | Statistic | Value |
  |---|---|
  | Median | 0.028 km/s |
  | 10–90 % | 0.020–0.144 km/s |
  | Max | 2.10 km/s |
  | > 0.1 km/s | 14 systems |
  | > 0.5 km/s | 2 systems |
  | > 1 km/s | 1 system |

  For all 118 the median is 0.029 km/s, and 4 systems have σ_γ > 0.5 km/s. All 4 have EBLM IV σ_add of 3–9 km/s
  (broad lines). The pre-registration F rule (median > 1 km/s) is not triggered.
- **Fitted K vs published K:** median |ΔK/K| = 0.03 % (90th percentile 0.3 %).
  - One system (in the 79) has K_fit ≈ 0.51 K_pub with both ephemerides.
  - EBLM IV's own f(m) is consistent with its published K, and the system has no Gaia orbit.
  - The amendment rules out a residual-based diagnosis. Per pre-registration A, the system is kept.
- **γ − simple average:** median |·| 1.7 km/s, max 10.3 km/s. This measures phase coverage, as expected.

**D5. Variants (γ_variant − γ_primary).**
- **`tab:params` P and T_peri:** median +0.001 km/s, range −0.07 to +0.15 km/s.
- **No drift subtracted:** |Δ| > 0.1 km/s for 4 of the 79 (max 2.3 km/s).
- **Drift t_ref = BJD 2,455,000:** up to 3.3 × 10³ km/s. This confirms D2 option (b) is not the published convention.

**D6. Cross-checks.**
- **Gaia DR3 two-body orbits:**
  - 68/118 systems have a solution; 65 are SB1 with a centre-of-mass velocity, all at the EBLM IV period (48 of
    the 79).
  - γ_fit − v_com: median −0.08 km/s (79: −0.13), robust SD 0.73 km/s; 18/65 differ by more than 1 km/s.
  - Pulls: robust SD 1.6, with 3 beyond 3σ. The Gaia NSS errors look underestimated by ~1.6×, consistent with the
    extra NSS scatter EBLM XIV reports.
  - K_fit − K_NSS: 10–90 % range −1.6 to +1.4 km/s. The scatter is Gaia-dominated.
- **Published γ:**
  - The compilation has 28 rows for 16 systems. Six rows were dropped as not comparable:
    - 4 secondary-star CCF velocities;
    - 2 ESPRESSO values derived by the compiling agent as γ + offset.
  - One name match was false. WASP-30 (Dec −10° 07′) is not EBLM IV J2353-10 (Dec −10° 53′); both its rows are
    marked EXCLUDED.
  - That leaves **15 systems, 11 of them in the 79**, with one reference value each (CORALIE first).
  - γ_fit − γ_pub: median +0.007 km/s, robust SD 0.14 km/s, range −0.67 to +0.15 km/s. None differs by more than
    1 km/s.
  - CORALIE references only (12 systems): median 0.000 km/s, robust SD 0.11 km/s.
  - Pulls: 5 beyond 3σ. Some published errors are implausibly small (10⁻⁵–10⁻⁴ km/s), and published fits with a
    drift refer to their own epochs.

### User decisions at STOP D1 (2026-11-02, before any EBLM kinematics was computed)
**D7. Drift epoch: option (a) is primary** (D2): t_ref = the midpoint of each system's RV span. Chosen at D1, before
any UVW, dispersion or Monte Carlo of the EBLMs.
- **Caveat:** for drift and triple systems, the fitted γ is the RV level at the mid-epoch, or the inner-binary
  centre of mass with the published outer orbit removed. It can differ from the system's (triple's) centre-of-mass
  velocity by up to a few km/s.
- **Coverage:** pre-registered test G6 (excluding the 9 drift and 2 triple systems among the 79) covers this.

**D8. The J1500-33 K discrepancy (D4): file check, without residuals.**
- **What was read:** only the `object`, `alpha` and `delta` columns of the keyword file
  `150057S332621_coralie_extract.rdb`. No RV, residual or other column.
- **Result:**
  - All 26 exposures are object `1SWASPJ150057.38-332620.7`, with header position 15h00m57.0 −33:26:21.
  - The object name's coordinates coincide with EBLM IV's J1500-33 position (separation 0.00″).
  - 07's header-based match (`rv_star_matches.csv`) also assigns this file to this system (4.7″ from Gaia).
  - It is the right star, so **J1500-33 stays in the primary analysis** (pre-registration A).
- **Added sensitivity test G11 (user, after D1 and before any kinematics; not in the frozen pre-registration):** the
  primary analysis without J1500-33. It is labelled as post-D1 wherever it is reported.

### D2 implementation (`src/eblm_kinematics.py`; recorded before the first run)
**D9. Field.**
- The Phase 1 field: same cuts, MS and UVW, with the Phase 1 host IDs removed as before. **The 118 EBLM IV primaries
  are also removed**, mirroring Phase 1, which removes the stars it tests from the field. The pre-registration is
  silent on this.
- Per-system EBLM tables (UVW from γ) are kept in `data/private/eblm_d2/`. Monte Carlo checkpoints there hold only
  per-iteration dispersions.
- **Fig. 3 for the EBLMs:** drawn with the ellipses only and no per-star velocity points, because the figure is
  committed and per-system values stay private (C7).

**D10. Details the pre-registration leaves open.**
- **B:**
  - The CI uses 40,000 bootstrap resamples, as stated.
  - The sensitivity tests (G) each use 5,000 Monte Carlo iterations (as stated) and 5,000 bootstrap resamples.
  - BCa is computed with jackknife acceleration.
- **C-2:**
  - Distance quintile edges come from the 79 EBLM distances; HJ hosts outside the EBLM distance range are not
    eligible.
  - Colour match: |Δ(BP−RP)₀| ≤ 0.05, with replacement.
  - HJ sample: the Phase 1 primary sample (297).
- **G2:** 5 quantile bins in z × 5 quantile bins in distance.
- **G3/G4:** |z| in 50 pc fixed bins.
- **G9:** F-type EBLMs vs the field (B), plus C-2 against F-type HJ hosts.
- **G10:** RV = the simple average (γ − (γ − mean)), through the same pipeline.
- **I:**
  - Colour bins start at the bluest EBLM colour.
  - HJ hosts outside the EBLM colour range are removed.
  - The "current archive" variant reruns the host list with `disc_year_max: null`.

### STOP D2 results (2026-11-02; `data/processed/eblm/d2/results.json`; report: EBLM_REPORT.md)
**D11. Run.**
- `python -m src.eblm_kinematics`: about 3.5 min, exit 0. Sample sizes:
  - **EBLMs:** 79 (F 54, G 24, K 1), as pre-registered.
  - **Field:** 376,786 (Phase 1's 376,791 minus 5 EBLM primaries, D9).
  - **HJ hosts:** 297.
- **Monte Carlo coverage:** no EBLM was ever unmatched, and none fell in an empty bin.
- **B, the primary result:**
  - σ_EBLM = 42.27 km/s; matched-field median 41.77 km/s (SD 2.84).
  - **Δσ = +0.50 km/s, 95 % CI [−4.23, +4.71]**: consistent with the field.
  - Bootstrap bias −0.27 km/s, below 0.2 × the half-width (0.89), so no BCa interval is needed.
- **C-1:** 22,864 of 40,000 controls are as cold or colder (0.572), a Gaussian-equivalent −0.18σ.
- **C-2 (EBLM − matched HJ hosts):** median +7.90 km/s, 95 % interval [+2.07, +13.74]. 2 EBLMs have no HJ host in
  their colour and distance cell, so about 2 are dropped per iteration.
- **H, which outcome:** the first case. The CI includes 0, its lower bound (−4.23) excludes an HJ-sized deficit
  (−6.6), and C-2 has the EBLMs hotter than the HJ hosts.
- **G sensitivity:** every variant gives a CI that includes 0.
  - Δσ ranges from −0.91 (G8, RMS statistic) to +1.23 (G6, without drift and triple systems).
  - G9, F stars only: Δσ −0.85, CI [−6.52, 4.04]; C-2 against F-type HJ hosts +6.27, CI [−0.08, 12.54].
  - G11 (post-D1, without J1500-33): Δσ −0.19, CI [−4.72, 3.86].

**D12. E ([Fe/H]): dropped, as the pre-registration requires.**
- **Count:** spectroscopic [Fe/H] was found for **27 of the 79** primaries, below the required 40 (50 %).
- **Search:** a background agent read the 27 cached EBLM/BEBOP arXiv sources
  (`data/processed/eblm/feh_feasibility.csv`). Contributions: EBLM XI (11), PAWS (6), EBLM XIX (4), EBLM V (3),
  Gill+2019 (2) and EBLM I (1).
- **Not searched:** sources that are not cached, such as Freckelton+2024 itself and GALAH DR4 for the EBLMs. The agent
  judged that they would not add the 13 systems needed. That is an assessment, not a measured count.
- Five systems have discrepant values between sources (up to 0.56 dex), which also argues against homogeneity.

**D13. F (errors).** The fitted median σ_γ is 0.028 km/s, below the 1 km/s trigger, so no deconvolved value is
reported.
