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

### After STOP D2: collaborator follow-ups (2026-11-02; exploratory, not pre-registered)
**D14. Fig. 4 redrawn.**
- Both panels (EBLM primaries; HJ hosts from Phase 1) now share one horizontal range, 30–58 km/s, with the same
  0.1 km/s bins (`figures.fig_mc(shared_x=True)`).
- The panel labels were shortened so the legend no longer overlaps the text block. The matching scheme is stated in
  the report caption instead.

**D15. Why the EBLM control distribution is wider** (`src/eblm_width.py`; EBLM_REPORT.md §11).
- **Result:** sample size. 200 random 79-host HJ subsamples give a control SD of 2.86 km/s (16–84 %: 2.80–2.93), and
  the per-star spreads of the matched pools agree (EBLM 25.0, HJ 25.1 km/s).
- **Distance matching** adds only +0.03–0.10 km/s to the width.
- **Design:**
  - (a) uses 2,000 control samples per setup;
  - (b) uses 200 subsamples × 1,000 iterations;
  - fixed seeds;
  - the D2 field (without the EBLM primaries) for the EBLMs, and the Phase 1 field for the HJ hosts.
- **Runtime:** 12 min. Each subsample rebuilds its control sampler.

## Stage E: Hu et al. (2026) challenge to HS19 (2026-11-04)

**E1. The paper.** Hu, Zhu, Huang & Zhang, "Are Hot Jupiters Tidally Disrupted During Stellar Main Sequence?",
arXiv:2608.12790 (v1, 2026-08-13; resubmitted to AJ after minor revisions), https://arxiv.org/abs/2608.12790.
Read from the arXiv HTML; no copy of the PDF is kept in this repo.
- **Claim:** HS19's cold HJ-host result comes mostly from the heterogeneity (mixed surveys, different depths and
  sky positions) and anisotropy of the HS19 sample, not from tidal destruction.
- **Method differences from HS19:**
  - Gaia DR3; extinction from GSP-Phot, else Bayestar19.
  - The vertical dispersion σ_W = sqrt((1/N) Σ (W_i − W̄)²) (their eq. 2, the standard RMS with 1/N) as the age
    indicator, instead of a total dispersion.
  - Their Table 1 matching per host: angular separation ≤ 10°; d_HJ − 100 pc ≤ d_FS ≤ d_HJ + 50 pc;
    G_HJ − 1 ≤ G_FS ≤ G_HJ + 0.1; |Δ(BP−RP)₀| ≤ 0.2. About 1,000 matches per host; one random match per host per
    iteration, with replacement; 10,000 iterations.
  - HS19-sample cuts: 0.5 ≤ (BP−RP)₀ ≤ 1.5 and, for subsamples, d < 600 pc. The stricter 0.8 ≤ (BP−RP)₀ ≤ 1.5 is used
    for the ultra-hot-Jupiter sample to remove F stars that may be evolved. Field stars with v_tot > 180 km/s are
    removed as halo stars. HJ hosts get only the astrometric cuts 1–3 of their list.
- **Results:**
  - Updated HS19 sample (282 hosts): σ_W 15.8 km/s, about 3σ below the matched field. With d < 600 pc, about 2σ.
  - Ground-based transit hosts with d < 600 pc, without HATSouth: 16.7 vs 18.1 km/s, about 1σ.
  - Per survey (d < 600 pc): WASP 102 hosts, σ_W 17.7 km/s, < 1σ from the field; HATNet 40 hosts, < 1σ;
    HATSouth 29 hosts, about 2σ low. WASP-South and SuperWASP-North are not separated.
  - TESS HJs (100 within 600 pc) and ultra-hot Jupiters (31, P < 2 d) match the field within 1σ.
  - Numerical field medians are given only for some subsamples, and no total-dispersion numbers are given for
    the subsamples.
- **What it independently confirms of ours:**
  - **HS19 eq. 1 is not a typo.** Hu et al. footnote 1: "we confirm that it is not due to typographical errors
    (Hamer & Schlaufman, private comm.)". We found that eq. 1 as typeset (mean residual speed, root inside the sum)
    reproduces HS19's Fig. 4 host value (36.18 vs 36.11 km/s; D14). Logged 2026-10-30 in the working session; first
    committed 2026-11-02 in d7dd82c (the repository's first commit).
  - **Cut 4 is missing a parenthesis.** Hu et al. footnote 7: "A parentheses is missing in Cut 4 in the appendix
    of [HS19]". We found the same against Lindegren et al. (2018) eq. C.1 (D1). Logged 2026-10-28; first committed
    2026-11-02 in d7dd82c.
  - Git dates both findings to 2026-11-02 only because the repository was created then. The 2026-10-28 and 10-30
    dates come from the working session's log.
- **How it affects this project:**
  - **Phase 1 is unaffected.** It reproduces HS19's method as published, and Hu et al. confirm that method.
  - **The HJ vs EBLM contrast (C-2, +7.90 km/s [+2.07, +13.74]) is at risk.** It compares a single-survey sample
    (EBLM IV, WASP-South/CORALIE) with the mixed-survey HS19 Table 1 hosts. Hu et al. find that WASP hosts alone
    look field-like. If so, C-2 measures survey mix, not a planet-specific effect. The EBLM-vs-field result
    (Δσ +0.50, CI [−4.23, +4.71]) is not affected.
  - **The pre-registered H interpretation of D2 (case 1: "supports the planet-specific interpretation") needs
    this caveat.** No frozen file is changed. A WASP-only test is drafted (E5).
  - Our statistic is the total dispersion (eq. 1). σ_W has not been computed for any sample in this project.

**E2. EBLM primaries vs Hu et al.'s cuts** (`python -m src.hu2026_checks eblm`; `notes/hu2026_eblm_cuts.md`).
Positions, photometry and BJ18 distances only. The recomputed unevolved set is the pre-registered 79 (same
source_ids, same distances and colours).
- d < 600 pc: 71 of the 79.
- 0.5 ≤ (BP−RP)₀ ≤ 1.5: 76.
- 0.8 ≤ (BP−RP)₀ ≤ 1.5: **12**.
- Both d < 600 pc and 0.5–1.5: 69.
- Both d < 600 pc and 0.8–1.5: 12.
- Evolved primaries (28 of the 107 after cuts 1–4) are counted in the note. Hu et al.'s main-sequence cut is the
  same as ours (1 mag above the Pleiades MS), so they would remove them too.
- Many of the 79 sit near the edges: 9 within 0.05 mag of 0.5, 17 within 0.05 of 0.8, 9 within 50 pc of 600 pc.
  DR3 distances or GSP-Phot extinction would move some of them across.
- The 0.8 cut removes the F stars, which are 68 % of the EBLMs, so an EBLM test with it has only 12 stars.

**E3. WASP-discovered hosts in our HS19 sample** (`python -m src.hu2026_checks wasp`;
`notes/hu2026_wasp_counts.md`). Counts only.
- **Label:** NASA Exoplanet Archive `disc_facility` of the host's hot Jupiter(s), queried 2026-11-04
  (`data/raw/exoarchive/pscomppars_disc_2026-11-04.csv`, public).
- **Phase 1 primary (297 hosts):** 98 have a WASP facility, 97 of them within 600 pc, 94 with d < 600 pc and
  0.5 ≤ (BP−RP)₀ ≤ 1.5. Hu et al. have 102 within 600 pc (DR3; RVs from Gaia DR3 where the archive has none).
- **Hemisphere:** the archive rarely says which. Of the 98: 15 "WASP-South"/"SuperWASP-South", 1 "SuperWASP-North",
  82 plain "SuperWASP". The 82 by declination: 16 at Dec > +15°, 24 at −15° to +15°, 42 at Dec < −15°.
  Declination does not identify the discovering instrument near the equator, where both cameras observed.
- **Name vs facility:** a name starting "WASP-" gives 100. The two extra are WASP-53 and WASP-81, whose archive
  facility is "Multiple Facilities" (Triaud et al. 2017).
- All 338 Table 1 hosts pass cuts 1–4; at that stage 105 have a WASP facility (103 within 600 pc).

**E4. Power of a WASP-only test** (`python -m src.hu2026_checks power` and `contrast`; `notes/hu2026_power.md`).
Field velocities only: no host or EBLM velocity is used, and no dispersion of a WASP subsample is computed.
- **Method:**
  - The control distribution X comes from field stars matched to the WASP hosts' z and colour (HS19-style), or
    to their position, distance, G and colour (Hu et al. Table 1); 10,000 samples.
  - "HS19-sized" means every host velocity scaled by 0.845 = 36.18/42.82 (Phase 1). The host statistic is then
    distributed as 0.845 X, because both statistics scale linearly with the velocities.
  - Assumption: σ_W shrinks by the same ratio as the total dispersion. That is −2.9 km/s on a ~19 km/s σ_W.
  - Two-sided 5 % test against X.
- **Results, 97 WASP hosts within 600 pc:**

  | Arm | Control median | Control SD | Power, HS19-sized | Offset at 80 % power |
  |---|---|---|---|---|
  | (a) total, HS19 matching | 43.0 | 2.56 | 0.80 | 6.7 km/s |
  | (b) σ_W, HS19 matching | 19.0 | 2.13 | 0.38 | 4.9 km/s |
  | (c) σ_W, Hu matching | not usable with DR2 (below) | | | |

  98 hosts (all distances) and 94 (d < 600 pc and 0.5–1.5) give the same to within 0.04 in power. A Dec < 0
  subsample (70) gives 0.68 and 0.28.
- **Arm (c) needs DR3.**
  - Our field needs a DR2 RV, so it stops at G ≈ 13. The G window G_HJ − 1 ≤ G_FS ≤ G_HJ + 0.1 then admits
    almost no field star fainter than the host.
  - 21 of 97 hosts have no match, and the rest have a median of 9. Hu et al. report about 1,000 with DR3.
  - The partial estimate (power 0.46, 76 hosts, drawn with replacement from ~10 stars) underestimates the
    control SD, so it overstates the power.
- **WASP vs EBLM contrast** (normal approximation; each sample's sampling SD taken as its field-control SD, the
  EBLMs' from z + distance matching on the D2 field: total 41.76 ± 2.79, σ_W 18.19 ± 2.29 km/s):
  - total: power 0.42 (smallest detectable difference at 80 % power, 10.6 km/s);
  - σ_W: power 0.16 (8.8 km/s).
- **Reading:**
  - A WASP-vs-field test in the total dispersion has about 80 % power for an HS19-sized deficit; in σ_W it has
    about 40 %.
  - A null in σ_W, or in the contrast, would therefore not distinguish the two papers.
  - The +7.9 km/s C-2 contrast is larger than the contrast's 80 %-power threshold for the total dispersion
    (10.6 km/s) only within its CI. A WASP-only contrast of that size would be detected with about 55 % power.
    This is a rough number from the normal approximation.

**E5. Draft pre-registration** `docs/DRAFT_PREREG_wasp_hj_vs_eblm_2026-11-04.txt` (draft for user review; not
frozen, no hash). It contains:
- both papers' predictions;
- arms (a)–(c) as specified by the user;
- the expected precision from E4;
- what counts as support for each paper and what a null means;
- open judgment calls Q1–Q10, left to the user.
- An identical copy is in `prereg_drafts/` (user request, 2026-11-04); the `docs/` file is kept.

### User decisions on the WASP-only draft (2026-11-04; no WASP-only dispersion computed)
**E6. Where pre-registrations live (convention).**
- Frozen pre-registrations are at the repository root, as `PREREG_<topic>_<date>.txt` with a `.sha256` (here
  `PREREG_eblm_kinematics_2026-11-02.txt`; project 07 does the same).
- Drafts live in `prereg_drafts/`. A draft moves to the root only when it is frozen and hashed.
- The duplicate `docs/DRAFT_PREREG_wasp_hj_vs_eblm_2026-11-04.txt` was deleted after checking it was identical
  to the `prereg_drafts/` copy (same sha256, `dd6f31ae…146e`). So E5's `docs/` path is superseded.
- The older `docs/DRAFT_07_PREREG_amendment_rv_gamma_2026-10-30.txt` (project 07's amendment, frozen in 07) was
  left where it is.

**E7. Design principle.** The primary arm differs from the frozen HS19 pipeline in exactly one way: only
WASP-discovered hosts are kept. This separates Hu et al.'s sample-heterogeneity argument from their metric
argument. Every other change goes in a sensitivity or exploratory arm.

**E8 (Q1). Gaia DR2**, as the frozen pipeline. Arm (c) (Hu et al.'s Table 1 matching) is removed: the DR2 RV field
gives too few matches (E4). It moves to Phase 2 (DR3).

**E9 (Q2). Membership:** the archive `disc_facility` (98 hosts) is primary. Sensitivity: add the hosts named WASP
without a WASP facility label (WASP-53, WASP-81; 100 hosts).

**E10 (Q3). Hemispheres:** the primary includes all WASP hosts, South and North. Exploratory: WASP-South only, with
hemisphere labels from the discovery papers, only if every host can be assigned. No declination proxy is used.
The labels have not been compiled yet.

**E11 (Q4). Distance:** no cut in the primary. Sensitivity: d < 600 pc, applied to both the WASP-vs-field test and
the WASP-vs-EBLM contrast.

**E12 (Q5). Colour:** no new cut in the primary. Sensitivity: 0.5 ≤ (BP−RP)₀ ≤ 1.5 (95 WASP hosts, 76 EBLMs). The
0.8–1.5 cut is not used (12 EBLMs).

**E13 (Q6). Halo cut** (field stars with v_tot > 180 km/s): sensitivity only. The draft notes that it applies to
field stars only and trims the field's velocity tail.

**E14 (Q7). Moot**, because arm (c) is removed. All other settings match the frozen pipeline: 40,000 control
samples, colour matching without reuse, seeds as D19.

**E15 (Q8). Primary = arm (a):** WASP hosts vs matched field in the total dispersion (eq. 1), power ≈ 0.80
(0.81 for 98 hosts). **Secondary:** the WASP-vs-EBLM contrast in the total dispersion, with power stated in
advance as ≈ 0.42.

**E16 (Q9). EBLM σ_W** is reported on its own with a 95 % CI, as exploratory. All σ_W results are exploratory,
with power stated in advance: about 0.38 for WASP vs field (0.36 for 98 hosts) and about 0.16 for the contrast.
A null in σ_W is declared uninformative.

**E17 (Q10). Hosts with close stellar companions or in clusters:** removing them is a sensitivity test only. The
draft notes that EBLMs are binaries by definition, so a companion-removal rule for HJ hosts needs justification.

### Draft revision 3 (user decisions, 2026-11-04; no WASP-only dispersion computed)
**E18. Decision rule for the primary test, and δ_HS19.**
- **δ_HS19 = −6.64 km/s:** the Phase 1 primary host dispersion (36.1825) minus its matched-field median
  (42.8234), from `data/processed/phase1_dr2/summary.json` (`mc.hj_table1`). The frozen EBLM pre-registration
  quotes it rounded to −6.6.
- **Why revision 2 said −6.7:** it applied the ratio 0.845 to the WASP hosts' expected control median (43.10,
  from E4), giving −6.68. That value depends on this test's own controls, so the fixed −6.64 is used instead.
  Power changes by < 0.01.
- **The rule:** six cases on the 95 % CI [L, U] relative to 0 and δ, exhaustive and mutually exclusive.
  1. Supports HS19.
  2. Larger cold offset than HS19's.
  3. Real but smaller deficit.
  4. Supports Hu et al.
  5. Inconclusive.
  6. Hotter than the field.
- **σ_W:** the HS19-sized deficit stays ratio-based (no frozen σ_W result exists).

**E19. Contrast matching.**
- Each half keeps its frozen matching: z only for the WASP hosts, z + distance for the EBLMs. The draft states
  this mismatch as a limitation.
- New sensitivity test D2: the WASP hosts matched on z and distance (5 quantile bins of their own distances).

**E20. Companion and cluster removal:** dropped from the sensitivity tests (supersedes E17), because Hu et al.'s
list of these hosts is not published. Moves to Phase 2.

**E21. "Without J1500-33"** is kept for the contrast. Correction to the description in the request: this test is
G11, added by the user after STOP D1 and before any EBLM kinematics (D8). It is not in the frozen EBLM
pre-registration file, and the draft says so.

**E22. Hemisphere labels from the discovery papers** (`notes/wasp_hemisphere_labels.csv` and `.md`). Metadata and
paper text only.
- **Method:** the archive's `disc_refname` bibcode was resolved through the ADS link gateway to arXiv, and the
  LaTeX source was read. For WASP-18 and WASP-24, which are not on arXiv, the publisher page was read.
- **Result:**
  - 59 WASP-South only;
  - 19 SuperWASP-North only;
  - 16 both;
  - **4 not assigned:** WASP-35 and WASP-48, whose paper describes both instruments for its systems jointly, and
    WASP-57 and WASP-127, whose camera is not stated.
- **Outcome:** not every host can be assigned, so the WASP-South-only exploratory test is removed (E10). No
  declination proxy is used.
- **Caveats:**
  - WASP-1/2 are labelled N by inference (2004 season, before WASP-South).
  - WASP-18 is labelled S from the publisher page's acknowledgements only.
  - For WASP-28 the archive's reference is a follow-up paper (Petrucci et al. 2015); the discovery paper
    (Anderson et al., arXiv:1402.1482) was used.
- Downloaded paper sources stayed in the session scratchpad. None is in the repository.

**E23. Draft revision 3:** `prereg_drafts/DRAFT_PREREG_wasp_hj_vs_eblm_2026-11-04.txt` with E18–E22. Not
frozen, no hash.

### Draft revision 4 (user decisions, 2026-11-04; no WASP-only dispersion computed)
**E24. Counts labelled as computed.** The WASP-host counts in the draft (98, 100, 97, 95) are marked as computed
by us from the Phase 1 sample (E3). The draft notes that Hu et al.'s 102 WASP hosts within 600 pc differ because
they use Gaia DR3 (distances and RVs) and their own selection: their updated HS19 list, exclusions and cuts.

**E25. σ_W offset:** the ratio rule is kept, −(1 − 0.845) × median(σ_W control). The draft adds that it differs
from the fixed −6.64 km/s for the total dispersion because no frozen σ_W result exists.

**E26. Distance bins in sensitivity D2:** each sample uses its own distance quintiles. Each is then matched to its
own field the same way: the frozen EBLM primary bins by the EBLMs' quintiles, so the WASP hosts are binned by
theirs.

**E27. Interval type:** the percentile bootstrap interval, the same as the frozen EBLM pre-registration
(`PREREG_eblm_kinematics_2026-11-02.txt`, section B: "the percentile interval stays primary"). BCa is reported
alongside only if |bias| > 0.2 × the half-width, and it does not change the decision case.

**E28. Cases 1 and 2** stay separate. A line is added that both count against the survey-mixing explanation.

### Freeze (user, 2026-11-04)
**E29. WASP-only pre-registration frozen** as `PREREG_wasp_kinematics_2026-11-04.txt`.
- **Hash:** sha256 `fece9b6608624c8213d1fb5fa96a6223bc07f5328adf820bd133db4a3490ecc4`, recorded in
  `PREREG_wasp_kinematics_2026-11-04.sha256`.
- **Format:** the 07 and 08 convention.
  - The revision 4 draft was moved from `prereg_drafts/` (`git mv`) and copied verbatim, except for its first five
    lines. The title "DRAFT pre-registration" became "Pre-registered", and the four draft status lines were
    replaced by a freeze header.
  - A diff against the committed revision 4 confirms that the body is unchanged.
- **File name:** the user's name is used. The draft's own text gives `PREREG_wasp_hj_vs_eblm_<date>.txt` (E6); that
  line is inside the frozen body and is left as written.
- **Nothing computed beforehand:** no dispersion (total or σ_W) of any WASP-only subsample, and no σ_W of any host
  or EBLM sample, was computed before freezing.
  - Stage E computed only counts and field-only control distributions.
  - The control distributions come from the power estimate (E4), which included field σ_W distributions matched to
    the WASP hosts and to the EBLMs.
- **Next:** the analysis has not been run.

### Running the frozen WASP-only analysis (2026-11-04)
**E30. Implementation details the frozen text leaves open** (`src/wasp_kinematics.py`, recorded before the first
run; none changes what the pre-registration specifies).
- **Hash check:** SHA-256 of `PREREG_wasp_kinematics_2026-11-04.txt` recomputed and matched `fece9b66…ecc4` and the
  `.sha256` file before the run. The script asserts it again.
- **Samples:**
  - **WASP hosts:** `hosts_hj_table1_final.csv` (Phase 1, 297) with the E3 facility labels (98). The name-based
    variant (D1) has 100.
  - **EBLMs:** the Stage D2 table, `data/private/eblm_d2/eblm_final_uvw.csv` (79).
  - **Fields:**
    - the Phase 1 field (`field_ms.csv.gz`, 376,791) for the WASP hosts;
    - the D2 field (376,786, D9) for the EBLMs.
- **B:**
  - `pipeline.monte_carlo` with `match_vars: [z]` and 40,000 iterations (chunked streams as D19; label
    `wasp_primary`).
  - The interval comes from `eblm_kinematics.bootstrap`, which is the frozen EBLM implementation: a 40,000-resample
    percentile interval, with bias reported and BCa only if |bias| > 0.2 × half-width.
  - The case comes from section B's rule, coded in `case()` and unit-tested for all six cases and the closed ends.
- **C:**
  - Δ_EBLM uses the frozen EBLM primary Monte Carlo (label `eblm_primary`, z + distance, 40,000; cached chunks
    reused). Its control median is checked to be 41.77.
  - The joint bootstrap draws the WASP and EBLM resamples independently: 40,000 each, with each control median
    fixed.
- **D (5,000 + 5,000):**
  - Only the half a test changes is rerun at 5,000. The unchanged half keeps its 40,000-sample reference median:
    - D1 and D2 change the WASP half only;
    - D6 changes the EBLM half only;
    - D3–D5 rerun both halves.
  - **D2:** `match_vars: [z, dist_pc]` with 5 quantile bins of the WASP hosts' own distances (config
    `match_var_bins`), as for the EBLMs.
  - **D3 and D4:** the EBLM subsets keep z + distance matching, with quintiles of their own distances.
  - **D5:** v_tot uses the field's LSR-frame UVW (same solar motion as Hu et al.). The cut is applied to both fields.
- **E:**
  - Control σ_W uses the same ControlSampler matching as A (z for the WASP hosts, z + distance for the EBLMs), with
    40,000 iterations (section A's number) on chunk streams as D19.
  - Δ is the full-sample σ_W minus the control median. The percentile bootstrap uses 40,000 resamples, as for B.
  - HS19-sized σ_W deficit: −(1 − 36.1825/42.8234) × the control median.
- **Outputs:**
  - `data/processed/wasp/results.json` (aggregates);
  - `mc_wasp_primary.npz`;
  - WASP Monte Carlo chunks in `data/processed/wasp/mc_checkpoints/`;
  - EBLM chunks in `data/private/eblm_d2/mc_checkpoints/`.

**E31. Results** (`python -m src.wasp_kinematics`, 2026-11-04, about 4 min, exit 0;
`data/processed/wasp/results.json`).
- **Deviations:** none. Every arm ran as specified, and no host was ever unmatched. The EBLM half reproduced the
  frozen control median (41.77 km/s, Δ_EBLM +0.50).
- **B, primary (98 WASP hosts, total dispersion):**
  - σ_WASP = 37.23 km/s; matched-field median 43.04 km/s (SD 2.58).
  - **Δ_W = −5.81 km/s, 95 % percentile CI [−9.69, −2.09].**
  - Bootstrap bias −0.20, below 0.2 × half-width (0.76), so no BCa interval.
  - HS19-style fraction: 285 of 40,000 controls as cold or colder (0.0071), a Gaussian-equivalent 2.25σ.
  - **Case 1** (U < 0 and L ≤ −6.64 ≤ U): "Supports HS19: a deficit consistent with HS19's size." Cases 1 and 2
    count against the survey-mixing explanation.
- **C, secondary (contrast Δ_EBLM − Δ_W):**
  - **+6.31 km/s, 95 % CI [+0.58, +12.02]** (40,000 joint resamples).
  - The whole CI is > 0, which reads as "supports an HJ-specific deficit among WASP hosts".
  - The power stated in advance was 0.42. The matching mismatch limitation in section C applies.
- **D, sensitivity** (B: Δ and CI, with case; C: contrast and CI):

  | Test | B | C |
  |---|---|---|
  | D1 by name (100) | −5.60 [−9.53, −1.82], case 1 | +6.09 [+0.25, +11.75] |
  | D2 WASP z + distance | −6.68 [−10.57, −2.92], case 1 | +7.18 [+1.27, +12.84] |
  | D3 d < 600 pc (97 / 71) | −6.23 [−9.98, −2.72], case 1 | +8.05 [+2.01, +13.54] |
  | D4 colour 0.5–1.5 (95 / 76) | −6.04 [−10.09, −2.38], case 1 | +6.88 [+0.91, +12.51] |
  | D5 halo cut (788 field stars removed) | −5.41 [−9.25, −1.69], case 1 | +6.28 [+0.55, +12.10] |
  | D6 without J1500-33 (C only) | — | +5.59 [−0.11, +11.20] |

  - Every B variant is case 1.
  - All C variants but D6 have the whole CI > 0. D6's lower bound is −0.11, so that CI includes 0.
- **E, EXPLORATORY σ_W** (a null is uninformative by pre-registration):
  - **(i) WASP vs field:** σ_W 18.31 vs 18.95 km/s; Δ −0.64, CI [−4.95, +3.30]. The CI includes 0, and the
    HS19-sized σ_W deficit would be −2.94. This null is uninformative (power ≈ 0.36).
  - **(ii) Contrast in σ_W:** −2.42, CI [−7.32, +2.53]. It includes 0, so it is uninformative.
  - **(iii) EBLMs:** σ_W 15.14 km/s (CI [12.21, 17.84]) vs matched field 18.19. Δ −3.06, CI [−5.99, −0.36]. The CI
    excludes 0, so this is reported as an exploratory finding: in σ_W alone, the EBLM primaries look colder than
    their field. That is not seen in their total dispersion (B of the frozen EBLM pre-registration).

## WASP membership join fault (2026-11-09; metadata only, no velocity or dispersion computed)

**E32. Audit of name-based joins and the fix of the WASP join.**
- **Found:** while writing `notes/cuts_and_choices.md` (898b81a).
  - `hu2026_checks.wasp_hosts` mapped the archive `disc_facility` onto the Table 1 hosts by `hostname`.
  - For Table 1 rows, `hostname` is HS19's **SimbadName** (`hosts.table1_hosts`), but the facility table is keyed by
    the **archive** host name.
- **Audit of every join by host or star name in `src/`** (row counts from the Phase 1 tables; Table 1: 338 after
  cuts 1–4, 297 final):

  | Join | Key | SIMBAD-vs-archive risk | Rows affected |
  |---|---|---|---|
  | SIMBAD IDs → host list (`hosts.py:138`) | the queried name itself | none: the cache's `hostname` is the query string (337/337 Table 1 rows identical) | 0 |
  | Table 1 → archive host columns (RV, `gaia_dr2_id`, ra/dec) (`hosts.py:194-195`) | HS19 ArchiveName vs today's archive `hostname` | yes, if the archive renamed a host | **1:** "CoRoTID 223977153" is now CoRoT-32, so it got no archive columns. Today's CoRoT-32 has no `st_radv` and no DR2 RV, so its outcome (no RV, not in the 297) is unchanged. All 297 final hosts have their archive row. |
  | Long-period minus HJ hosts (`hosts.py:66-68`, variant only) | archive `hostname` on both sides | none | 0 |
  | **WASP facility (`hu2026_checks.py`, old lines 127/145)** | SimbadName vs archive `hostname` | **yes** | **10 after cuts 1–4, 9 on the MS, 9 in the final 297** (WASP-153 fails the MS cut) |
  | WASP name test `name_wasp` (old line 146) | SimbadName | yes | the same 9 |
  | WASP labels → final table (`wasp_kinematics.py:143`) | `hostname`, both sides from the same CSV, `validate="one_to_one"` | none | 0 |
  | Hemisphere labels (`notes/wasp_hemisphere_labels.csv`, E22; built outside `src/`) | archive name, from the 98-host list | the 9 are **missing** (not mislabelled) | 9; not used by any test (WASP-South-only dropped, E22) |
  | DR3 RV fill (`hosts.py:200-238`), EBLM DR2→DR3 | `source_id` (`gaiadr3.dr2_neighbourhood`) | none | 0 |
  | Table 1 ∩ archive list (`phase1.py:100-103`) | `source_id` | none (HATS-67's wrong star is by design, D6b) | 0 |
  | EBLM IV tables, γ fits, `eblm_kinematics`, `eblm_compare`, `gamma_summary` (internal) | EBLM IV names, all parsed from the same paper; CORALIE files by coordinates (30″) | none | 0 |
  | Published γ compilation (`gamma_summary.py:119-120`) | names assigned by the compilation | a false name match is possible; not SIMBAD vs archive | 1 known (WASP-30 ≠ J2353-10), already EXCLUDED (D6) |
  | `pipeline.monte_carlo` host names | labels in the unmatched-host report only | none | 0 |

- **Fix:**
  - `wasp_hosts(cfg, join="source_id")` is now the default. Facilities are grouped by the archive's `gaia_dr2_id`
    (`facility_map`) and matched to the host's Gaia DR2 `source_id`. The name test uses HS19's ArchiveName.
  - Every one of the 297 final hosts has its source_id in `pscomppars.gaia_dr2_id`. For the 9, it equals the
    archive id.
  - `join="hostname_legacy"` reproduces the old labels (98 / 100). `wasp_kinematics.py` is pinned to it, so the
    frozen run (E31) stays reproducible as executed.
  - Tests: `tests/test_wasp_membership.py`.
- **Corrected counts** (final 297):
  - **107 WASP hosts** (98 + 9; none lost; no existing label changed). Facility labels: hemisphere not given 84,
    South 22, North 1.
  - Name variant 109.
  - d < 600 pc: 106.
  - 0.5 ≤ (BP−RP)₀ ≤ 1.5: 103 (WASP-172 has 0.46).
  - After cuts 1–4: 115 (old 105). On the MS: 109 (old 100).
- **Hosts added** (archive name / SimbadName / disc_facility / disc_year): WASP-70 A / TYC 5782-970-1 / SuperWASP /
  2014; WASP-147 / TYC 6411-312-1 / SuperWASP-South / 2018; WASP-151 / EPIC 246441449 / SuperWASP-South / 2017;
  WASP-156 / TYC 38-617-1 / SuperWASP / 2017; WASP-158 / TYC 5264-1048-1 / SuperWASP-South / 2018; WASP-161 /
  TYC 5435-1109-1 / SuperWASP-South / 2018; WASP-168 / TYC 8103-143-1 / SuperWASP-South / 2018; WASP-172 /
  2MASS J13174412-4714152 / SuperWASP-South / 2018; WASP-173 A / TYC 7518-468-1 / SuperWASP-South / 2018.
  - Why: each has a WASP `disc_facility`, but its HS19 SimbadName is a catalogue designation.
- **`disc_year ≤ 2019`:**
  - Not in the frozen pre-registration text (section A names only `disc_facility` and the 2026-11-04 query), nor
    in E3 or E30. It entered through the code: `wasp_hosts` reuses the Phase 1 archive HJ selection
    (`hosts.select_planets`, which applies `hosts.disc_year_max`).
  - It excludes **none** of the 9 (years 2014–2018). The corrected count is 107 with the year limit, without it, or
    with every planet of the host.
- **Not changed:**
  - The frozen pre-registration and its results (E31).
  - `notes/hu2026_wasp_counts.md` and `notes/hu2026_power.md` (made with the faulty list; cited by the frozen text).
  - `notes/wasp_hemisphere_labels.csv`.
  - Rerunning `python -m src.hu2026_checks wasp|power|contrast` now uses the corrected join and would overwrite
    those two notes. This was not run.

**E33. Draft amendment** `prereg_drafts/DRAFT_AMENDMENT_wasp_membership_2026-11-09.txt` (not frozen, no hash). It:
- describes the fault;
- discloses that all 98-host results have been seen;
- changes only the join (107 hosts);
- reruns every pre-registered test (B, C, D1–D6, E), reusing the membership-independent EBLM parts;
- keeps the original results as reported.

Open choices for the user: Q1 the disc_year statement (no effect); Q2 which run is the headline, and what to say if
the cases differ; Q3 new Monte Carlo labels; Q4 ArchiveName for D1; Q5 recompute field-only power for 107; Q6 Stage E
notes and hemisphere labels; Q7 scope. **Stopped here:** nothing in the amendment has been run.

### Amendment 1 finalized as a draft (user decisions, 2026-11-09; not frozen, nothing in the rerun run)
**E34. User answers Q1–Q7, written into the draft as revision 2**
(`prereg_drafts/DRAFT_AMENDMENT_wasp_membership_2026-11-09.txt`):
- **Q1:** the planet selection behind membership is stated explicitly, including `disc_year ≤ 2019`. It excludes
  none of the 107 hosts.
- **Q2:** the 107-host run is the headline, and the 98-host run is reported alongside as the original. If the two
  runs land in different decision cases, the 107 case is reported as the result, with a plain statement that the
  conclusion depends on the 9 added hosts.
- **Q3:** new Monte Carlo and bootstrap labels (`wasp_a1_*`), listed in the draft (section 5) with the D19 seed
  scheme.
  - The EBLM-half references are reused from the frozen run (`eblm_primary`, `eblm_wD3`–`wD6`, `wasp_E_E`,
    `wasp_Eiii`).
  - The EBLM bootstrap inside each contrast is redrawn under the new `_E` labels.
  - No two labels in the project share a character-code sum (checked).
- **Q4:** the name-based variant D1 uses HS19's ArchiveName: 109 hosts.
- **Q5:** field-only power recomputed for 107 hosts (E35); the numbers are in draft section 7.
- **Q6:** the Stage E notes and hemisphere labels are unchanged apart from a header, and `hu2026_checks` cannot
  overwrite them (E35).
- **Q7:** the amendment changes only the join and what follows from it (draft sections 3–9).
- Results of the rerun will go to new files (`data/processed/wasp_amendment1/`).

**E35. Field-only power for 107 hosts, and protection of the Stage E outputs.**
- **New command** `python -m src.hu2026_checks power_a1` (`power_amendment1`). Field velocities only, the E4
  method: HS19-style z + colour matching, 10,000 controls, ratio 0.845, two-sided 5 %; contrast by the normal
  approximation.
  - New labels: `a1_pow_hs19_all`, `a1_pow_hs19_d600`, `a1_contrast_eblm`.
  - Outputs: `notes/wasp_amendment1_power.md`, `data/processed/hu2026/power_amendment1.json`,
    `logs/wasp_amendment1_power.log`.
  - Run 2026-11-09, exit 0.
- **Reproduction check:** the faulty 98-host list with the original E4 label gives exactly the frozen numbers:
  power 0.81 (control SD 2.55), σ_W 0.36, contrast 0.42.
- **107 hosts:**

  | Arm | Power | Detail |
  |---|---|---|
  | Total, WASP vs field | **0.85** | control median 43.11, SD 2.46; smallest offset at 80 % power 6.25 km/s |
  | Contrast, total | **0.43** | 10.53 km/s |
  | σ_W | **0.38** | median 19.07, SD 2.07 |
  | σ_W contrast | **0.16** | – |

  - d < 600 pc (106): 0.84 / 0.42 / 0.37 / 0.16.
  - EBLM field-control SD under the new label: 2.85 (total), 2.32 (σ_W), against 2.79 / 2.29 in E4. This is Monte
    Carlo noise.
- **Protection:**
  - `hu2026_checks._write` refuses `notes/hu2026_wasp_counts.md`, `notes/hu2026_power.md`,
    `data/processed/hu2026/wasp_counts.json` and `power.json` (`PROTECTED`).
  - The `wasp`, `power` and `contrast` commands exit with a message. The default command is now `eblm` only.
  - Test: `test_protected_stage_e_outputs_cannot_be_overwritten`.
- **Headers** (dated, pointing to amendment 1, saying the faulty list was used) were added to
  `notes/hu2026_wasp_counts.md`, `notes/hu2026_power.md` and `notes/wasp_hemisphere_labels.md`. The `.csv` got a
  leading `#` comment line; read it with `comment="#"`; still 98 rows. No other content changed.
- **Stopped:** the amendment is not frozen, and none of its rerun has been run.

### Amendment 1 review items a–h (2026-11-10; draft revision 3; not frozen, nothing in the rerun run)
**E36. Review items written into `prereg_drafts/DRAFT_AMENDMENT_wasp_membership_2026-11-09.txt`.** No host or EBLM
velocity, RV, or dispersion was computed. Step 0 checks passed: 4c922a9 present, both frozen pre-registrations match
their `.sha256`, 49/49 tests. The DR3 filing (user's step 1) was **not done**: the request gave no path to the files.
- **(a) 109 vs 107** (`python -m src.amendment1_join_audit` → `notes/amendment1_join_audit.md`; reads only names, ids,
  positions, G, parallax).
  - Name set (D1, ArchiveName) minus source_id set: WASP-53 and WASP-81. Both routes give the same DR2 source_id
    (5136559128264765184, 4244243354278985472), so the separation is 0″. G 12.35 and 12.33; parallaxes 4.938 and
    2.482 mas. Their archive facility is "Multiple Facilities".
  - Reverse: none.
  - *Judgment:* classified "other" (a definitional difference, not wrong-source, missing id or duplicate). The
    107-host headline sample is therefore not affected, so the work did not stop.
  - Also reported: 98 → 107 (the 9 have the same id by both routes, 0″); 296/297 ids agree, the exception being
    HATS-67 (D6b, not WASP); no duplicated source_id.
- **(b) Seeds.**
  - New `src/seeds.py` holds a registry of every label in the project with its stream tails (96 labels: 72 earlier,
    24 `wasp_a1_*`), `stream_entropy`, and `label_code`.
  - Every `sum(map(ord, label))` in `src/` is now `seeds.label_code(label)`. It returns the same integer, so no
    stored result changes (tested).
  - *Finding:* the revision-2 statement "no two labels share a character-code sum" (E34) was **wrong**: eblm_G2 and
    w_e_zd share 632, and eblm_G5 and w_h_zd share 635. Their streams still differ (tails 0/99 vs 11), so no stream
    was shared, and nothing reported changes.
  - *Judgment:* the labels were not renamed, because that would change reported (exploratory) results. The pairs are
    listed in `seeds.KNOWN_CODE_SHARES`, and the draft's statement is corrected.
  - numpy note: SeedSequence ignores trailing zeros, so [seed, code] = [seed, code, 0]. The test compares normalized
    entropy and generated states.
  - Tests (`tests/test_seeds.py`):
    - no two labels give the same stream;
    - `wasp_a1_*` codes are unique project-wide;
    - the draft's `wasp_a1_*` list equals the registry;
    - checkpointed labels on disk are registered.
  - The draft gives the derivation in one sentence (section 5).
- **(c) EBLM resampling.**
  - *Judgment, interpretation:* "the Monte Carlo SD of the EBLM resampling" is read as the Monte Carlo precision of the
    EBLM field-control SDs that differ between the two power columns (2.79/2.29 vs 2.85/2.32).
  - The draft states them with the normal-approximation Monte Carlo SD s/√(2(N−1)): 0.020 and 0.016 km/s. The
    draws were not stored, so nothing was recomputed.
  - The differences are 2.0σ and 1.4σ. *Judgment:* "is Monte Carlo noise" is softened to these numbers.
  - The same-draw comparison (98 → 0.42, 107 → 0.43) is added, plus the fact that E4's frozen contrast row is the
    d < 600 pc set (97).
  - **Diagnostic pre-registered, not run:** the 98-host list (D1 100, D3 97, D4 95) through every WASP-half
    computation under the **same** `wasp_a1_*` labels. *Judgment:* the same labels, not new ones, so that "diagnostic
    vs 107" isolates the join (exactly for the EBLM `_E` bootstrap, approximately for the WASP half). Labelled
    DIAGNOSTIC; no decision rule applies. Output dir `data/processed/wasp_amendment1/diagnostic_98/`.
- **(d) Power.** 6.25 km/s = smallest offset with ≥ 80 % power (two-sided 5 %, field-only, ratio grid 0.001). The
  draft notes that |−5.81| is smaller. *Judgment:* power at −5.81 was not computed (not requested; post-hoc).
- **(e) Hashes and the CSV readers.**
  - No manifest or frozen file records a hash of the Stage E notes or the hemisphere CSV. They are searched by full
    hash and 8-character prefix, in the working tree and all of git history. The frozen text cites them by path.
  - *Judgment:* their SHA-256 at the freeze commit 694eb97 (before the 4c922a9 headers) is written into draft
    section 8, so the cited content is pinned.
  - No code reads the CSV (`src/`, `tests/`, `notebooks/` (empty), and repos 01–07). The only readers were ad hoc.
  - Test `test_hemisphere_csv_readers_use_comment_hash` makes any future `read_csv` of it pass `comment=`. Without
    it, pandas silently returns 99 × 1.
- **(f) Crossmatches.**
  - EBLM IV → Gaia DR2 is matched by position (10″ cone, `src/eblm.py`).
  - CORALIE files → EBLM IV systems are matched by file-name coordinates (30″, `gamma_fit.rv_files`).
  - γ → targets are matched by EBLM IV name, internal to one paper.
  - No source_id or cross-catalogue name join is involved.
  - The draft did not say explicitly that the matched field is regenerated. Section 4 now states that the WASP-half
    field controls are regenerated for the corrected hosts under the original procedure.
- **(g) Guard.**
  - No `wasp_a1_*` rerun command exists yet. *Judgment:* the guard sits in `seeds.label_code`, so **any** code path
    that starts a `wasp_a1_*` stream refuses, including future commands (`pipeline._run_chunked`,
    `wasp_kinematics.boot_draws`/`mc_sigma_w`, `eblm_kinematics.bootstrap`).
  - It requires exactly one `PREREG_wasp_kinematics_amendment1_*.txt` at the root, matching its `.sha256`, both
    committed and clean. The file must cite the frozen pre-registration's hash, and that file must still match.
  - *Judgment:* the git-committed condition is included so that the hash can't be produced together with an
    uncommitted edit.
  - Tests use a temporary git repo; the entry points refuse today.
- **(h) Privacy scan.**
  - Checked against 791 extended-sample Gaia ids, 1,025 name-like strings from 07's private extended-sample files,
    and 13,178 RV tokens from 390 CORALIE `.rdb` files. A positive control on the private files detects them.
  - Result: the draft, the new notes, src and tests have **0** matches.
  - `data/processed/wasp_amendment1/` does not exist (nothing has been run).
  - The scanner lives in the session scratchpad, not in the repo.

### Amendment 1, user decisions on E36 and further read-only items (2026-11-10; draft revision 4; not frozen)
**E37. User decisions and results.** No velocity, RV, or dispersion was computed. Nothing was frozen, pushed, or rerun.
- **DR3 filing deferred (user).** The collaborator's DR3 files are on the Tufts HPC, and the user is handling them
  separately. Step 1 of the 2026-11-10 request was not done in this repository: nothing was filed, hashed, joined or
  read.
- **User decisions on E36:**
  - (c) The normal-approximation Monte Carlo SD is accepted; no exact recomputation. The 98-host diagnostic under the
    `wasp_a1_*` labels is the actual join comparison (the draft now says so).
  - The code-sum pairs eblm_G2/w_e_zd and eblm_G5/w_h_zd stay unchanged, with the correction documented (E36).
  - The committed-in-git requirement of the guard is kept.
  - The diagnostic covers all WASP-half tests.
- **Guard extended.**
  - `seeds.require_frozen_amendment` now runs `git fetch origin victoria-barriola-08`. It then requires that the
    freeze commit is an ancestor of `origin/victoria-barriola-08`.
  - It refuses if the fetch fails.
  - *Judgment:* "freeze commit" = the last commit that touched the frozen `.txt` or its `.sha256`. The guard already
    requires both to be clean in HEAD, so this is the commit holding the content in use.
  - *Consequence:* the rerun can only run after the user has pushed the freeze commit to that branch.
  - Tests use a temporary bare remote and cover four cases: branch missing; branch present without the freeze
    commit ("not an ancestor"); freeze commit pushed (allowed); fetch failure (refused). No network is used.
- **(a) WASP-53 and WASP-81.**
  - The draft now calls them WASP-discovered hosts excluded by the frozen facility rule (archive: "Multiple
    Facilities"). It no longer calls them "other". The membership rule is unchanged.
  - `notes/amendment1_join_audit.md` uses the same wording.
- **(a) Other survey names** (archive snapshot 2026-10-28, `disc_facility` 2026-11-04; metadata only):
  - Within the 297, all 107 hosts with a WASP facility have a "WASP-" ArchiveName.
  - One has another survey's catalogue ID as its HS19 SimbadName: WASP-151 = EPIC 246441449. It is one of the 9
    added in E32.
  - No SIMBAD main identifier of the 107 is another survey's name.
  - In the whole archive, no planet with a WASP facility has a non-"WASP-" host name.
  - Three WASP-named hot-Jupiter hosts in the archive have a non-WASP facility: WASP-85 A (K2), WASP-87 (AAT) and
    WASP-167 (Multiple Observatories). None is in HS19 Table 1, so they cannot enter the sample.
- **(b) EBLM crossmatch audit** (`data/processed/eblm/eblm4_dr2_match.csv`, `gamma_fit.rv_files` file names only).
  - **Cones:** 41 of the 118 10″ cones held more than one DR2 source (31 among the 79). The match is the nearest
    source after propagating Gaia to J2000 with the DR2 proper motion (C11).
  - **Distances:** every match is ≤ 1.84″ away, and every next source ≥ 1.31″. Every next source within 5″ is
    ≥ 3.2 mag fainter.
  - **J0040+01:** its next source is the only one < 3 mag fainter (8.0″ away, ΔG = 1.5; match at 0.08″). It is the
    "confused" flag already inspected and kept in C11.
  - **Separations:** no two EBLM IV systems are within 30″ (closest pair 1,397″). Every CORALIE file matches exactly
    one system (at most two files per system: the CORALIE07 and CORALIE14 copies under one coordinate stem). No DR2
    source_id is used twice.
  - *Judgment:* nothing is ambiguous, so the work did not stop.
- **(c) Privacy.**
  - The full diff `origin/victoria-barriola-08..HEAD` (20cc253..8d663b0, local remote-tracking ref, not fetched) has
    0 extended-sample ids, 0 extended-sample name strings, 0 CORALIE RV values, and no `data/` or private paths.
  - **Five untracked logs** (`eblm4_dr3_nss`, `eblm_d2`, `eblm_width`, `hu2026_power`, `wasp_kinematics`): 0 ids, 0
    names.
    - One numeric match, `0.326` in `hu2026_power.log`, is a power value that coincides with an RV value. It is not
      CORALIE data.
    - Read by eye: aggregates only (counts, Δ, CIs, control widths). The one system name is EBLM IV's public
      J1500-33.
    - The five are added to `.gitignore` by name and stay local.

### Amendment 1, draft revision 5 (user edits, 2026-11-10; not frozen, nothing run)
**E38. Revision 5 of `prereg_drafts/DRAFT_AMENDMENT_wasp_membership_2026-11-09.txt`**, exactly the user's eight edits.
No velocity, RV, or dispersion was computed.
1. **Section 4, first line:** B, C, D1–D6 and E(i)–(ii) are rerun; E(iii) is reported from the frozen run.
2. **Section 5, seed sentence completed:** the tail is c, 99, or a test-specific tail (2 and 3 in eblm_kinematics,
   11 and 12 in eblm_width, none in hu2026_checks and eblm_compare); `src/seeds.py` is the complete registry.
   - *Judgment:* `eblm_compare` has two streams without a label: SeedSequence([20191021, 7]) and
     `astrometric_error`'s default_rng(1), which both callers use. They are added to `seeds.UNLABELLED`, so the
     statement "complete registry" is true. The collision test now includes them and passes. No result changes.
3. **Section 7:** the frozen contrast cell 0.42 is footnoted as E4's d < 600 pc row (97 hosts).
4. **Section 2, publication:**
   - The push date comes from the local reflog of `origin/victoria-barriola-08`: 20cc253, pushed 2026-11-04 14:05
     −0500. That commit carries E31's numbers and the commit message; `results.json` is gitignored.
   - *Judgment:* this is stated precisely as what was pushed.
   - **The Slack date is left as the placeholder `[DATE]`**, because it is not known here. It must be filled
     before the freeze.
5. **Section 3, mass condition:** the exact rule, from `hosts.select_planets` and `config.yaml`:
   - pl_orbper < 10 (non-null);
   - pl_bmassj > 0.1;
   - pl_bmassprov ∈ {"Mass", "Msini", "Msin(i)/sin(i)"};
   - pl_bmassjlim not used (`drop_mass_upper_limits: false`);
   - disc_year ≤ 2019.
6. **Section 4, contrast references:** a table of the WASP and EBLM reference for C, D1–D6 C and E(ii). It is taken
   from `wasp_kinematics.py` as frozen (E30). D1 C and D2 C use eblm_primary (40,000). D3–D6 use
   eblm_wD3–eblm_wD6 (5,000). E(ii) uses wasp_E_E.
7. **Section 4, sample sizes:** "5,000 control samples and 5,000 bootstrap resamples for each D test" (40,000 each
   for B, C, E).
8. **Section 8:**
   - The repository scope now reads "repositories cloned on this machine (AST_192_Fall26 01–07)".
   - The nine hosts are absent from the 98-row `wasp_hemisphere_labels.csv`.
   - Section 3's hemisphere counts (84/22/1) come from the archive `disc_facility`.
- **Noted, not changed (outside the eight edits):** the section 7 "Check" bullet says the 98-host check reproduces
  "contrast 0.42". The check row is the all-distance 98 list with the a1 EBLM draws (0.4165), while the frozen 0.42
  is E4's d < 600 pc row (0.4212). They agree only after rounding.
- **Final revision-5 changes (user, 2026-11-10):**
  - The Slack date placeholder now reads `[YYYY-MM-DD]`.
  - The section 7 "Check" line now states only what was reproduced. The check reproduced, exactly, the E4 power for
    the all-distance 98 list under the E4 label `hu_pow_hs19_WASP, all distances`: B 0.8135 and σ_W 0.3558, with
    identical control summaries. The E4 d < 600 pc contrast row (0.4212) was **not** recomputed with the E4 label,
    so it is not claimed. The rounded "contrast 0.42" is removed.
  - Section 4 now reads "the EBLM reference of every contrast (table below)".
  - The `seeds.UNLABELLED` addition is kept (user).
  - Revision 5 is approved by the user for freezing.

### Freeze of amendment 1 (user, 2026-11-10)
**E39. Amendment 1 frozen** as `PREREG_wasp_kinematics_amendment1_2026-11-10.txt`.
- **Hash:** sha256 `533a6437d3bd32348a71bf7436647113ffe761c5164c59076053cde456be4e21`, recorded in
  `PREREG_wasp_kinematics_amendment1_2026-11-10.sha256`.
- **Format (E6, as in E29):** revision 5 (commit bc69542) was moved from `prereg_drafts/` with `git mv`.
  - Only its first five lines changed. The title "DRAFT amendment 1" became "Amendment 1", and the status lines
    became "FROZEN 2026-11-10", listing revisions 2–5 and the hash file.
  - A diff against bc69542 shows only these lines and the rename.
- **Verification:** both earlier frozen pre-registrations still match their hashes. The frozen amendment cites the
  WASP pre-registration's hash.
- **Frozen as written:**
  - Section 2 contains the placeholder `[YYYY-MM-DD]` for the Slack date, as the user set it.
  - Section 2 of the frozen file is not edited. The actual date, if recorded, goes in DECISIONS.md or a later
    amendment.
- **Test:** `test_amendment1_label_list_matches_the_draft` now reads the frozen file (the draft path no longer
  exists).
- **Stale paths, not changed:** the headers of the four Stage E files (E35) still point to the old draft path
  `prereg_drafts/DRAFT_AMENDMENT_wasp_membership_2026-11-09.txt`. They were left as written, because the amendment
  describes them as unchanged apart from that header.
- **Not run:** the rerun. The guard requires the freeze commit to be on `origin/victoria-barriola-08`, so the user
  must push first.
- **Test adjusted:** `test_wasp_a1_streams_refuse_without_a_valid_frozen_amendment` checked the real repository with "frozen ⇒ allowed". It now checks the real repository only before a freeze, so tests never fetch from GitHub; the temporary-repository cases are unchanged.

### Rerun under amendment 1 (2026-11-10)
**E40. Pre-run checks and implementation** (`src/wasp_amendment1.py`, recorded and committed before the run).
- **Pre-run checks:**
  - The user pushed; `origin/victoria-barriola-08` = 33fc3bd (the freeze commit).
  - `src/seeds.py` is unchanged since 809f764 apart from `UNLABELLED` (E38), and the guard logic is identical. All 8
    seed tests pass, including the four temporary-remote cases.
  - `seeds.require_frozen_amendment()` passes on this repository after `git fetch`. It was not bypassed; the script
    calls it first, and every `wasp_a1_*` stream calls it again.
- **Implementation:** the frozen functions of `src/wasp_kinematics.py` (`case`, `contrast`, `boot_draws`,
  `mc_sigma_w`, `sigma_w`) and `eblm_kinematics.bootstrap` / `mc` are used unchanged.
  - *Judgment:* `vs_field` is repeated with one difference: its checkpoint directory is
    `data/processed/wasp_amendment1/mc_checkpoints/`, so frozen outputs are untouched.
  - Samples come from `hu2026_checks.wasp_hosts`:
    - `join="source_id"` for the headline (107, D1 109);
    - `join="hostname_legacy"` for the diagnostic (98, D1 100).
    - *Judgment:* both are merged onto the Phase 1 table on source_id (one_to_one).
  - Fields, EBLM table, matching, settings and labels are as in amendment sections 4–5.
- **Reused from the frozen run (section 4):**
  - `eblm_primary` and `eblm_wD3`–`wD6` control medians come from their cached frozen chunks. The script asserts
    that no new EBLM checkpoint directory appears, so nothing is recomputed, and that eblm_primary = 41.7745 as in
    `results.json`.
  - The `wasp_E_E` σ_W control median and E(iii) are read from `data/processed/wasp/results.json`.
- **Diagnostic:** the same code and labels, with the decision case, C reading and interpretation fields removed
  from its output and logs.
- **Outputs (aggregates):**
  - `data/processed/wasp_amendment1/results.json`;
  - `diagnostic_98/results.json`;
  - `logs/wasp_amendment1.log`.

**E41. Results of the amendment 1 rerun** (`python -m src.wasp_amendment1`, 2026-11-10, about 7.5 min, exit 0; the guard
passed; commit 648e393; log `logs/wasp_amendment1.log`; table in `notes/wasp_amendment1_results.md`).
- **Headline, 107 hosts:**
  - B: Δ_W = −5.65 km/s [−9.37, −2.04], **case 1**.
  - C: +6.15 [+0.40, +11.60], whole CI > 0.
  - D1–D5: B variants all case 1, and every C has the whole CI > 0.
  - D6 C: +5.43 [−0.05, +10.85], CI includes 0.
  - E(i): −1.18 [−5.21, +2.59]; E(ii): −1.88 [−6.58, +2.98]. Both include 0, and both are exploratory.
  - E(iii): reused from the frozen run, −3.06 [−5.99, −0.36].
  - BCa was not triggered.
- **Frozen 98-host run (E31), reported alongside:** B −5.81 [−9.69, −2.09], case 1; C +6.31 [+0.58, +12.02].
- **Diagnostic (98 hosts, wasp_a1_* labels; no reading applies):** B −5.78 [−9.66, −2.03]; C +6.27 [+0.57, +11.87].
  - Stream change: frozen → diagnostic = +0.03 km/s in B.
  - Join: diagnostic → headline = +0.13 km/s in B.
- **Section 6 rule:** every decision case and C reading of the 107-host run equals the frozen run's. The conclusion
  does not depend on the 9 added hosts.
- **Reuse checks:** the EBLM references were reused from cached frozen chunks. No new EBLM checkpoint directory
  appeared, and eblm_primary = 41.7745 was asserted.
- **Unchanged:** the frozen outputs in `data/processed/wasp/`.

**E42. Checks before the Slack correction (user request, 2026-11-10; nothing recomputed).**
- **Copied function: passed.** `wasp_amendment1.vs_field` was compared with the frozen `wasp_kinematics.mc_total` +
  `vs_field`. The calls are the same, in the same order. The only differences:
  - `mc_total` is inlined;
  - `CKPT` points to `data/processed/wasp_amendment1/mc_checkpoints`;
  - `case` and `CASES` are referenced as `wk.case` / `wk.CASES`;
  - the unused raw `res` is not returned.
- **Since 20cc253, `wasp_kinematics.py` itself changed only in two places:** the seed code (`seeds.label_code`,
  same integer, 8d663b0), and the pin of `main()` to the legacy join (c521506), which the rerun does not call.
- **EBLM references: passed.**
  - `eblm_primary` = 41.7745 (8 cached chunks, 2026-11-02 21:03:20–21:03:54). It matches the written record:
    E30/E31 41.77, the frozen `results.json` 41.77448979222274, and EBLM_REPORT 41.77.
  - The `wasp_E_E` σ_W reference, 18.1950 from the frozen `results.json`, matches E31's 18.19.
  - The 20cc253 commit message contains no EBLM reference values.
- **D3–D6 references have no committed record:** `eblm_wD3` 41.6571, `eblm_wD4` 41.9615, `eblm_wD5` 41.3987,
  `eblm_wD6` 41.8675. They were verified in three ways:
  - the cache file times (2026-11-04 13:56:44, 13:56:57, 13:57:10, 13:57:14) match the frozen run's log
    (`logs/wasp_kinematics.log`, lines 16, 18, 20, 21);
  - there is one cache folder per label, and no EBLM cache file is newer than the freeze;
  - the EBLM term of each contrast (EBLM dispersion − reference = contrast + Δ_W) is identical to the last digit
    in the frozen, headline and diagnostic runs, for C, D1–D6 C and E(ii).

**E43. Housekeeping: Stage E headers now point to the frozen amendment** (user, 2026-11-10; resolves the stale paths
noted in E39).
- In the header line of `notes/hu2026_wasp_counts.md`, `notes/hu2026_power.md`, `notes/wasp_hemisphere_labels.md`
  and `notes/wasp_hemisphere_labels.csv`, only the path changed:
  `prereg_drafts/DRAFT_AMENDMENT_wasp_membership_2026-11-09.txt` → `PREREG_wasp_kinematics_amendment1_2026-11-10.txt`.
- The pre-header SHA-256 values in amendment section 8 still verify, with the header removed (2 lines for the
  `.md` files, 1 for the `.csv`).
- Tests 58/58; all three frozen files match their hashes.

**E44. Figure for the Slack correction** (presentation only; no new statistics). `python -m src.fig_wasp_amendment1` →
`notes/figures/fig_wasp_amendment1_forest.png`.
- **Inputs:** values are read from `data/processed/wasp_amendment1/results.json` (107) and the frozen
  `data/processed/wasp/results.json` (98). Nothing is retyped. The script asserts that both runs cite the same
  frozen pre-registration hash and δ_HS19 = −6.64.
- **Panels:**
  - left: B, D1–D5 B and E(i) σ_W, with lines at 0 and δ_HS19 (dashed);
  - right: C, D1–D6 C and E(ii) σ_W, with a line at 0.
- **Markers:** 107-host values as filled blue markers, frozen 98-host values as open gray markers, each with its 95 %
  percentile interval. The diagnostic is not shown.
- *Judgment:* the σ_W rows are labelled "(exploratory)" in gray italics and set off by a dotted line. Each row
  label gives N (107 run vs 98 run); contrasts give N_WASP/N_EBLM.
- *Judgment:* the style follows the project palette (`src/figures.py`). The reference-line labels sit above the
  panels so that the lines don't cross them (checked by viewing the PNG).
- **Privacy scan:** the script and the PNG have 0 matches. The PNG metadata holds only Matplotlib's Software tag
  and the dpi. The figure shows only aggregates already in E31/E41.

## DR3 data from a collaborator (filed 2026-11-10; no analysis)
**E45. Provenance and inventory of Kayla's Gaia DR3 files.** Filing only. No WASP host was selected, nothing was
crossmatched, and no velocity, dispersion or other kinematic quantity was computed or summarized.
- **Sender:** Kayla. Originally received **2026-11-06 on the Tufts cluster**.
- **This copy:**
  - Location: `/Users/victoriabarriola/Documents/Astro2026/AST_192_Fall26/hu2026/`.
  - Contents: 271 files, 3.6 GB. File times run from 2026-11-04 10:15 to 14:18 (they predate the 2026-11-06
    receipt, so they are presumably the original creation times).
  - **No code is included**, and no CLAUDE.md or other text file: only `.csv`, `.adql` and one `.key`. All files
    were treated as data.
- **Location / git:**
  - The folder is outside the 08 repository.
  - **It is visible to the parent repository `AST_192_Fall26`** (remote github.com/vivikuss/AST_192_Fall26), where
    it shows as untracked (`?? hu2026/`), not ignored. Nothing is tracked.
  - Nothing was changed there; the fix is the user's decision.
  - Permissions: directory 755, files 644 (readable by other local accounts).
- **Checksums:**
  - SHA-256 of every file (`find . -type f`, sorted by path) is saved in `~/dr3_Kayla_sha256_local.txt` (271
    lines). That list's own SHA-256 is `09414836e5b01ef6083dd44762a3d39023f84da450e0ca6a36628d9abb895f16`.
  - `~/dr3_Kayla_sha256_cluster.txt` does not exist. **This copy is unverified** against the cluster original.
- **Inventory** (headers and counts only; row counts from newlines agree with a CSV parse):

  | Files | Rows | Columns |
  |---|---|---|
  | 101 `gaia_dr3_field_*.csv` (leaves) | 15,683,031 total (6,075–655,555 per leaf) | source_id, ra, dec, parallax, pmra, pmdec, radial_velocity, phot_g_mean_mag, bp_rp, ag_gspphot, ebpminrp_gspphot |
  | `field_main_sequence_uvw.csv` (1.39 GB) | 8,512,817 | source_id, ra, dec, d, phot_g_mean_mag, color0, ext_source, U, V, W |
  | `gaia_dr3_hosts_00/01/02.csv` | 500 / 500 / 122 | the field columns + parallax_over_error, visibility_periods_used, astrometric_chi2_al, astrometric_n_good_obs_al |
  | `hosts_hs19.csv`, `hosts_tess.csv`, `hosts_wcj.csv` | 282, 105, 315 | 39: archive planet columns, dr2, dr3, label, survey, Gaia columns, rv, d, A_G, E_BPRP, ext_source, color0, mg0, zams_offset, U, V, W, v_tot |
  | `hosts_uhj.csv` | 34 | the same 39 + `origin` (40) |
  | `nea_pscomppars.csv` | 6,451 | pl_name, hostname, pl_orbper, pl_bmassj, discoverymethod, disc_facility, disc_locale, disc_year, st_radv, sy_dist, gaia_dr2_id, gaia_dr3_id |
  | `sample_differences.csv` | 41 | hostname, dr3, label, status, sample, our_cut |
  | `dr2_to_dr3.csv` | 1 | dr2_source_id, dr3_source_id, angular_distance, magnitude_difference |

- **`field_main_sequence_uvw.csv.key`:**
  - 40 bytes: 40 lowercase hex characters with no newline, the length of a SHA-1 digest.
  - It is **not** the SHA-1 or BLAKE2b-160 of `field_main_sequence_uvw.csv`, nor the SHA-1 of any other file in
    the folder.
  - Its purpose is unknown, so it is treated as possibly sensitive. Its contents were not copied anywhere.
- **Field query split tree** (from the `.adql` text):
  - Shape: 154 queries, 48 top-level ranges `00`–`47` of equal width 2^57 in source_id. The ranges are contiguous
    from 0 to 6,917,529,027,641,081,856 = 12·4^12·2^35, the upper bound of Gaia source_ids.
  - Splits:
    - 53 split nodes, each with both halves `a`/`b`; the halves tile the parent range exactly, with no gap or
      overlap;
    - leaves at depths 0–7;
    - every leaf has a `.csv`, no split node has one, and there are no orphans.
  - All 154 queries share one WHERE clause apart from the source_id range.
  - **Truncation:** no leaf's row count equals a common Gaia archive limit (2,000 … 5,000,000), and none is a
    multiple of 1,000. Unsplit leaves hold up to 655,555 rows while other ranges were split down to 6,075, so the
    splits look like timeouts rather than row caps. This is an inference; the split nodes have no output to
    check.
  - Every row lies inside its query's source_id range. **No source_id appears in more than one leaf**, nor twice
    within a leaf.
- **Field cuts in the ADQL** (`gaiadr3.gaia_source`):
  - `radial_velocity IS NOT NULL`;
  - parallax > 0.4 mas; parallax_over_error > 10; visibility_periods_used > 8;
  - the HS19/Lindegren UWE cut: √(χ²/(n−5)) < 1.2·exp(−0.2(G−19.5)) for G < 19.5, otherwise < 1.2;
  - astrometric_gof_al < 3; ruwe < 1.4; rv_nb_transits > 5;
  - 1.0 + **0.0015**·bp_rp² < phot_bp_rp_excess_factor < 1.3 + 0.06·bp_rp²;
  - BP and RP flux_over_error > 10.
  - *Noted, not judged:* the excess-factor coefficient is HS19's printed 0.0015. Our Phase 1 primary uses
    Lindegren's 0.015 (D2), with 0.0015 only as a variant. The ADQL also has a RUWE cut and no
    astrometric_excess_noise_sig cut (our cut 6).
- **Host queries:** `gaiadr3.gaia_source` with an IN-list of source_ids only (no quality cuts):
  - 500 / 500 / 122 IDs, 1,122 unique, with no overlap between files;
  - every ID returned exactly once.
- **`dr2_to_dr3.csv`:** one query of `gaiadr3.dr2_neighbourhood` for a single DR2 source_id, returning 1 row.
  Nothing in the folder records how the other DR2→DR3 IDs (the `dr2`/`dr3` columns of the host files) were
  obtained.
- **Not done:** no code was run on these files beyond checksums, header lines, row counts, reading `source_id`
  for the duplicate and range checks, and ADQL parsing. No inventory script or data file was added to git.

**E46. Kayla's DR3 files: moved, permissions restricted, open questions for Phase 2** (user decisions on E45,
2026-11-10). Still no analysis.
- **New location:** `~/data/hu2026_dr3_Kayla/`, moved from `AST_192_Fall26/hu2026/` (supersedes the location in E45).
  - All 271 SHA-256 values recomputed after the move are identical to `~/dr3_Kayla_sha256_local.txt`.
  - **No git repository sees it.** Neither the folder nor `~` is inside a repo. `git status --ignored` in
    `AST_192_Fall26` and in `08_HS19_Kinematics` lists no `hu2026` path; 08's own ignored
    `logs/hu2026_power.log` is unrelated.
  - The parent repo's `.gitignore` was not touched.
- **Permissions:** `chmod -R go-rwx`, so the directory is 700 and all files 600 (owner only).
- **Checksum comparison with the cluster: still pending.** The user reported copying the cluster list to
  `~/dr3_Kayla_sha256_cluster.txt`, but no such file exists (checked 2026-11-10, also searched `~` to depth 3), so
  nothing was compared.
  - When it arrives, the agreed rule applies: compare by file name and hash, ignoring the `hu2026/` prefix and
    the cluster's CLAUDE.md line (added on the cluster on 2026-11-10, not part of Kayla's data).
  - The copy stays **unverified** until then.
- **`field_main_sequence_uvw.csv.key`:** left in place, not copied or displayed. The user will ask Kayla.
- **Open questions for the Phase 2 (DR3) pre-registration; no decision yet:**
  1. *Excess-factor cut:* Kayla's field query uses 1.0 + 0.0015·bp_rp² (HS19's printed coefficient). Our Phase 1
     primary uses Lindegren's 0.015 (D2).
  2. *RUWE:* her query adds ruwe < 1.4, which is not in HS19 or our Phase 1 cuts.
  3. *Excess noise:* her query has no astrometric_excess_noise_sig cut (our cut 6).
  4. *Which cuts can be reapplied:* the field chunks keep only 11 columns (source_id, ra, dec, parallax, pmra,
     pmdec, radial_velocity, phot_g_mean_mag, bp_rp, ag_gspphot, ebpminrp_gspphot).
     - Changing or adding any quality cut therefore needs a new query: parallax_over_error,
       visibility_periods_used, UWE, gof, RUWE, excess factor, flux errors, rv_nb_transits, excess noise.
     - Only the parallax, magnitude and colour cuts could be recomputed from these files.
     - The host files add parallax_over_error, visibility_periods_used, astrometric_chi2_al and
       astrometric_n_good_obs_al (UWE), but no gof, RUWE, excess factor, flux errors, rv_nb_transits or excess
       noise.
  5. *DR2→DR3 crossmatch:* `dr2_to_dr3.csv` holds a single-star test query (`gaiadr3.dr2_neighbourhood` for one
     DR2 source_id, 1 row). **No documented crossmatch of the hosts exists** in the files. How the `dr2`/`dr3`
     columns of the host files were filled is unknown.

**E47. Kayla's answers about her DR3 and DR2 files** (received via the user, logged 2026-11-11). **Information
only: no decision is taken here**, nothing was checked against the files, and no analysis was run.
- **`field_main_sequence_uvw.csv.key`:** a SHA-1 cache fingerprint of the build settings, not data. **Resolved**;
  the file stays in place. This answers the "purpose unknown" in E45 and E46. Its contents were still not
  displayed or copied.
- **Her DR3 field (`field_main_sequence_uvw.csv`):**
  - Quality cuts: Hu et al. (2026) Appendix A cuts 1–9, plus parallax > 0.4 mas.
  - The excess-factor lower bound 1.0 + **0.0015**·bp_rp² is hard-coded in cut 7. The files have **no
    excess-factor column**, so that cut cannot be redone from them (consistent with E46 question 4).
  - Distances: 1/parallax.
  - Extinction: GSP-Phot; otherwise Bayestar19 × 0.884; otherwise zero. Coefficients R_G, R_BP, R_RP = 2.740,
    3.374, 2.035.
  - Main sequence: less than 1 mag above the Pleiades ZAMS, and 0.5 ≤ (BP−RP)₀ ≤ 1.5.
  - Solar motion (U, V, W)☉ = (11.1, 12.24, 7.25) km/s.
  - Halo cut: v_tot > 180 km/s.
  - **Planet hosts are NOT removed** from the field file.
- **Her DR2 pipeline:** HS19 cuts 1–10, Bailer-Jones distances, STILISM extinction, no colour limit, colour
  matching within 0.025 mag. Her DR3 pipeline matches colour within 0.2 mag.
- **Host IDs:**
  - DR3 IDs: taken from the Hu et al. tables and matched to the NASA Exoplanet Archive by name.
  - HS19 DR2 IDs: from SIMBAD.
  - The `gaiadr3.dr2_neighbourhood` fallback (`dr2_to_dr3.csv`, E45/E46 question 5) was used only for nu Oct A
    (parallax 126 mas, bright binary).
  - TOI-5261: her ID differs from the NASA Exoplanet Archive's (no further detail given).
- **Her two pipelines differ in several factors at once:** cuts, distances, extinction, colour limit and
  colour-matching width. This is an **open question for the Phase 2 pre-registration**, alongside E46 questions
  1–5.
- **Not checked:** whether the ADQL WHERE clause recorded in E45 matches Hu Appendix A cuts 1–9; whether the one
  row of `dr2_to_dr3.csv` is nu Oct A; the TOI-5261 discrepancy. Checking any of these would mean selecting or
  joining hosts in DR3, which is not authorized.
- **Cluster checksum list (status, unrelated to Kayla's answers):** `~/dr3_Kayla_sha256_cluster.txt` now exists
  (125 bytes, two lines of text). It holds the two shell commands meant to create it (a `pbpaste` redirect and a `wc -l`
  check), not the checksum list. The clipboard apparently held the commands rather than the list. The comparison
  is still pending and the copy is still **unverified**.

**E48. Gaia DR2 → DR3 crossmatch of the 107 WASP hosts** (2026-11-11; Phase 2 groundwork, user-requested).
The user authorized this crossmatch, so the "no joining hosts to DR3" limit noted in E45–E47 no longer applies to
it. The ban on velocities stays: **no velocity, RV or kinematic quantity was read or computed for any host.**
Script: `src/dr3_crossmatch.py`. Output (gitignored): `data/processed/dr3_crossmatch/`.
- **Membership:** the frozen amendment 1 sample (source_id join, E36/E39).
  - It is rebuilt with column-restricted reads, not with `wasp_amendment1.samples()`, which loads the whole final
    host table including `radial_velocity`, `rv`, `U`, `V`, `W`.
  - Columns read: `hosts_hj_table1_final.csv` → archive_name_2019, hostname, source_id, ra, dec,
    phot_g_mean_mag; the 2026-10-28 archive snapshot → every column except velocity-like ones (a regex guard
    also refuses any velocity-like column in every table the script loads).
  - Facility label: as in E32/E36 (`hu2026_checks.facility_map` on the DR2 source_id).
  - **Check:** 297 hosts → 107 (22 WASP-South, 1 SuperWASP-North, 84 hemisphere not given). The D1 set
    (WASP name or facility) has 109, and 109 − 107 = WASP-53 and WASP-81, as recorded in E36.
- **Privacy:** before querying the Gaia archive, the IDs sent were checked against the private files, per the
  privacy rule. Only the 107 DR2 source_ids were sent.
- **Route 1 (primary):** host DR2 source_id → every row of Kayla's `nea_pscomppars.csv` with that `gaia_dr2_id`
  → its `gaia_dr3_id`. Columns read: pl_name, hostname, gaia_dr2_id, gaia_dr3_id only (not st_radv or sy_dist).
- **Route 2 (check):** `gaiadr3.dr2_neighbourhood` for all 107 DR2 source_ids, **keeping every candidate**, joined
  to `gaiadr3.gaia_source` for the candidate's ra, dec and phot_g_mean_mag only.
  - A synchronous TAP call (pyvo) got HTTP 503, so the query ran as an async job through `src/gaia_tap.py`.
  - The ADQL and raw answer are cached in `data/raw/gaia/dr3_crossmatch_wasp107_neighbourhood.*` and
    `data/processed/dr3_crossmatch/neighbourhood_raw.csv` / `neighbourhood_query.json`.
- **Agreement rule:** the routes agree if the NEA DR3 id equals the closest neighbourhood candidate (by
  angular_distance). Otherwise the host is flagged as "NEA id is a candidate but not the closest", "NEA id not
  among candidates", "no neighbourhood candidate", "NEA route gives no single DR3 id", or "no DR3 match".
- **Results:**
  - **Route 1:** all 107 hosts have exactly one NEA DR3 id. 102 hosts have one NEA planet row, 4 have two and 1
    has three; the rows of each host agree.
  - **Route 2:** every host has **exactly one** neighbourhood candidate (107 rows), and every
    proper_motion_propagation flag is true.
  - **The routes agree for all 107 hosts.** No host has more than one candidate, and no host lacks a DR3 match.
  - Angular distance (mas): min 0.024, median 0.173, max 0.467 (WASP-11).
  - magnitude_difference: −0.013 to +0.009 mag, |median| 0.006.
  - **DR3 id ≠ DR2 id for 6 hosts** (all matched by both routes). Their DR3 source_ids are recorded only in the
    local table `wasp107_dr2_dr3_summary.csv` (SHA-256 below):

    | Host | DR2 source_id | ang. dist. (mas) | Δmag |
    |---|---|---|---|
    | WASP-29 | 6534414719318224512 | 0.179 | +0.002 |
    | WASP-34 | 3537110833333561728 | 0.133 | +0.001 |
    | WASP-120 | 4837286842167819776 | 0.208 | −0.000 |
    | WASP-94 A | 6780546169633475456 | 0.175 | +0.005 |
    | WASP-100 | 4675352109658261376 | 0.248 | +0.001 |
    | WASP-140 | 5094154336330994688 | 0.155 | +0.005 |

  - The other 101 keep their DR2 source_id in DR3.
- **Files** (SHA-256):
  - `wasp107_dr2_dr3_candidates.csv` (one row per host × candidate) `b32f383b128a33a81c8cac6951a5e7caf43e9e66b7a16de8f3d15c77a860d363`;
  - `wasp107_dr2_dr3_summary.csv` (one row per host) `9298bacfd22cb0fdb0522e3b1d98c8d61bc5c8b187d69f66ec246623d37b12df`;
  - `neighbourhood_raw.csv` `c1d9d6339c41167fb800dec548e04731ac7285994d13d2c6cc7e1f61657567ce`.
- **Noted, not judged:**
  - Kayla's NEA file has no recorded download date. The membership's facility labels come from our 2026-10-28
    archive snapshot.
  - WASP-94 A has a known stellar companion (WASP-94 B), but the neighbourhood table lists only one candidate
    for its DR2 id.

**E49. Blocking privacy hooks** (2026-11-11; user decision). The privacy scan now runs as git hooks that abort on
any match, instead of as a manual step.
- **Why: the manual scan was bypassed by a masked exit status.**
  - Earlier on 2026-11-11, a local commit went ahead without its manual pre-commit scan being enforced. The
    scan script only *printed* counts and always exited 0, and the command chain (`scan | tail … && git
    commit …`) did not act on the printout. Even a failing status would have been masked by the `tail` in a
    pipeline run without `pipefail`.
  - That commit (4d2dbb7) was local and never pushed. It was amended under the privacy rule, and the amended
    branch scans clean.
  - 4d2dbb7 is no longer on any branch but stays in the local object store until git prunes it.
- **Files** (`scripts/hooks/`, committed; install with `scripts/hooks/install.sh`, which sets the repository-local
  `core.hooksPath = scripts/hooks`):
  - `privacy_scan.py`: standard library only, any `python3`. It never prints a private value, only location
    and kind. Exit 1 on any match.
    - **It fails closed (exit 2)** if the reference files are missing or empty.
    - Reference directory: `$AST192_PRIVATE_DIR`, default 07's `data/private` (extended_sample_all.csv,
      extended_sample_pool_k.csv, `rv/coralie*/dbase/*_coralie.rdb`).
  - `pre-commit`: scans the lines added in the index and the paths touched.
  - `commit-msg`: scans the commit message. This hook was added beyond the request, because pre-commit cannot see
    the message.
  - `pre-push`: scans **every commit being pushed**: `rev-list <pushed tip> --not
    origin/victoria-barriola-08`, i.e. origin/victoria-barriola-08..tip. For each commit it scans the message,
    the added lines and the paths. Without the base ref it scans the tip's full history.
  - All wrappers run under `set -euo pipefail`.
- **Match rules** (unchanged from the manual scan):
  - an extended-sample numeric ID;
  - an extended-sample name-like string (letters and digits, 7–40 characters; CSV headers excluded);
  - a CORALIE vrad value as written or rounded to 3 or 4 decimals;
  - a private path: `data/private/`, `data/eblm/`, `*.rdb`, or a path containing "coralie" or "extended_sample".
- **False-positive rate, measured before choosing "abort on any match":**
  - Replayed over the 28 commits then on the branch, the rules would have blocked 2. Both were 3-decimal
    coincidences with a rounded vrad value: one is the power value already explained in E37, the other a
    statistic in another committed file. No raw-precision or 4-decimal match occurred.
  - Only added lines are scanned, so existing lines block only if they are re-added or edited.
  - **On a coincidental match, the remedy is to change the number's formatting or wording.** It is not
    `--no-verify` and not an allowlist: a committed allowlist would itself disclose rounded RV values.
- **Tests** (`tests/test_privacy_hooks.py`, with an invented ID, name and RV value; no real private data):
  - In a throwaway repository with the committed hooks installed, pre-commit blocks a planted ID, a name, an RV
    value as written, an RV value rounded to 3 decimals, and a private path. commit-msg blocks a planted ID.
  - pre-push blocks a planted ID that skipped the commit hooks (`--no-verify`) one commit below the tip, and an RV
    value in a commit message. The remote is unchanged and the value is not echoed.
  - Clean commits and pushes pass.
  - Both hooks fail closed when the reference directory is missing.
  - Every wrapper uses `set -euo pipefail`.
- **Checked on the real repository:**
  - The current unpushed range (4 commits) is clean.
  - The scanner was also run on the discarded commit 4d2dbb7, as a check that it blocks without printing any
    value.
- **Limits:**
  - Hooks are local to each clone, so another clone (e.g. on the HPC) needs `scripts/hooks/install.sh`. There, a
    missing reference directory blocks every commit and push until `AST192_PRIVATE_DIR` points at a copy.
  - `git commit/push --no-verify` skips the hooks.
  - Binary files are checked by path only.

**E50. Kayla's two pipeline scripts: provenance; not to be executed before the Phase 2 freeze** (filed
2026-11-11; information only, no decision).
- **Sender:** Kayla. **Received 2026-11-11 via Slack**: the macOS quarantine attribute names Slack, with
  download times 13:09:26 and 13:09:31 local, the same as the files' creation times.
- **Where they arrived vs where they are now:**
  - They were expected at `~/data/hu2026_dr3_Kayla/code/`, but that folder did not exist. Both files were at
    the **root of the 08 repository, untracked** (`??`), never committed, never imported (no `.pyc`).
  - They were moved (`mv`, not copied) to `~/data/hu2026_dr3_Kayla/code/`, which is outside every repository.
    Permissions: directory 700, files 600. SHA-256 is identical before and after the move, and the repo
    root holds no `.py` file afterwards.

  | File | SHA-256 | Bytes | Lines |
  |---|---|---|---|
  | `hu2026_replication.py` | `6f7bd10a539e13de70b3b8a73f21720aaafbb2ee12c3550f7b58d81057057887` | 60,452 | 1,312 |
  | `velocity_dispersion_replication.py` | `7829cee33923f81b654e3354a943272127fa6a16db385a0300852d319260e119` | 56,757 | 1,209 |

- **They must not be run, imported or executed before the Phase 2 pre-registration is frozen.**
  - `hu2026_replication.py` `main()` computes σ_W of Gaia DR3 host samples, including the WASP hosts.
  - `velocity_dispersion_replication.py` has **no `main()` guard**: its whole analysis (field build, Monte
    Carlo, figure) runs at module level, so **importing it runs it**.
  - Both also create a `data/` folder next to themselves at import.
- **What was done:** both files were read as text only. They are treated as data, not instructions.
- **Contents** (from reading):
  - `velocity_dispersion_replication.py` is her HS19 (DR2) replication. Its docstring says its choices were
    aligned with our `notes/cuts_and_choices.md` at our commit 55aa85d. That commit exists in this repository
    ("Figure of the frozen WASP-only results").
  - `hu2026_replication.py` is her Hu et al. (2026) DR3 replication. As received, `USE_PUBLISHED_SAMPLES =
    True`, so it analyses the hosts of Hu Tables 2 and 3.
  - The code agrees with the field ADQL recorded in E45 and with Kayla's answers in E47 (the `.key` files are
    the SHA-1 of the build settings; the DR2→DR3 fallback; colour tolerances 0.025 and 0.2).
- **Differences from our frozen DR2 pipeline:** `notes/phase2_pipeline_differences.md`, 29 rows covering cuts,
  distances, extinction, MS selection, halo cut, field exclusions, matching, statistic, inference, RNG, host
  list, WASP membership, exclusions, host RVs, IDs and NEA snapshot.
  - K-DR3 differs in almost every block at once.
  - K-DR2 differs mainly in:
    - the cut-8 coefficient (0.0015);
    - the host-RV filters (`st_radvlim` missing → rejected; a broader Gaia-DR3 reference filter);
    - extra field exclusions (warm-Jupiter and Hu Table 3 TESS hosts);
    - the inference reported (P and a Gaussian equivalent; no CI);
    - its RNG.
  - These are open questions for the Phase 2 pre-registration, alongside E46 and E47. No effect size was
    computed.

**E51. Kayla's DR3 field file: its cache key reproduces from `hu2026_replication.py`** (2026-11-11; user
request; information only). **MATCH.**
- **Method:**
  - `cache_key()` and the inputs `get_field()` passes to it were copied verbatim from the source text of
    `hu2026_replication.py` (sha256 `6f7bd10a…`, E50) into a separate script
    (`~/data/hu2026_dr3_Kayla/checks/cache_key_check.py`, owner-only (dir 700, file 600), sha256
    `4bed043a816ba3db3b448cb4907df6098d23c9ea4918f2bb66c53edb6f18fe0e`; not committed, since it contains
    Kayla's code). It was written in the session scratchpad and moved there on 2026-11-11, with the SHA-256
    checked before and after the move.
    - What was copied: `field_queries()` with `field_edges()` and `field_query()`, `FIELD_CUTS_ADQL`, `_U`,
      `GAIA_COLUMNS`, `N_FIELD_CHUNKS`, `FIELD_MAX_DIST_PC`, and the constants in the call: `EXTINCTION_MAP`,
      `LOS_STEP_PC`, `BAYESTAR_MODE`, `R_G`/`R_BP`/`R_RP`, `BAYESTAR19_TO_EBV`, `ZAMS_COEFFS`, `MS_DELTA_MAG`,
      `COLOR_RANGE`, `SUN_UVW`, `HALO_VTOT`.
  - **Kayla's script was not imported or executed.** It was read as text only.
  - Before computing, the check script asserts four things:
    - each of its 10 copied blocks, and the `get_field()` call, occurs exactly once verbatim in the source;
    - each constant is assigned exactly once at module level;
    - there is no `global` statement;
    - the run uses Python 3.10.21 (`-I`).
- **Result:**
  - The recomputed SHA-1 `40c665fcd587794f4cac4e60e5cc14154436a602` **equals** the content of
    `field_main_sequence_uvw.csv.key`.
  - The `.key` file is unchanged since E45 (its SHA-256 equals the entry in `~/dr3_Kayla_sha256_local.txt`).
  - Both it and `field_main_sequence_uvw.csv` are dated 2026-11-04 14:18:32.
- **What the match confirms:** `field_main_sequence_uvw.csv` was written by a version of the code whose
  `get_field()` key inputs equal these. That covers:
  - the field query text of the 48 top-level source_id ranges: RV required, parallax > 0.4 mas, App. A cuts
    1–9 with excess coefficient 0.0015;
  - Bayestar19 (median) as the extinction fallback, × 0.884;
  - the CV18 coefficients;
  - the HS19 ZAMS with the 1 mag cut, and 0.5 ≤ (BP−RP)₀ ≤ 1.5;
  - solar motion (11.1, 12.24, 7.25);
  - the halo cut at 180 km/s.

  It **agrees with Kayla's answers in E47**.
- **What it does not confirm:**
  - The key hashes settings, not data or code. It does not cover the contents or version of
    `bayestar2019.h5`, nor the code of `add_photometry` / `add_kinematics` / `fetch_field_range`.
  - It does not cover the split sub-range queries. They share the WHERE clause (E45), but are not in the key.
  - It does not cover the contents of the CSV itself.
  - Planet-host removal happens after the cached file is read, so it is not in the key. This is consistent with
    E47: hosts are not removed in the file.

## Phase 2: DR3 WASP analysis, EXPLORATORY (2026-11-11)
**E52. The DR3 WASP analysis is run as EXPLORATORY, not pre-registered** (user decision, 2026-11-11). Decided
before any WASP DR3 kinematics were computed. **Results are preliminary; this analysis cannot later be
pre-registered as confirmatory.**
- **Already known before any computation:**
  - Hu et al. (2026) published WASP hosts with d < 600 pc: σ_W = 17.7 km/s vs field 18.4 km/s (102 hosts).
  - **Kayla's WASP DR3 output: the user has not seen her WASP DR3 numbers, tables or figures.** But her Slack
    message of 2026-11-11 states that her host σ_W matches Hu et al. to ≤ 0.1 km/s in every subsample,
    which includes WASP d < 600 (Hu: 17.7 vs 18.4, 102 hosts).
  - The same message reported, for non-WASP samples:
    - HS19 full sample σ_W 15.8 vs field 18.8 (this sample contains WASP hosts);
    - TESS HJs 18.1 vs 17.7;
    - CLS warm/cold Jupiters 19.9 vs 17.6.
  - The agent has not seen any WASP DR3 output of Kayla's. It read only the column names and row counts of her
    host files (E45). Her folder contains no results file.
- No WASP host velocity, dispersion or result has been computed in this project from DR3.

**E53. DR3 EXPLORATORY run, Stage 1 (2026-11-11): plan, code, seeds; Method A stopped (Bayestar19 missing);
Method B field download started.** EXPLORATORY: not pre-registered (E52). Full plan, with every choice and
adaptation: `notes/dr3_exploratory_plan.md`. **No WASP host velocity, dispersion or result was computed, and
no Gaia query was made for any WASP host.**
- **Seeds:** `dr3x_A_val_hs19`, `dr3x_A_val_tess`, `dr3x_A_wasp` (no tail); `dr3x_B_wasp` (chunks 0–7, 99).
  The registry tests pass, and no new label shares a code with any label.
- **Method A** (`src/dr3x_method_a.py`):
  - Our own implementation of Kayla's `hu2026_replication.py` (sha256 `6f7bd10a…`, E50), written from reading
    it as text. Her script was not imported, run or copied.
  - A missing Bayestar19 file raises an error (no zero fallback). The WASP run refuses without `--stage2`.
  - Unit-tested on synthetic data only: the HEALPix NESTED hierarchy and equal-area pixels, the Bayestar19
    rules on a fake map, inclusive pool bounds, the statistics, the one-sided σ.
  - **STOPPED: no Bayestar19 file** anywhere under `~`, and `dustmaps` is not installed. Per the user's rule,
    neither the host query nor the validation was run.
  - Hu et al. Tables 2/3 (the validation input) are also not available locally.
  - Four of the 107 are on her close-companion list (WASP-70 A, WASP-94 A, WASP-94 B, WASP-173 A). This is not
    applied; it is a pending decision.
- **Method B** (`src/dr3x_field_b.py`):
  - Cut port: cuts 1, 3, 5–7, 9, 10 port exactly. Cut 4 ports exactly for 5-parameter solutions only. Cut 8 is
    exact on the raw factor with 0.015.
  - Cut 2 (`mean_varpi_factor_al`) has no DR3 column. The 5-parameter restriction, the corrected C* and the
    Katz RV corrections are open. Distances: BJ21 `r_med_geo` is proposed.
  - Server-side query: RV required; cuts 1, 3, 5, 6, 7, 9, 10; cut 8 with the looser 0.0015. Everything else
    is applied locally, and all its inputs are downloaded. So no pending choice is fixed by the download.
  - **Service:**
    - ESA aborts anonymous async jobs after 60 s. A smoke test on 1/196,608 of the source_id range timed out,
      split three times, then got HTTP 429 (rate limit). No further ESA attempts were made.
    - ARI Heidelberg (as A6) answered the same slice in 17 s: 50 rows, all with a BJ21 distance, 3 with
      6-parameter solutions, 180 B/row gzipped. Extrapolation: about 10M rows, about 1.8 GB (an order of
      magnitude only).
  - **Download started** 2026-11-11, detached (pid 47363). Sequential, resumable, stops below 4 GB free.
    Output and log in `~/data/dr3x_field_b/` (outside every repository, owner-only).
- **Decisions needed:** listed in the plan, §5. No host kinematics are computed until the user has reviewed
  them.

**E54. DR3 EXPLORATORY run: the user's decisions on the Stage 1 questions** (2026-11-11). They supersede the
open items of E53 and `notes/dr3_exploratory_plan.md` §5.
1. **Bayestar19:** the full map is not downloaded.
   - Hosts without GSP-Phot extinction only are queried through the dustmaps **web** query (version
     bayestar2019, median mode, × 0.884).
   - Positions with no map coverage get zero, as in Kayla's method.
   - **If the web service fails: STOP, no silent fallback.**
2. **Validation input:** only the name and DR3-ID columns of Kayla's `hosts_hs19.csv` and `hosts_tess.csv`
   are read, no other column. **Hu et al. Tables 2/3 were not available locally.**
3. **Exclusions:**
   - Method A: Kayla's close-companion exclusions are applied for the **main** Method A result; Method A
     without them is reported as a labelled extra line.
   - Method B: no exclusions, as frozen.
4. **5-parameter solutions:** Method B requires `astrometric_params_solved = 31` in field and hosts, and any
   WASP host it removes is reported. Method A has no such requirement, as in Kayla's field.
5. **Cut 2** (`mean_varpi_factor_al`) is dropped, because it has no DR3 equivalent.
6. **Cut 8:** the raw `phot_bp_rp_excess_factor` with 0.015 (exact port). The Riello et al. (2021) corrected
   excess is noted as a possible later variant and is not used.
7. **No Katz et al. (2023) RV corrections** in either method. This is a known limitation.
8. **Distances (Method B):** Bailer-Jones et al. (2021) `r_med_geo`.
9. **Six-case rule (Method B):** the case is reported descriptively only, labelled "for comparison only".
10. **Archive RVs citing Gaia DR3** are still rejected (D6), and the Gaia DR3 RV is then used directly.
- Also: in E52, "[date not given]" was replaced by "[YYYY-MM-DD]" (the date of Kayla's Slack message, to be
  filled in).

**Implementation notes, 2026-11-11:**
- **Environment change, reverted:**
  - Installing `dustmaps` 1.0.14 with pip pulled in `healpy` 1.20.1, and pip **replaced numpy 1.26.4 with
    2.2.6**. The project's environment had numpy 1.26.4, per earlier sessions; `fleck` requires numpy < 2.
  - This was reverted within minutes: numpy 1.26.4 was reinstalled, with `healpy` 1.19.0 and `dustmaps`
    1.0.14 resolved against it. `pip check` reports no broken requirements, and all 68 tests pass.
  - No analysis ran under numpy 2.2.6. The field download (pid 47363) had started at 14:35 under 1.26.4, before
    the change at 14:45.
- **Bayestar19 web service failed:**
  - `argonaut.skymaps.info` resolves (to `argonaut.rc.fas.harvard.edu`), but the dustmaps query timed out
    connecting. Plain HTTP and HTTPS requests also got no answer (25 s each), while ESA answered normally.
  - Tested with arbitrary sky positions only; no host position was sent.
  - **Per decision 1: STOP for everything that needs Bayestar19.** The Method A extinction step for hosts
    without GSP-Phot and the validation (HS19-all, TESS d<600) were **not run**. There is no fallback.

**E55. DR3 EXPLORATORY Stage 1: host step counts for the 107 WASP hosts** (2026-11-11; `src/dr3x_hosts.py`;
output `data/processed/dr3_exploratory/host_steps_wasp107.json`, headed "EXPLORATORY — not pre-registered").
**No velocity, dispersion or other kinematic quantity was computed. RVs were used only for whether an RV
exists and where it comes from.**
- **Inputs:**
  - Gaia DR3 rows of the 107 DR3 IDs (E48), queried at ARI and joined to BJ21 `gaiaedr3_complements.distances`:
    107 of 107 returned.
  - 101 of the 107 have a GSP-Phot extinction pair (`ag_gspphot`, `ebpminrp_gspphot`). The other 6 are exactly
    the 6 PENDING hosts below.
- **Method A, main result (with Kayla's close-companion exclusions):** 107 → **99**.
  - Close companion: −4 (WASP-70 A, WASP-94 A, WASP-94 B, WASP-173 A).
  - No DR3 row: −0.
  - App. A cuts 1–3: −1 (WASP-11).
  - No RV: −0.
  - MS + 0.5 ≤ (BP−RP)₀ ≤ 1.5: −2 (WASP-90, WASP-172), of the hosts with GSP-Phot.
  - d = 1000/parallax < 600 pc: −1 (WASP-142).
  - **6 hosts are PENDING at the MS step:** WASP-32, WASP-66, WASP-68, WASP-80, WASP-93, WASP-107. They have
    no GSP-Phot, so they need Bayestar19, whose web service failed (E54). They are neither kept nor removed,
    so the final count is between 93 and 99.
- **Method A, extra line (without exclusions):** 107 → **103**, with the same steps and the same 6 pending
  (between 97 and 103).
- **Method B:** 107 → **101**.
  - No DR3 row: −0.
  - **`astrometric_params_solved = 31`: −5** (WASP-29, WASP-34, WASP-120, WASP-94 A, WASP-100). These are
    5 of the 6 hosts whose DR3 source_id differs from their DR2 one (E48); the sixth, WASP-140, passes.
  - Host cuts 1, 3, 4: −1 (WASP-11).
  - No RV (archive under D6, else Gaia DR3): −0.
  - No BJ21 `r_med_geo`: −0.
  - MS (STILISM along `r_med_geo`, no colour limit): −0.
- **Not run:** the Method A extinction for the 6 pending hosts, and the validation (HS19-all, TESS d<600).
  Both need Bayestar19, and its web service failed: STOP (E54).
- **Field download** (Method B, pid 47363): still running at the time of this entry. It is left to finish.

**E56. DR3 EXPLORATORY follow-up instructions** (user, 2026-11-11).
1. **E52:** the date of Kayla's Slack message is **2026-11-11**. It replaces "[YYYY-MM-DD]" (the user gave the
   date as "11/11/26").
2. **Bayestar19: wait.**
   - When the Method B field download finishes, the web query is retried **once**.
   - If it is still down, the size of the full bayestar2019 file on Harvard Dataverse and the free disk space
     are checked and reported **before any download**.
   - **The validation and Stage 2 are not run while the 6 hosts are pending** (E55).
3. **New exploratory comparison line (Stage 2, not run yet):** the frozen DR2 primary method (test B) on DR2
   data for exactly the 101 hosts that survive Method B, so that DR2 vs DR3 compares the same stars.
   - Seed label `dr3x_B_dr2match` (chunks 0–7, 99); it shares no code with any label.
   - Output labelled EXPLORATORY. Plan §3.6.
4. **Software environment:** future package installs go in a **separate virtual environment**, not the main
   `ast192-hj` environment. This follows the E54 incident, in which installing `dustmaps` replaced numpy
   1.26.4 (since reverted).

**E57. Method B field downloader: retries, a completeness check, and `caffeinate`** (user instructions,
2026-11-11). EXPLORATORY run (E52).
- **Why:** the download stopped twice after 15 of 48 ranges (00–14 complete, 686 MB).
  - 15:44: a DNS failure (`gaia.ari.uni-heidelberg.de` did not resolve).
  - 15:50, two minutes after a manual restart: ARI closed the connection without a response.
  - The downloader split ranges on ARI time limits but had no retry for network errors, so one dropped request
    ended the run.
- **Retries** (`with_retries`, as Kayla's `with_retries` does, E50, but with increasing waits):
  - What is retried: dropped connections, DNS failures, timeouts, and HTTP 429/5xx. The error's chain of
    causes is checked, because pyvo wraps some network errors.
  - Schedule: up to **10 tries per range**, waiting **1, 2, 4, 8, 16, 30, 30, 30, 30 min** between them.
  - **Every retry is logged** (range, try number, error type, wait). After the 10th failure the run stops.
  - **Errors that indicate a bad query are not retried and stop the run:** a rejected submission, a job that
    ends in ERROR for a reason other than time, or `DALQueryError`.
  - Time limits and the row cap still split the range, as before.
- **Completeness check** (`verify`): a range counts as complete only if all of the following hold.
  - The whole gzip stream reads; a truncated stream fails.
  - The data rows are below the row limit (10,000,000).
  - The last line ends with a newline and parses to the header's number of fields.
  - The written row count equals the job's.

  Files are written to `.tmp`, verified, then renamed. A `.ok` marker records the verified size, so unchanged
  files are not re-read on every restart. **A range that fails the check, or an interrupted `.tmp`, is deleted
  and rerun from scratch.** The finished ranges 00–14 are verified on the next run.
- **Tests** (`tests/test_dr3x_field_b.py`, offline, with a fake ARI job):
  - A simulated DNS failure (socket.gaierror inside requests.ConnectionError) and a dropped connection
    (RemoteDisconnected), each raw and wrapped by pyvo, are retried with waits of 1, 2 and 4 min and logged.
  - 10 failures stop the run after 9 waits.
  - A bad query (`BadQuery`, `DALQueryError`) is not retried.
  - Truncated, newline-less and wrong-field files, and the row limit, all fail `verify`.
  - A broken leftover and an interrupted `.tmp` are rerun; a verified range is not.
  - Splitting still works.
  - All 75 tests pass.
- **Restart:** under `caffeinate -i`, so the Mac cannot idle-sleep during the download. The finished ranges
  00–14 are kept.

**E58. Download restarted; Bayestar19 web query retried once: still down** (2026-11-11).
- **Field download:** restarted at about 17:52 under `caffeinate -i` (python pid 52003) with the E57 code.
  - **All 15 finished ranges (00–14) passed `verify`.** Their row counts equal those logged when they were
    written; for example, range 14 has 738,961 rows. They are kept.
  - Range 15 onward is running.
  - Log: `~/data/dr3x_field_b/download.log`.
- **Bayestar19 web query, one retry** (user instruction: independent of the download, one arbitrary sky
  position):
  - Query at 17:52 for (l, b) = (90°, 10°), d = 500 pc, `BayestarWebQuery(version="bayestar2019")`, median mode.
    No host position was sent.
  - **It failed again:** ConnectTimeout to `argonaut.skymaps.info:80` after 75 s, as on the first attempt (E54).
  - Nothing else was run. The validation and Stage 2 stay on hold while the 6 hosts are pending (E55, E56).
  - The check of the full map's size and the free disk space is still due when the download finishes (E56
    item 2).

**E59. Bayestar19 from the local file instead of the web query** (2026-11-12; user decision).
- **Why:** the dustmaps web service (`argonaut.skymaps.info`) has been down since about 14:50 on 2026-11-11.
  The retries at 17:52 and 20:22 both ended in ConnectTimeout (E54, E58). This replaces E54 item 1's route.
- **What changes:** only the route. Method A reads `bayestar2019.h5` from Harvard Dataverse
  (doi:10.7910/DVN/2EJ9TX; 0.73 GB, md5 ab815d2fd3068d1b81a1bd61fb18a722 as listed by the Dataverse API).
  Stored outside the repo at `~/data/bayestar19/`, owner-only permissions.
- **Unchanged:** median mode, × 0.884, CV18 coefficients. Positions with no map coverage get zero, as in
  Kayla's method.
- **Consistency check before use** (outside the repo, like `cache_key_check.py`): Kayla's `bayestar_ebv()`
  reimplemented from her source text (her script not imported or run), compared with dustmaps'
  `BayestarQuery` (median) on about 1,000 arbitrary positions and distances (no host positions). If the largest
  E(B-V) difference exceeds 0.001 mag, stop and report.
- **Then, only if they agree:** resolve the 6 pending hosts and run the validation (HS19-all and TESS
  d < 600 pc; no WASP-only numbers). No Stage 2.
- **Result (2026-11-12): the check FAILED the 0.001 mag limit, so the run stopped** before the hosts and the
  validation.
  - **File:** downloaded 05:51 to `~/data/bayestar19/bayestar2019.h5` (727,902,836 B, mode 600, directory
    700). md5 matches Dataverse. sha256 f68f14df…e32ef99.
  - **Check script:** `~/data/hu2026_dr3_Kayla/checks/bayestar_check.py`, sha256 85a69ccf…990aea, mode 600.
    1,000 positions uniform on the sphere, d log-uniform 20–3000 pc, seed 20261112 (a check, not a result
    stream).
  - **Pixel lookup and coverage agree.** Her HEALPix function matches healpy at all five nsides (0 / 1000
    mismatches). The same 758 positions are covered; 242 have no coverage (zero in both).
  - **Kayla's method vs dustmaps median: max |ΔE(B-V)| = 0.0219 mag.** All 24 positions above 0.001 mag lie
    between distance bins (24 of 587). Closer than the first bin and no coverage: 0 difference.
  - **Kayla's method vs the repository's `src.dr3x_method_a.bayestar_ebv`: 0 difference** at all 1,000.
  - **Cause (verified by a scratchpad diagnostic):** the order of median and interpolation. Kayla takes the
    median of the map's 5 samples at each distance bin, then interpolates in distance modulus. dustmaps
    interpolates each sample, then takes the median. Recomputing interpolate-then-median from the file
    reproduces dustmaps to 5e-8 mag. Over the 758 covered positions, the two orders differ by > 0.001 mag at
    24 (3 %); 95th percentile 0.0000, 99th percentile 0.0044, max 0.0219 mag.
  - Not run: the 6 pending hosts, the validation, Stage 2. They wait for the user's decision on which median
    order Method A uses.

**E60. Method A keeps Kayla's Bayestar19 order** (2026-11-12; user decision after E59).
- **Decision:** Method A takes the **median of the 5 map samples at each tabulated distance, then interpolates**
  in distance modulus (Kayla's order), via the repository's `src.dr3x_method_a.bayestar_ebv`. E59 showed it
  equals her method at all 1,000 test positions. dustmaps' order (interpolate each sample, then median) is not
  used.
- **Reason: host/field consistency.** Her field file (`field_main_sequence_uvw.csv`, E51) was built with this
  order. Host and field (BP-RP)_0 and M_G0 must come from the same extinction calculation, because the pools
  match on colour and the main-sequence cut uses both.
- **Logged check, no switch:** for the 6 hosts pending since E55, E(B-V) is computed in both orders, and any
  change in their main-sequence or colour-window outcome is reported. Extinction values only, no kinematics.
- **Then:** resolve the 6 hosts in Kayla's order, update the Method A host counts, and run the validation
  (HS19-all and TESS d < 600 pc). No WASP-only numbers, no Stage 2.
- **Six-host check, both orders** (`~/data/hu2026_dr3_Kayla/checks/bayestar_order_6hosts.py`, outside the
  repo). **max |ΔE(B-V)| < 0.00005 mag; no host's main-sequence or colour-window outcome differs.**
  - Bayestar19 values: WASP-32 0.0177, WASP-68 0.0619, WASP-93 0.0265, WASP-107 0.0000 (map value zero at its
    distance).
  - Zero by design: WASP-66 (no map coverage), WASP-80 (closer than the first distance bin).
- **The 6 hosts resolved (Kayla's order).** WASP-32, -66, -68 and -93 are kept. WASP-80 ((BP-RP)_0 = 1.79)
  and WASP-107 (1.54) fail the 0.5–1.5 colour window. Both are on the main sequence.
  `src/dr3x_hosts.py` no longer has a pending state.
- **Final Method A host counts (107 WASP hosts):**
  - **Main, with close-companion exclusions: 97.** 107 → −4 companions (WASP-70 A, -94 A, -94 B, -173 A)
    → 103 → −1 cuts 1–3 (WASP-11) → 102 → RV 102 → −4 main sequence + colour (WASP-90, -80, -107, -172) → 98
    → −1 d ≥ 600 pc (WASP-142) → 97.
  - **Without exclusions: 101** (same removals, without the companion step).
  - **Method B: 101** (unchanged, E55).
- **Validation** (`python -m src.dr3x_method_a validate` → `data/processed/dr3_exploratory/method_a_validation.json`):
  - **Inputs:** Hu et al. samples from the name and `dr3` columns of Kayla's `hosts_hs19.csv` (282) and
    `hosts_tess.csv` (105) only (E54). Gaia DR3 rows are from ARI (ESA anonymous jobs fail, E54).
  - **Steps:** RV = NEA `st_radv` by host name (her host-table rule), else Gaia DR3 RV. Photometry as for
    the hosts, UVW with her SUN_UVW (11.1, 12.24, 7.25). As in her pinned samples, no cuts: only hosts with no
    astrometry or RV are dropped (none were).
  - Her field file; NEA hosts removed. Seeds `dr3x_A_val_hs19`, `dr3x_A_val_tess`; 10,000 draws.
  - **HS19-all: N = 279** (282 published; 51 Peg, tau Boo and HD 217107 dropped with empty pools, under the
    same rule as her `run_analysis`). **σ_W hosts 15.77 (paper 15.8), field median 18.83 (paper 18.7)**,
    1σ 17.70–20.01, −2.89σ. Extinction: GSP-Phot 267, Bayestar19 9, zero 6.
  - **TESS d < 600 pc: N = 100 (paper 100;** 5 of 105 removed by the distance cut). **σ_W hosts 18.07 (paper
    18.1), field median 17.65 (paper 17.5)**, 1σ 15.87–19.65, +0.23σ. Extinction: GSP-Phot 98, Bayestar19 1,
    zero 1.
  - Host σ_W matches the paper to ≤ 0.05 km/s. The field medians are 0.13 and 0.15 km/s above the paper's.
  - Also computed on the same draws (HS19-all, hosts vs field median): σ_tot 40.25 vs 48.06 and eq. 1 36.05
    vs 41.89, both with no control below the hosts (≤ −3.72σ, the 1/N floor). TESS: σ_tot 43.56 vs 45.73
    (−0.79σ), eq. 1 38.16 vs 39.73 (−0.68σ).
  - **Code:** `validate` added to `src/dr3x_method_a.py`. Host rows now come from ARI (`query_hosts`). New
    tests cover `add_uvw` (solar motion, axes) and the RV-by-name rule; 77 tests pass.
- **Not run:** Stage 2 (no WASP σ_W, no Method B kinematics).
- **Note on the field medians (user, 2026-11-12, after the push).** The validation field medians are 0.13–0.15
  km/s above Hu et al.'s. They match the values Kayla reported (18.8 for HS19-all, 17.7 for TESS d < 600), so
  the offset comes from the shared field file (her `field_main_sequence_uvw.csv`), not from our implementation.
  Kayla's values are as relayed by the user; her outputs were not read here.

**E61. Stage 2 (EXPLORATORY): run as the plan specifies; implementation readings written before the run**
(2026-11-12).
- **User instruction:** run exactly what `notes/dr3_exploratory_plan.md` specifies, with no changes; if an input
  is missing or anything forces a change, stop.
- **Rows:**
  - Method A main (97) and without exclusions (101): σ_W (primary), σ_tot, HS19 eq. 1.
  - Method B (101, DR3): Δ with the frozen test B 95 % interval.
  - Method B on DR2 for the same 101 (`dr3x_B_dr2match`).
  - For comparison: the frozen DR2 107-host headline and Hu et al.'s WASP d < 600.
  - The frozen six-case rule for the Method B rows, "for comparison only" (E54 item 9).
- **Implementation readings.** None changes a rule of the plan; each fixes how a stated rule is carried out.
  1. **Field processing order.** The 48 verified leaves are processed one at a time (8 GB RAM; 15,968,381
    rows). The local cuts (§3.1: 1, 3, 4, 5, 6, 7, 8 at 0.015, 9, 10; `astrometric_params_solved = 31`; RV
    present) and the BJ21 `r_med_geo` requirement come first. STILISM extinction, the MS flag, UVW and z then
    run on the survivors only. Extinction does not enter any cut, so the final field is the same as
    extinction-first.
  2. **Planet-host removal (D17 adaptation, §3.3).** The lists are the three frozen Phase 1 host lists
    (`phase1_dr2/hosts_{hj_table1,hj,longp}_after_cuts.csv`: 338, 457, 510 rows), the set the frozen run removed.
    Each host's DR3 ID comes from the list's NEA `gaia_dr3_id` (our 2026-10-28 archive snapshot); 1 of 338
    HS19 Table 1 rows has none. It then takes the closest (`angular_distance`) `gaiadr3.dr2_neighbourhood`
    match of its DR2 ID, queried at ARI.
  3. **Method A seeds.** The plan registers one label, `dr3x_A_wasp`, so the main and no-exclusion rows both
    draw from it. Separate samples; no new label.
  4. **Method A WASP host RV and UVW.** NEA `st_radv` matched on the DR3 ID (plan A3), else Gaia DR3 RV. UVW
    as in the validation (`add_uvw`).
  5. **Method B hosts.** The 101 survivors of the E55 steps (`dr3x_hosts.method_b`), prepared like the field
    (r_med_geo, STILISM, MS, UVW, z). The DR2 line takes the same 101 by DR2 source_id (E48) from the frozen
    Phase 1 rows (`hosts_hj_table1_final.csv`) and the Phase 1 field.
  6. **Method B Monte Carlo** as `wasp_amendment1.vs_field` (pipeline.monte_carlo, eq. 1, z matching, 40,000
    controls; eblm_kinematics.bootstrap, 40,000 resamples). Chunk checkpoints are in
    `data/processed/dr3_exploratory/mc_checkpoints/`.
- **Run record (2026-11-12).**
  - **Method B field** (`python -m src.dr3x_stage2 field`, 07:28–08:00, `field_b_log.json`):
    - 15,968,381 rows → 15,835,037 after `astrometric_params_solved = 31`. Cuts 1, 3–10 (8 at 0.015) and the BJ21
      requirement remove no further star; cut 4 also removed none in the frozen DR2 field.
    - −138 planet hosts (941 unique DR3 IDs; CoRoT-32 via `dr2_neighbourhood`).
    - **7,867,685 main-sequence stars** (50 %; the DR2 field had 27 %).
  - **Two code fixes, no change of method:**
    - `write_csv` wrote plain text to a `.gz` path. The field file was gzipped in place (same content), and
      `write_csv` now gzips `.gz` paths (new test; 78 pass).
    - The first `run` stopped at Method B on a missing `source_id` column, needed for the checkpoint key. The
      host DR3 ID now fills it. The rerun recomputed every row from scratch.
- **Results (`data/processed/dr3_exploratory/stage2_results.json`, EXPLORATORY):**
  - **Method A**, her field (8,510,700 rows; 2,117 NEA hosts removed); no host dropped for an empty pool. Pool
    median 413, minimum 9. RV: NEA `st_radv` for every host. Extinction: GSP-Phot 93 / 97 (WASP-32, -68, -93
    Bayestar19; WASP-66 zero, outside the map).

    | Row | N | σ_W host / field median, 1σ, 2σ, one-sided σ | σ_tot | eq. 1 |
    |---|---|---|---|---|
    | main | 97 | 18.01 / 18.57, [16.69, 20.67], [15.06, 22.95], −0.29σ | 41.36 / 47.24, −2.12σ | 36.89 / 41.34, −1.95σ |
    | no exclusions | 101 | 17.78 / 18.64, [16.82, 20.67], [15.18, 22.85], −0.46σ | 42.00 / 47.26, −1.96σ | 37.44 / 41.37, −1.76σ |

  - **Method B, DR3** (101; `dr3x_B_wasp`):
    - σ = 37.75, control median 45.83: **Δ = −8.08 [−11.95, −4.36]**.
    - Fraction of controls ≤ host 0.000625. No unmatched host or empty z bin. Usable field 6,381,441.
    - Six-case: case 1 (for comparison only). RV: archive `st_radv` for all 101.
  - **Method B, DR2, same 101 hosts** (`dr3x_B_dr2match`, Phase 1 field 376,791):
    - σ = 37.86, control median 43.27: **Δ = −5.41 [−9.26, −1.71]**.
    - Fraction 0.011. Six-case: case 1 (for comparison only).
  - **Comparison:**
    - Frozen DR2 107-host headline Δ = −5.65 [−9.37, −2.04] (case 1; read from
      `data/processed/wasp_amendment1/results.json`).
    - Hu et al. WASP d < 600: 102 hosts, σ_W 17.7 vs 18.4 (published values as given by the user).
  - Between the two Method B lines the host σ barely changes (37.86 → 37.75). The control median rises
    (43.27 → 45.83) with the DR3 field. Recorded, not interpreted.

**E62. POST-HOC DIAGNOSTIC: why the Method B DR3 control median (45.83) exceeds DR2's (43.27)** (2026-11-12;
user request after seeing the Stage 2 results; written before any diagnostic run).
- **Status:** post-hoc diagnostic of the EXPLORATORY Stage 2. It changes nothing in the E61 results. Outputs go
  to `data/processed/dr3_exploratory/diag_*`, headed "EXPLORATORY — not pre-registered" and labelled POST-HOC
  DIAGNOSTIC.
- **1. Matched control stars, DR2 vs DR3.** The control stars actually drawn are reproduced, not resampled anew.
  - Chunk 0 (5,000 iterations × 101 hosts) of each E61 run is replayed: `dr3x_B_dr2match` on the Phase 1 field,
    and `dr3x_B_wasp` on the DR3 field.
  - Same stream, same `ControlSampler` calls, so these are the same picks as in those runs.
  - Properties of the picked stars, counted with multiplicity: distance (BJ18 `r_est` / BJ21 `r_med_geo`), G,
    `grvs_mag` (DR3 only; DR2 has none), `radial_velocity_error`, `parallax_over_error` (DR2: parallax /
    parallax_error) and |z|. Read from the DR2 prefilter cache and the DR3 leaves by source_id.
  - Also the per-component RMS of U, V, W of the picked stars.
- **2. DR3 field restricted to stars also in the DR2 field.**
  - "The DR2 field" is the frozen Phase 1 matching field (`field_ms.csv.gz`, 376,791 main-sequence stars).
  - Crossmatch: `gaiadr3.dr2_neighbourhood` at ARI, for DR2 stars passing the frozen DR2 server WHERE, in
    HEALPix level-2 source_id ranges. Each DR2 field star takes its closest (`angular_distance`) DR3 candidate,
    as E48.
  - The DR3 Method B field is kept where its source_id is one of those matches. Its own DR3 cuts, MS flag and
    values are unchanged.
  - Same 101 hosts, 40,000 controls + 40,000 bootstrap as in E61. New label **`dr3x_B_diag_dr2field`** (chunks
    0–7, 99).
- **3. D5 halo cut**, as the frozen sensitivity test D5: field stars with √(U²+V²+W²) > 180 km/s removed. The
  frozen test used 5,000 controls + 5,000 bootstrap, and so does this. New label **`dr3x_B_diag_halo`** (chunk
  0, 99). The user named only the label for item 2; this second label is our choice, needed for a separate
  stream.
- **Results** (2026-11-12; `diag_results.json`, `diag_run.log`; POST-HOC DIAGNOSTIC, EXPLORATORY):
  - Run: the ARI crossmatch took about 10 h, 192 level-2 ranges. Status-poll errors (DALServiceError, including
    a network outage around 12:38) were retried by `ari_tap`; no range failed.
  - **1. Replay:** both chunk-0 replays match the saved checkpoints exactly (max difference 0.0). There are
    505,000 picks each: 218,517 unique stars in DR2 and 476,563 in DR3.

    | matched controls, median [16–84 %] | DR2 | DR3 |
    |---|---|---|
    | distance (pc) | 423 [254, 638] | 752 [462, 1103] |
    | G (mag) | 12.64 [12.18, 13.17] | 14.16 [13.48, 14.63] |
    | radial_velocity_error (km/s) | 1.83 [0.99, 4.38] | 4.53 [2.52, 7.65] |
    | parallax_over_error | 78.5 [47.4, 147.1] | 73.5 [44.0, 131.8] |
    | \|z\| (pc) | 149 [59, 271] | 153 [63, 275] |
    | RMS U / V / W (km/s) | 38.41 / 26.20 / 19.50 | 41.43 / 28.27 / 19.99 |
    | eq. 1 of all picks (km/s) | 43.57 | 46.16 |
    | fraction v_tot > 180 km/s | 0.21 % | 0.38 % |

    `grvs_mag` (DR3 only): median 13.33.
  - **2. DR3 field restricted to DR2-field stars:**
    - All 376,791 DR2 field stars have a DR3 candidate (at most 4 each; the closest gives 376,791 unique DR3 ids).
    - Of these, 124,142 are in the DR3 Method B field.
    - Result: control median **43.25**, Δ = −5.49 [−9.35, −1.78] (40,000 + 40,000). Case 1, for comparison only.
  - **3. D5 halo cut:** 41,051 field stars removed. Control median **45.09**, Δ = −7.34 [−11.19, −3.63]
    (5,000 + 5,000). Case 1, for comparison only.
  - **Reading** (descriptive; post-hoc):
    - Matching is on z only, and the |z| distributions agree.
    - The DR3 picks are about 1.8× farther and about 1.5 mag fainter, with about 2.5× larger RV errors. This is
      the deeper DR3 RVS sample.
    - On the same stars that make up the DR2 field, the DR3 median (43.25) equals the DR2 median (43.27).
    - The halo cut accounts for only about 0.7 of the 2.6 km/s rise.
    - So the higher DR3 control median comes from the DR3 field's additional, more distant and fainter stars.
      It does not come from DR3 changing the values of the stars common to both fields.
    - This diagnostic does not separate how much of the rise is a real population difference and how much is
      from larger velocity errors.
