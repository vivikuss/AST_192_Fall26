# WASP-only analysis: amendment 1 rerun (107 hosts), frozen run (98 hosts) and diagnostic

Generated from `data/processed/wasp_amendment1/results.json`, `.../diagnostic_98/results.json` and the frozen
`data/processed/wasp/results.json` (DECISIONS E31, E41). Amendment: `PREREG_wasp_kinematics_amendment1_2026-11-10.txt`
(sha256 533a6437…). Δ in km/s; 95 % percentile bootstrap intervals (BCa was not triggered in any run).

- **Headline (107):** the result (amendment section 6). Case = section B six-case rule (δ_HS19 = −6.64 km/s);
  C reading = whole CI > 0 / includes 0 / < 0.
- **Frozen (98):** the original run, with the membership join fault (amendment 1), reported alongside, unchanged.
- **Diagnostic (98, new labels):** DIAGNOSTIC, not a result; no decision rule or reading applies. Frozen vs
  diagnostic = Monte Carlo stream change; diagnostic vs headline = effect of the join.
- E(iii) is not rerun; it is reported from the frozen run.

| Test | N (107 run) | Headline Δ [95 % CI] | Case / reading | N (frozen) | Frozen Δ [95 % CI] | Case / reading | N (diag.) | Diagnostic Δ [95 % CI] |
|---|---|---|---|---|---|---|---|---|
| B (primary) | 107 | -5.65 [-9.37, -2.04] | **1** | 98 | -5.81 [-9.69, -2.09] | 1 | 98 | -5.78 [-9.66, -2.03] |
| C (contrast) | 107/79 | +6.15 [+0.40, +11.60] | **whole CI > 0** | 98/79 | +6.31 [+0.58, +12.02] | whole CI > 0 | 98/79 | +6.27 [+0.57, +11.87] |
| D1 (name) B | 109 | -5.52 [-9.27, -1.93] | **1** | 100 | -5.60 [-9.53, -1.82] | 1 | 100 | -5.60 [-9.43, -1.84] |
| D1 (name) C | 109/79 | +6.02 [+0.35, +11.47] | **whole CI > 0** | 100/79 | +6.09 [+0.25, +11.75] | whole CI > 0 | 100/79 | +6.10 [+0.20, +11.61] |
| D2 (z + dist) B | 107 | -6.49 [-10.23, -2.84] | **1** | 98 | -6.68 [-10.57, -2.92] | 1 | 98 | -6.70 [-10.62, -2.97] |
| D2 (z + dist) C | 107/79 | +6.99 [+1.17, +12.34] | **whole CI > 0** | 98/79 | +7.18 [+1.27, +12.84] | whole CI > 0 | 98/79 | +7.20 [+1.15, +12.92] |
| D3 (d < 600 pc) B | 106 | -6.01 [-9.77, -2.52] | **1** | 97 | -6.23 [-9.98, -2.72] | 1 | 97 | -6.23 [-9.99, -2.49] |
| D3 (d < 600 pc) C | 106/71 | +7.83 [+1.63, +13.15] | **whole CI > 0** | 97/71 | +8.05 [+2.01, +13.54] | whole CI > 0 | 97/71 | +8.04 [+2.02, +13.80] |
| D4 (colour) B | 103 | -5.95 [-9.74, -2.30] | **1** | 95 | -6.04 [-10.09, -2.38] | 1 | 95 | -6.02 [-10.00, -2.12] |
| D4 (colour) C | 103/76 | +6.79 [+0.79, +12.40] | **whole CI > 0** | 95/76 | +6.88 [+0.91, +12.51] | whole CI > 0 | 95/76 | +6.86 [+0.95, +12.65] |
| D5 (halo cut) B | 107 | -5.29 [-9.07, -1.84] | **1** | 98 | -5.41 [-9.25, -1.69] | 1 | 98 | -5.49 [-9.37, -1.86] |
| D5 (halo cut) C | 107/79 | +6.16 [+0.78, +11.53] | **whole CI > 0** | 98/79 | +6.28 [+0.55, +12.10] | whole CI > 0 | 98/79 | +6.36 [+0.47, +11.88] |
| D6 C (no J1500-33) | 107/78 | +5.43 [-0.05, +10.85] | **CI includes 0** | 98/78 | +5.59 [-0.11, +11.20] | CI includes 0 | 98/78 | +5.56 [-0.00, +11.23] |
| E(i) σ_W WASP vs field | 107 | -1.18 [-5.21, +2.59] | includes 0 | 98 | -0.64 [-4.95, +3.30] | includes 0 | 98 | -0.66 [-4.98, +3.29] |
| E(ii) σ_W contrast | 107/79 | -1.88 [-6.58, +2.98] | includes 0 | 98/79 | -2.42 [-7.32, +2.53] | includes 0 | 98/79 | -2.40 [-7.36, +2.63] |
| E(iii) σ_W EBLMs (frozen, reused) | 79 | -3.06 [-5.99, -0.36] | excludes 0 | 79 | same | same | – | – |

Notes: the diagnostic's D6 lower limit is −0.005 km/s (shown as −0.00). The E(i) HS19-sized σ_W offset for the
107 hosts is −2.95 km/s. B details, 107 / 98 frozen / 98 diagnostic: host σ 37.49 / 37.23 / 37.23 km/s; control
median 43.14 / 43.04 / 43.01 km/s; fraction of controls ≤ host 0.0062 / 0.0071 / 0.0074.

**Reading (amendment section 6):** the 107-host run gives case 1 for B, which supports HS19: a deficit
consistent with HS19's size. Every D variant is also case 1. The C reading is "whole CI > 0" for C and D1–D5,
and "CI includes 0" for D6. E(i) and E(ii) include 0; E(iii) is unchanged. Every case and reading equals the
frozen run's, so the conclusion does not depend on the 9 added hosts.
