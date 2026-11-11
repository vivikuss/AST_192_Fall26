> **Header added 2026-11-09 (DECISIONS E35).** This file was made with the faulty WASP host list: the hostname join missed 9 WASP hosts (98 instead of 107; DECISIONS E32). It is kept unchanged as the record cited by the frozen pre-registration. See `PREREG_wasp_kinematics_amendment1_2026-11-10.txt` (amendment 1); the 9 added hosts have no label here.

# Discovery-camera labels for the 98 WASP hosts (DECISIONS E22)

`wasp_hemisphere_labels.csv` lists one row per host:
- the archive's discovery facility and bibcode;
- the paper that was read;
- the label: **S** (WASP-South only), **N** (SuperWASP-North only), **Both**, or **Not assigned**;
- the basis, taken from the paper's statement about which instrument supplied the discovery photometry.

Only metadata and paper text were used; no velocities. Source: the discovery paper listed in the NASA Exoplanet
Archive (`disc_refname`, 2026-11-04). Each bibcode was resolved to arXiv through the ADS link gateway, and the
paper's LaTeX source was read. Two papers are not on arXiv (WASP-18 in Nature, WASP-24 in ApJ); for those the
publisher page was read.

| Label | Hosts |
|---|---|
| S | 59 |
| N | 19 |
| Both | 16 |
| Not assigned | 4 (WASP-35, WASP-48, WASP-57, WASP-127) |

Caveats:
- **WASP-1 and WASP-2:** labelled N by inference. Their data come from the 2004 SuperWASP season, and only the
  northern instrument operated then.
- **WASP-18:** labelled S from the publisher page's acknowledgements and author contributions; the methods text
  is paywalled.
- **WASP-28:** the archive lists a follow-up paper (Petrucci et al. 2015), so the discovery paper (Anderson et
  al., arXiv:1402.1482) was used instead.
- **WASP-68:** the label comes from a paper-level statement (the title).
