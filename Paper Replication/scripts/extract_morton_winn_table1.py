"""One-off transcription of Morton & Winn (2014) Table 1 from the PDF into
data/morton_winn_table1.csv. Needs pypdf (not in ast192; run in a separate env):
    python scripts/extract_morton_winn_table1.py
Asymmetric radius errors are printed in the PDF as +hi then -lo; the source
footnote letters (a-e) are kept as provenance columns."""
import sys
from pathlib import Path
import pypdf

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT.parents[1] / "Papers" / "Morton_2014_ApJ_796_47.pdf"

import re, csv
txt = "\n".join(p.extract_text() for p in pypdf.PdfReader(PDF).pages)
# Row: K00003 a0.77±0.03\n0.02\nb29.47 ± 0.13 a2.26 ± 0.50 1 89.4
pat = re.compile(r"K(\d{5}) ([a-e])([\d.]+)±([\d.]+)\s*\n([\d.]+)\s*\n([a-e])([\d.]+) ± ([\d.]+) ([a-e])([\d.]+) ± ([\d.]+) (\d) ([\d.]+)")
SRC = {'a':'Buchhave+2012','b':'McQuillan+2014','c':'Hirano+2012b','d':'Hirano+2014','e':'CFOP'}
rows = []
for m in pat.finditer(txt):
    k,rs,r,rhi,rlo,ps,p,pe,vs,v,ve,n,ul = m.groups()
    rows.append(dict(koi=int(k), R_star=r, R_star_err_hi=rhi, R_star_err_lo=rlo, P_rot=p, P_rot_err=pe,
        vsini=v, vsini_err=ve, n_cands=n, inc_95UL_deg=ul,
        R_star_source=SRC[rs], P_rot_source=SRC[ps], vsini_source=SRC[vs]))
print(len(rows), 'rows;', sum(r['n_cands']=='1' for r in rows), 'singles')
w = csv.DictWriter(open(ROOT / 'data' / 'morton_winn_table1.csv', 'w', newline=''), fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
