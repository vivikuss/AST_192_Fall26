"""Read HS19's Fig. 4 values from the vector figure in the arXiv source (A9).

    python -m src.paper_fig4      # writes data/processed/paper_fig4_values.json

combined_hist.pdf (matplotlib 3.0.2) draws the host value as an orange vertical
line (1 0.647 0 RG) and the 2/16/50/84/98 percentiles as black lines (0 G).
Data coordinates come from the x tick marks (36 and 38 km/s), not the labels.
"""
import json
import re
import tarfile
import zlib

from .config import RAW, ROOT

SRC = RAW / "arxiv" / "1908.06998_source.tar.gz"


def _streams(pdf):
    out = []
    for m in re.finditer(rb"stream\r?\n(.*?)\r?\nendstream", pdf, re.S):
        try:
            out.append(zlib.decompress(m.group(1)).decode("latin-1"))
        except zlib.error:
            pass
    return "\n".join(out)


def extract():
    with tarfile.open(SRC) as t:
        s = _streams(t.extractfile("combined_hist.pdf").read())
    # tick mark (short vertical segment from the axis) immediately before each x tick label
    ticks = {int(lab): float(x) for x, lab in
             re.findall(r"([\d.]+) [\d.]+ m\s+[\d.]+ [\d.]+ l\s+B\s+[\d.]+ [\d.]+ m\s+[\d.]+ [\d.]+ l\s+B\s+BT\s+"
                        r"/F1 [\d.]+ Tf\s+[\d.]+ [\d.]+ Td\s+\((\d\d)\) Tj", s)}
    x36, x38 = ticks[36], ticks[38]
    to_kms = lambda x: 36 + (x - x36) * 2 / (x38 - x36)
    col, lines = None, {"top": [], "bottom": []}
    pat = r"([\d.]+) ([\d.]+) ([\d.]+) RG|(?<![\d.])([\d.]+) G\b|([-\d.]+) ([-\d.]+) m\s+([-\d.]+) ([-\d.]+) l"
    for tok in re.finditer(pat, s):
        if tok.group(1):
            col = "orange" if float(tok.group(2)) > 0.5 and float(tok.group(1)) > 0.9 else "other"
        elif tok.group(4):
            col = "black" if float(tok.group(4)) == 0 else "other"
        else:
            x0, y0, x1, y1 = (float(tok.group(i)) for i in (5, 6, 7, 8))
            if x0 == x1 and abs(y1 - y0) > 20 and x36 < x0 < 480:       # full-height lines inside the axes
                lines["top" if y0 > 250 else "bottom"].append((col, round(to_kms(x0), 2)))
    res = {"note": "values read from the vector figure; HS19 give no numbers in the text",
           "pt_per_kms": (x38 - x36) / 2}
    for panel, name in (("top", "hj"), ("bottom", "longp")):
        host = [v for c, v in lines[panel] if c == "orange"]
        pct = sorted(v for c, v in lines[panel] if c == "black")
        res[name] = {"host_dispersion": host[0], "control_pct_2_16_50_84_98": pct}
    return res


if __name__ == "__main__":
    r = extract()
    (ROOT / "data" / "processed" / "paper_fig4_values.json").write_text(json.dumps(r, indent=1))
    print(json.dumps(r, indent=1))
