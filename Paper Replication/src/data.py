"""Loading the Morton & Winn (2014) Table 1 sample, with guards against the two
bookkeeping errors the spec warns about."""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MW14_TABLE = ROOT / "data" / "morton_winn_table1.csv"

ALPHA_SUN = 0.23  # MW14 Sec. 4.2 differential-rotation shear


def load_mw14_table(*, p_rot_err_if_zero, path=MW14_TABLE):
    """Table 1 of Morton & Winn (2014), 70 KOIs.

    p_rot_err_if_zero : float (days), required.
        Replacement for rotation-period errors printed as 0.00 (KOI-63 only).
        Every other P_rot_err is used as printed.

    The 'vsini' column is MW14's "f V sin I": already divided by (1 - α/2) for
    differential rotation (Buchhave values) or taken from Hirano et al., which
    were already corrected. It must be used as is. The returned frame carries
    vsini_diffrot_corrected = True so downstream code can assert on it.
    """
    df = pd.read_csv(path)
    assert len(df) == 70, f"expected 70 KOIs, found {len(df)}"
    assert df.koi.is_unique

    # multiplicity: n_cands = 1 -> single, n_cands >= 2 -> multi
    assert (df.n_cands >= 1).all(), "n_cands must be >= 1"
    df["group"] = np.where(df.n_cands == 1, "single", "multi")
    assert ((df.group == "single") == (df.n_cands == 1)).all()
    assert ((df.group == "multi") == (df.n_cands >= 2)).all()
    assert (df.group == "single").sum() == 45 and (df.group == "multi").sum() == 25

    # rotation-period error floor, applied only where the table prints 0.00
    zero = df.P_rot_err <= 0
    assert set(df.koi[zero]) <= {63}, f"unexpected zero P_rot_err: {list(df.koi[zero])}"
    df["P_rot_err_raw"] = df.P_rot_err
    df.loc[zero, "P_rot_err"] = p_rot_err_if_zero
    df["P_rot_err_source"] = np.where(zero, "floor (printed 0.00)", "table")

    # v sin i is already corrected for differential rotation: never divide again.
    df["vsini_diffrot_corrected"] = True
    check_vsini_not_double_corrected(df)
    return df


def check_vsini_not_double_corrected(df, alpha=ALPHA_SUN):
    """KOI-975: MW14 Sec. 4.1 adopts V sin I = 7.5 km/s (CFOP); Table 1 lists
    7.5 / (1 - α/2) = 8.47. If the column has been divided a second time, it
    would read 8.47 / 0.885 = 9.57."""
    row = df.loc[df.koi == 975].iloc[0]
    assert row.vsini_diffrot_corrected
    assert abs(7.5 / (1 - alpha / 2) - 8.47) < 0.01
    assert abs(row.vsini - 8.47) < 0.005, (
        f"KOI-975 vsini = {row.vsini}: expected 8.47 (corrected once)")
    assert abs(row.vsini * (1 - alpha / 2) - 7.5) < 0.01
