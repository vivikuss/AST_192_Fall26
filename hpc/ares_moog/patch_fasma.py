"""Python 3 / numpy >= 1.24 fixes for FASMA (DanielAndreasen/FASMA @ 1d40427, last updated 2019).

The equivalent-width code path (aresDriver, ewDriver, minimization, interpolation, utils, loggf_update) has three
Python-2-only constructs; nothing in the method changes. Each replacement must apply exactly as listed (or already
be applied), otherwise this script stops, so a different FASMA version cannot be patched silently.

    python hpc/ares_moog/patch_fasma.py $ARESMOOG_HOME/FASMA
"""
import os
import sys

PATCHES = [
    # numpy removed the np.float alias in 1.24
    ('aresDriver.py', 'return np.float(rounded)', 'return float(rounded)', 1),
    # dict.itervalues() is Python 2 only (outlier removal: 'allIter' / 'allOnce')
    ('ewDriver.py', 'for wavelength in outliers.itervalues():', 'for wavelength in outliers.values():', 2),
    # filter() returns an iterator in Python 3 (Readmoog.elements; not used by the EW method, fixed for safety)
    ('utils.py', "line = filter(None, line.split('abundance =')[1].split(' '))",
     "line = list(filter(None, line.split('abundance =')[1].split(' ')))", 1),
    ('utils.py', "line = filter(None, line.split(' '))\n                element.append",
     "line = list(filter(None, line.split(' ')))\n                element.append", 1),
]


def main(fasma):
    for fname, old, new, count in PATCHES:
        path = os.path.join(fasma, fname)
        with open(path) as f:
            src = f.read()
        n_old, n_new = src.count(old), src.count(new)
        if n_old == count:
            src = src.replace(old, new)
            with open(path, 'w') as f:
                f.write(src)
            print(f'patched  {fname}: {old.splitlines()[0][:60]}')
        elif n_old == 0 and n_new == count:
            print(f'ok       {fname}: already patched')
        else:
            sys.exit(f'{fname}: expected {count} x {old!r}, found {n_old} (and {n_new} patched); '
                     'is this FASMA commit 1d40427?')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ['ARESMOOG_HOME'], 'FASMA'))
