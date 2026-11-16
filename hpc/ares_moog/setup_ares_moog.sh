#!/bin/bash
# One-time install of the ARES+MOOG equivalent-width method (the SWEET-Cat recipe) on the cluster, via FASMA.
#
#   export ARESMOOG_HOME=/cluster/tufts/<lab>/<user>/ares_moog     # ~1 GB; not $HOME
#   bash hpc/ares_moog/setup_ares_moog.sh                          # from the project root, on a login node
#
# What it does (re-running is safe; each step is skipped or redone cleanly):
#   1. conda env at $ARESMOOG_HOME/env (conda-forge): python 3.10, numpy<2, pandas, scipy, astropy, statsmodels,
#      plus the C and Fortran compilers and the gsl + cfitsio libraries ARES needs (no cluster modules required).
#      The env AND conda's package cache live under $ARESMOOG_HOME, not in $HOME (the home quota is too small).
#   2. FASMA (DanielAndreasen/FASMA, pinned): python driver for ARES+MOOG (Andreasen et al. 2017), with
#        - MOOG (Sneden 1973; the 2014 version FASMA ships), built non-interactively as MOOGSILENT
#        - Kurucz ATLAS9 model atmospheres ("kurucz95"), the grid SWEET-Cat uses
#        - line lists: Sousa et al. 2008 (Sousa2007_opt_kurucz.lst) and Tsantaki et al. 2013 (cool stars)
#      patched for Python 3 by hpc/ares_moog/patch_fasma.py (3 one-line fixes, no method changes)
#   3. ARES v2 (sousasag/ARES, pinned), compiled with OpenMP; checked against its own solar test output
#   4. $ARESMOOG_HOME/bin/{ARES,MOOGSILENT}; then a solar test of the whole chain (hpc/ew_pipeline.py sun-test)
#
# Options: EW_ENV (env path, default $ARESMOOG_HOME/env), MODULES (e.g. "anaconda/2023.07" if conda comes from a
#          module), MOOG_DATA_LINK (see step 2), SKIP_SUN_TEST=1.
set -euo pipefail
: "${ARESMOOG_HOME:?set ARESMOOG_HOME to the install directory (project space, ~3 GB with the conda env)}"
PROJECT=$(cd "$(dirname "$0")/../.." && pwd)
ENV=${EW_ENV:-$ARESMOOG_HOME/env}
export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS:-$ARESMOOG_HOME/conda_pkgs}   # keep conda's download cache out of $HOME
FASMA_COMMIT=1d4042705c5180a6e609d08e7455d3d893e86392   # 2019-11-14, last FASMA (EW) commit
ARES_COMMIT=f85c9d73ec6fadedd4242349714b277b5bdba05d    # 2026-10-24
log() { echo "$(date '+%F %T')  $*"; }
log "setup_ares_moog.sh v2: running from $PROJECT"
log "  conda env -> $ENV"
log "  package cache -> $CONDA_PKGS_DIRS"
case "$(realpath -m "$ENV")/ $(realpath -m "$CONDA_PKGS_DIRS")/" in
    *"$(realpath "$HOME")/"*) echo "refusing to install into \$HOME (quota); set ARESMOOG_HOME to project space"; exit 1 ;;
esac

if [ -n "${MODULES:-}" ]; then
    module purge
    # shellcheck disable=SC2086
    module load $MODULES
fi
source "$(conda info --base)/etc/profile.d/conda.sh"

# ---------------------------------------------------------------- 1. conda environment
command -v git > /dev/null || { echo "git not found (load a git module)"; exit 1; }
env_ok() {   # an env left half-built by an interrupted or failed create must not be reused
    [ -x "$ENV/bin/python" ] && "$ENV/bin/python" -c 'import numpy, scipy, pandas, astropy, yaml, statsmodels' 2>/dev/null \
        && [ -e "$ENV/include/fitsio.h" ] && [ -e "$ENV/include/gsl/gsl_vector.h" ] && ls "$ENV"/bin/*-gfortran > /dev/null 2>&1
}
if ! env_ok; then
    [ -e "$ENV" ] && { log "removing incomplete env $ENV"; rm -rf "$ENV"; }
    log "creating conda env $ENV (package cache: $CONDA_PKGS_DIRS)"
    mkdir -p "$CONDA_PKGS_DIRS"
    conda create -y -p "$ENV" -c conda-forge --override-channels \
        python=3.10 "numpy>=1.24,<2" scipy "pandas>=1.5,<2.3" astropy pyyaml statsmodels patsy \
        gsl cfitsio c-compiler fortran-compiler make
    env_ok || { echo "conda env $ENV is incomplete after creation"; exit 1; }
    conda clean -y --tarballs > /dev/null   # the extracted packages stay (hard-linked); the archives go
fi
set +u; conda activate "$ENV"; set -u   # conda-forge compiler activation scripts are not set -u clean
export OMP_NUM_THREADS=1
log "env $ENV: $(python --version), CC=${CC:-cc}, FC=${FC:-gfortran}"

mkdir -p "$ARESMOOG_HOME/bin"
cd "$ARESMOOG_HOME"

clone_at() {   # url dir commit
    [ -d "$2/.git" ] || git clone -q "$1" "$2"
    git -C "$2" fetch -q origin
    git -C "$2" checkout -q -f "$3"
    log "$2 at $(git -C "$2" log -1 --format='%h %ad' --date=short)"
}

# ---------------------------------------------------------------- 2. FASMA + MOOG + models + line lists
clone_at https://github.com/DanielAndreasen/FASMA.git FASMA "$FASMA_COMMIT"   # checkout -f also resets old patches
python "$PROJECT/hpc/ares_moog/patch_fasma.py" "$ARESMOOG_HOME/FASMA"

if [ ! -d FASMA/models/kurucz95 ]; then
    log "extracting the Kurucz (ATLAS9) model atmospheres"
    tar -xzf FASMA/models/kurucz95.tar.gz -C FASMA/models
fi
[ -d FASMA/models/kurucz95 ] || { echo "kurucz95 models missing after extraction"; exit 1; }

# MOOG keeps the path to its data files (Barklem damping constants) in a CHARACTER*60 variable, so the path must
# be short: point it at a symlink in $HOME if the real path is too long.
MOOG_DIR="$ARESMOOG_HOME/FASMA/MOOG/"
if [ ${#MOOG_DIR} -gt 58 ]; then
    LINK=${MOOG_DATA_LINK:-$HOME/.moog_fasma}
    ln -sfn "$ARESMOOG_HOME/FASMA/MOOG" "$LINK"
    MOOG_DIR="$LINK/"
    [ ${#MOOG_DIR} -le 58 ] || { echo "MOOG data path $MOOG_DIR too long (>58 chars); set MOOG_DATA_LINK"; exit 1; }
fi
log "building MOOGSILENT (data path $MOOG_DIR)"
(
    cd FASMA/MOOG
    python - "$MOOG_DIR" <<'EOF'
import re, sys
src = open('Moogsilent.f').read()
new, n = re.subn(r"(moogpath =\n\s+\.\s+)'[^']*'", lambda m: m.group(1) + "'" + sys.argv[1] + "'", src)
assert n == 1, 'moogpath assignment not found in Moogsilent.f'
# fixed-form Fortran: code ends at column 72 (MOOG keeps comments beyond it, so the line length must stay standard)
line = next(l for l in new.splitlines() if sys.argv[1] in l)
assert len(line) <= 72, f'moogpath line is {len(line)} > 72 columns; use a shorter MOOG_DATA_LINK'
open('Moogsilent.f', 'w').write(new)
EOF
    make -s -f Makefile.rh64silent clean > /dev/null
    # flags as in FASMA's own Makefile.rh64silent (gfortran -std=legacy -fsecond-underscore), warnings silenced
    make -s -f Makefile.rh64silent \
        FC="${FC:-gfortran} -g -w -std=legacy -fsecond-underscore" > make.log 2>&1 \
        || { tail -30 make.log; echo "MOOG build failed (full log: FASMA/MOOG/make.log)"; exit 1; }
    [ -x MOOGSILENT ] || { echo "MOOGSILENT was not built"; exit 1; }
)
ln -sf "$ARESMOOG_HOME/FASMA/MOOG/MOOGSILENT" bin/MOOGSILENT

# ---------------------------------------------------------------- 3. ARES
clone_at https://github.com/sousasag/ARES.git ARES "$ARES_COMMIT"
log "building ARES"
(
    cd ARES
    ${CC:-gcc} -O2 -o ARES ARES_v2.c -I"$CONDA_PREFIX/include" -L"$CONDA_PREFIX/lib" -Wl,-rpath,"$CONDA_PREFIX/lib" \
        -lcfitsio -lgsl -lgslcblas -lm -fopenmp
    # ARES's own test: measure the solar (Ganymede) HARPS spectrum and compare with the shipped result
    rm -rf _test && mkdir _test && cp mine.opt linelist.dat sun_harps_ganymede.fits _test/
    (cd _test && ../ARES > ares_test.log 2>&1) || { tail _test/ares_test.log; echo "ARES test run failed"; exit 1; }
    python - <<'EOF'
import numpy as np
a, b = np.loadtxt('_test/test.ares'), np.loadtxt('test.ares_oric')
n = min(len(a), len(b))
d = np.abs(a[:n, 4] - b[:n, 4])
print(f'ARES test: {len(a)} lines measured ({len(b)} expected); max |dEW| = {d.max():.2f} mA, median {np.median(d):.3f} mA')
if len(a) != len(b) or d.max() > 1.0:
    print('WARNING: ARES output differs from its reference result; check the build before trusting EWs')
EOF
    rm -rf _test
)
ln -sf "$ARESMOOG_HOME/ARES/ARES" bin/ARES

# ---------------------------------------------------------------- 4. end-to-end solar test
export PATH="$ARESMOOG_HOME/bin:$PATH"
if [ "${SKIP_SUN_TEST:-0}" != 1 ]; then
    log "solar test of ARES -> MOOG -> FASMA (a few minutes)"
    METAL_ROOT=${METAL_ROOT:-$ARESMOOG_HOME/test_root} python "$PROJECT/hpc/ew_pipeline.py" sun-test
fi
log "done. Before running the EW pipeline: export ARESMOOG_HOME=$ARESMOOG_HOME"
