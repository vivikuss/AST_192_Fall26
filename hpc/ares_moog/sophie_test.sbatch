#!/bin/bash
# SOPHIE HE diagnostic (hpc/ares_moog/sophie_test.py): select -> degrade -> run (10 variants x 14 stars) -> report.
# Submit from the project root, with the main metallicity directory as METAL_ROOT:
#
#   export METAL_ROOT=/cluster/tufts/martinlab/ndahle01/metallicity
#   export ARESMOOG_HOME=/cluster/tufts/martinlab/ndahle01/ares_moog
#   sbatch hpc/ares_moog/sophie_test.sbatch
#
# Writes only to $METAL_ROOT/sophie_test/. Options: STAGES (default "select degrade run report"), VARIANTS (for the
# run stage, e.g. VARIANTS=he_pedcorr), EW_ENV, MODULES.
#   e.g. STAGES="pedestal run report lines" VARIANTS=he_pedcorr sbatch hpc/ares_moog/sophie_test.sbatch
#SBATCH --job-name=sophie-test
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=sophie_test_%j.log

set -euo pipefail
: "${METAL_ROOT:?set METAL_ROOT (export it before sbatch)}"
: "${ARESMOOG_HOME:?set ARESMOOG_HOME (export it before sbatch)}"
cd "${SLURM_SUBMIT_DIR:-.}"
[ -f hpc/ares_moog/sophie_test.py ] || { echo "submit from the project root"; exit 1; }
export METAL_ROOT ARESMOOG_HOME OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -n "${MODULES:-}" ]; then
    module purge
    # shellcheck disable=SC2086
    module load $MODULES
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u; conda activate "${EW_ENV:-$ARESMOOG_HOME/env}"; set -u
export PATH="$ARESMOOG_HOME/bin:$PATH"

for stage in ${STAGES:-select degrade run report}; do
    echo "$(date '+%F %T')  == $stage"
    extra=()
    [ "$stage" = run ] && [ -n "${VARIANTS:-}" ] && extra=(--variants "$VARIANTS")
    python -u hpc/ares_moog/sophie_test.py "$stage" "${extra[@]+"${extra[@]}"}"
done
echo "$(date '+%F %T')  done; results in $METAL_ROOT/sophie_test/"
