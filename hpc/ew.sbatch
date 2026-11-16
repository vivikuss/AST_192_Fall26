#!/bin/bash
# ARES+MOOG [Fe/H] (hpc/ew_pipeline.py: coadd -> fit -> collect) as ONE Slurm job. Submit from the project root:
#
#   export METAL_ROOT=/cluster/tufts/<lab>/<user>/ew_validation   # or the main metallicity directory
#   export ARESMOOG_HOME=/cluster/tufts/<lab>/<user>/ares_moog     # where setup_ares_moog.sh installed
#   COMPARE=hpc/ares_moog/validation_targets.csv sbatch hpc/ew.sbatch
#
# Run the manifest (custom targets only) and download stages first, on a login node (they need internet):
#   python hpc/ew_pipeline.py manifest --targets hpc/ares_moog/validation_targets.csv
#   python hpc/metallicity_pipeline.py download --task 0 --ntasks 1
#
# Options (environment variables, or override the #SBATCH lines on the command line, e.g. `sbatch -c 32`):
#   STAGES      stages to run (default "coadd fit collect"); e.g. STAGES="fit collect" to refit only;
#               after the main run: STAGES="rescue collect adopt" (fixed-logg re-fits, then one [Fe/H] per star)
#   ONLY        comma-separated target_ids to restrict coadd/fit to, e.g. TOI-426
#   LINELIST    force one line list for every star (sousa | tsantaki); default: per-target column, else sousa
#   COMPARE     CSV with published values for the collect stage (e.g. hpc/ares_moog/validation_targets.csv)
#   COADD_ENV   (ast192), EW_ENV ($ARESMOOG_HOME/env), MODULES (e.g. "anaconda/2023.07", loaded with `module load`)
#SBATCH --job-name=ew-feh
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=ew_%j.log

set -euo pipefail
: "${METAL_ROOT:?set METAL_ROOT (export it before sbatch)}"
: "${ARESMOOG_HOME:?set ARESMOOG_HOME (export it before sbatch)}"
cd "${SLURM_SUBMIT_DIR:-.}"
[ -f hpc/ew_pipeline.py ] || { echo "submit from the project root"; exit 1; }
export METAL_ROOT ARESMOOG_HOME
STAGES=${STAGES:-coadd fit collect}

if [ -n "${MODULES:-}" ]; then
    module purge
    # shellcheck disable=SC2086
    module load $MODULES
fi
source "$(conda info --base)/etc/profile.d/conda.sh"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR="${TMPDIR:-/tmp}/ew_${SLURM_JOB_ID:-local}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

P="python -u hpc/ew_pipeline.py"
only=(${ONLY:+--only "$ONLY"})
ll=(${LINELIST:+--linelist "$LINELIST"})
log() { echo "$(date '+%F %T')  $*"; }
log "host $(hostname), ${SLURM_CPUS_PER_TASK:-?} CPUs, stages: $STAGES, METAL_ROOT=$METAL_ROOT"

for stage in $STAGES; do
    log "== $stage"
    case $stage in
        coadd)   set +u; conda activate "${COADD_ENV:-ast192}"; set -u
                 $P coadd --task 0 --ntasks 1 ${only[@]+"${only[@]}"} ;;
        fit)     set +u; conda activate "${EW_ENV:-$ARESMOOG_HOME/env}"; set -u
                 export PATH="$ARESMOOG_HOME/bin:$PATH"
                 $P fit --task 0 --ntasks 1 ${only[@]+"${only[@]}"} ${ll[@]+"${ll[@]}"} ;;
        collect) set +u; conda activate "${EW_ENV:-$ARESMOOG_HOME/env}"; set -u
                 $P collect ${COMPARE:+--compare "$COMPARE"} ;;
        rescue)  set +u; conda activate "${EW_ENV:-$ARESMOOG_HOME/env}"; set -u
                 export PATH="$ARESMOOG_HOME/bin:$PATH"
                 $P rescue --task 0 --ntasks 1 ${only[@]+"${only[@]}"} ;;
        adopt)   set +u; conda activate "${EW_ENV:-$ARESMOOG_HOME/env}"; set -u
                 $P adopt ;;
        *)       echo "unknown stage: $stage"; exit 1 ;;
    esac
done
log "done; results in $METAL_ROOT/results"
