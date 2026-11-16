#!/bin/bash
# The whole metallicity pipeline (manifest -> download -> coadd -> fit -> combine) as ONE Slurm job.
# Submit from the project root:
#
#   export METAL_ROOT=/cluster/tufts/<lab>/<user>/metallicity     # ~60 GB free; not $HOME
#   sbatch hpc/metallicity.sbatch
#
# Options (environment variables, or override the #SBATCH lines on the command line, e.g. `sbatch -p preempt -c 32`):
#   STAGES      stages to run (default "manifest download coadd fit combine"); e.g. STAGES="fit combine" to refit
#   ONLY        comma-separated target_ids to restrict every stage to (for a quick test), e.g. WASP-43
#   N_DOWNLOAD  parallel download streams, 4 connections each (default 8)
#   CONDA_ENV   (ast192), MODULES (e.g. "anaconda/2023.07", loaded with `module load`)
# The manifest and download stages need internet on the compute node. If it has none, run those two on a login
# node first (see hpc/README.md) and submit with STAGES="coadd fit combine".
# For the multi-job version (job arrays, one dependency chain) use hpc/submit_all.sh instead.
#SBATCH --job-name=metallicity
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=metallicity_%j.log

set -euo pipefail
: "${METAL_ROOT:?set METAL_ROOT (export it before sbatch)}"
cd "${SLURM_SUBMIT_DIR:-.}"
[ -f hpc/metallicity_pipeline.py ] || { echo "submit from the project root"; exit 1; }
export METAL_ROOT
STAGES=${STAGES:-manifest download coadd fit combine}
N_DOWNLOAD=${N_DOWNLOAD:-8}

if [ -n "${MODULES:-}" ]; then
    module purge
    # shellcheck disable=SC2086
    module load $MODULES
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-ast192}"

# BLAS single-threaded: the parallelism comes from the process pool (one worker per CPU)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR="${TMPDIR:-/tmp}/metal_${SLURM_JOB_ID:-local}"
mkdir -p "$TMPDIR" "$METAL_ROOT/logs"
trap 'rm -rf "$TMPDIR"' EXIT

P="python -u hpc/metallicity_pipeline.py"
only=(${ONLY:+--only "$ONLY"})
log() { echo "$(date '+%F %T')  $*"; }
log "host $(hostname), ${SLURM_CPUS_PER_TASK:-?} CPUs, stages: $STAGES, output in $METAL_ROOT"

for stage in $STAGES; do
    log "== $stage"
    case $stage in
        download)
            # N independent streams; a failed stream doesn't stop the job (failures are listed in
            # catalogs/download_status_*.csv, and re-running the stage only fetches what is missing)
            pids=()
            for i in $(seq 0 $((N_DOWNLOAD - 1))); do
                $P download --task "$i" --ntasks "$N_DOWNLOAD" ${only[@]+"${only[@]}"} > "$METAL_ROOT/logs/download_$i.log" 2>&1 &
                pids+=($!)
            done
            failed=0
            for p in "${pids[@]}"; do wait "$p" || failed=$((failed + 1)); done
            tail -qn1 "$METAL_ROOT"/logs/download_*.log
            [ $failed -eq 0 ] || log "warning: $failed download streams failed (see $METAL_ROOT/logs/download_*.log)"
            ;;
        manifest|coadd|fit|combine)
            # coadd and fit use every CPU of the job through a process pool
            if [ "$stage" = combine ]; then $P combine; else $P "$stage" --task 0 --ntasks 1 ${only[@]+"${only[@]}"}; fi
            ;;
        *)
            echo "unknown stage: $stage"; exit 1 ;;
    esac
done
log "done; results in $METAL_ROOT/results"
