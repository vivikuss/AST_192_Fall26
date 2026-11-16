#!/bin/bash
# The PAWS re-analysis (prepare -> run -> collect) as ONE Slurm job. Needs the spectra from metallicity.sbatch
# (at least its download stage) and a one-time `bash hpc/setup_paws_env.sh`. Submit from the project root:
#
#   export METAL_ROOT=...  PAWS_HOME=...      # the same METAL_ROOT as the metallicity job
#   sbatch hpc/paws.sbatch
#   sbatch -d afterok:<metallicity job id> hpc/paws.sbatch     # or queue it behind the metallicity job
#
# Options (environment variables):
#   PAWS_TARGETS  published (default: the EBLMs with a published PAWS [Fe/H], ~87 groups) or all (~505 groups)
#   PAWS_ARGS     extra arguments for the run stage, e.g. "--fix-mh" (fix [M/H] in synthesis, as in the paper)
#   STAGES        (default "prepare run collect")
#   ONLY          comma-separated target_ids to restrict every stage to (for a quick test)
#   CONDA_ENV (ast192, for prepare/collect), PAWS_ENV (paws, for run), MODULES
# iSpec is single-threaded, so the job runs one PAWS group per CPU; each group takes ~0.5-3 h. With 32 CPUs the
# published set takes ~10 h. For PAWS_TARGETS=all ask for more CPUs or time (e.g. `sbatch -c 64 -t 3-00:00:00`),
# or use the job-array version, hpc/submit_paws.sh. Finished groups are skipped when the job is re-submitted.
#SBATCH --job-name=paws
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=paws_%j.log

set -euo pipefail
: "${METAL_ROOT:?set METAL_ROOT (export it before sbatch)}"
: "${PAWS_HOME:?set PAWS_HOME (see hpc/setup_paws_env.sh)}"
cd "${SLURM_SUBMIT_DIR:-.}"
[ -f hpc/paws_pipeline.py ] || { echo "submit from the project root"; exit 1; }
export METAL_ROOT PAWS_HOME
STAGES=${STAGES:-prepare run collect}
NPROC=${SLURM_CPUS_PER_TASK:-1}

if [ -n "${MODULES:-}" ]; then
    module purge
    # shellcheck disable=SC2086
    module load $MODULES
fi
source "$(conda info --base)/etc/profile.d/conda.sh"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR="${TMPDIR:-/tmp}/paws_${SLURM_JOB_ID:-local}"
mkdir -p "$TMPDIR" "$METAL_ROOT/logs"
trap 'rm -rf "$TMPDIR"' EXIT

P="python -u hpc/paws_pipeline.py"
only=(${ONLY:+--only "$ONLY"})
log() { echo "$(date '+%F %T')  $*"; }
log "host $(hostname), $NPROC CPUs, stages: $STAGES, targets: ${PAWS_TARGETS:-published}"

for stage in $STAGES; do
    log "== $stage"
    case $stage in
        prepare)
            conda activate "${CONDA_ENV:-ast192}"
            $P prepare --targets "${PAWS_TARGETS:-published}" ${only[@]+"${only[@]}"}
            ;;
        run)
            conda activate "${PAWS_ENV:-paws}"
            # one worker per CPU, each taking every NPROC-th group; a failed group is recorded in its
            # paws_result.json and doesn't stop the others
            pids=()
            for i in $(seq 0 $((NPROC - 1))); do
                # shellcheck disable=SC2086
                $P run --task "$i" --ntasks "$NPROC" ${only[@]+"${only[@]}"} ${PAWS_ARGS:-} > "$METAL_ROOT/logs/paws_run_$i.log" 2>&1 &
                pids+=($!)
            done
            failed=0
            for p in "${pids[@]}"; do wait "$p" || failed=$((failed + 1)); done
            [ $failed -eq 0 ] || log "warning: $failed PAWS workers failed (see $METAL_ROOT/logs/paws_run_*.log)"
            ;;
        collect)
            conda activate "${CONDA_ENV:-ast192}"
            $P collect
            ;;
        *)
            echo "unknown stage: $stage"; exit 1 ;;
    esac
done
log "done; results in $METAL_ROOT/results/paws_results.csv"
