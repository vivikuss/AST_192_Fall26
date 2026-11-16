#!/bin/bash
# One stage of a pipeline, as a single job or one task of an array job. Submitted by submit_all.sh / submit_paws.sh:
#   sbatch [--array=0-N] hpc/run_stage.sbatch <script.py> <stage> [extra args]
# Environment: METAL_ROOT (output directory), CONDA_ENV (default ast192), MODULES (optional `module load` list).
#SBATCH --job-name=metal
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -euo pipefail
script=$1; shift

if [ -n "${MODULES:-}" ]; then
    module purge
    # shellcheck disable=SC2086
    module load $MODULES
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-ast192}"

# keep BLAS single-threaded: parallelism comes from the process pool (one worker per CPU)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR="${TMPDIR:-/tmp}/metal_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR"

cd "${SLURM_SUBMIT_DIR:-.}"
echo "$(date)  host $(hostname)  task ${SLURM_ARRAY_TASK_ID:-0}/${SLURM_ARRAY_TASK_COUNT:-1}  $script $*"
python -u "hpc/$script" "$@"
echo "$(date)  done"
rm -rf "$TMPDIR"
