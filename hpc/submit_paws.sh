#!/bin/bash
# Submit the PAWS re-analysis (after metallicity_pipeline.py has downloaded the spectra) as a chain of Slurm jobs.
#
#   export METAL_ROOT=...  PAWS_HOME=...       # same METAL_ROOT as submit_all.sh; PAWS_HOME from setup_paws_env.sh
#   bash hpc/submit_paws.sh                     # EBLMs with a published PAWS [Fe/H] (reproduce BEBOP V / EBLM XIX)
#   PAWS_TARGETS=all bash hpc/submit_paws.sh    # every SOPHIE + HARPS group (method-vs-method comparison)
#
# Optional: PAWS_ARGS (extra args for the run stage, e.g. "--fix-mh"), N_PAWS (64), PARTITION (batch),
# PAWS_ENV (paws), CONDA_ENV (ast192, for the prepare/collect stages), AFTER (a job id to wait for, e.g. the
# metallicity download array), MODULES.
set -euo pipefail
: "${METAL_ROOT:?set METAL_ROOT}"
: "${PAWS_HOME:?set PAWS_HOME (see hpc/setup_paws_env.sh)}"
cd "$(dirname "$0")/.."
export METAL_ROOT PAWS_HOME MODULES="${MODULES:-}"

PARTITION=${PARTITION:-batch}
N_PAWS=${N_PAWS:-64}
LOGS="$METAL_ROOT/logs"
mkdir -p "$LOGS"
S="hpc/run_stage.sbatch paws_pipeline.py"
dep=()
if [ -n "${AFTER:-}" ]; then dep=(-d "afterany:$AFTER"); fi

j1=$(sbatch --parsable --export=ALL,CONDA_ENV="${CONDA_ENV:-ast192}" -p "$PARTITION" -J paws-prepare -c 1 --mem=8G \
     -t 04:00:00 "${dep[@]}" -o "$LOGS/paws_prepare_%j.log" $S prepare --targets "${PAWS_TARGETS:-published}")
# one group at a time per task (iSpec is single-threaded); a SOPHIE co-add + EW + synthesis takes ~0.5-3 h
j2=$(sbatch --parsable --export=ALL,CONDA_ENV="${PAWS_ENV:-paws}" -p "$PARTITION" -J paws-run -c 1 --mem=12G \
     -t 2-00:00:00 --array=0-$((N_PAWS - 1)) -d afterok:$j1 -o "$LOGS/paws_run_%A_%a.log" $S run ${PAWS_ARGS:-})
j3=$(sbatch --parsable --export=ALL,CONDA_ENV="${CONDA_ENV:-ast192}" -p "$PARTITION" -J paws-collect -c 1 --mem=4G \
     -t 00:30:00 -d afterany:$j2 -o "$LOGS/paws_collect_%j.log" $S collect)
echo "submitted: prepare $j1 -> run $j2 -> collect $j3"
