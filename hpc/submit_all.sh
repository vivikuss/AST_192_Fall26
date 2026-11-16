#!/bin/bash
# Submit the whole metallicity pipeline as a chain of Slurm jobs (each stage waits for the previous one to succeed).
#
#   export METAL_ROOT=/cluster/tufts/<lab>/<user>/metallicity     # somewhere with ~60 GB free
#   bash hpc/submit_all.sh                                         # from the project root
#
# Optional overrides: PARTITION (batch), DL_PARTITION (partition whose nodes have internet; default $PARTITION),
# N_DOWNLOAD (8), N_COADD (32), N_FIT (32), CPUS (8), CONDA_ENV (ast192), MODULES (e.g. "anaconda/2023.07").
# Re-running is safe: downloads, co-adds and the cached catalogues are reused; delete $METAL_ROOT/fits to refit.
set -euo pipefail
: "${METAL_ROOT:?set METAL_ROOT to the output directory (cluster scratch or project space)}"
cd "$(dirname "$0")/.."
export METAL_ROOT CONDA_ENV="${CONDA_ENV:-ast192}" MODULES="${MODULES:-}"

PARTITION=${PARTITION:-batch}
DL_PARTITION=${DL_PARTITION:-$PARTITION}
N_DOWNLOAD=${N_DOWNLOAD:-8}      # parallel streams per archive is the limit here, not CPUs; keep this modest
N_COADD=${N_COADD:-32}
N_FIT=${N_FIT:-32}
CPUS=${CPUS:-8}
LOGS="$METAL_ROOT/logs"
mkdir -p "$LOGS"

S="hpc/run_stage.sbatch metallicity_pipeline.py"
common=(--parsable --export=ALL)

j1=$(sbatch "${common[@]}" -p "$DL_PARTITION" -J metal-manifest -c 2 --mem=8G -t 02:00:00 \
     -o "$LOGS/manifest_%j.log" $S manifest)
j2=$(sbatch "${common[@]}" -p "$DL_PARTITION" -J metal-download -c 4 --mem=4G -t 24:00:00 \
     --array=0-$((N_DOWNLOAD - 1)) -d afterok:$j1 -o "$LOGS/download_%A_%a.log" $S download)
# afterany: a few failed downloads should not block the co-adds (failures are listed in catalogs/download_status_*)
j3=$(sbatch "${common[@]}" -p "$PARTITION" -J metal-coadd -c "$CPUS" --mem=$((CPUS * 2))G -t 06:00:00 \
     --array=0-$((N_COADD - 1)) -d afterany:$j2 -o "$LOGS/coadd_%A_%a.log" $S coadd)
j4=$(sbatch "${common[@]}" -p "$PARTITION" -J metal-fit -c "$CPUS" --mem=$((CPUS * 2))G -t 06:00:00 \
     --array=0-$((N_FIT - 1)) -d afterok:$j3 -o "$LOGS/fit_%A_%a.log" $S fit)
j5=$(sbatch "${common[@]}" -p "$PARTITION" -J metal-combine -c 1 --mem=8G -t 01:00:00 \
     -d afterok:$j4 -o "$LOGS/combine_%j.log" $S combine)

echo "submitted: manifest $j1 -> download $j2 -> coadd $j3 -> fit $j4 -> combine $j5"
echo "logs in $LOGS; results in $METAL_ROOT/results; watch with: squeue -u $USER"
