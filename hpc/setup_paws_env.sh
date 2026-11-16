#!/bin/bash
# One-time set-up of iSpec v2023.08.04 + PAWS (pinned commit) in a dedicated conda env, on a Linux x86-64 cluster.
#
#   export PAWS_HOME=/cluster/tufts/<lab>/<user>/paws_software     # ~4 GB once unpacked
#   bash hpc/setup_paws_env.sh
#
# Run it on a node with internet and a C compiler (a login node is usually fine; it takes ~15 min).
# PAWS uses WIDTH (shipped with iSpec as Linux binaries) for the EW stage and SPECTRUM (compiled here) for synthesis.
set -euo pipefail
: "${PAWS_HOME:?set PAWS_HOME to where iSpec and PAWS should be installed}"
ENV=${PAWS_ENV:-paws}
ISPEC_URL=https://www.cfa.harvard.edu/~sblancoc/iSpec/iSpec_v20230804.tar.gz
PAWS_COMMIT=8ec055c960cb856cd477ead71d7444772b011899

if [ -n "${MODULES:-}" ]; then module load $MODULES; fi
source "$(conda info --base)/etc/profile.d/conda.sh"

# iSpec v2023.08.04 predates numpy 2 / Cython 3; pin to versions of that era
if ! conda env list | grep -q "^$ENV "; then
    conda create -y -n "$ENV" -c conda-forge python=3.10 "numpy=1.23" "scipy=1.10" "astropy=5.3" "pandas=1.5" \
        "matplotlib=3.7" "cython=0.29" statsmodels lockfile dill requests make c-compiler fortran-compiler
fi
conda activate "$ENV"

mkdir -p "$PAWS_HOME"
cd "$PAWS_HOME"
if [ ! -d iSpec_v20230804 ]; then
    curl -L --retry 5 -o iSpec_v20230804.tar.gz "$ISPEC_URL"     # 1.6 GB: code plus input/ (atmospheres, line lists)
    tar xzf iSpec_v20230804.tar.gz
    rm iSpec_v20230804.tar.gz
fi
# SPECTRUM is K&R-era C: GCC >= 15 defaults to C23, where `double f();` means no arguments, and GCC >= 14 turns
# implicit declarations / int conversions into errors. Compile it as gnu17 with those relaxed back to warnings.
# The flags go in CC on the make command line: the synthesizer/spectrum sub-Makefile sets its own CFLAGS, and
# command-line variables override that and are passed down to the sub-make.
(cd iSpec_v20230804 && CFLAGS="${CFLAGS:-} -std=gnu17 -fpermissive" make spectrum CC="${CC:-gcc} -std=gnu17 -fpermissive")
chmod +x iSpec_v20230804/synthesizer/atmos/bin.amd64/*.exe

if [ ! -d PAWS ]; then
    git clone https://github.com/alixviolet/PAWS.git
fi
(cd PAWS && git fetch -q origin && git checkout -q "$PAWS_COMMIT")

python - <<EOF
import sys
sys.path.insert(0, "$PAWS_HOME/iSpec_v20230804")
import ispec
sp = ispec.read_spectrum("$PAWS_HOME/iSpec_v20230804/input/spectra/templates/NARVAL.Sun.370_1048nm/template.txt.gz")
print("iSpec OK:", len(sp), "template pixels")
EOF
echo "done. Use: CONDA_ENV=$ENV PAWS_HOME=$PAWS_HOME bash hpc/submit_paws.sh"
