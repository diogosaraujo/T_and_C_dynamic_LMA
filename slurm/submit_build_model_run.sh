#!/bin/bash
## Build the model_run tree from the T&C template + the preprocessed inputs.
##
##     sbatch slurm/submit_build_model_run.sh --dry-run     # what is ready, writes nothing
##     sbatch slurm/submit_build_model_run.sh               # the two example stations
##     sbatch slurm/submit_build_model_run.sh --stations all
##
## Writes OUTSIDE the repo, alongside input_data:
##     $TC_ROOT/model_run/Code/                shared T&C source (145 files, copied once)
##     $TC_ROOT/model_run/GRAPH_MOD.m          shared
##     $TC_ROOT/model_run/<STATION>/era5_land/{fixed_lma,dyn_lma}/
##     $TC_ROOT/model_run/<STATION>/{hist_gcm,ssp126,ssp245,ssp370,ssp585}/   empty for now
##
## A station is only built when every input is present. Anything missing is
## reported per station rather than filled with a template default -- an
## unsubstituted MOD_PARAM would run with US_xRM's soil and canopy and give
## plausible, wrong answers.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J build_run
#SBATCH -t 01:00:00
#SBATCH -o slurm/logs/build_run_%j.out
#SBATCH -e slurm/logs/build_run_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

MODEL_RUN="${MODEL_RUN:-$(dirname "$TC_INPUT_DATA")/model_run}"
TEMPLATE="${TEMPLATE:-$REPO_ROOT/T&C/Thanos_US_xRM}"
METEO_DIR="${METEO_DIR:-$TC_INPUT_DATA/meteo}"

mkdir -p "$REPO_ROOT/slurm/logs"

echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "model_run  : $MODEL_RUN"
echo "template   : $TEMPLATE"
echo "meteo      : $METEO_DIR$([ -d "$METEO_DIR" ] || echo '   <-- NOT FOUND')"
echo

if [ ! -d "$TEMPLATE/Code" ]; then
    echo "ERROR: template Code/ not found at $TEMPLATE/Code" >&2
    exit 1
fi
# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"

cd "$REPO_ROOT/preprocessing"
python build_model_run.py \
    --root "$MODEL_RUN" \
    --template "$TEMPLATE" \
    --meteo "$METEO_DIR" "$@"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status  (1 = a generated file is untrustworthy, or nothing built;"
echo "                        stations blocked for missing inputs are reported, not fatal)"
du -sh "$MODEL_RUN" 2>/dev/null || true
exit $status
