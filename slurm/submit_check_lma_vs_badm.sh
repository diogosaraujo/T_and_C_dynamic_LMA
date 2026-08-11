#!/bin/bash
## Compare the PLSR LMA driving the model against measured LMA in AmeriFlux BADM.
##
##     sbatch slurm/submit_check_lma_vs_badm.sh
##
## LMA is the treatment variable, so a bias here matters more than a bias in any
## parameter. At US-Ha2 the PLSR series (mean 102.6, range 78-119 g/m2) sits 22%
## below the BADM measurement of 131, which lies outside the whole 37-year series;
## substituting the measurement reproduces the reported LAI of 4.4 almost exactly.
## This asks how many stations can be checked and whether that offset is general.
##
## Reads only what is already on disk (BADM + the generated LMA_<ST>.mat files).
## Writes lma_vs_badm.csv next to the other coverage reports.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J lma_badm
#SBATCH -t 00:30:00
#SBATCH -o slurm/logs/lma_badm_%j.out
#SBATCH -e slurm/logs/lma_badm_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
MODEL_RUN="${MODEL_RUN:-$(dirname "$TC_INPUT_DATA")/model_run}"

BADM_DIR="${BADM_DIR:-$TC_INPUT_DATA/ameriflux}"
OUT_CSV="${OUT_CSV:-$TC_INPUT_DATA/lma_vs_badm.csv}"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "badm dir   : $BADM_DIR$([ -d "$BADM_DIR" ] || echo '   <-- NOT FOUND')"
echo "model_run  : $MODEL_RUN"
echo "output     : $OUT_CSV"
echo

[ -d "$BADM_DIR" ] || { echo "ERROR: $BADM_DIR not found -- has download_ameriflux.py run?" >&2; exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"

cd "$REPO_ROOT/preprocessing"
python check_lma_vs_badm.py \
    --badm "$BADM_DIR" --model-run "$MODEL_RUN" --csv "$OUT_CSV" "$@"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
exit $status
