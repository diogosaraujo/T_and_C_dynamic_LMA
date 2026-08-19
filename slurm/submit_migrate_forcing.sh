#!/bin/bash
## ONE-TIME: move already-built forcing out of input_data and into model_run.
##
##     sbatch slurm/submit_migrate_forcing.sh --dry-run   # what would move
##     sbatch slurm/submit_migrate_forcing.sh             # do it
##
## input_data holds what was DOWNLOADED; a T&C forcing file is something we built,
## so it belongs with the runs that read it:
##
##     model_run/<ST>/<scenario>/<GCM>/Meteo_<ST>_<GCM>_<scen>_<years>.mat
##     model_run/<ST>/era5_land/Meteo_<ST>_<years>.mat
##
## one copy per (station, scenario, GCM), directly above the fixed_lma/dyn_lma pair
## that shares it. That makes model_run self-contained: "rsync -a model_run/" is
## the whole transfer to Amarel, with no absolute symlinks into input_data to
## dangle on the far side and no duplicated forcing.
##
## New runs need none of this -- build_gcm_meteo.py and build_meteo_input.py stamp
## a dest_dir into each raw file and finish_meteo.m writes straight to it. This is
## only for what was built before that change.
##
## The move is os.rename inside one filesystem: instant, and it does NOT need a
## second copy of the 54 GB. Run it before submit_gcm_model_run.sh, because the
## run-tree builders now look for the forcing in its new home.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J migrate_fx
#SBATCH -t 02:00:00
#SBATCH -o slurm/logs/migrate_fx_%j.out
#SBATCH -e slurm/logs/migrate_fx_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "input_data : $TC_INPUT_DATA"
echo "model_run  : $MODEL_RUN$([ -d "$MODEL_RUN" ] || echo '   <-- NOT FOUND')"
echo

[ -d "$MODEL_RUN" ] || { echo "ERROR: no $MODEL_RUN -- run submit_build_model_run.sh first" >&2; exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python migrate_forcing.py --root "$MODEL_RUN" \
    --meteo "$TC_INPUT_DATA/gcm_meteo" --era5-meteo "$TC_INPUT_DATA/meteo" "$@"
rc=$?

echo
echo "left in input_data:"
find "$TC_INPUT_DATA/gcm_meteo" "$TC_INPUT_DATA/meteo" -name 'Meteo_*.mat' 2>/dev/null | wc -l
echo "now in model_run  :"
find "$MODEL_RUN" -mindepth 2 -maxdepth 4 -name 'Meteo_*.mat' 2>/dev/null | wc -l
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc  (1 = something was left behind; read the unmatched list above)"
exit $rc
