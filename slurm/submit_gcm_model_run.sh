#!/bin/bash
## Build the GCM arms of the model_run tree.
##
##     sbatch slurm/submit_gcm_model_run.sh --dry-run    # what is ready, writes nothing
##     sbatch slurm/submit_gcm_model_run.sh              # build everything available
##     sbatch slurm/submit_gcm_model_run.sh --gcm GFDL-ESM4 --scenario ssp585
##
##     model_run/<STATION>/<scenario>/<GCM>/{fixed_lma,dyn_lma}/
##
## PREREQUISITES, in order:
##   1. submit_build_model_run.sh --stations all   the era5_land arms, whose
##                                                 MOD_PARAM this reuses
##   2. submit_ssp_co2.sh                          CO2 series
##   3. submit_gcm_extract.sh  + submit_gcm_meteo.sh   the GCM forcing
##
## MOD_PARAM is a property of the SITE, not the forcing, so this copies the
## already-verified era5_land file and patches ONLY the Sl_H line rather than
## re-running the soil/root/canopy substitution. Two copies of that logic would
## eventually disagree, and the failure mode is a plausible, wrong MOD_PARAM.
##
## Writes $MODEL_RUN/run_list_gcm.txt as "<station> <scenario> <GCM> <arm>",
## which submit_gcm_tc_run.sh indexes.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -p SOE_main
#SBATCH -J gcm_build
#SBATCH -t 02:00:00
#SBATCH -o slurm/logs/gcm_build_%j.out
#SBATCH -e slurm/logs/gcm_build_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
MODEL_RUN="${MODEL_RUN:-$(dirname "$TC_INPUT_DATA")/model_run}"

export ECOREGION_ROOT="${ECOREGION_ROOT:-/vol_efthymios/NFS07/dd1136/ecoregions}"
PLSR_ROOT="${PLSR_ROOT:-$ECOREGION_ROOT/GCMs/PLSR_future_trait_predictions_GCM_clim_DOY_python/LMA}"
GCM_METEO="${GCM_METEO:-$TC_INPUT_DATA/gcm_meteo}"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "model_run  : $MODEL_RUN$([ -d "$MODEL_RUN" ] || echo '   <-- NOT FOUND')"
echo "plsr       : $PLSR_ROOT$([ -d "$PLSR_ROOT" ] || echo '   <-- NOT FOUND')"
echo "gcm meteo  : $GCM_METEO$([ -d "$GCM_METEO" ] || echo '   <-- NOT FOUND')"
echo

[ -d "$MODEL_RUN" ] || { echo "ERROR: run submit_build_model_run.sh first" >&2; exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python build_gcm_model_run.py \
    --root "$MODEL_RUN" --plsr-root "$PLSR_ROOT" --meteo "$GCM_METEO" "$@"
status=$?

echo
if [ -f "$MODEL_RUN/run_list_gcm.txt" ]; then
    N=$(wc -l < "$MODEL_RUN/run_list_gcm.txt")
    echo "run list   : $MODEL_RUN/run_list_gcm.txt  ($N arms)"
    echo "next       : sbatch --array=1-$N%NN slurm/submit_gcm_tc_run.sh"
    echo "             (MaxArraySize is 1001 on Amarel -- chunk if N exceeds it)"
fi
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
exit $status
