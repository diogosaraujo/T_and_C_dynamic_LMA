#!/bin/bash
## Write a harvested initial state into run directories that already exist.
##
##   step 5 -- historical arms, both from the one spin-up baseline:
##     sbatch -p main slurm/submit_apply_state.sh \
##         --from '*/historical/*/*_lma' --ic-key 'historical/{gcm}/spinup' \
##         --run-list run_list_gcm_historical.txt
##
##   step 6 -- ssp arms, each continuing from its OWN historical arm:
##     sbatch -p main slurm/submit_apply_state.sh \
##         --from '*/ssp*/*/*_lma' --ic-key 'historical/{gcm}/{arm}' \
##         --run-list run_list_gcm_ssp.txt
##
## USE THIS RATHER THAN REBUILDING. Between a pre-spin-up build and a post-spin-up
## one, four lines change: LAI_H, B_H, PHE_S_H, AgeL_H. Everything else in the arm
## directory is already correct -- the soil/root/canopy substitution, Sl_H (the
## GCM's own 1985-2014 mean, independent of the initial state), GO's ms and its
## '../Meteo_*.mat' load, the MAIN_FRAME/MAIN_FRAME_SLA choice, and LMA_<ST>.mat
## with the yearly SLA series.
##
## submit_gcm_model_run.sh would regenerate all of it identically, but needs two
## things that are not part of a RUN to do so: the era5_land MOD_PARAM it patches
## as a template, and the PLSR projection CSVs it rebuilds the LMA series from.
## Requiring those on the run cluster breaks the property model_run was
## restructured to have. Job 60692281 is the demonstration: 92 stations blocked on
## "no era5_land MOD_PARAM", for a change of four lines per file.
##
## Needs nothing but model_run and its initial_state.csv, so it runs wherever the
## runs run.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J apply_ic
#SBATCH -t 01:00:00
#SBATCH -o slurm/logs/apply_ic_%j.out
#SBATCH -e slurm/logs/apply_ic_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "model_run  : $MODEL_RUN$([ -d "$MODEL_RUN" ] || echo '   <-- NOT FOUND')"
echo

[ -d "$MODEL_RUN" ] || { echo "ERROR: no $MODEL_RUN" >&2; exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python apply_state.py --root "$MODEL_RUN" "$@"
rc=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc  (1 = at least one arm was refused; read the list above)"
exit $rc
