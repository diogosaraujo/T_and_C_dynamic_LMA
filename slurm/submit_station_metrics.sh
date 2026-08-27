#!/bin/bash
## Per-station fixed-vs-dynamic metrics and climate sensitivities.
##
##     source slurm/config.sh
##     sbatch -p SOE_legacy -A efthymios slurm/submit_station_metrics.sh
##
## Writes two tidy tables to $TC_RESULTS, both keyed per station AND PER GCM --
## the five members are never averaged here, so any pooling downstream is
## visible in the code that does it:
##
##   station_metrics.csv       SDs, SD ratios, and means of (dyn-fixed)/|fixed|
##                             by freq (annual/monthly/seasonal), by period, and
##                             by subset (all/drought/normal)
##   station_sensitivity.csv   yearly flux regressed on Ta, SPEI-12 and LMA, per
##                             arm, with delta_slope = dyn - fixed
##
## DELTA_SLOPE IS THE RESULT, the per-arm slopes are context. Both arms see
## identical forcing, so the moisture confounding that Jung et al. (2017) and
## Humphrey et al. (2018) raise about apparent temperature sensitivity sits in
## both slopes and cancels in the difference. Same logic as the C4MIP gamma
## feedback, which is also read off paired factorial runs.
##
## NEEDS FIRST: drought_periods_{era5,gcm}.csv from submit_drought_periods.sh,
## and period tables that contain Ta -- re-run submit_period_effect.sh if they
## predate that. Missing inputs are named and skipped, exit 1.
##
## No minimum-n guard, by request: a cell with one drought year yields NaN and
## its n is written, so filtering happens afterwards on real counts.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH -p SOE_main
#SBATCH -J st_metrics
#SBATCH -t 06:00:00
#SBATCH -o slurm/logs/st_metrics_%j.out
#SBATCH -e slurm/logs/st_metrics_%j.err

set -uo pipefail
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
tc_check_args "$@" || exit 1
echo "results    : $TC_RESULTS"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python -u station_metrics.py "$@"
rc=$?
echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc  (1 = an input was missing; read the list above)"
exit $rc
