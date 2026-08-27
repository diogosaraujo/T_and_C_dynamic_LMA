#!/bin/bash
## The three hourly figures, from hourly_stats.csv.
##
##     source slurm/config.sh
##     sbatch -p SOE_legacy -A efthymios slurm/submit_hourly_figures.sh
##
## Writes to $TC_FIGURES: hourly_skill_maps.png, hourly_taylor.png,
## hourly_errorbars.png. Needs submit_hourly_stats.sh to have run first.
##
## HOURLY ONLY, and that is the point. Annual and seasonal skill were never
## defensible: the model record is 1985-2020 and many towers start after 2015,
## so at US-HBK the annual overlap is TWO YEARS and a correlation from two
## points is +/-1 by construction. Hourly gives tens of thousands of matched
## steps at every station.
##
## The CONUS outline is searched for automatically -- the maps shipped without
## one for several rounds purely because --basemap was never passed and nothing
## said so. The log now names the shapefile it found, or says it found none.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -p SOE_main
#SBATCH -J hrly_fig
#SBATCH -t 02:00:00
#SBATCH -o slurm/logs/hrly_fig_%j.out
#SBATCH -e slurm/logs/hrly_fig_%j.err

set -uo pipefail
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
tc_check_args "$@" || exit 1
mkdir -p "$TC_FIGURES" || { echo "ERROR: cannot create $TC_FIGURES" >&2; exit 1; }
echo "results    : $TC_RESULTS"
echo "figures    : $TC_FIGURES"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python -u figure_hourly.py "$@"
rc=$?
echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc"
exit $rc
