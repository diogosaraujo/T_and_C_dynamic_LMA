#!/bin/bash
## Fixed-vs-dynamic LMA per MONTH or SEASON, with the year kept.
##
##     source slurm/config.sh
##     sbatch slurm/submit_period_effect.sh --pair 'era5_land:*_ic' \
##         --freq monthly  --out era5_monthly.csv
##     sbatch slurm/submit_period_effect.sh --pair 'era5_land:*_ic' \
##         --freq seasonal --out era5_seasonal.csv
##     sbatch slurm/submit_period_effect.sh --pair 'ssp585/*:*' \
##         --freq monthly  --out gcm_monthly_ssp585.csv
##
## THE DIFFERENCE FROM submit_daily_effect.sh. That one writes a day-of-year
## CLIMATOLOGY -- every row a multi-year mean, no year column -- which answers
## "when in the year does the treatment bite" and cannot answer "how did 2012
## differ from 2011". This keeps the year, so the output is a TIME SERIES and
## trends, single drought years and event composites are all group-bys on it.
##
## Daily-with-year was rejected on size: ~50M rows for ERA5 alone and 15x that
## across the GCMs. Monthly divides that by ~30 and still resolves the seasonal
## cycle and the onset of a drought.
##
## NO --drought HERE. The year is in the output, so labels join on
## (station, gcm, scenario, year) afterwards. Baking the classes in would triple
## the rows and pin the table to one threshold -- the fixed cut and the
## percentile cut could not both be applied to the same file.
##
## SEASONS. DJF is filed under the year of its January, so December 2001 is DJF
## 2002. The first December and last Jan-Feb of a record land in short seasons;
## n_days shows it, and those are the ones to drop before fitting a trend.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH -p SOE_main
#SBATCH -J period_fx
#SBATCH -t 08:00:00
#SBATCH -o slurm/logs/period_fx_%j.out
#SBATCH -e slurm/logs/period_fx_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
tc_check_args "$@" || exit 1
mkdir -p "$TC_RESULTS" || { echo "ERROR: cannot create $TC_RESULTS" >&2; exit 1; }
echo "results    : $TC_RESULTS"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "model_run  : $MODEL_RUN$([ -d "$MODEL_RUN" ] || echo '   <-- NOT FOUND')"
echo

[ -d "$MODEL_RUN" ] || { echo "ERROR: no $MODEL_RUN" >&2; exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python -u analyze_period_effect.py --root "$MODEL_RUN" "$@"
rc=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc  (1 = at least one pair was skipped; read the list above)"
exit $rc
