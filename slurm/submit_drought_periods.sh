#!/bin/bash
## Per-period drought labels, on the accumulation window that matches the step.
##
##     source slurm/config.sh
##     sbatch -p SOE_legacy -A efthymios slurm/submit_drought_periods.sh --source era5
##     sbatch -p SOE_legacy -A efthymios slurm/submit_drought_periods.sh --source gcm
##
## Writes drought_periods_<source>.csv to $TC_RESULTS:
##     source,gcm,scenario,station,freq,year,period,spei,class
##
## THE WINDOW MATCHES THE STEP:
##     annual    SPEI-12 at SEPTEMBER, the water-year end in that calendar year
##     monthly   SPEI-3 at that month
##     seasonal  SPEI-3 at the season's LAST month (DJF->Feb, JJA->Aug, ...)
##
## This REPLACES the drought_years*.csv that classify_drought.py wrote. Those
## held the ANNUAL MEAN of monthly SPEI-12 -- an average over twelve overlapping
## 12-month windows, which is not a quantity anyone wants. Verified on
## ACCESS-CM2/US-Bar/2000: September SPEI-12 is -0.7268 while the annual mean is
## -0.6143, and this table stores -0.7268.
##
## One table, one definition, joined by everything downstream, so the figures and
## the metrics cannot drift apart the way they already had.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH -p SOE_main
#SBATCH -J drought_p
#SBATCH -t 04:00:00
#SBATCH -o slurm/logs/drought_p_%j.out
#SBATCH -e slurm/logs/drought_p_%j.err

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
echo

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python -u drought_labels.py "$@"
rc=$?
echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc  (1 = a station or stack was skipped; read the list)"
exit $rc
