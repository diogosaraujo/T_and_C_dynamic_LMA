#!/bin/bash
## Read every effect table and write one summary of what the experiment found.
##
##     source slurm/config.sh
##     sbatch -p SOE_legacy -A efthymios slurm/submit_summarise.sh
##     sbatch -p SOE_legacy -A efthymios slurm/submit_summarise.sh --freqs monthly,seasonal
##
## Walks daily/monthly/seasonal/annual for ERA5 and all three GCM scenarios in
## $TC_RESULTS and writes effect_summary.md plus effect_summary_per_station.csv
## beside them. Reports effect size, direction and station agreement, where in
## the year it peaks, the drought contrast and the deciduous/evergreen split.
##
## Sizes are the MEDIAN ACROSS STATIONS of rel_ann_pct. Median because one
## pathological station once turned a fleet doing ~1.2% into a printed 129.54%;
## rel_ann because rel_pct divides by a denominator that vanishes out of season
## and once reported leakage at 50465967%, which was a February denominator
## rather than an effect.
##
## Tables without a rel_ann_pct column are REFUSED and named, not summarised --
## they predate that fix and would produce numbers that look fine and are not.
## Exit 1 means at least one table was skipped; read the list.
##
## Memory: the SSP monthly tables are ~5.9M rows each and several are held at
## once, hence 32G rather than the 16G the other analysis jobs use.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -p SOE_main
#SBATCH -J summarise
#SBATCH -t 04:00:00
#SBATCH -o slurm/logs/summarise_%j.out
#SBATCH -e slurm/logs/summarise_%j.err

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

python -u summarise_effect.py "$@"
rc=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc  (1 = a table was skipped; read the list above)"
exit $rc
