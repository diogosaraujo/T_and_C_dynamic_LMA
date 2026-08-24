#!/bin/bash
## Label each station-year drought or normal, from SPEI.
##
##     sbatch slurm/submit_classify_drought.sh --out $TC_INPUT_DATA/drought_years.csv
##     sbatch slurm/submit_classify_drought.sh --out ... --months 5,6,7,8,9
##     sbatch slurm/submit_classify_drought.sh --out ... --percentile 20
##
## Produces the station,year,class CSV that submit_daily_effect.sh --drought
## consumes, with the SPEI value kept alongside so a composite can be re-cut at a
## different threshold without re-reading the stacks.
##
## SPEI, not SPI: SPEI is precipitation minus potential evapotranspiration, so a
## hot year with normal rainfall counts as a drought. That is exactly where the
## treatment is expected to act -- dynamic LMA sheds leaf area the fixed arm
## cannot, so the arms should diverge most when evaporative demand is high, and
## SPI cannot see those years at all.
##
## The classification is deliberately INDEPENDENT of the model. Deriving drought
## from the runs' own precipitation would make "the effect is larger in dry years"
## circular; an external index makes it a real statement.
##
## WHICH CUT. --threshold -1.0 is comparable across stations but leaves wet sites
## contributing no drought years at all, so they drop out of the composite.
## --percentile 20 guarantees every station contributes, at the cost of "drought"
## meaning something different at each. Run both -- they answer different
## questions and share one schema, so a figure can be redone either way.
##
## Note that annual means are smoother than monthly values: P(monthly SPEI < -1)
## is about 16% by construction, but averaging twelve months shrinks the
## variance, so -1.0 on an annual mean is a rarer and more severe event -- 9-12%
## of years at the stations checked so far.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -p SOE_main
#SBATCH -J drought_c
#SBATCH -t 02:00:00
#SBATCH -o slurm/logs/drought_c_%j.out
#SBATCH -e slurm/logs/drought_c_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python -u classify_drought.py "$@"
rc=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc  (1 = a station was skipped; read the list above)"
exit $rc
