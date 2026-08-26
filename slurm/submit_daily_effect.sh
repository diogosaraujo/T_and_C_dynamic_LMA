#!/bin/bash
## Seasonal signature of the LMA treatment: dyn - fixed by day of year.
##
## Output goes to $TC_RESULTS -- <TC_ROOT>/result_summary, a sibling of
## model_run -- and a bare --out filename resolves there automatically. It used
## to resolve against the working directory, i.e. preprocessing/, which put a
## 108 MB table inside the repo; "git add -A" swept it into a commit and GitHub
## declined the push at its 100 MB hard limit.
##
##     sbatch slurm/submit_daily_effect.sh --pair 'era5_land:*_ic' \
##         --out era5_daily.csv --drought drought_years.csv
##     sbatch slurm/submit_daily_effect.sh --pair 'historical/*:*' \
##         --out hist_daily.csv --drought drought_years_gcm.csv
##
## --drought is an INPUT and is read from $TC_RESULTS too, since that is where
## submit_classify_drought.sh writes it. Pass an absolute path to override.
##
## Complements analyze_lma_effect.py rather than replacing it. That one reduces
## each run to ANNUAL series, which is right for the fleet-wide headline and the
## future trends. This one keeps the day of year, because an annual mean hides
## compensating changes -- +8% in spring against -4% in summer averages to +1%
## and looks like nothing happened.
##
## It is worth having because the mechanism is PHENOLOGICAL: LMA enters T&C only
## as Sl in LAI = Sl*B(1) and never touches Vmax, so the treatment acts on leaf
## area, and the effect should concentrate around leaf-out, peak LAI and
## senescence -- differently for deciduous and evergreen canopies.
##
## RUN IT ON THE ERA5 PAIRS FIRST. Those are the runs that can be compared with
## tower observations: a GCM does not reproduce actual weather, so 2003 in
## GFDL-ESM4 is not the real 2003 and the GCM runs are comparable to towers only
## climatologically, never year by year.
##
## Reads every hourly array in every RES it touches, so it is I/O heavy: budget
## roughly a minute per pair for the 36-year ERA5 runs and three times that for
## an 86-year ssp pair.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -p SOE_main
#SBATCH -J daily_fx
#SBATCH -t 08:00:00
#SBATCH -o slurm/logs/daily_fx_%j.out
#SBATCH -e slurm/logs/daily_fx_%j.err

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

python -u analyze_daily_effect.py --root "$MODEL_RUN" "$@"
rc=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc  (1 = at least one pair was skipped; read the list above)"
exit $rc
