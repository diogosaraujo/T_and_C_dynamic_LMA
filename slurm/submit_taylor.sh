#!/bin/bash
## Taylor diagrams: how each arm stands against the towers, station by station.
##
##     source slurm/config.sh
##     sbatch slurm/submit_taylor.sh --step all
##     sbatch slurm/submit_taylor.sh --step JJA --gpp DT
##
## Six 4x2 portrait figures -- annual, monthly, DJF, MAM, JJA, SON. Rows are
## GPP/ET/LE/H; the left column uses every step and the right only drought steps.
## Each panel is a Taylor diagram and each station contributes two markers:
##
##     shape   circle = deciduous, triangle = evergreen
##     fill    FILLED = fixed LMA, OUTLINED = dynamic LMA
##
## so the question is whether the outlined marker sits closer to REF than the
## filled one. Radius is sd(model)/sd(tower) and angle is the correlation, which
## makes distance from REF the centred RMS difference.
##
## Normalised by the observed sigma so REF is at 1.0 in every panel -- otherwise
## one panel would mix stations whose fluxes differ by an order of magnitude and
## the spread would be about site productivity rather than model skill.
##
## RUN IT TWICE, --gpp NT and --gpp DT. Tower GPP is not measured -- it is
## partitioned from net exchange by nighttime extrapolation or a daytime
## light-response fit, and the two methods disagree. Their spread IS the
## observational uncertainty on that row, and a GPP result that flips sign
## between them is not a result.
##
## Figures go to $TC_FIGURES (<TC_ROOT>/figures, a sibling of model_run), the
## tables come from $TC_RESULTS, and the tower archives from wherever
## submit_fluxnet_download.sh put them.
##
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -p SOE_main
#SBATCH -J taylor_fig
#SBATCH -t 02:00:00
#SBATCH -o slurm/logs/taylor_fig_%j.out
#SBATCH -e slurm/logs/taylor_fig_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

TOWER_DIR="${TOWER_DIR:-$REPO_ROOT/preprocessing/fluxnet}"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
tc_check_args "$@" || exit 1
mkdir -p "$TC_FIGURES" || { echo "ERROR: cannot create $TC_FIGURES" >&2; exit 1; }
echo "results    : $TC_RESULTS"
echo "figures    : $TC_FIGURES"
echo "towers     : $TOWER_DIR$([ -d "$TOWER_DIR" ] || echo '   <-- NOT FOUND')"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo

[ -d "$TOWER_DIR" ] || {
    echo "ERROR: no tower archives at $TOWER_DIR" >&2
    echo "       Run slurm/submit_fluxnet_download.sh --agree-policy first." >&2
    exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python -u figure_taylor.py \
    --model-dir "$TC_RESULTS" --tower-dir "$TOWER_DIR" "$@"
rc=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc"
exit $rc
