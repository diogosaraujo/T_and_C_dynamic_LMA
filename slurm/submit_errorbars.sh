#!/bin/bash
## RMSE and skill score by forest type: 4x4 error-bar panels, ERA5-Land only.
##
##     source slurm/config.sh
##     sbatch -p SOE_legacy -A efthymios slurm/submit_errorbars.sh --step all --gpp NT
##
## Rows GPP/ET/H/LE; columns evergreen-all, evergreen-drought, deciduous-all,
## deciduous-drought. Three error bars per panel: RMSE dynamic, RMSE fixed, and
## the skill score, as mean +/- 1 SD across that type's stations with the
## individual stations drawn behind.
##
## RMSE is on the LEFT axis in the flux's own units, SS on the RIGHT axis and
## dimensionless. One scale per row so the four columns are comparable.
##
## ERA5-Land only: a GCM does not reproduce the weather a tower measured, so a
## station-by-station RMSE against observations means nothing for those runs.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -p SOE_main
#SBATCH -J errbars
#SBATCH -t 02:00:00
#SBATCH -o slurm/logs/errbars_%j.out
#SBATCH -e slurm/logs/errbars_%j.err

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
echo "figures    : $TC_FIGURES"
echo "towers     : $TOWER_DIR$([ -d "$TOWER_DIR" ] || echo '   <-- NOT FOUND')"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo
[ -d "$TOWER_DIR" ] || { echo "ERROR: no tower archives at $TOWER_DIR" >&2; exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python -u figure_errorbars.py \
    --model-dir "$TC_RESULTS" --tower-dir "$TOWER_DIR" "$@"
rc=$?
echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc"
exit $rc
