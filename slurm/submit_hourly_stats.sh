#!/bin/bash
## Model-versus-tower statistics at HOURLY resolution, per station and arm.
##
##     source slurm/config.sh
##     sbatch -p SOE_legacy -A efthymios slurm/submit_hourly_stats.sh --gpp NT
##
## Writes hourly_stats.csv to $TC_RESULTS:
##   station,pft,variable,arm,n,mean_obs,sd_obs,mean_mod,sd_mod,
##   rmse,rsr,bias,r,skill_score
##
## WHY HOURLY. The annual comparison rests on the overlap between the model
## record (1985-2020) and the tower record, and at many sites that is nearly
## nothing: US-HBK runs 2016-2024, so the annual overlap is TWO YEARS and a
## correlation from two points is +/-1 by construction. At n=3 no correlation is
## significant, at n=10 the threshold is |r|=0.63, and only past n~30 does
## |r|=0.36 become distinguishable from zero. Hourly gives ~40,000 matched steps
## at US-HBK, so the sampling problem disappears.
##
## TIME ZONE IS HANDLED. AmeriFlux is LOCAL STANDARD TIME, T&C runs UTC. The
## per-site offset comes from the FLUXNET BIF, not from longitude. Verified at
## US-HBK: UTC_OFFSET=-5 and peak GPP lands at UTC 16, i.e. 11:00 local, which
## is where a deciduous canopy peaks.
##
## CARBON UNITS ARE HANDLED. At HH the tower reports umol CO2 m-2 s-1, not the
## gC m-2 d-1 of its own daily files, and T&C's An_H + Rdark_H is the same unit
## (VEG_DYN_RES forms GPP = 1.0368*(An+Rdark) and 1.0368 = 12.0*86400*1e-6).
## No conversion. The DAILY GPP used elsewhere is a different quantity.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH -p SOE_main
#SBATCH -J hrly_stats
#SBATCH -t 08:00:00
#SBATCH -o slurm/logs/hrly_stats_%j.out
#SBATCH -e slurm/logs/hrly_stats_%j.err

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
echo "results    : $TC_RESULTS"
echo "towers     : $TOWER_DIR$([ -d "$TOWER_DIR" ] || echo '   <-- NOT FOUND')"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo
[ -d "$TOWER_DIR" ] || { echo "ERROR: no tower archives at $TOWER_DIR" >&2; exit 1; }
[ -d "$MODEL_RUN" ] || { echo "ERROR: no $MODEL_RUN" >&2; exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python -u hourly_stats.py --root "$MODEL_RUN" --tower-dir "$TOWER_DIR" "$@"
rc=$?
echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc  (1 = a station was skipped; read the list above)"
exit $rc
