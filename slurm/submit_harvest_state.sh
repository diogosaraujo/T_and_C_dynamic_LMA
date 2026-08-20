#!/bin/bash
## Take the final vegetation state out of finished runs, to restart from.
##
##     sbatch slurm/submit_harvest_state.sh --from '*/era5_land/fixed_lma' --drift
##     sbatch slurm/submit_harvest_state.sh --from '*/historical/*/spinup'
##
## Writes/updates $MODEL_RUN/initial_state.csv, one row per (station, key), where
## key is the harvested arm directory relative to the station. Rounds merge: a
## later harvest replaces its own rows and leaves the others alone.
##
## THE CHAIN (Dr. Paschalis, 2026-08-19: the evergreen sites start with too little
## leaf carbon, so LAI ramps up inside the analysis window)
##
##   1. harvest  '*/era5_land/fixed_lma'      -> key era5_land/fixed_lma
##   2. build+run the GCM spin-up pass, seeded from that state:
##        build_gcm_model_run.py --arms spinup --scenario historical \
##            --ic initial_state.csv --ic-key era5_land/fixed_lma \
##            --require-era5-state
##      then harvest '*/historical/*/spinup'  -> key historical/<GCM>/spinup
##   3. build the real GCM arms with --ic-key 'historical/{gcm}/spinup'
##   4. the ERA5 restart arms, from step 1's state:
##        build_model_run.py --arms fixed_lma_ic,dyn_lma_ic \
##            --ic initial_state.csv --ic-key era5_land/fixed_lma
##
## Step 2 exists because an ERA5-equilibrated state dropped into a GCM run spends
## its first years re-equilibrating to that model's own climate, and those years
## land inside 1985-2014 -- the validation window. The spin-up pass absorbs that
## transient in output nobody analyses. It costs one arm per (station, GCM),
## ~133 core-hours for the whole fleet.
##
## --drift reports how far the live pools (B1-B4) moved over the last decade. A
## few percent means equilibrated; tens of percent means the wood pools are still
## moving and another cycle is warranted before trusting the restart.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -p SOE_main
#SBATCH -J harvest
#SBATCH -t 04:00:00
#SBATCH -o slurm/logs/harvest_%j.out
#SBATCH -e slurm/logs/harvest_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "model_run  : $MODEL_RUN$([ -d "$MODEL_RUN" ] || echo '   <-- NOT FOUND')"
echo

[ -d "$MODEL_RUN" ] || { echo "ERROR: no $MODEL_RUN" >&2; exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python harvest_state.py --root "$MODEL_RUN" "$@"
rc=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc  (1 = nothing harvested, or a state was unusable)"
exit $rc
