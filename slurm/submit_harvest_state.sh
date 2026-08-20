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
##   1. harvest  '*/era5_land/fixed_lma'       -> key era5_land/fixed_lma
##
##   2. GCM spin-up pass, seeded from it (disposable output):
##        build_gcm_model_run.py --arms spinup --scenario historical \
##            --ic initial_state.csv --ic-key era5_land/fixed_lma \
##            --require-era5-state
##      run it, then harvest '*/historical/*/spinup'
##                                              -> key historical/<GCM>/spinup
##
##   3. the real GCM HISTORICAL arms, both from that one baseline:
##        build_gcm_model_run.py --scenario historical \
##            --ic initial_state.csv --ic-key 'historical/{gcm}/spinup'
##      run them, then harvest '*/historical/*/*_lma'
##                                              -> key historical/<GCM>/<arm>
##
##   4. the SSP arms, each continuing from its OWN historical arm:
##        build_gcm_model_run.py --scenario ssp126 --scenario ssp585 \
##            --ic initial_state.csv --ic-key 'historical/{gcm}/{arm}'
##
##   5. the ERA5 restart arms, from step 1's state -- no spin-up pass, because
##      the forcing does not change and the first run already equilibrated:
##        build_model_run.py --arms fixed_lma_ic,dyn_lma_ic \
##            --ic initial_state.csv --ic-key era5_land/fixed_lma
##
## Step 2 exists because an ERA5-equilibrated state dropped into a GCM run spends
## its first years re-equilibrating to that model's own climate, and those years
## land inside 1985-2014 -- the validation window. The spin-up pass absorbs that
## transient in output nobody analyses. It costs one arm per (station, GCM),
## ~133 core-hours for the whole fleet.
##
## WHY STEP 3 IS ONE ROUND AND NOT TWO. The scenarios are identical over
## 1985-2014 and only diverge after 2015, so there is ONE historical run per
## (station, GCM, arm). ssp126 and ssp585 both read that same state -- the key in
## step 4 says 'historical/...' whatever scenario is being built -- and the
## historical period is never simulated twice.
##
## WHY STEP 3 USES 'spinup' FOR BOTH ARMS AND STEP 4 USES '{arm}'. The
## fixed-vs-dynamic contrast requires a COMMON baseline at the start of the
## experiment, which is 1985: both historical arms therefore begin from the same
## spin-up state, and the SLA treatment is the only difference between them. By
## 2015 the treatment has legitimately been running for 30 years, so each SSP arm
## continues from its own 2014 state and the historical+future pair is one
## continuous simulation rather than two with a jump at the join.
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
