#!/bin/bash
## Is the fixed-vs-dynamic LMA treatment actually doing anything?
##
##     sbatch slurm/submit_check_treatment.sh US-Ha2 US-HBK
##     sbatch slurm/submit_check_treatment.sh --all
##     sbatch slurm/submit_check_treatment.sh US-Wrc --pair 'ssp585/*'
##
## Run this BEFORE any analysis. The first full array produced dyn_lma runs that
## were bit-identical to fixed_lma at every station -- MAIN_FRAME_SLA updated Sl_H
## yearly, but VEGETATION_DYNAMIC reads Sl from VegH_Param_Dyn, which
## Restating_parameters filled once before the time loop. Both arms looked healthy
## in every single-run figure; the failure was only visible in the difference.
##
## EVERY arm pair is checked -- <station>/era5_land AND <station>/<scenario>/<GCM>,
## 16 per fully-run station. Jobs 37691/37692 are why: the checker only knew the
## era5_land path, so it vetted the ERA5 pair, printed "the treatment is live", and
## never opened the 30 GCM arms it had been submitted to vet. --pair takes a glob
## over the pair label when 16 x 101 is more than you want to read.
##
## Exits 1 if any PAIR's arms are identical, or if a named station has no pair at
## all -- finding nothing is not a pass.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -p SOE_main
#SBATCH -J tc_treat
#SBATCH -t 00:30:00
#SBATCH -o slurm/logs/tc_treat_%j.out
#SBATCH -e slurm/logs/tc_treat_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
export MODEL_RUN="${MODEL_RUN:-$(dirname "$TC_INPUT_DATA")/model_run}"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "model_run  : $MODEL_RUN"
echo

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"

cd "$REPO_ROOT/preprocessing"
python check_treatment_effect.py --root "$MODEL_RUN" "$@"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status  (1 = at least one station has NO treatment effect)"
exit $status
