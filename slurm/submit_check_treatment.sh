#!/bin/bash
## Is the fixed-vs-dynamic LMA treatment actually doing anything?
##
##     sbatch slurm/submit_check_treatment.sh US-Ha2 US-HBK
##     sbatch slurm/submit_check_treatment.sh --all
##
## Run this BEFORE any analysis. The first full array produced dyn_lma runs that
## were bit-identical to fixed_lma at every station -- MAIN_FRAME_SLA updated Sl_H
## yearly, but VEGETATION_DYNAMIC reads Sl from VegH_Param_Dyn, which
## Restating_parameters filled once before the time loop. Both arms looked healthy
## in every single-run figure; the failure was only visible in the difference.
##
## Exits 1 if any station's arms are identical.

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
