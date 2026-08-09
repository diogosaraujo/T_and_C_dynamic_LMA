#!/bin/bash
## Parse and lint the MATLAB that build_model_run.py generates, before running it.
##
##     sbatch slurm/submit_check_generated.sh
##
## Two failures so far were generated MATLAB that could not run -- an unterminated
## character vector (job 35691) and a variable referenced ten lines before it was
## defined (35696). Both cost a full MATLAB startup to discover, roughly three to
## six minutes each. This finds them in seconds.
##
## Reads : $TC_ROOT/model_run/*/era5_land/*/{GO_*.m,MOD_PARAM_*.m}
## Writes: nothing. Exit 1 if any file fails to parse or uses an undefined name.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J check_gen
#SBATCH -t 00:30:00
#SBATCH -o slurm/logs/check_gen_%j.out
#SBATCH -e slurm/logs/check_gen_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
MODEL_RUN="${MODEL_RUN:-$(dirname "$TC_INPUT_DATA")/model_run}"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "model_run  : $MODEL_RUN"
echo

if ! command -v module >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    [ -f /opt/apps/lmod/lmod/init/profile ] && source /opt/apps/lmod/lmod/init/profile
fi
module load Matlab/2025a 2>/dev/null || ml Matlab/2025a 2>/dev/null || true
command -v matlab >/dev/null 2>&1 || { echo "ERROR: matlab not on PATH" >&2; exit 1; }

cd "$REPO_ROOT/preprocessing"
matlab -nodisplay -nosplash -batch "exit(check_generated('$MODEL_RUN'))"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status  (1 = a generated file will not run)"
exit $status
