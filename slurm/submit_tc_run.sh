#!/bin/bash
## Run T&C for one station/arm, or a job array over many.
##
##     sbatch slurm/submit_tc_run.sh US-Ha2 dyn_lma
##     sbatch slurm/submit_tc_run.sh US-HBK fixed_lma
##     sbatch --array=1-4 slurm/submit_tc_run.sh          # every arm in run_list.txt
##
## The array form reads $MODEL_RUN/run_list.txt, one "<station> <arm>" per line,
## which submit_build_model_run.sh writes. Build it by hand for a subset.
##
## Each job cds into $MODEL_RUN/<station>/era5_land/<arm> and runs GO_<ST>.m, so
## every relative path in the generated launcher (Code/, the .mat inputs, RES_*)
## resolves the same way it would interactively.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -p SOE_main
#SBATCH -J tc_run
#SBATCH -t 3-00:00:00
#SBATCH -o slurm/logs/tc_run_%A_%a.out
#SBATCH -e slurm/logs/tc_run_%A_%a.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
MODEL_RUN="${MODEL_RUN:-$(dirname "$TC_INPUT_DATA")/model_run}"

STATION="${1:-}"
ARM="${2:-}"
if [ -z "$STATION" ]; then
    LIST="$MODEL_RUN/run_list.txt"
    IDX="${SLURM_ARRAY_TASK_ID:-1}"
    if [ ! -f "$LIST" ]; then
        echo "ERROR: no station given and no $LIST for the array form" >&2
        exit 1
    fi
    read -r STATION ARM < <(sed -n "${IDX}p" "$LIST")
fi
if [ -z "${STATION:-}" ] || [ -z "${ARM:-}" ]; then
    echo "ERROR: usage: $0 <station> <fixed_lma|dyn_lma>" >&2
    exit 1
fi

RUNDIR="$MODEL_RUN/$STATION/era5_land/$ARM"
MNAME="${STATION//-/_}"          # US-Ha2 -> US_Ha2, the MATLAB-safe name

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}  task ${SLURM_ARRAY_TASK_ID:-none}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "station    : $STATION   arm: $ARM"
echo "rundir     : $RUNDIR"
echo

if [ ! -d "$RUNDIR" ]; then
    echo "ERROR: run directory not found -- has submit_build_model_run.sh run?" >&2
    exit 1
fi
for f in "GO_$MNAME.m" "MOD_PARAM_$MNAME.m" "LMA_$MNAME.mat"; do
    [ -f "$RUNDIR/$f" ] || { echo "ERROR: missing $RUNDIR/$f" >&2; exit 1; }
done
if ! ls "$RUNDIR"/Meteo_*.mat >/dev/null 2>&1; then
    echo "ERROR: no Meteo_*.mat in $RUNDIR -- the forcing builder has not run" >&2
    exit 1
fi

# LMOD is only initialised for shells that read .bashrc, so source it defensively.
if ! command -v module >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    [ -f /opt/apps/lmod/lmod/init/profile ] && source /opt/apps/lmod/lmod/init/profile
fi
module load Matlab/2025a 2>/dev/null || ml Matlab/2025a 2>/dev/null || true
command -v matlab >/dev/null 2>&1 || { echo "ERROR: matlab not on PATH" >&2; exit 1; }

cd "$RUNDIR"
echo "matlab     : $(command -v matlab)"
echo
matlab -nodisplay -nosplash -batch "GO_$MNAME"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
ls -lh "$RUNDIR"/RES_*.mat 2>/dev/null || echo "  (no RES_*.mat produced)"
exit $status
