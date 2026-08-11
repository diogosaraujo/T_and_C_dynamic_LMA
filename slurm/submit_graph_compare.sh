#!/bin/bash
## Fixed-LMA vs dynamic-LMA comparison figures for one station.
##
##     sbatch slurm/submit_graph_compare.sh US-HBK
##     sbatch --array=1-98%6 slurm/submit_graph_compare.sh    # every station in run_list
##
## GRAPH_MOD plots one run at a time, so the treatment -- which lives entirely in
## the DIFFERENCE between the arms -- is invisible in it. This draws the same
## quantities with both arms overlaid or side by side, plus a table of the
## parameters the run used.
##
## Writes to $MODEL_RUN/<STATION>/figures_compare/ (11 x 6.5 in PNGs).
## Needs BOTH arms to have a RES_*.mat; stations missing either are skipped.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH -p SOE_main
#SBATCH -J tc_cmp
#SBATCH -t 02:00:00
#SBATCH -o slurm/logs/tc_cmp_%A_%a.out
#SBATCH -e slurm/logs/tc_cmp_%A_%a.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
MODEL_RUN="${MODEL_RUN:-$(dirname "$TC_INPUT_DATA")/model_run}"

STATION="${1:-}"
if [ -z "$STATION" ]; then
    LIST="$MODEL_RUN/${RUN_LIST:-run_list.txt}"
    IDX="${SLURM_ARRAY_TASK_ID:-1}"
    [ -f "$LIST" ] || { echo "ERROR: no station given and no $LIST" >&2; exit 1; }
    # One entry per STATION here, not per arm: the comparison needs both arms.
    mapfile -t STATIONS < <(awk '{print $1}' "$LIST" | sort -u)
    if [ "$IDX" -gt "${#STATIONS[@]}" ]; then
        echo "task $IDX is beyond the ${#STATIONS[@]} stations -- nothing to do"; exit 0
    fi
    STATION="${STATIONS[$((IDX-1))]}"
fi

SDIR="$MODEL_RUN/$STATION"
mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}  task ${SLURM_ARRAY_TASK_ID:-none}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "station    : $STATION"
echo

for arm in fixed_lma dyn_lma; do
    if ! ls "$SDIR/era5_land/$arm"/RES_*.mat >/dev/null 2>&1; then
        echo "no RES_*.mat for $arm -- both arms are required, nothing to compare"
        exit 0
    fi
done

# LMOD is only initialised for shells that read .bashrc, so source it defensively.
if ! command -v module >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    [ -f /opt/apps/lmod/lmod/init/profile ] && source /opt/apps/lmod/lmod/init/profile
fi
module load Matlab/2025a 2>/dev/null || ml Matlab/2025a 2>/dev/null || true
command -v matlab >/dev/null 2>&1 || { echo "ERROR: matlab not on PATH" >&2; exit 1; }

# Figure creation is node-dependent here (works on soeepyc05, dies at the first
# figure() on soeepyc06/07). QT_PLATFORM lets that be set without editing this
# file; see slurm/submit_figures.sh for the evidence.
if [ -n "${QT_PLATFORM:-}" ]; then
    export QT_QPA_PLATFORM="$QT_PLATFORM"
    export MW_QT_PLATFORM="$QT_PLATFORM"
    echo "Qt platform : $QT_PLATFORM"
fi

OUT="$SDIR/figures_compare"
before=$(ls "$OUT"/*.png 2>/dev/null | wc -l)

matlab -nodisplay -nosplash -batch \
    "addpath('$REPO_ROOT/preprocessing'); exit(GRAPH_COMPARE('$SDIR'))"
status=$?

after=$(ls "$OUT"/*.png 2>/dev/null | wc -l)
# Judge on what this run produced, never on what happens to be in the directory.
if [ "$status" != "0" ] && [ "$after" -gt "$before" ]; then
    echo "NOTE: matlab exited $status but the PNG count rose $before -> $after;"
    echo "      treating as success."
    status=0
fi

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
echo "png files   : $after  (was $before)"
exit $status
