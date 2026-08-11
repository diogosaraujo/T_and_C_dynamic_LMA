#!/bin/bash
## Draw GRAPH_MOD figures from runs that have already finished -- no re-simulation.
##
##     sbatch slurm/submit_figures.sh US-Ha2 dyn_lma          # one arm
##     sbatch --array=1-196%10 slurm/submit_figures.sh        # every arm in run_list.txt
##     FORCE=1 sbatch --array=1-196%10 slurm/submit_figures.sh   # redraw existing figures
##
## GRAPH_MOD's inputdlg calls threw under 'matlab -batch', so every completed run
## saved its RES_*.mat but no figures. The results are fine; only the plotting has
## to be repeated, and RES_<ST>.mat holds the entire saved workspace, so the
## figures can be drawn from it directly at ~2 min per arm instead of repeating a
## ~0.5 h simulation.
##
## Arms with no RES_*.mat (still failing, never run) exit 0 without doing
## anything, so the full 1-196 range can be submitted without filtering first.
## Arms that already have PNGs are skipped unless FORCE=1.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -p SOE_main
#SBATCH -J tc_figs
#SBATCH -t 02:00:00
#SBATCH -o slurm/logs/tc_figs_%A_%a.out
#SBATCH -e slurm/logs/tc_figs_%A_%a.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
MODEL_RUN="${MODEL_RUN:-$(dirname "$TC_INPUT_DATA")/model_run}"

STATION="${1:-}"
ARM="${2:-}"
if [ -z "$STATION" ]; then
    LIST="$MODEL_RUN/${RUN_LIST:-run_list.txt}"
    IDX="${SLURM_ARRAY_TASK_ID:-1}"
    if [ ! -f "$LIST" ]; then
        echo "ERROR: no station given and no $LIST for the array form" >&2
        exit 1
    fi
    NLINES=$(wc -l < "$LIST")
    if [ "$IDX" -gt "$NLINES" ]; then
        echo "task $IDX is beyond the $NLINES entries in $LIST -- nothing to do"
        exit 0
    fi
    read -r STATION ARM < <(sed -n "${IDX}p" "$LIST")
fi
if [ -z "${STATION:-}" ] || [ -z "${ARM:-}" ]; then
    echo "ERROR: usage: $0 <station> <fixed_lma|dyn_lma>" >&2
    exit 1
fi

RUNDIR="$MODEL_RUN/$STATION/era5_land/$ARM"

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
if ! ls "$RUNDIR"/RES_*.mat >/dev/null 2>&1; then
    echo "no RES_*.mat -- this arm has not produced results yet, nothing to plot"
    exit 0
fi

# model_run/GRAPH_MOD.m is a COPY of the template, refreshed by
# submit_build_model_run.sh. Editing the repo and going straight to the figures
# runs the old script and fails identically to before the fix (jobs 36111-36114),
# which reads as "the fix did not work" rather than "the fix was not deployed".
# Compare instead of trusting: this is cheap and the failure mode is expensive.
TEMPLATE="${TEMPLATE:-$REPO_ROOT/T&C/Thanos_US_xRM}"
if [ -f "$TEMPLATE/GRAPH_MOD.m" ] && [ -f "$MODEL_RUN/GRAPH_MOD.m" ]; then
    if ! cmp -s "$TEMPLATE/GRAPH_MOD.m" "$MODEL_RUN/GRAPH_MOD.m"; then
        echo "ERROR: $MODEL_RUN/GRAPH_MOD.m differs from the repo template." >&2
        echo "       The run tree has a stale copy; figures would be drawn by the" >&2
        echo "       old script (or fail the way it used to). Refresh it first:" >&2
        echo "         sbatch slurm/submit_build_model_run.sh --stations all" >&2
        exit 1
    fi
elif [ ! -f "$MODEL_RUN/GRAPH_MOD.m" ]; then
    echo "ERROR: no $MODEL_RUN/GRAPH_MOD.m -- run submit_build_model_run.sh first" >&2
    exit 1
fi

# LMOD is only initialised for shells that read .bashrc, so source it defensively.
if ! command -v module >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    [ -f /opt/apps/lmod/lmod/init/profile ] && source /opt/apps/lmod/lmod/init/profile
fi
module load Matlab/2025a 2>/dev/null || ml Matlab/2025a 2>/dev/null || true
command -v matlab >/dev/null 2>&1 || { echo "ERROR: matlab not on PATH" >&2; exit 1; }

FORCE_ARG="false"
[ "${FORCE:-0}" = "1" ] && FORCE_ARG="true"

# Whether MATLAB can create a figure at all is NODE-DEPENDENT here: it works on
# soeepyc05 (jobs 36119-36122 wrote 14 PNGs in ~96 s) and dies at the first
# figure() on soeepyc06/07 with "no Qt platform plugin could be initialized"
# (36295/36296, 0 PNGs). QT_PLATFORM lets that be tested instead of guessed --
# the assertion itself lists linuxfb, minimal, offscreen, vnc, xcb.
#   QT_PLATFORM=offscreen  renders correctly but is very slow (>3.7 h in job 36261)
#   QT_PLATFORM=minimal    untested, lighter than offscreen
#   unset                  crashes on nodes without a usable plugin
if [ -n "${QT_PLATFORM:-}" ]; then
    export QT_QPA_PLATFORM="$QT_PLATFORM"
    export MW_QT_PLATFORM="$QT_PLATFORM"
    echo "Qt platform : $QT_PLATFORM"
fi

# Count BEFORE, so success can be judged on what this run produced.
npng_before=$(ls "$RUNDIR"/figures/*.png 2>/dev/null | wc -l)
echo "png before  : $npng_before"
echo

matlab -nodisplay -nosplash -batch \
    "addpath('$REPO_ROOT/preprocessing'); exit(make_figures('$RUNDIR', $FORCE_ARG))"
status=$?

npng=$(ls "$RUNDIR"/figures/*.png 2>/dev/null | wc -l)
# Judge on what THIS run produced. Keying the check on "PNGs exist" let 14
# leftovers from a run on another node mask a MATLAB that wrote nothing, and the
# job reported success (36240/36241/36292). Only an increase proves work happened.
if [ "$status" != "0" ] && [ "$status" != "2" ] && [ "$npng" -gt "$npng_before" ]; then
    echo "NOTE: matlab exited $status but the PNG count rose $npng_before -> $npng,"
    echo "      so figures were written; treating as success."
    status=0
elif [ "$status" != "0" ] && [ "$status" != "2" ]; then
    echo "FAILED: matlab exited $status and the PNG count did not rise"
    echo "        ($npng_before -> $npng). No figures were produced by this run."
fi

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status  (0 figures written, 2 nothing to do, 1 GRAPH_MOD failed)"
echo "png files   : $npng"
[ "$status" = "2" ] && exit 0
exit $status
