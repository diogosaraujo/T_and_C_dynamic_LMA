#!/bin/bash
## Run T&C for one GCM station/scenario/arm, or a job array over many.
##
##   SOE:
##     sbatch slurm/submit_gcm_tc_run.sh US-Wrc ssp585 GFDL-ESM4 dyn_lma
##     sbatch --array=1-500 slurm/submit_gcm_tc_run.sh         # from run_list_gcm.txt
##
##   AMAREL -- the partition MUST be given on the command line, because sbatch
##   parses "#SBATCH -p" before the script body runs and it cannot read config.sh:
##     sbatch -p main --array=1-500 slurm/submit_gcm_tc_run.sh
##     OFFSET=500  sbatch -p main --array=1-500 slurm/submit_gcm_tc_run.sh
##     OFFSET=1000 sbatch -p main --array=1-500 slurm/submit_gcm_tc_run.sh
##     ...   2,760 entries = 6 chunks of 500 (MaxSubmitPU on 'main' is 500)
##   The script prints the partition it landed on and warns if it looks wrong for
##   the host, so a forgotten -p shows up in the first lines of the log.
##
## Separate from submit_tc_run.sh because the run list carries four fields
## ("<station> <scenario> <GCM> <arm>") against that script's two, and the
## directory is <STATION>/<scenario>/<GCM>/<arm> rather than <STATION>/era5_land/
## <arm>. The ERA5 script is left untouched: it produced every completed run and
## there is nothing to gain from making it polymorphic.
##
## SIZING. At the measured 0.0099 core-hours per simulated year, a 35-year
## historical run is ~0.35 h and an 86-year SSP run ~0.85 h. The default 8 h wall
## clock is ample; memory is the thing to watch, because an SSP run holds ~2.4x
## the hourly arrays of the 36-year ERA5 runs the 32 G request was validated on.
## Check a completed one with
##     sacct -j <jobid> --format=JobID,MaxRSS,Elapsed
## and raise --mem if it is close.
##
## MaxArraySize is 1001 on Amarel and MaxSubmitPU on the 'main' QOS is 500, so
## chunk the list into arrays of at most 500 and submit the next as one drains.
##
## STORAGE. Full hourly RES is kept deliberately, so the 2,760 runs come to about
## 1.59 TB against a 2 TB /scratch hard limit. That fits only if results are
## drained to the SOE HPC DURING the campaign rather than after it -- run the
## chunks in order and rsync each one off before launching the next.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH -p SOE_main
#SBATCH -J gcm_tc
#SBATCH -t 08:00:00
#SBATCH -o slurm/logs/gcm_tc_%A_%a.out
#SBATCH -e slurm/logs/gcm_tc_%A_%a.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
MODEL_RUN="${MODEL_RUN:-$(dirname "$TC_INPUT_DATA")/model_run}"

STATION="${1:-}"; SCEN="${2:-}"; GCM="${3:-}"; ARM="${4:-}"
if [ -z "$STATION" ]; then
    LIST="$MODEL_RUN/${RUN_LIST:-run_list_gcm.txt}"
    IDX="${SLURM_ARRAY_TASK_ID:-1}"
    # OFFSET lets a >1001-entry list be run in chunks despite MaxArraySize:
    #   OFFSET=1000 sbatch --array=1-1000 ...   runs entries 1001-2000
    IDX=$(( IDX + ${OFFSET:-0} ))
    [ -f "$LIST" ] || { echo "ERROR: no $LIST -- run submit_gcm_model_run.sh" >&2; exit 1; }
    NLINES=$(wc -l < "$LIST")
    if [ "$IDX" -gt "$NLINES" ]; then
        echo "task $IDX is beyond the $NLINES entries in $LIST -- nothing to do"
        exit 0
    fi
    read -r STATION SCEN GCM ARM < <(sed -n "${IDX}p" "$LIST")
fi
if [ -z "${STATION:-}" ] || [ -z "${SCEN:-}" ] || [ -z "${GCM:-}" ] || [ -z "${ARM:-}" ]; then
    echo "ERROR: usage: $0 <station> <scenario> <GCM> <fixed_lma|dyn_lma>" >&2
    exit 1
fi

RUNDIR="$MODEL_RUN/$STATION/$SCEN/$GCM/$ARM"
MNAME="${STATION//-/_}"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}  task ${SLURM_ARRAY_TASK_ID:-none}"
echo "node       : $(hostname)"
tc_check_partition
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "station    : $STATION   scenario: $SCEN   gcm: $GCM   arm: $ARM"
echo "rundir     : $RUNDIR"
echo

[ -d "$RUNDIR" ] || { echo "ERROR: $RUNDIR not found -- run submit_gcm_model_run.sh" >&2; exit 1; }
for f in "GO_$MNAME.m" "MOD_PARAM_$MNAME.m" "LMA_$MNAME.mat"; do
    [ -f "$RUNDIR/$f" ] || { echo "ERROR: missing $RUNDIR/$f" >&2; exit 1; }
done
# The forcing is one level UP, shared by fixed_lma and dyn_lma, and GO loads it as
# '../Meteo_*.mat'. It is not in the arm directory and must not be looked for there.
ls "$(dirname "$RUNDIR")"/Meteo_*.mat >/dev/null 2>&1 || {
    echo "ERROR: no Meteo_*.mat in $(dirname "$RUNDIR") -- the GCM forcing builder" >&2
    echo "       has not run, or the tree predates the move and needs" >&2
    echo "       preprocessing/migrate_forcing.py" >&2
    exit 1; }

tc_load_matlab || exit 1

cd "$RUNDIR" || exit 1
echo "matlab     : $(command -v matlab)"
echo

matlab -nodisplay -nosplash -batch "GO_$MNAME"
status=$?

echo
if ls "$RUNDIR"/RES_*.mat >/dev/null 2>&1; then
    echo "result     : $(ls -lh "$RUNDIR"/RES_*.mat | awk '{print $5, $9}')"
    # RES on disk is the real success marker: MATLAB can exit non-zero on a Qt
    # teardown after the science is already saved.
    status=0
else
    echo "result     : NO RES_*.mat -- the run did not finish"
    [ $status -eq 0 ] && status=1
fi
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
exit $status
