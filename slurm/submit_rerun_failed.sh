#!/bin/bash
## Write $MODEL_RUN/rerun_list.txt -- every station/arm that has no RES_*.mat --
## so a fix can be retried without re-running the arms that already finished.
##
##     sbatch slurm/submit_rerun_failed.sh
##     # then, with the count it prints:
##     sbatch --array=1-N%10 --export=ALL,RUN_LIST=rerun_list.txt slurm/submit_tc_run.sh
##
## A finished run writes RES_*.mat as its last act, so its absence is the success
## marker. A run that died in MOD_PARAM or inside ode45 leaves the directory
## otherwise complete, which is why the check looks for the output rather than
## for the inputs.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -p SOE_main
#SBATCH -J rerun_list
#SBATCH -t 00:10:00
#SBATCH -o slurm/logs/rerun_list_%j.out
#SBATCH -e slurm/logs/rerun_list_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
MODEL_RUN="${MODEL_RUN:-$(dirname "$TC_INPUT_DATA")/model_run}"

mkdir -p "$REPO_ROOT/slurm/logs"
SRC="$MODEL_RUN/run_list.txt"
OUT="$MODEL_RUN/rerun_list.txt"

if [ ! -f "$SRC" ]; then
    echo "ERROR: $SRC not found -- has submit_build_model_run.sh run?" >&2
    exit 1
fi

echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "full list  : $SRC"
echo

: > "$OUT"
done_n=0
while read -r STATION ARM; do
    [ -n "${STATION:-}" ] || continue
    RUNDIR="$MODEL_RUN/$STATION/era5_land/$ARM"
    if ls "$RUNDIR"/RES_*.mat >/dev/null 2>&1; then
        done_n=$((done_n + 1))
    else
        printf '%s %s\n' "$STATION" "$ARM" >> "$OUT"
    fi
done < "$SRC"

todo_n=$(wc -l < "$OUT")
total=$((done_n + todo_n))

echo "$done_n of $total arm(s) already have a RES_*.mat"
echo "$todo_n arm(s) written to $OUT"
echo
if [ "$todo_n" -gt 0 ]; then
    echo "  by station:"
    awk '{print $1}' "$OUT" | sort | uniq -c | sort -rn | head -20
    echo
    echo "next:"
    echo "  sbatch --array=1-$todo_n%10 --export=ALL,RUN_LIST=rerun_list.txt slurm/submit_tc_run.sh"
else
    echo "nothing to re-run"
fi
echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
