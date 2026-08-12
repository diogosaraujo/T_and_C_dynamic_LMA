#!/bin/bash
## Fleet-wide fixed-vs-dynamic LMA assessment: read every station's RES files and
## write a report we can read off the cluster.
##
##     sbatch --array=1-101%8 slurm/submit_lma_effect.sh          # extract, parallel
##     sbatch slurm/submit_lma_effect.sh --report                 # then synthesise
##
##     sbatch slurm/submit_lma_effect.sh --all                    # extract serially
##     sbatch slurm/submit_lma_effect.sh --station US-Ha2         # one station
##
## Two stages on purpose. Each RES_*.mat is ~300 MB and there are ~200 of them, so
## EXTRACT reduces each station to a small JSON of annual series in
## $TC_INPUT_DATA/lma_effect/ (~40 kB per station). REPORT then works entirely off
## those, which means the 60 GB read happens once and every later re-analysis --
## a new metric, a different grouping -- is seconds rather than hours.
##
## The array form is over STATIONS, not arms; each task reads both arms. A range
## wider than the station list is fine, extra tasks exit 0. Throttle with %N: the
## per-user limit is ~15 concurrent jobs, and these are I/O bound on beegfs, so
## %8 is plenty -- more just queues behind itself at the filesystem.
##
## Already-cached stations are skipped, so re-running after a partial failure
## costs nothing. Use --force to re-extract (needed after any re-run of the model).
##
## Outputs, all under $TC_INPUT_DATA/lma_effect/:
##     lma_effect_report.md            the document
##     lma_effect_metrics.csv          per station x flux: effect size, slope, r
##     lma_effect_input_quality.csv    per station: is the LMA series signal or noise
##     lma_effect_annual.csv           every year, both arms -- for re-plotting

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -p SOE_main
#SBATCH -J lma_effect
#SBATCH -t 04:00:00
#SBATCH -o slurm/logs/lma_effect_%A_%a.out
#SBATCH -e slurm/logs/lma_effect_%A_%a.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
MODEL_RUN="${MODEL_RUN:-$(dirname "$TC_INPUT_DATA")/model_run}"
CACHE="${CACHE:-$TC_INPUT_DATA/lma_effect}"

mkdir -p "$REPO_ROOT/slurm/logs" "$CACHE"

echo "job        : ${SLURM_JOB_ID:-interactive}${SLURM_ARRAY_TASK_ID:+  task $SLURM_ARRAY_TASK_ID}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "model_run  : $MODEL_RUN$([ -d "$MODEL_RUN" ] || echo '   <-- NOT FOUND')"
echo "cache      : $CACHE"
echo

[ -d "$MODEL_RUN" ] || { echo "ERROR: $MODEL_RUN not found" >&2; exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing -- run setup_env.sh" >&2; exit 1; }

cd "$REPO_ROOT/preprocessing" || exit 1

# No arguments + inside an array => extract this task's station. Anything else is
# passed straight through, so --report / --all / --station work unchanged.
if [ "$#" -eq 0 ] && [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    set -- --index "$SLURM_ARRAY_TASK_ID"
elif [ "$#" -eq 0 ]; then
    echo "no arguments and not an array job -- extracting every station serially"
    echo "(this reads ~60 GB; the array form is faster)"
    echo
    set -- --all
fi

python analyze_lma_effect.py --model-run "$MODEL_RUN" --cache "$CACHE" "$@"
status=$?

echo
echo "cache holds: $(ls -1 "$CACHE"/*.json 2>/dev/null | wc -l) station file(s)"
if [ -f "$CACHE/lma_effect_report.md" ]; then
    echo
    echo "report     : $CACHE/lma_effect_report.md"
    echo "  copy it down with:"
    echo "    scp -P 22 $USER@soenfs1.hpc.rutgers.edu:$CACHE/lma_effect_*.{md,csv} ."
fi
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
exit $status
