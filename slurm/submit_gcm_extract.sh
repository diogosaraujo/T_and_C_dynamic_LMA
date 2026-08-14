#!/bin/bash
## Stage 0: pull the station pixels out of the NEX-GDDP-CMIP6 store.
##
##     sbatch --array=1-105%20 slurm/submit_gcm_extract.sh    # 5 GCM x 3 scen x 7 var
##     sbatch slurm/submit_gcm_extract.sh --report            # what is on disk
##     sbatch slurm/submit_gcm_extract.sh --gcm GFDL-ESM4     # one model
##
## One array task = one (GCM, scenario, variable). Each opens every yearly global
## file ONCE and takes all ~92 station pixels from it. Doing it the other way --
## looping stations outside files -- would reopen each file 92 times and turn
## 7,245 reads into 640,000 on a shared filesystem.
##
## Reads : $NEXGDDP_ROOT/<GCM>/<scenario>/<var>/<var>_day_*.nc   (read-only, not ours)
## Writes: $TC_INPUT_DATA/gcm_stations/<GCM>/<scenario>/<var>.npz
##
## Already-extracted tasks are skipped, so a partial failure costs nothing to
## resume; --force re-extracts.
##
## Throttle with %N. This is I/O bound against /vol_efthymios, so more than ~20
## concurrent readers will not go faster and may slow everyone else down.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J gcm_extract
#SBATCH -t 24:00:00
#SBATCH -o slurm/logs/gcm_extract_%A_%a.out
#SBATCH -e slurm/logs/gcm_extract_%A_%a.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

export NEXGDDP_ROOT="${NEXGDDP_ROOT:-/vol_efthymios/NFS07/Data/CMIP6/NEXGDDP}"
OUT_DIR="${OUT_DIR:-$TC_INPUT_DATA/gcm_stations}"

mkdir -p "$REPO_ROOT/slurm/logs" "$OUT_DIR"
echo "job        : ${SLURM_JOB_ID:-interactive}${SLURM_ARRAY_TASK_ID:+  task $SLURM_ARRAY_TASK_ID}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "nexgddp    : $NEXGDDP_ROOT$([ -d "$NEXGDDP_ROOT" ] || echo '   <-- NOT FOUND')"
echo "output     : $OUT_DIR"
echo

[ -d "$NEXGDDP_ROOT" ] || { echo "ERROR: $NEXGDDP_ROOT not found" >&2; exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
python - <<'PY' || { echo "ERROR: netCDF4 not installed -- pip install netCDF4" >&2; exit 1; }
import netCDF4  # noqa: F401
PY

cd "$REPO_ROOT/preprocessing" || exit 1

# Add the array index unless the caller already chose what to work on. Passing a
# NON-selector flag (--force, --dry-run) must not suppress it: job 36991 ran
# "sbatch --array=1-105 ... --force", $# was 1 rather than 0, --index was never
# added, and all 105 tasks died in argparse one second in.
have_selector=0
for _a in "$@"; do
    case "$_a" in --index|--all|--report|--gcm|--scenario|--variable) have_selector=1 ;; esac
done
if [ "$have_selector" -eq 0 ] && [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    set -- "$@" --index "$SLURM_ARRAY_TASK_ID"
elif [ "$have_selector" -eq 0 ]; then
    echo "no selector and not an array job -- extracting everything serially"
    set -- "$@" --all
fi
echo "args       : $*"

python extract_gcm_stations.py --out "$OUT_DIR" "$@"
status=$?

echo
echo "extracted  : $(find "$OUT_DIR" -name '*.npz' 2>/dev/null | wc -l) of 105 npz files"
du -sh "$OUT_DIR" 2>/dev/null || true
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
exit $status
