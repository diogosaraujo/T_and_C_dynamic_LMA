#!/bin/bash
## Download hourly ERA5-Land forcing for all AmeriFlux stations, 1985-2021.
##
## Submit from the repo root:
##     sbatch slurm/submit_era5_download.sh
##
## Monitor:
##     squeue -u $USER
##     tail -f slurm/logs/era5_dl_<jobid>.out
##
## This is a NETWORK-bound job, not a compute-bound one: 80 requests that queue
## server-side at ECMWF. One node, a handful of threads, and a long walltime is the right
## shape. Do not scale it up with more CPUs -- the CDS throttles per user, not per core.

#SBATCH -N 1                          # single node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4             # one per concurrent CDS request
#SBATCH --mem=8G
#SBATCH -p SOE_main                   # new Epyc nodes; SOE_legacy also fine for this
#SBATCH -J era5_dl
#SBATCH -t 7-00:00:00                 # default is 3 days; CDS queueing can exceed that
#SBATCH -o slurm/logs/era5_dl_%j.out
#SBATCH -e slurm/logs/era5_dl_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
JOBS="${JOBS:-4}"
OUT_DIR="${OUT_DIR:-$TC_INPUT_DATA/era5_land}"

mkdir -p "$REPO_ROOT/slurm/logs" "$OUT_DIR"

echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "repo       : $REPO_ROOT"
echo "output     : $OUT_DIR"
echo

if [ ! -d "$TC_VENV" ]; then
    echo "ERROR: venv not found at $TC_VENV -- run 'bash slurm/setup_env.sh' first" >&2
    exit 1
fi
if [ ! -f "$HOME/.cdsapirc" ]; then
    echo "ERROR: ~/.cdsapirc not found -- see slurm/setup_env.sh" >&2
    exit 1
fi

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"

# Output goes straight to the shared filesystem rather than node-local /tmp. The usual
# advice is the opposite, but this job is network-bound and resumable: writing in place
# means a walltime kill loses only the file in flight, and a resubmit picks up where it
# stopped. Staging to /tmp would throw away days of downloads if the job were killed
# before the copy-back step.
cd "$REPO_ROOT/preprocessing"

python download_era5_land.py \
    --out "$OUT_DIR" \
    --jobs "$JOBS" \
    "$@"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
if [ $status -ne 0 ]; then
    echo "Some cells failed. Resubmit the same script -- completed files are skipped."
fi
du -sh "$OUT_DIR" 2>/dev/null || true
exit $status
