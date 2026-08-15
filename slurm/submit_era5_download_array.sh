#!/bin/bash
## ERA5-Land download as a SLURM job array -- for the GRIDDED fallback mode only.
##
##     sbatch slurm/submit_era5_download_array.sh --mode gridded
##
## Do NOT use this for the default timeseries mode. That is only 80 requests; the CDS
## throttles per user, so 8 tasks would just sit in the same ECMWF queue while holding 8
## SLURM allocations. Use submit_era5_download.sh instead.
##
## The gridded fallback is ~2,960 requests, where spreading tasks across the array does
## help absorb per-request overhead and lets each task finish inside its walltime.
##
## Each task takes a round-robin slice of the GRID CELLS (never split mid-cell, so the
## shared-cell copy always finds its source). Tasks write manifest.shard<NNN>.csv; merge
## them when the array finishes -- see slurm/README.md.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH -p SOE_main
#SBATCH -J era5_dl_arr
#SBATCH -t 7-00:00:00
#SBATCH --array=0-7%4                 # 8 shards, at most 4 running at once
#SBATCH -o slurm/logs/era5_dl_%A_%a.out
#SBATCH -e slurm/logs/era5_dl_%A_%a.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
OUT_DIR="${OUT_DIR:-$TC_INPUT_DATA/era5_land}"
JOBS="${JOBS:-2}"

SHARD="${SLURM_ARRAY_TASK_ID:-0}"
# Derive the shard count from the array bounds so it stays in sync with --array above.
NUM_SHARDS="${NUM_SHARDS:-$(( ${SLURM_ARRAY_TASK_MAX:-0} - ${SLURM_ARRAY_TASK_MIN:-0} + 1 ))}"

mkdir -p "$REPO_ROOT/slurm/logs" "$OUT_DIR"

echo "array task : $SHARD of $NUM_SHARDS"
echo "node       : $(hostname)"
tc_check_partition
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo

if [ ! -d "$TC_VENV" ]; then
    echo "ERROR: venv not found at $TC_VENV -- run 'bash slurm/setup_env.sh' first" >&2
    exit 1
fi

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"
cd "$REPO_ROOT/preprocessing"

python download_era5_land.py \
    --out "$OUT_DIR" \
    --jobs "$JOBS" \
    --shard "$SHARD" \
    --num-shards "$NUM_SHARDS" \
    "$@"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
exit $status
