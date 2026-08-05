#!/bin/bash
## Verify the downloaded ERA5-Land station files.
##
## Submit from the repo root, after the download job finishes:
##     sbatch slurm/submit_era5_verify.sh
##
## Or just run it directly on the login node -- it is a few minutes of reading, no
## network:
##     source ~/envs/tc-preproc/bin/activate
##     python preprocessing/verify_era5_land.py
##
## Exits non-zero if any station fails, so the log tail tells you whether the data is
## safe to feed the forcing builder.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J era5_verify
#SBATCH -t 04:00:00
#SBATCH -o slurm/logs/era5_verify_%j.out
#SBATCH -e slurm/logs/era5_verify_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
DATA_DIR="${DATA_DIR:-$TC_INPUT_DATA/era5_land}"
REPORT="${REPORT:-$TC_INPUT_DATA/era5_land_verification.json}"

mkdir -p "$REPO_ROOT/slurm/logs"

echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "data       : $DATA_DIR"
echo

if [ ! -d "$TC_VENV" ]; then
    echo "ERROR: venv not found at $TC_VENV -- run 'bash slurm/setup_env.sh' first" >&2
    exit 1
fi

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"
cd "$REPO_ROOT/preprocessing"

python verify_era5_land.py --dir "$DATA_DIR" --report "$REPORT" "$@"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
if [ $status -ne 0 ]; then
    echo "At least one station FAILED -- do not feed this data to the forcing builder yet."
fi
exit $status
