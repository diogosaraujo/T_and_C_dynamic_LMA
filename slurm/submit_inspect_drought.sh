#!/bin/bash
## What the SPEI/SPI stacks contain at our stations -- a look before anything
## classifies a year as drought. Writes nothing.
##
##     sbatch slurm/submit_inspect_drought.sh
##     sbatch slurm/submit_inspect_drought.sh --stations US-Ha2,US-NR1 --show-years
##     sbatch slurm/submit_inspect_drought.sh --index SPEI6_ts --threshold -0.8
##
## Reads through era5_predictors.Era5Monthly, which already encodes the parts
## that are easy to get wrong and were derived from the MATLAB source: which
## stack is indexed [lat, lon] on the shifted -180..180 grid and which is
## [lon, lat] on the native 0..360 one, the 1800-column longitude split, and
## h5py handing back MATLAB's dimensions reversed.
##
## The statistics prove a number was read. What proves the pixel is the RIGHT
## one is whether the driest years look like the droughts that actually
## happened -- 2012 across the Midwest and Front Range, 2002 and 2012 in
## Colorado, 2002/2007/2012 in the Southwest. Check that before thresholding.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -p SOE_main
#SBATCH -J drought_q
#SBATCH -t 01:00:00
#SBATCH -o slurm/logs/drought_q_%j.out
#SBATCH -e slurm/logs/drought_q_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
tc_check_args "$@" || exit 1
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python -u inspect_drought.py "$@"
rc=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc  (1 = a stack is missing, or a station read all-NaN)"
exit $rc
