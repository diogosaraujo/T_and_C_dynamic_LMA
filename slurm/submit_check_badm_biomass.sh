#!/bin/bash
## How many stations report biomass in AmeriFlux BADM, and in what units?
##
##     sbatch slurm/submit_check_badm_biomass.sh
##
## Answers whether the T&C heartwood pool B(6) can be set from site data instead
## of transplanted from US_xRM. Vegetation_Structural_Attributes.m computes
##     TBio = 0.02*(B(1)+B(2)+B(3)+B(4)+B(6))    [ton DM/ha]
## and TBio drives Allocation_Coefficients, so with B(6)=0 a mature forest is
## presented to the model as a ~21 t DM/ha sapling. If BADM reports aboveground
## biomass we can set B(6) = max(0, 50*TBio_target - (B1+B2+B3+B4)) per station.
##
## Reads only what download_ameriflux.py already fetched -- no network access.
## Writes badm_biomass.csv next to the other coverage reports.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J badm_bio
#SBATCH -t 00:30:00
#SBATCH -o slurm/logs/badm_bio_%j.out
#SBATCH -e slurm/logs/badm_bio_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

BADM_DIR="${BADM_DIR:-$TC_INPUT_DATA/ameriflux}"
OUT_CSV="${OUT_CSV:-$TC_INPUT_DATA/badm_biomass.csv}"

mkdir -p "$REPO_ROOT/slurm/logs"

echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "badm dir   : $BADM_DIR$([ -d "$BADM_DIR" ] || echo '   <-- NOT FOUND')"
echo "output     : $OUT_CSV"
echo

if [ ! -d "$BADM_DIR" ]; then
    echo "ERROR: $BADM_DIR not found -- has download_ameriflux.py run?" >&2
    exit 1
fi

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"

cd "$REPO_ROOT/preprocessing"
python check_badm_biomass.py --dir "$BADM_DIR" --csv "$OUT_CSV" "$@"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
exit $status
