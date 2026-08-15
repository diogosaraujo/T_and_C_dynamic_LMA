#!/bin/bash
## Verify the ERA5-Land predictor port before any fallback value is trusted.
##
##     sbatch slurm/submit_verify_era5_predictors.sh              # eco1
##     sbatch slurm/submit_verify_era5_predictors.sh 7 49 75      # several ecoregions
##
## build_lma_input.py --reconstruct fills years the preprocessed table never
## emitted (it only keeps forest pixel-years) by recomputing the predictors from
## the ERA5-Land monthly stacks. That recomputation is a port of
## build_climate_predictor_row from PLSR_PREPROCESS_PIXEL_CLIM_DOY_CORE.m, and a
## port can be wrong in ways that still produce plausible numbers -- a transposed
## grid, a shifted time axis, an unflipped PET sign.
##
## This job recomputes predictors for rows that ARE in the table and compares them
## against the stored values. A pass means the port reproduces the MATLAB on data
## the MATLAB itself produced; a failure names the predictor families involved.
##
## Reads:
##     $PREDICTOR_ROOT/LMA_ecoregion_no<ii>.csv     (reference values)
##     $ERA5_MONTHLY/                               (the stacks to recompute from)
## Writes: nothing.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -p SOE_main
#SBATCH -J era5_pred_verify
#SBATCH -t 01:00:00
#SBATCH -o slurm/logs/era5_pred_verify_%j.out
#SBATCH -e slurm/logs/era5_pred_verify_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

ECOREGION_ROOT="${ECOREGION_ROOT:-/vol_efthymios/NFS07/dd1136/ecoregions}"
PREDICTOR_ROOT="${PREDICTOR_ROOT:-$ECOREGION_ROOT/PLSR_inputs_pixel_climatology_DOY/LMA}"
ERA5_MONTHLY="${ERA5_MONTHLY:-/vol_efthymios/NFS07/Data/ERA5_Land/monthly}"
ROWS="${ROWS:-25}"
TOL="${TOL:-1e-6}"

ECOS=("$@")
if [ ${#ECOS[@]} -eq 0 ]; then ECOS=(1); fi

mkdir -p "$REPO_ROOT/slurm/logs"

echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "pred root  : $PREDICTOR_ROOT"
echo "era5       : $ERA5_MONTHLY"
echo "ecoregions : ${ECOS[*]}   (rows=$ROWS tol=$TOL)"
echo

if [ ! -d "$ERA5_MONTHLY" ]; then
    echo "ERROR: not a directory: $ERA5_MONTHLY" >&2
    echo "       override with ERA5_MONTHLY=... before sbatch" >&2
    exit 1
fi
ls -la "$ERA5_MONTHLY" | head -20
echo

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"
cd "$REPO_ROOT/preprocessing"

# The name port is checkable without touching the stacks -- do it first, since a
# column-order mismatch would invalidate every numeric comparison below.
python era5_predictors.py names --check "$PREDICTOR_ROOT/LMA_ecoregion_no${ECOS[0]}.csv"
name_status=$?
echo

status=0
for ii in "${ECOS[@]}"; do
    tab="$PREDICTOR_ROOT/LMA_ecoregion_no${ii}.csv"
    echo "=== ecoregion $ii ==="
    if [ ! -f "$tab" ]; then
        echo "  missing: $tab" >&2
        status=1
        continue
    fi
    python era5_predictors.py verify "$tab" \
        --era5-root "$ERA5_MONTHLY" --rows "$ROWS" --tol "$TOL" || status=1
    echo
done

echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "name check : $name_status (0 = ok)"
echo "exit status: $status"
exit $(( status || name_status ))
