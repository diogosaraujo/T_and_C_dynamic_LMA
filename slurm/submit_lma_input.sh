#!/bin/bash
## Build the per-station time-varying LMA series from the PLSR temporal-CV output.
##
##     sbatch slurm/submit_lma_input.sh
##
## Targeted runs -- flags pass straight through:
##     sbatch slurm/submit_lma_input.sh --stations US-HBK,US-Ha2
##     sbatch slurm/submit_lma_input.sh --dry-run
##     sbatch slurm/submit_lma_input.sh --fill-gaps
##
## Reads (already on the cluster, nothing is downloaded):
##     $PLSR_ROOT/eco<ii>/time/PLSR_predictions_eco<ii>_<forest>_oofcv.mat
##     $ECOREGION_ROOT/ecoregion_no<ii>.csv
##
## Writes one subfolder per station, outside the repo:
##     $TC_INPUT_DATA/lma/<station>/<station>_LMA_observed.csv
##     $TC_INPUT_DATA/lma/<station>/<station>_LMA_modelled.csv
##     $TC_INPUT_DATA/lma/lma_manifest.csv
##
## LMA is written in g/m2. The SLA conversion (Sl = 1/(LMA*f_C)) is applied when
## the T&C .mat inputs are built, so changing f_C never means re-running this.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J lma_input
#SBATCH -t 02:00:00
#SBATCH -o slurm/logs/lma_input_%j.out
#SBATCH -e slurm/logs/lma_input_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

ECOREGION_ROOT="${ECOREGION_ROOT:-/vol_efthymios/NFS07/dd1136/ecoregions}"
PLSR_ROOT="${PLSR_ROOT:-$ECOREGION_ROOT/PLSR_temporal_cv_pixel_climatology_DOY/LMA}"
OUT_DIR="${OUT_DIR:-$TC_INPUT_DATA/lma}"

mkdir -p "$REPO_ROOT/slurm/logs" "$OUT_DIR"

echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "plsr root  : $PLSR_ROOT"
echo "eco root   : $ECOREGION_ROOT"
echo "output     : $OUT_DIR"
echo

if [ ! -d "$TC_VENV" ]; then
    echo "ERROR: venv not found at $TC_VENV -- run 'sbatch slurm/submit_setup_env.sh' first" >&2
    exit 1
fi

# Fail here rather than after reading 28 prediction files.
for d in "$PLSR_ROOT" "$ECOREGION_ROOT"; do
    if [ ! -d "$d" ]; then
        echo "ERROR: not a directory: $d" >&2
        echo "       override with PLSR_ROOT=... / ECOREGION_ROOT=... before sbatch" >&2
        exit 1
    fi
done

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"

# The PLSR outputs are MATLAB -v7.3, i.e. HDF5; scipy.io.loadmat cannot read them.
if ! python -c "import h5py" 2>/dev/null; then
    echo "ERROR: h5py is missing -- re-run 'sbatch slurm/submit_setup_env.sh'" >&2
    exit 1
fi

cd "$REPO_ROOT/preprocessing"

python build_lma_input.py \
    --plsr-root "$PLSR_ROOT" \
    --ecoregion-root "$ECOREGION_ROOT" \
    --out "$OUT_DIR" "$@"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
du -sh "$OUT_DIR" 2>/dev/null || true
exit $status
