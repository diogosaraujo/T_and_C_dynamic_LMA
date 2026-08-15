#!/bin/bash
## Fetch GEDI/Landsat canopy height (Potapov et al. 2021, UMD GLAD 30 m) for every
## station in the site lists.
##
## Submit from the repo root:
##     sbatch slurm/submit_canopy_height.sh
##
## Targeted runs -- every flag passes straight through to the Python script:
##     sbatch slurm/submit_canopy_height.sh --stations US-HBK,US-Ha2
##     sbatch slurm/submit_canopy_height.sh --radius 150
##     sbatch slurm/submit_canopy_height.sh --dry-run
##
## Output lands OUTSIDE the repo, in its own folder alongside the other input data:
##     $TC_INPUT_DATA/canopy_height/
##
## No merging with AmeriFlux happens here -- this just downloads. The selection rule
## (BADM HEIGHTC where present, GEDI as fallback) is applied later, when the .mat
## forcing/parameter files are built.
##
## Network-bound and short: the multi-GB mosaic is read over /vsicurl/, so only the
## scanlines covering each station are transferred. ~5 minutes for 118 stations.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J canopy_ht
#SBATCH -t 04:00:00
#SBATCH -o slurm/logs/canopy_height_%j.out
#SBATCH -e slurm/logs/canopy_height_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
OUT_DIR="${OUT_DIR:-$TC_INPUT_DATA/canopy_height}"
AMF_DIR="${AMF_DIR:-$TC_INPUT_DATA/ameriflux}"

mkdir -p "$REPO_ROOT/slurm/logs" "$OUT_DIR"

echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "output     : $OUT_DIR"
echo

if [ ! -d "$TC_VENV" ]; then
    echo "ERROR: venv not found at $TC_VENV -- run 'sbatch slurm/submit_setup_env.sh' first" >&2
    exit 1
fi

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"
if ! python -c "import rasterio" 2>/dev/null; then
    echo "ERROR: rasterio is not installed in $TC_VENV." >&2
    echo "       It was added to requirements.txt for this step -- re-run" >&2
    echo "       'sbatch slurm/submit_setup_env.sh' to install it." >&2
    exit 1
fi

cd "$REPO_ROOT/preprocessing"

# Compare against measured BADM HEIGHTC when the AmeriFlux download is present. This
# only annotates the output with a validation column; it does not merge the values.
extra=()
if [ -d "$AMF_DIR" ] && [ -f "$AMF_DIR/badm_values.csv" ]; then
    echo "cross-checking against BADM HEIGHTC in $AMF_DIR"
    extra+=(--ameriflux-dir "$AMF_DIR")
else
    echo "note: no BADM at $AMF_DIR - skipping the validation column"
fi
echo

# ${extra[@]+...} guards against "unbound variable" when the array is empty on
# bash < 4.4, which would fail at run time on the compute node rather than here.
python fetch_canopy_height.py --out "$OUT_DIR" ${extra[@]+"${extra[@]}"} "$@"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
du -sh "$OUT_DIR" 2>/dev/null || true
exit $status
