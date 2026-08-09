#!/bin/bash
## Build the T&C meteorological forcing from the downloaded ERA5-Land netCDFs.
##
##     sbatch slurm/submit_meteo.sh --dry-run
##     sbatch slurm/submit_meteo.sh                    # the two example stations
##     sbatch slurm/submit_meteo.sh --stations all
##
## Two stages, in one job:
##   1. build_meteo_input.py   netCDF -> Meteo_<ST>_raw.mat   (unit conversions)
##   2. finish_meteo.m         adds SAB/SAD/PAR/N via C_Automatic_Radiation_Partition
##                             -> Meteo_<ST>_1985_2020.mat
##
## Stage 2 reuses the existing MATLAB partition rather than porting several hundred
## lines of solar geometry and cloud physics to Python.
##
## Reads : $TC_INPUT_DATA/era5_land/<ST>/*.nc          the forcing
##         $TC_INPUT_DATA/ameriflux/badm_values.csv     Zbas (site elevation)
##         T&C/Diogo/Ca_Data.mat                        hourly CO2
## Writes: $TC_INPUT_DATA/meteo/Meteo_<ST>_raw.mat and Meteo_<ST>_1985_2020.mat
##
## Units are taken from the working Meteo_US_xRM_1985_2020.mat, not assumed. The
## one that matters: Pre is in MILLIBAR, and ERA5 sp is in Pa.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH -p SOE_main
#SBATCH -J meteo
#SBATCH -t 12:00:00
#SBATCH -o slurm/logs/meteo_%j.out
#SBATCH -e slurm/logs/meteo_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

ERA5_DIR="${ERA5_DIR:-$TC_INPUT_DATA/era5_land}"
AMERIFLUX_DIR="${AMERIFLUX_DIR:-$TC_INPUT_DATA/ameriflux}"
METEO_DIR="${METEO_DIR:-$TC_INPUT_DATA/meteo}"
PARTITION_DIR="${PARTITION_DIR:-$REPO_ROOT/T&C/Diogo}"
YEAR_TAG="${YEAR_TAG:-1985_2020}"

mkdir -p "$REPO_ROOT/slurm/logs" "$METEO_DIR"

echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "era5       : $ERA5_DIR$([ -d "$ERA5_DIR" ] || echo '   <-- NOT FOUND')"
echo "ameriflux  : $AMERIFLUX_DIR$([ -d "$AMERIFLUX_DIR" ] || echo '   <-- NOT FOUND (Zbas)')"
echo "meteo out  : $METEO_DIR"
echo "partition  : $PARTITION_DIR"
echo

[ -d "$ERA5_DIR" ] || { echo "ERROR: no ERA5-Land directory" >&2; exit 1; }
[ -f "$PARTITION_DIR/C_Automatic_Radiation_Partition.m" ] || {
    echo "ERROR: C_Automatic_Radiation_Partition.m not in $PARTITION_DIR" >&2; exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"

echo "=== stage 1: netCDF -> raw .mat ==="
cd "$REPO_ROOT/preprocessing"
python build_meteo_input.py --era5 "$ERA5_DIR" --out "$METEO_DIR" \
    --ameriflux "$AMERIFLUX_DIR" \
    --start-year "${YEAR_TAG%%_*}" --end-year "${YEAR_TAG##*_}" "$@"
s1=$?
echo "stage 1 exit: $s1"

# --dry-run writes nothing, so there is nothing for stage 2 to finish.
if [[ " $* " == *" --dry-run "* ]]; then
    echo "dry run -- stage 2 skipped"
    exit $s1
fi

echo
echo "=== stage 2: radiation partition (MATLAB) ==="
if ! command -v module >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    [ -f /opt/apps/lmod/lmod/init/profile ] && source /opt/apps/lmod/lmod/init/profile
fi
module load Matlab/2025a 2>/dev/null || ml Matlab/2025a 2>/dev/null || true
command -v matlab >/dev/null 2>&1 || { echo "ERROR: matlab not on PATH" >&2; exit 1; }

matlab -nodisplay -nosplash -batch \
    "finish_meteo('$METEO_DIR','$METEO_DIR','$PARTITION_DIR','$YEAR_TAG')"
s2=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: stage1=$s1 stage2=$s2"
ls -lh "$METEO_DIR"/Meteo_*_"$YEAR_TAG".mat 2>/dev/null | head
du -sh "$METEO_DIR" 2>/dev/null || true
exit $(( s1 || s2 ))
