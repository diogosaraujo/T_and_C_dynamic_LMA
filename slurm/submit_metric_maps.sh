#!/bin/bash
## Metric maps: where in CONUS each LMA metric is large.
##
##     source slurm/config.sh
##     sbatch -p SOE_legacy -A efthymios slurm/submit_metric_maps.sh
##     sbatch -p SOE_legacy -A efthymios slurm/submit_metric_maps.sh --metrics lma,flux
##
## Five figures -- map_{lma,flux,sd,spei,ta}.png -- each a 6x4 grid: rows are
## GPP/LAI/TR/ET/LE/H, columns are ERA5-Land, GCM historical, ssp126, ssp585.
## PuOr diverging, white at zero, one colour scale per row so the four columns
## can be read against each other. Marker size is the station's mean annual flux
## in GLOBAL quartiles -- edges from all four datasets pooled, so a large marker
## means the same flux in every column.
##
## Also writes map_stations.csv, the per-point table behind the figures
## (table_map_stations). It is emitted by the figure script rather than by
## figure_tables.py so the table and the figures cannot disagree.
##
## ONLY THE 82 STATIONS COMPLETE IN ALL FOUR DATASETS ARE DRAWN. The datasets
## cover 92/85/82/82, and on a map an unequal sample reads as a geographic
## pattern: the ERA5 column would carry ten dots the ssp columns lack. Pass
## --fleet all to keep them, which prints the per-dataset counts.
##
## Reads station_metrics.csv and station_sensitivity.csv from $TC_RESULTS, so
## submit_station_metrics.sh must have run first; the flux metric additionally
## reads the annual effect tables.
##
## --basemap is optional and takes the us_eco_l3 shapefile already downloaded by
## submit_verify_pairing.sh. Without it the panels plot but have no CONUS
## outline, which leaves the stations floating and is much harder to read --
## pass it if the shapefile is on disk.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -p SOE_legacy
#SBATCH -J mapfig
#SBATCH -t 02:00:00
#SBATCH -o slurm/logs/mapfig_%j.out
#SBATCH -e slurm/logs/mapfig_%j.err

set -uo pipefail
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
tc_check_args "$@" || exit 1
mkdir -p "$TC_FIGURES" || { echo "ERROR: cannot create $TC_FIGURES" >&2; exit 1; }
echo "results    : $TC_RESULTS"
echo "figures    : $TC_FIGURES"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python -u figure_metric_maps.py "$@"
rc=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc"
exit $rc
