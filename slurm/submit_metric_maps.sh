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
## THE CONUS OUTLINE IS AUTOMATIC. The wrapper looks for the us_eco_l3 shapefile
## that submit_verify_pairing.sh downloads to $TC_INPUT_DATA/ecoregions and
## passes it as --basemap; override with an explicit --basemap, or set ECO_DIR.
## If the file is missing the job still runs and says so loudly on stderr.
##
## It is drawn with pyshp + pyproj, both in requirements.txt. NOT geopandas,
## which requirements.txt deliberately excludes -- the earlier geopandas version
## of this could not draw an outline under any arguments, since a valid
## --basemap raised ImportError and an absent one printed a note.

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

# THE OUTLINE IS ON BY DEFAULT. It used to need an explicit --basemap, and a run
# without one produced perfectly good figures with no CONUS behind them and a
# single "note:" line to say so -- easy to miss in a log, and the figures look
# finished. Auto-detect the shapefile submit_verify_pairing.sh already downloads,
# and only fall back to no outline if it genuinely is not there.
# SEARCH SEVERAL LOCATIONS, not just the download target. The shapefile is not
# necessarily where submit_verify_pairing.sh would put it -- on SOE the copy in
# use lives in the shared /vol_efthymios/NFS07/Data tree, so an auto-detect that
# only checked $TC_INPUT_DATA/ecoregions reported NOT FOUND and the maps were
# drawn without a CONUS outline twice before anyone noticed. $ECO_SHP overrides
# the search outright.
ECO_DIR="${ECO_DIR:-$TC_INPUT_DATA/ecoregions}"
ECO_CANDIDATES="
${ECO_SHP:-}
$ECO_DIR/us_eco_l3.shp
/vol_efthymios/NFS07/Data/vegetation_indices/us_eco_l3.shp
$TC_INPUT_DATA/vegetation_indices/us_eco_l3.shp
"
case " $* " in
    *" --basemap "*) ;;                     # caller chose; leave it alone
    *)  SHP=""
        for c in $ECO_CANDIDATES; do
            [ -n "$c" ] && [ -f "$c" ] && { SHP="$c"; break; }
        done
        if [ -n "$SHP" ]; then
            echo "basemap    : $SHP"
            set -- "$@" --basemap "$SHP"
        else
            echo "basemap    : NOT FOUND -- panels will have no CONUS outline." >&2
            echo "             looked in:" >&2
            for c in $ECO_CANDIDATES; do
                [ -n "$c" ] && echo "               $c" >&2
            done
            echo "             pass --basemap <path>, set ECO_SHP, or run" >&2
            echo "             slurm/submit_verify_pairing.sh to download it." >&2
        fi ;;
esac

python -u figure_metric_maps.py "$@"
rc=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc"
exit $rc
