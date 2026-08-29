#!/bin/bash
## The study's station table as a Word document, ready to paste.
##
##     source slurm/config.sh
##     sbatch -p SOE_legacy -A efthymios slurm/submit_station_table.sh
##     sbatch -p SOE_legacy -A efthymios slurm/submit_station_table.sh --layout single
##
## TWO STEPS IN ONE JOB. Step 1 (dump_station_inputs.py) extracts everything the
## tables need from sources that exist only here -- MOD_PARAM, the forcing .mat
## files, the uncapped root-depth table, the FLUXNET archives -- into three small
## CSVs. Step 2 builds the tables from those. Copy the CSVs down and step 2 runs
## anywhere, so the Word formatting can be iterated without a SLURM round trip.
##
## Three 6.5-inch portrait tables at 10 pt: Table 1 is site identity, coordinates,
## vegetation type, elevation and the AmeriFlux DOIs; Table 2 is the
## station-specific model parameters; Table 3 is the parameters prescribed per
## vegetation type, with deciduous and evergreen columns. --layout single puts
## Tables 1-2 in one landscape table. --font sets the point size (10-12).
##
## Station markers: * not carried through the GCM set; a rooting depth capped at
## the soil column depth; b tower measurements entered the flux comparison.
##
## STEP 1 MUST RUN ON THE CLUSTER. The parameter columns come from each station's
## MOD_PARAM under $MODEL_RUN -- what the model actually used, which is not the
## same as the fetched inputs: ZR95_H is capped at the soil column depth, and
## hc_H prefers BADM HEIGHTC over the GEDI/Potapov raster. Elevation is read from
## the forcing .mat rather than the AmeriFlux registry, because those two can
## disagree and only one of them was in the run.
##
## DOIs come from preprocessing/ameriflux_dois.csv, which is cached IN THE REPO
## (~90 scraped pages, and the citation list has to be reproducible later). It is
## already populated; refresh it only when sites publish new versions:
##
##     python preprocessing/fetch_ameriflux_dois.py --site-lists --refresh
##
## Writes station_table.docx and station_table.csv to $TC_RESULTS. The CSV holds
## every column including the per-value source notes and the full citations, so
## the reference list can be built from it without re-scraping.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p SOE_legacy
#SBATCH -J sttable
#SBATCH -t 00:30:00
#SBATCH -o slurm/logs/sttable_%j.out
#SBATCH -e slurm/logs/sttable_%j.err

set -uo pipefail
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

TOWER_DIR="${TOWER_DIR:-$REPO_ROOT/preprocessing/fluxnet}"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
tc_check_args "$@" || exit 1
echo "results    : $TC_RESULTS"
echo "model run  : $MODEL_RUN"
echo "towers     : $TOWER_DIR$([ -d "$TOWER_DIR" ] || echo '   <-- NOT FOUND, no b markers')"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo

[ -d "$MODEL_RUN" ] || {
    echo "ERROR: MODEL_RUN='$MODEL_RUN' is not a directory; every parameter" >&2
    echo "       column would be an em dash." >&2
    exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python -c "import docx" 2>/dev/null || {
    echo "ERROR: python-docx is not installed in $TC_VENV." >&2
    echo "       pip install python-docx    (it is in requirements.txt)" >&2
    exit 1; }

export MODEL_RUN

# STEP 1 -- the summary. Reads MOD_PARAM, the forcing .mat files, the uncapped
# root-depth table and the FLUXNET archives, all of which exist only here, and
# writes three small CSVs. Copy those down and build_station_table.py will run
# anywhere, so the Word formatting can be iterated off-cluster:
#
#   scp -P 222 dd1136@soemaster2.hpc.rutgers.edu:$TC_RESULTS/'{station_inputs,mod_param_values,tower_overlap}.csv' result_summary/
#
echo "--- step 1: extracting the inputs summary ---"
python -u dump_station_inputs.py \
    --tower-dir "$TOWER_DIR" \
    --root-depth "$TC_INPUT_DATA/root_depth/root_depth_schenk_jackson.csv"
rc=$?
[ $rc -eq 0 ] || { echo "ERROR: the summary step failed; not building tables" >&2; exit $rc; }

# STEP 2 -- the tables, from what step 1 just wrote.
echo
echo "--- step 2: building the tables ---"
python -u build_station_table.py --csv station_table.csv "$@"
rc=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc"
exit $rc
