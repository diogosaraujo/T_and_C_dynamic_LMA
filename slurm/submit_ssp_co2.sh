#!/bin/bash
## Build the annual CO2 series T&C needs for Ca, one per scenario.
##
##     sbatch slurm/submit_ssp_co2.sh                 # convert whatever is in co2/raw/
##     sbatch slurm/submit_ssp_co2.sh --report        # what is already built
##     sbatch slurm/submit_ssp_co2.sh --from-netcdf /some/other/dir/*.nc
##
## NEX-GDDP-CMIP6 carries no CO2 -- it is nine surface weather fields and nothing
## else -- so Ca has to come from outside. CMIP6 ScenarioMIP runs are
## CONCENTRATION-driven: every model is prescribed the same harmonised
## greenhouse-gas concentrations (Meinshausen et al. 2020, GMD 13, 3571,
## distributed through input4MIPs), so ONE series per scenario serves all five
## GCMs and there is nothing model-specific to fetch.
##
## With no arguments it DOWNLOADS the three input4MIPs files from ESGF and then
## converts them; already-downloaded files are reused. If the download fails it
## falls back to converting whatever is already in $TC_INPUT_DATA/co2/raw/.
##
## The datasets, should you ever need to fetch them by hand:
##
##   historical 1980-2014  input4MIPs.CMIP6.CMIP.UoM.UoM-CMIP-1-2-0
##   ssp126     2015-2100  input4MIPs.CMIP6.ScenarioMIP.UoM.UoM-IMAGE-ssp126-1-2-1
##   ssp585     2015-2100  input4MIPs.CMIP6.ScenarioMIP.UoM.UoM-REMIND-MAGPIE-ssp585-1-2-1
##
## from https://esgf-node.llnl.gov/search/input4mips/ (ESGF search, needs a free
## account for some nodes) or http://greenhousegases.science.unimelb.edu.au
## (the same data, no account). The scenario is read from each filename.
##
## The ESGF index is queried rather than any URL being hard-coded, because ESGF is
## mid-migration to ESGF-NG and the endpoints move: esgf-node.llnl.gov now 302s to
## an ORNL bridge, and three more index nodes are tried after it. Two quirks are
## handled and both return an EMPTY result rather than an error, so neither is
## self-announcing: the CMIP historical collection indexes the variable with
## HYPHENS while ScenarioMIP uses UNDERSCORES, and the year-0 start of the
## historical file needs has_year_zero on the calendar.
##
## Writes $TC_INPUT_DATA/co2/co2_<scenario>.csv, which build_gcm_meteo.py
## interpolates onto the hourly stamp.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -p SOE_main
#SBATCH -J ssp_co2
#SBATCH -t 00:20:00
#SBATCH -o slurm/logs/ssp_co2_%j.out
#SBATCH -e slurm/logs/ssp_co2_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

CO2_DIR="${CO2_DIR:-$TC_INPUT_DATA/co2}"
RAW_DIR="${RAW_DIR:-$CO2_DIR/raw}"

mkdir -p "$REPO_ROOT/slurm/logs" "$CO2_DIR" "$RAW_DIR"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "raw inputs : $RAW_DIR"
echo "output     : $CO2_DIR"
echo

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

# --report and an explicit --from-netcdf/--from-csv pass straight through.
if [ "$#" -gt 0 ]; then
    python fetch_ssp_co2.py --out "$CO2_DIR" --raw-dir "$RAW_DIR" "$@"
    status=$?
else
    python fetch_ssp_co2.py --out "$CO2_DIR" --raw-dir "$RAW_DIR" --download
    status=$?
    if [ $status -ne 0 ]; then
        shopt -s nullglob
        RAW=( "$RAW_DIR"/*.nc "$RAW_DIR"/*.csv )
        shopt -u nullglob
        if [ ${#RAW[@]} -gt 0 ]; then
            echo
            echo "download failed; converting the ${#RAW[@]} file(s) already in $RAW_DIR"
            python fetch_ssp_co2.py --out "$CO2_DIR" --from-netcdf "${RAW[@]}"
            status=$?
        fi
    fi
fi

echo
echo "series built:"
ls -l "$CO2_DIR"/co2_*.csv 2>/dev/null | awk '{print "  ", $5, $9}' || echo "  none"
echo
echo "sanity check -- historical should end near 397 ppm, ssp126 near 446 and"
echo "ssp585 near 1135 by 2100. A factor of 1e6 out means the mole-fraction"
echo "to ppm conversion did not fire."
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
exit $status
