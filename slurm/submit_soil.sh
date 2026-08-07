#!/bin/bash
## Build depth-resolved soil profiles and column depths for the AmeriFlux stations.
##
##     sbatch slurm/submit_soil.sh --probe        # coverage check first, writes nothing
##     sbatch slurm/submit_soil.sh                # build the profiles, then report
##     sbatch slurm/submit_soil.sh --report       # re-summarise an existing run
##     sbatch slurm/submit_soil.sh --stations US-Ho1,US-xRM
##
## Source order, first hit wins, recorded per station:
##     AmeriFlux BADM  ->  SSURGO (USDA Soil Data Access)  ->  POLARIS  ->  SoilGrids
##
## Only SSURGO reports a named restriction (lithic/paralithic bedrock, densic,
## fragipan) and its horizons terminate there, so texture and column depth come from
## one internally consistent profile. POLARIS is itself disaggregated from SSURGO and
## serves as the CONUS gap filler; SoilGrids is global and properties-only.
##
## Reads (network):
##     https://sdmdataaccess.sc.egov.usda.gov/Tabular/post.rest   (no auth)
##     http://hydrology.cee.duke.edu/POLARIS/...                  (/vsicurl/)
##     https://rest.isric.org/soilgrids/v2.0/...
## Reads (local):
##     $TC_INPUT_DATA/ameriflux/<SITE>/AMF_*_BIF_*.xlsx           (in-situ texture)
##     $TC_INPUT_DATA/root_depth/root_depth.csv                   (ZR95, so the column
##                                                                 is never shallower)
## Writes:
##     $TC_INPUT_DATA/soil/soil_sites.csv          one row per station + provenance
##     $TC_INPUT_DATA/soil/soil_profiles.csv       one row per Zs layer
##     $TC_INPUT_DATA/soil/soil_horizons_raw.csv   the horizons as fetched
##     $TC_INPUT_DATA/soil/soil_provenance.json
##
## Porg is written as a FRACTION, which is what Soil_parameters expects. Each source
## states organic matter differently (SSURGO percent, POLARIS log10-percent,
## SoilGrids SOC in dg/kg needing x1.72) and all are normalised on the way in.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J soil
#SBATCH -t 04:00:00
#SBATCH -o slurm/logs/soil_%j.out
#SBATCH -e slurm/logs/soil_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

OUT_DIR="${OUT_DIR:-$TC_INPUT_DATA/soil}"
BADM_DIR="${BADM_DIR:-$TC_INPUT_DATA/ameriflux}"
ROOT_DEPTH="${ROOT_DEPTH:-$TC_INPUT_DATA/root_depth/root_depth.csv}"
SOURCES="${SOURCES:-badm,ssurgo,polaris,soilgrids}"

mkdir -p "$REPO_ROOT/slurm/logs" "$OUT_DIR"

echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "sources    : $SOURCES"
echo "badm dir   : $BADM_DIR"
echo "root depth : $ROOT_DEPTH"
echo "output     : $OUT_DIR"
echo

if [ ! -d "$TC_VENV" ]; then
    echo "ERROR: venv not found at $TC_VENV -- run 'sbatch slurm/submit_setup_env.sh' first" >&2
    exit 1
fi
# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"

# Reachability, checked with Python because curl is not installed on these nodes --
# its absence looks exactly like a firewall block and sends you after the wrong bug.
python - <<'PY'
import json, sys, urllib.request
url = "https://sdmdataaccess.sc.egov.usda.gov/Tabular/post.rest"
body = json.dumps({"query": "SELECT TOP 1 mukey FROM mapunit",
                   "format": "JSON+COLUMNNAME"}).encode()
try:
    req = urllib.request.Request(url, data=body,
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        json.loads(r.read().decode())
    print("SDA reachable")
except Exception as exc:
    print(f"SDA NOT reachable: {type(exc).__name__}: {exc}", file=sys.stderr)
    sys.exit(1)
PY
if [ $? -ne 0 ]; then
    echo "ERROR: Soil Data Access is unreachable from this node" >&2
    exit 1
fi
echo

cd "$REPO_ROOT/preprocessing"
python fetch_soil.py \
    --out "$OUT_DIR" \
    --badm-dir "$BADM_DIR" \
    --root-depth "$ROOT_DEPTH" \
    --sources "$SOURCES" "$@"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
du -sh "$OUT_DIR" 2>/dev/null || true
exit $status
