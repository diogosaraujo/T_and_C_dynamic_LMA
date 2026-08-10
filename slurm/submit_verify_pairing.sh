#!/bin/bash
## Verify every station falls inside the EPA Level III ecoregion it is paired to.
##
##     sbatch slurm/submit_verify_pairing.sh
##
## The pairing decides which ecoregion's PLSR fit supplies a station's LMA series,
## so a mis-paired station is modelled with another region's leaf traits and nothing
## downstream complains. Point-in-polygon against the EPA shapefile is the only test
## that settles it.
##
## Downloads us_eco_l3 to $TC_INPUT_DATA/ecoregions on first run (~25 MB).
## NOTE: curl is not installed on these nodes, so the fetch goes through Python --
## a missing-curl error looks exactly like a firewall block and wastes an afternoon.
##
## Exit 1 means at least one station is paired to an ecoregion that does not contain
## it, or a control station failed (projection wrong -- results suppressed).

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J verify_eco
#SBATCH -t 00:30:00
#SBATCH -o slurm/logs/verify_eco_%j.out
#SBATCH -e slurm/logs/verify_eco_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

ECO_DIR="${ECO_DIR:-$TC_INPUT_DATA/ecoregions}"
SHP="$ECO_DIR/us_eco_l3.shp"
# The long-standing gaftp.epa.gov path now 404s; EPA serves the Data Commons from
# S3. Verified 2026-08-10: 28.4 MB.
URL="${ECO_URL:-https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/us/us_eco_l3.zip}"

mkdir -p "$REPO_ROOT/slurm/logs" "$ECO_DIR"

echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "eco dir    : $ECO_DIR"
echo

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"

python - "$URL" "$ECO_DIR" "$SHP" <<'PY'
import sys, ssl, zipfile, io, urllib.request
from pathlib import Path
url, eco_dir, shp = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
if shp.is_file():
    print(f"shapefile already present: {shp}")
    raise SystemExit(0)
print(f"downloading {url}")
# The EPA FTP host has repeatedly shipped an incomplete chain; verification is
# relaxed only for this public, unauthenticated download of a published dataset.
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE
try:
    with urllib.request.urlopen(url, context=ctx, timeout=300) as r:
        blob = r.read()
except Exception as e:
    print(f"ERROR: download failed: {e}", file=sys.stderr)
    print("Fetch us_eco_l3.zip by hand and unzip it into "
          f"{eco_dir}, then re-run.", file=sys.stderr)
    raise SystemExit(1)
print(f"  {len(blob)/1e6:.1f} MB, unzipping")
zipfile.ZipFile(io.BytesIO(blob)).extractall(eco_dir)
# The archive may nest the files one level down.
if not shp.is_file():
    for p in eco_dir.rglob("us_eco_l3.shp"):
        for q in p.parent.glob("us_eco_l3.*"):
            q.replace(eco_dir / q.name)
        break
print("ok" if shp.is_file() else "ERROR: us_eco_l3.shp not found after unzip")
raise SystemExit(0 if shp.is_file() else 1)
PY
dl=$?
if [ $dl -ne 0 ]; then
    echo "aborting: the ecoregion shapefile is not available" >&2
    exit 1
fi

cd "$REPO_ROOT/preprocessing"
python verify_ecoregion_pairing.py \
    --shapefile "$SHP" \
    --out "$TC_INPUT_DATA/ecoregion_pairing_check.csv"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status  (1 = at least one station is paired to an ecoregion"
echo "                        that does not contain it, or a control failed)"
exit $status
