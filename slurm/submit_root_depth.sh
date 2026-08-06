#!/bin/bash
## Fetch rooting depth (Schenk & Jackson, ISLSCP II, 1 degree) for every station.
##
## Needs a free NASA Earthdata Login (https://urs.earthdata.nasa.gov). Provide it by
## either of these BEFORE submitting -- SLURM copies the submitting environment into
## the job:
##
##     # option A: ~/.netrc  (preferred; nothing to re-export each session)
##     printf 'machine urs.earthdata.nasa.gov login <user> password <pass>\n' >> ~/.netrc
##     chmod 600 ~/.netrc
##
##     # option B: environment variables
##     export EARTHDATA_USER=<user> EARTHDATA_PASS=<pass>
##
##     sbatch slurm/submit_root_depth.sh
##
## Targeted runs -- flags pass straight through:
##     sbatch slurm/submit_root_depth.sh --stations US-HBK,US-Ha2
##     sbatch slurm/submit_root_depth.sh --dry-run
##
## Output goes outside the repo, in its own folder:
##     $TC_INPUT_DATA/root_depth/
##
## Downloads only. The selection rule (AmeriFlux first, gridded fallback) is applied
## later, when the .mat files are built. Note BADM contains NO rooting-depth variable at
## all, so unlike canopy height this product is the only source for every station.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -p SOE_main
#SBATCH -J root_depth
#SBATCH -t 02:00:00
#SBATCH -o slurm/logs/root_depth_%j.out
#SBATCH -e slurm/logs/root_depth_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
OUT_DIR="${OUT_DIR:-$TC_INPUT_DATA/root_depth}"

mkdir -p "$REPO_ROOT/slurm/logs" "$OUT_DIR"

echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "output     : $OUT_DIR"
echo

if [ ! -d "$TC_VENV" ]; then
    echo "ERROR: venv not found at $TC_VENV -- run 'sbatch slurm/submit_setup_env.sh' first" >&2
    exit 1
fi

# Fail here with a clear message rather than after the job has been queued and started.
if [ -z "${EARTHDATA_TOKEN:-}" ] && [ -z "${EARTHDATA_USER:-}" ] \
   && ! grep -qs "urs.earthdata.nasa.gov" "$HOME/.netrc"; then
    echo "ERROR: no NASA Earthdata credentials visible to this job." >&2
    echo "       Register free at https://urs.earthdata.nasa.gov, then either" >&2
    echo "         add to ~/.netrc:  machine urs.earthdata.nasa.gov login <u> password <p>" >&2
    echo "         or export EARTHDATA_USER / EARTHDATA_PASS before sbatch" >&2
    exit 1
fi

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"
cd "$REPO_ROOT/preprocessing"

python fetch_root_depth.py --out "$OUT_DIR" "$@"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
du -sh "$OUT_DIR" 2>/dev/null || true
exit $status
