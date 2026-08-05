#!/bin/bash
## Download AmeriFlux BASE measurements + BADM metadata, then report parameter coverage.
##
## Credentials must be exported BEFORE sbatch -- SLURM copies the submitting environment
## into the job, so this is how they reach the compute node without ever being written
## into the repo:
##
##     export AMF_USER_ID=<your ameriflux username>
##     export AMF_USER_EMAIL=<your email>
##     sbatch slurm/submit_ameriflux_download.sh --agree-policy
##
## Test on two sites first, and mark it as a test so site teams are not emailed:
##     sbatch slurm/submit_ameriflux_download.sh --agree-policy --is-test \
##            --stations US-HBK,US-Ha2
##
## Public site metadata needs no account at all:
##     python preprocessing/download_ameriflux.py --metadata-only

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J amf_dl
#SBATCH -t 12:00:00
#SBATCH -o slurm/logs/amf_dl_%j.out
#SBATCH -e slurm/logs/amf_dl_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"
OUT_DIR="${OUT_DIR:-$TC_INPUT_DATA/ameriflux}"

mkdir -p "$REPO_ROOT/slurm/logs" "$OUT_DIR"

echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "output     : $OUT_DIR"
echo "amf user   : ${AMF_USER_ID:-<unset>}"
echo

if [ ! -d "$TC_VENV" ]; then
    echo "ERROR: venv not found at $TC_VENV -- run 'bash slurm/setup_env.sh' first" >&2
    exit 1
fi
# --metadata-only uses the public site registry, which needs no account.
metadata_only=0
for arg in "$@"; do
    [ "$arg" = "--metadata-only" ] && metadata_only=1
done
if [ "$metadata_only" -eq 0 ] && { [ -z "${AMF_USER_ID:-}" ] || [ -z "${AMF_USER_EMAIL:-}" ]; }; then
    echo "ERROR: AMF_USER_ID / AMF_USER_EMAIL are not set in the job environment." >&2
    echo "       export them in your shell BEFORE running sbatch, or pass" >&2
    echo "       --metadata-only to fetch just the public site registry." >&2
    exit 1
fi

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"
cd "$REPO_ROOT/preprocessing"

python download_ameriflux.py --out "$OUT_DIR" "$@"
status=$?

# Nothing to inspect after a metadata-only run -- no BADM has been downloaded.
if [ $status -eq 0 ] && [ "$metadata_only" -eq 0 ]; then
    echo
    echo "=== BADM parameter coverage ==="
    python inspect_ameriflux_badm.py --dir "$OUT_DIR"
    status=$?
fi

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
du -sh "$OUT_DIR" 2>/dev/null || true
exit $status
