#!/bin/bash
## Download the AmeriFlux FLUXNET (ONEFlux) product for the study's stations.
##
## WHY A SECOND DOWNLOAD. BASE carries FC, the NET CO2 flux, and no GPP column at
## all -- verified against AMF_US-Ha2_BASE_HH_16-5.csv, which has zero GPP or
## RECO fields. GPP is a partitioned quantity and exists only in this product. A
## tower comparison of modelled GPP therefore cannot be done from BASE; against
## NEE it would be a test against a net flux that includes soil respiration,
## which is a much weaker check on a mechanism acting through leaf area.
##
## CCBY4.0 ONLY. Requesting FLUXNET for a LEGACY site returns nothing, so the
## downloader drops those sites and names them. Checked 2026-08-26 against the
## public site_availability endpoint: 109 of our 118 stations qualify. The nine
## that do not -- US-Blk, US-CZ2, US-CZ3, US-CZ4, US-LPH, US-MRf, US-NR2,
## US-SB3, US-WBW -- keep their BASE record; only partitioned GPP is missing.
## Policy eligibility is necessary but not sufficient: a CCBY4.0 site still only
## has a FLUXNET product if ONEFlux processed it, and there is no working
## availability endpoint to pre-check that, so read the per-site report below.
##
## Credentials must be exported BEFORE sbatch -- SLURM copies the submitting
## environment into the job, which is how they reach the compute node without
## ever being written into the repo:
##
##     export AMF_USER_ID=<your ameriflux username>
##     export AMF_USER_EMAIL=<your email>
##
## Test two sites first, marked as a test so site teams are not emailed:
##     sbatch slurm/submit_fluxnet_download.sh --agree-policy --is-test \
##            --stations US-HBK,US-Ha2
##
## Then the full set:
##     sbatch slurm/submit_fluxnet_download.sh --agree-policy
##
## OUT_DIR defaults to preprocessing/fluxnet inside the repo, which is what was
## asked for. It is gitignored: these archives are tens of MB per site and a
## 108 MB file has already had one push declined at GitHub's 100 MB hard limit.
## Override with OUT_DIR=... to put them elsewhere.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p SOE_main
#SBATCH -J flx_dl
#SBATCH -t 12:00:00
#SBATCH -o slurm/logs/flx_dl_%j.out
#SBATCH -e slurm/logs/flx_dl_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

OUT_DIR="${OUT_DIR:-$REPO_ROOT/preprocessing/fluxnet}"

mkdir -p "$REPO_ROOT/slurm/logs"
echo "job        : ${SLURM_JOB_ID:-interactive}"
echo "node       : $(hostname)"
tc_check_partition
tc_check_args "$@" || exit 1
echo "out_dir    : $OUT_DIR"
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo

if [ -z "${AMF_USER_ID:-}" ] || [ -z "${AMF_USER_EMAIL:-}" ]; then
    echo "ERROR: AMF_USER_ID and AMF_USER_EMAIL must be exported before sbatch." >&2
    echo "       SLURM copies the submitting environment; it cannot prompt here." >&2
    exit 1
fi

mkdir -p "$OUT_DIR" || { echo "ERROR: cannot create $OUT_DIR" >&2; exit 1; }

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

python -u download_ameriflux.py --product FLUXNET --out "$OUT_DIR" "$@"
rc=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc"
exit $rc
