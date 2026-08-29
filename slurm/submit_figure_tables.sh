#!/bin/bash
## Set 2: violin + error bar figures of the LMA effect, model vs model.
##
##     source slurm/config.sh
##     sbatch -p SOE_legacy -A efthymios slurm/submit_effect_figures.sh
##
## Reads station_metrics.csv and station_sensitivity.csv from $TC_RESULTS, so
## submit_station_metrics.sh must have run first. Writes 20 PNGs to
## $TC_FIGURES: 5 metrics x {era5, gcm, gcm_median, gcm_by_model}.
##
## Metrics: variability, flux, sens_Ta, sens_SPEI12, sens_LMA. Subset with
## --metrics and --sets, e.g.
##     sbatch ... slurm/submit_effect_figures.sh --metrics flux --sets era5
##
## sens_LMA is NOT a percent change. The fixed arm holds LMA constant, so its
## regression slope is NaN by construction and a percent change against it is
## 0/0; that panel shows the dynamic arm's standardised slope instead.
##
## No silent skipping: any metric whose rows are absent is named on stderr and
## the job exits non-zero, rather than quietly producing fewer figures.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -p SOE_main
#SBATCH -J figtab
#SBATCH -t 00:30:00
#SBATCH -o slurm/logs/figtab_%j.out
#SBATCH -e slurm/logs/figtab_%j.err

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

# table_scatter.csv covers the zr95 and hc axes as well as phi, and those read
# ZR95_H/hc_H back out of the MOD_PARAM files the runs actually used. Without
# MODEL_RUN set, 10 of the 15 combinations skip with a note instead of failing,
# so the table would look complete and silently be a third of itself.
# config.sh already sets MODEL_RUN; this only has to export it. Written with
# ${TC_MODEL_RUN:-} rather than a bare $TC_MODEL_RUN because that variable does
# not exist in config.sh, and under "set -u" expanding it would abort the job.
export MODEL_RUN="${MODEL_RUN:-${TC_MODEL_RUN:-}}"
[ -d "$MODEL_RUN" ] || echo "  ! MODEL_RUN='$MODEL_RUN' is not a directory;" \
    "the zr95 and hc axes of table_scatter.csv will be skipped" >&2
python -u figure_tables.py "$@"
rc=$?
echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $rc"
exit $rc
