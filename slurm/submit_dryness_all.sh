#!/bin/bash
## Fan build_dryness.py out over 16 jobs, then merge.
##
##     bash slurm/submit_dryness_all.sh          # submits 16 + a merge job
##
## One job for ERA5 and one per (GCM x scenario), each writing its OWN csv so
## nothing races on a shared file -- the same reason submit_station_metrics is
## run with --out-prefix rather than --update when parallel.
##
## The merge waits on all 16 via --dependency=afterok and concatenates into
## station_dryness.csv. If any job fails the merge does not run, so a partial
## table is never written.
set -uo pipefail
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

GCMS="GFDL-ESM4 IPSL-CM6A-LR MPI-ESM1-2-HR MRI-ESM2-0 UKESM1-0-LL"
SUB="sbatch -p ${TC_PART:-SOE_legacy} -A ${TC_ACCT:-efthymios} --parsable"
ids=()

id=$($SUB "$REPO_ROOT/slurm/submit_dryness.sh" \
      --datasets era5 --out dryness_part_era5.csv)
echo "  era5                        -> job $id"; ids+=("$id")

for scen in historical ssp126 ssp585; do
  for g in $GCMS; do
    id=$($SUB "$REPO_ROOT/slurm/submit_dryness.sh" \
          --datasets "$scen" --gcms "$g" \
          --out "dryness_part_${scen}_${g}.csv")
    echo "  $scen $g -> job $id"; ids+=("$id")
  done
done

dep=$(IFS=:; echo "${ids[*]}")
mid=$(sbatch -p "${TC_PART:-SOE_legacy}" -A "${TC_ACCT:-efthymios}" --parsable \
      --dependency=afterok:"$dep" --job-name=drymerge \
      --output="$REPO_ROOT/slurm/logs/drymerge_%j.out" \
      --error="$REPO_ROOT/slurm/logs/drymerge_%j.err" \
      --wrap "source $REPO_ROOT/slurm/config.sh; \
              source \$TC_VENV/bin/activate; \
              python -u $REPO_ROOT/preprocessing/merge_dryness.py")
echo "  merge (after all 16)        -> job $mid"
