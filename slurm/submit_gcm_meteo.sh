#!/bin/bash
## Stages 1 and 2: daily GCM station series -> hourly T&C forcing.
##
##     sbatch --array=1-15 slurm/submit_gcm_meteo.sh        # 5 GCM x 3 scenarios
##     sbatch slurm/submit_gcm_meteo.sh --dry-run           # diagnostics, writes nothing
##     sbatch slurm/submit_gcm_meteo.sh --gcm GFDL-ESM4 --scenario historical
##
## One array task = one (GCM, scenario), all stations. Two stages in one job,
## mirroring submit_meteo.sh so both forcings are built the same way:
##
##   1. build_gcm_meteo.py  daily npz -> Meteo_<ST>_<GCM>_<scen>_raw.mat
##                          disaggregation, humidity, barometric pressure, CO2
##   2. finish_meteo.m      adds SAB1/SAB2/SAD1/SAD2, PARB/PARD and N via
##                          C_Automatic_Radiation_Partition
##                          -> Meteo_<ST>_<GCM>_<scen>_<years>.mat
##
## Stage 2 reuses the SAME tested MATLAB partition as the ERA5 path, called with
## t_bef/t_aft = 0/1 instead of 1/0. That is not a tuning knob: build_gcm_meteo.py
## constructs the hourly series with hour H covering (H, H+1], the opposite of
## de-accumulated ERA5-Land, and the two ends must agree or the solar geometry is
## displaced by an hour.
##
## Requires: gcm_stations/ (submit_gcm_extract.sh) and co2/ (fetch_ssp_co2.py).
## Without the CO2 series the build still runs but Ca comes out NaN, which T&C
## will not tolerate -- the log says so explicitly.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH -p SOE_main
#SBATCH -J gcm_meteo
#SBATCH -t 12:00:00
#SBATCH -o slurm/logs/gcm_meteo_%A_%a.out
#SBATCH -e slurm/logs/gcm_meteo_%A_%a.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

STATION_DIR="${STATION_DIR:-$TC_INPUT_DATA/gcm_stations}"
METEO_DIR="${METEO_DIR:-$TC_INPUT_DATA/gcm_meteo}"
CO2_DIR="${CO2_DIR:-$TC_INPUT_DATA/co2}"
PARTITION_DIR="${PARTITION_DIR:-$REPO_ROOT/T&C/Diogo}"
PRECIP_SCHEME="${PRECIP_SCHEME:-block}"

mkdir -p "$REPO_ROOT/slurm/logs" "$METEO_DIR"
echo "job        : ${SLURM_JOB_ID:-interactive}${SLURM_ARRAY_TASK_ID:+  task $SLURM_ARRAY_TASK_ID}"
echo "node       : $(hostname)"
tc_check_partition
echo "started    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "stations   : $STATION_DIR$([ -d "$STATION_DIR" ] || echo '   <-- NOT FOUND')"
echo "co2        : $CO2_DIR$([ -d "$CO2_DIR" ] || echo '   <-- NOT FOUND, Ca will be NaN')"
echo "output     : $METEO_DIR"
echo "precip     : $PRECIP_SCHEME"
echo

[ -d "$STATION_DIR" ] || { echo "ERROR: run submit_gcm_extract.sh first" >&2; exit 1; }
[ -f "$PARTITION_DIR/C_Automatic_Radiation_Partition.m" ] || {
    echo "ERROR: C_Automatic_Radiation_Partition.m not in $PARTITION_DIR" >&2; exit 1; }

# ---------------------------------------------------------------- stage 1
# shellcheck disable=SC1091
source "$TC_VENV/bin/activate" || { echo "ERROR: venv $TC_VENV missing" >&2; exit 1; }
cd "$REPO_ROOT/preprocessing" || exit 1

# Add the array index unless the caller already chose what to work on. Passing a
# NON-selector flag (--force, --dry-run) must not suppress it: job 36991 ran
# "sbatch --array=1-105 ... --force", $# was 1 rather than 0, --index was never
# added, and all 105 tasks died in argparse one second in.
have_selector=0
for _a in "$@"; do
    case "$_a" in --index|--all|--gcm|--scenario|--station) have_selector=1 ;; esac
done
ARGS=("$@")
if [ "$have_selector" -eq 0 ] && [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    ARGS+=(--index "$SLURM_ARRAY_TASK_ID")
elif [ "$have_selector" -eq 0 ]; then
    ARGS+=(--all)
fi
echo "args       : ${ARGS[*]}"

python build_gcm_meteo.py \
    --stations-root "$STATION_DIR" --out "$METEO_DIR" --co2-dir "$CO2_DIR" \
    --precip-scheme "$PRECIP_SCHEME" "${ARGS[@]}"
s1=$?

# A dry run writes no _raw.mat, so stage 2 would find nothing and error.
case " ${ARGS[*]} " in *" --dry-run "*)
    echo; echo "dry run -- stage 2 skipped"; exit $s1 ;;
esac
[ $s1 -eq 0 ] || { echo "stage 1 failed -- not running the partition" >&2; exit $s1; }

# ---------------------------------------------------------------- stage 2
tc_load_matlab || exit 1

# The year tag is taken from the raw files themselves rather than assumed, because
# historical (1980-2014) and the SSPs (2015-2100) differ.
# One directory per (scenario, GCM), and this task only touches the ones stage 1
# just wrote. finish_meteo.m globs a whole directory, so a shared one means every
# concurrent task re-partitions every file and several write the same output at
# once -- job 37232 corrupted files that way.
s2=0
for SUB in "$METEO_DIR"/*/*/; do
    tag=$(basename "$(dirname "$SUB")")
    ls "$SUB"/Meteo_*_raw.mat >/dev/null 2>&1 || continue
    case "$tag" in
        historical) YEARS=1980_2014 ;;
        ssp*)       YEARS=2015_2100 ;;
        *)          continue ;;
    esac
    echo "  partition: $tag/$(basename "$SUB") -> $YEARS  ($(ls "$SUB"/Meteo_*_raw.mat | wc -l) file(s))"
    matlab -nodisplay -nosplash -batch         "finish_meteo('$SUB','$SUB','$PARTITION_DIR','$YEARS',0,1)"
    s2=$(( s2 || $? ))
done

echo
echo "raw  files : $(find "$METEO_DIR" -name 'Meteo_*_raw.mat' 2>/dev/null | wc -l)"
echo "final files: $(find "$METEO_DIR" -name 'Meteo_*_raw.mat' -prune -o -name 'Meteo_*.mat' -print 2>/dev/null | wc -l)"
du -sh "$METEO_DIR" 2>/dev/null || true
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: stage1=$s1 stage2=$s2"
exit $(( s1 || s2 ))
