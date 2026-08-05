#!/bin/bash
# Preflight check: can a COMPUTE node reach the Climate Data Store?
#
# Clusters commonly firewall compute nodes while leaving login nodes routable. If that
# were the case here, a submitted download job would fail every request and burn its
# walltime -- and the failure would look like a CDS problem, not a network one.
#
#     bash slurm/check_cds_access.sh
#
# Uses Python for the reachability probes rather than curl, which is not installed on the
# SOE nodes. The authoritative test is the real CDS retrieval at the end: a reachability
# probe only proves a route exists, not that credentials and the licence are in order.

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TC_VENV="${TC_VENV:-$HOME/envs/tc-preproc}"
PARTITION="${PARTITION:-SOE_main}"
PYTHON_MODULE="${PYTHON_MODULE:-Python/3.13.7}"

if [ ! -d "$TC_VENV" ]; then
    echo "ERROR: venv not found at $TC_VENV -- run 'bash slurm/setup_env.sh' first" >&2
    exit 1
fi

if ! command -v module >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source /opt/apps/lmod/lmod/init/profile
    export MODULEPATH_ROOT=/opt/apps/lmod/lmod/modulefiles
    export MODULEPATH=$MODULEPATH_ROOT/Core
    export LMOD_PACKAGE_PATH=/opt/apps/lmod/lmod/libexec
fi
module purge >/dev/null 2>&1
module load "$PYTHON_MODULE" >/dev/null 2>&1

PY="$TC_VENV/bin/python"

# An HTTP error response still proves a route exists; only a transport failure does not.
read -r -d '' PROBE <<'PY'
import sys, urllib.error, urllib.request
try:
    urllib.request.urlopen("https://cds.climate.copernicus.eu/api", timeout=20)
    print("    reachable (HTTP 200)")
except urllib.error.HTTPError as e:
    print(f"    reachable (HTTP {e.code} -- a response is what matters here)")
except Exception as e:
    print(f"    NO route to the CDS: {e}")
    sys.exit(1)
PY

echo "==> probing CDS reachability from the login node"
"$PY" -c "$PROBE"

echo
echo "==> probing from a compute node on $PARTITION"
srun --partition="$PARTITION" --nodes=1 --cpus-per-task=1 --mem=2G --time=00:05:00 \
     bash -c "echo \"    running on \$(hostname)\"; '$PY' -c '$PROBE'"

echo
echo "==> real CDS request from a compute node (one variable, two days, one point)"
srun --partition="$PARTITION" --nodes=1 --cpus-per-task=1 --mem=4G --time=00:30:00 \
     bash -c "
        source '$TC_VENV/bin/activate'
        python - <<'PY'
import cdsapi, os, tempfile
c = cdsapi.Client()
out = os.path.join(tempfile.mkdtemp(), 'probe.nc')
c.retrieve('reanalysis-era5-land-timeseries', {
    'variable': ['2m_temperature'],
    'location': {'latitude': 45.2, 'longitude': -68.7},
    'date': ['2020-01-01/2020-01-02'],
    'data_format': 'netcdf',
}, out)
print('    OK, retrieved', os.path.getsize(out), 'bytes')
PY
     "
probe_status=$?

echo
if [ $probe_status -eq 0 ]; then
    cat <<'EOF'
==> PASS. Compute nodes can reach the CDS, credentials work, and the ERA5-Land licence
    is accepted. Submit the download with:

        sbatch slurm/submit_era5_download.sh
EOF
else
    cat <<'EOF'
==> FAIL. The real retrieval did not succeed from a compute node.

If the reachability probes passed but this did not, the problem is credentials or the
licence, not the network: check ~/.cdsapirc and accept the ERA5-Land licence at
https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-timeseries

If the compute-node probe also failed, compute nodes have no outbound internet. Options:
  1. Ask SOE support (help@soe.rutgers.edu) about an HTTPS proxy, then export
     https_proxy in the submit script.
  2. Run the download on the login node inside tmux. It is network-bound, not CPU-bound,
     so it is light on the node -- keep --jobs low (2) to stay polite.
  3. Download locally and rsync the ~1 GB of netCDF up to the cluster.
EOF
fi
exit $probe_status
