#!/bin/bash
# Preflight check: can a COMPUTE node reach the Climate Data Store?
#
# This is the single most likely thing to break the download on a cluster. Login nodes
# often have outbound internet while compute nodes sit behind a firewall with no route
# out. If that is the case here, a submitted download job will burn its walltime failing
# every request.
#
# Run this before submitting anything:
#     bash slurm/check_cds_access.sh
#
# It grabs a short interactive allocation on SOE_main and tests connectivity from there.

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TC_VENV="${TC_VENV:-$HOME/envs/tc-preproc}"
PARTITION="${PARTITION:-SOE_main}"

echo "==> testing CDS reachability from the login node"
if curl -sS -m 20 -o /dev/null -w 'login node HTTP %{http_code}\n' \
        https://cds.climate.copernicus.eu/api; then
    echo "    login node: OK"
else
    echo "    login node: NO route to the CDS -- check the RU VPN / campus network"
fi

echo
echo "==> requesting a short allocation on $PARTITION to test from a compute node"
srun --partition="$PARTITION" --nodes=1 --cpus-per-task=1 --mem=2G --time=00:05:00 \
     bash -c "
        echo \"    running on \$(hostname)\"
        curl -sS -m 20 -o /dev/null -w '    compute node HTTP %{http_code}\n' \
            https://cds.climate.copernicus.eu/api \
            || echo '    compute node: NO outbound route to the CDS'
     "

echo
echo "==> testing an actual CDS request (tiny: one variable, one day, one point)"
srun --partition="$PARTITION" --nodes=1 --cpus-per-task=1 --mem=4G --time=00:30:00 \
     bash -c "
        source '$TC_VENV/bin/activate'
        cd '$REPO_ROOT/preprocessing'
        python - <<'PY'
import cdsapi, tempfile, os
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

cat <<'EOF'

If the compute-node checks failed but the login node worked, the compute nodes have no
outbound internet. Options, in order of preference:
  1. Ask SOE support (help@soe.rutgers.edu) whether an HTTPS proxy is available, then set
     https_proxy in the submit script.
  2. Run the download on the login node inside tmux. It is network-bound, not CPU-bound,
     so it is light on the node -- keep --jobs low (2) to stay polite.
  3. Download to your local machine and rsync the ~1 GB of netCDF up to the cluster.
EOF
