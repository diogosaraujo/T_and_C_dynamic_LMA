#!/bin/bash
# One-time setup of the Python environment for the preprocessing step on the SOE cluster.
#
# Run this ONCE on the login node (it only pip-installs, no heavy compute):
#     bash slurm/setup_env.sh
#
# Creates a venv at $TC_VENV (default ~/envs/tc-preproc) with cdsapi installed, and
# checks that the CDS credentials file exists.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TC_VENV="${TC_VENV:-$HOME/envs/tc-preproc}"
PYTHON_MODULE="${PYTHON_MODULE:-Python/3.13.7}"

# LMOD is not initialised in non-interactive shells unless .bashrc ran, so do it here.
if ! command -v module >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source /opt/apps/lmod/lmod/init/profile
    export MODULEPATH_ROOT=/opt/apps/lmod/lmod/modulefiles
    export MODULEPATH=$MODULEPATH_ROOT/Core
    export LMOD_PACKAGE_PATH=/opt/apps/lmod/lmod/libexec
fi

echo "==> loading $PYTHON_MODULE"
module purge
module load "$PYTHON_MODULE"
python3 --version

if [ ! -d "$TC_VENV" ]; then
    echo "==> creating venv at $TC_VENV"
    python3 -m venv "$TC_VENV"
else
    echo "==> reusing existing venv at $TC_VENV"
fi

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"
python -m pip install --upgrade pip
python -m pip install -r "$REPO_ROOT/preprocessing/requirements.txt"

echo
echo "==> installed:"
python -c "import cdsapi; print('cdsapi', cdsapi.__version__)"

echo
if [ -f "$HOME/.cdsapirc" ]; then
    echo "==> found $HOME/.cdsapirc"
    chmod 600 "$HOME/.cdsapirc"
else
    cat <<'EOF'
==> MISSING ~/.cdsapirc

Create it with your key from https://cds.climate.copernicus.eu/profile :

    cat > ~/.cdsapirc <<'KEY'
    url: https://cds.climate.copernicus.eu/api
    key: <your-api-key>
    KEY
    chmod 600 ~/.cdsapirc

Also accept the ERA5-Land licence on the dataset page, or every request fails.
EOF
    exit 1
fi

echo
echo "Setup complete. Next: bash slurm/check_cds_access.sh"
