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
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

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
# The import must succeed (set -e makes this fatal). cdsapi does not expose
# __version__, so read the version from package metadata, and treat that as
# best-effort -- a missing version string is cosmetic, a failed import is not.
python -c "import cdsapi" && echo "    cdsapi imports OK"
python -c "import importlib.metadata as md; print('    cdsapi', md.version('cdsapi'))" \
    2>/dev/null || echo "    (version lookup unavailable)"

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
