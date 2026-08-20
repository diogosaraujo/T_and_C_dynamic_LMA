#!/bin/bash
# Sets up the Python environment for the preprocessing steps on the SOE cluster.
#
# Do not run this directly in the login shell -- submit it:
#     sbatch slurm/submit_setup_env.sh
#
# Creates a venv at $TC_VENV (default ~/envs/tc-preproc) from requirements.txt, and
# checks that the CDS credentials file exists. Re-run it whenever requirements.txt
# changes; the venv is reused, so it is quick.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO_ROOT/slurm/config.sh"

# LMOD is not initialised in non-interactive shells unless .bashrc ran, so do it
# here -- via config.sh's tc_init_lmod, which tries several init paths. Sourcing
# /opt/apps/lmod/... directly is the SOE location and does not exist on Amarel,
# so hardcoding it made this script SOE-only.
if ! tc_init_lmod; then
    echo "ERROR: could not initialise LMOD on $(hostname)" >&2
    echo "       'module avail' by hand will show whether modules work at all." >&2
    exit 1
fi

# An empty PYTHON_MODULE means "use the interpreter already on PATH" -- see the
# note in config.sh about Amarel's broken python/3.8.2.
if [ -n "${PYTHON_MODULE:-}" ]; then
    echo "==> loading $PYTHON_MODULE"
    module purge
    module load "$PYTHON_MODULE"
else
    echo "==> no python module configured; using $(command -v python3)"
fi
python3 --version

# Fail here, clearly, rather than 30 lines into pip's retry storm. A Python
# whose ssl module will not import cannot talk to PyPI at all, and the error it
# produces down in pip is indistinguishable from a network block.
if ! python3 -c "import ssl" 2>/dev/null; then
    echo "ERROR: this Python cannot import ssl, so pip cannot reach PyPI." >&2
    python3 -c "import ssl" || true
    echo "       Try a different interpreter, or install the wheels offline:" >&2
    echo "         pip install --no-index --find-links=<dir> -r <requirements>" >&2
    exit 1
fi

# A venv is BOUND to the interpreter that created it: bin/python is a link back
# to it and the whole tree hardcodes its version. So "the directory exists" is
# not the same as "the venv works". Job 60681690 reused a venv built from the
# broken python/3.8.2 module and died at exit 127 with
#     libpython3.8.so.1.0: cannot open shared object file
# because the module was no longer loaded. Test it and rebuild if it is dead.
if [ -z "${TC_VENV:-}" ]; then
    echo "ERROR: TC_VENV is empty -- refusing to touch it" >&2
    exit 1
fi
if [ -d "$TC_VENV" ] && ! "$TC_VENV/bin/python" -c 'import sys' >/dev/null 2>&1; then
    echo "==> existing venv at $TC_VENV is unusable (its interpreter will not"
    echo "    start -- most often built against a module that is gone). Recreating."
    "$TC_VENV/bin/python" -c 'import sys' || true
    rm -rf "${TC_VENV:?}"
fi
if [ ! -d "$TC_VENV" ]; then
    echo "==> creating venv at $TC_VENV with $(command -v python3)"
    python3 -m venv "$TC_VENV"
else
    echo "==> reusing existing venv at $TC_VENV"
fi

# shellcheck disable=SC1091
source "$TC_VENV/bin/activate"
python -m pip install --upgrade pip
# Which dependency set. The full preprocessing stack is the default; the run
# side (Amarel) needs only numpy/h5py/scipy and must not be made to build GDAL.
#     sbatch slurm/submit_setup_env.sh --run
REQ="$REPO_ROOT/preprocessing/requirements.txt"
for arg in "$@"; do
    case "$arg" in
        --run)  REQ="$REPO_ROOT/preprocessing/requirements-run.txt" ;;
        --requirements=*) REQ="${arg#--requirements=}" ;;
    esac
done
echo "==> installing from $REQ"
python -m pip install -r "$REQ"

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
