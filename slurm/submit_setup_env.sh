#!/bin/bash
## Build the Python environment on a compute node.
##
##     sbatch slurm/submit_setup_env.sh              # full preprocessing stack
##     sbatch -p main slurm/submit_setup_env.sh --run  # Amarel: numpy/h5py/scipy only
##     tail -f slurm/logs/setup_env_<jobid>.out
##
## Wraps setup_env.sh so no pip work happens in the login shell. Needs outbound HTTPS to
## PyPI, which the same compute nodes already have (confirmed by check_cds_access.sh).
##
## Re-run this whenever requirements.txt changes -- the venv is reused, so it is quick.

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH -p SOE_main
#SBATCH -J tc_setup
#SBATCH -t 00:45:00
#SBATCH -o slurm/logs/setup_env_%j.out
#SBATCH -e slurm/logs/setup_env_%j.err

set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
mkdir -p "$REPO_ROOT/slurm/logs"

echo "node    : $(hostname)"
echo "started : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo

bash "$REPO_ROOT/slurm/setup_env.sh" "$@"
status=$?

echo
echo "finished   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit status: $status"
if [ $status -ne 0 ]; then
    echo "Setup did not complete -- most often a missing ~/.cdsapirc. See the log above."
fi
exit $status
