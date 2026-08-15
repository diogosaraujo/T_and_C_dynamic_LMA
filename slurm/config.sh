#!/bin/bash
## Shared configuration, sourced by every submit script.
##
## Detects which cluster it is running on and sets the site-specific names, so
## the same scripts work on the SOE HPC and on Amarel without editing.
##
## ONE THING CANNOT BE AUTOMATED: the "#SBATCH -p ..." directive. sbatch parses
## those lines BEFORE the script body runs, so they cannot read a shell variable.
## The directives say SOE_main; on Amarel override on the command line, which
## takes precedence:
##
##     sbatch -p main --array=1-500 slurm/submit_gcm_tc_run.sh
##
## Every script prints the partition it actually landed on, and warns when that
## looks wrong for the host, so a forgotten -p is visible in the first lines of
## the log rather than after the job fails.

# ---------------------------------------------------------------- cluster detect
_host="$(hostname -f 2>/dev/null || hostname)"
case "$_host" in
    *amarel*|amarel*|slepner*|hal*|gpu[0-9]*|cuda[0-9]*|mem[0-9]*|node[0-9]*)
        TC_CLUSTER="${TC_CLUSTER:-amarel}" ;;
    *soe*|soeepyc*|soemaster*|soenfs*)
        TC_CLUSTER="${TC_CLUSTER:-soe}" ;;
    *)
        TC_CLUSTER="${TC_CLUSTER:-unknown}" ;;
esac

if [ "$TC_CLUSTER" = "amarel" ]; then
    # /scratch is 1 TB soft / 2 TB hard and purged after 90 days of inactivity.
    # Model output lives here and is drained to the SOE HPC afterwards.
    TC_ROOT_DEFAULT="/scratch/$USER/T_and_C"
    PARTITION="${PARTITION:-main}"
    # Amarel's MATLAB module naming differs from SOE's; several are tried in
    # order and the first that loads wins.
    # Confirmed from "module avail" on amarel1 (2026-08): R2024a is the default,
    # R2023a also present. Nothing else MATLAB-like exists, so the old guesses at
    # R2022b/R2021a are dropped rather than left to fail silently down the list.
    MATLAB_MODULES="${MATLAB_MODULES:-MATLAB/R2024a MATLAB/R2023a MATLAB}"
    # Amarel's newest is python/3.8.2 (python/2.7.12 is the only other). That is
    # old enough to matter: the preprocessing venv should be built against it, or
    # better, only MATLAB work should run here and the Python stages stay on SOE.
    PYTHON_MODULE="${PYTHON_MODULE:-python/3.8.2}"
    # main allows 3 days; MaxArraySize is 1001 and MaxSubmitPU on 'main' is 500,
    # so arrays must be chunked at 500 (submit_gcm_tc_run.sh takes OFFSET).
    MAX_ARRAY_CHUNK="${MAX_ARRAY_CHUNK:-500}"
else
    TC_ROOT_DEFAULT="/vol_efthymios/NFS07/$USER/T_and_C"
    PARTITION="${PARTITION:-SOE_main}"
    MATLAB_MODULES="${MATLAB_MODULES:-Matlab/2025a Matlab matlab MATLAB}"
    PYTHON_MODULE="${PYTHON_MODULE:-Python/3.13.7}"
    MAX_ARRAY_CHUNK="${MAX_ARRAY_CHUNK:-1000}"
fi

TC_INPUT_DATA="${TC_INPUT_DATA:-$TC_ROOT_DEFAULT/input_data}"
MODEL_RUN="${MODEL_RUN:-$(dirname "$TC_INPUT_DATA")/model_run}"
TC_VENV="${TC_VENV:-$HOME/envs/tc-preproc}"

export TC_CLUSTER TC_INPUT_DATA TC_VENV PARTITION PYTHON_MODULE MODEL_RUN
export MATLAB_MODULES MAX_ARRAY_CHUNK

# ------------------------------------------------------------------ LMOD helper
# LMOD is only initialised for shells that read .bashrc, so batch scripts have to
# source it defensively -- and the init path is not the same on both clusters.
tc_init_lmod() {
    command -v module >/dev/null 2>&1 && return 0
    for p in /opt/apps/lmod/lmod/init/profile \
             /usr/share/lmod/lmod/init/profile \
             /opt/sw/packages/lmod/lmod/init/profile \
             /etc/profile.d/lmod.sh /etc/profile.d/modules.sh; do
        # shellcheck disable=SC1090
        [ -f "$p" ] && . "$p" && command -v module >/dev/null 2>&1 && return 0
    done
    return 1
}

# Load the first MATLAB module that works, and say which. Returns 1 if none do.
tc_load_matlab() {
    tc_init_lmod || true
    if command -v matlab >/dev/null 2>&1; then
        echo "matlab     : $(command -v matlab)  (already on PATH)"
        return 0
    fi
    for m in $MATLAB_MODULES; do
        if module load "$m" >/dev/null 2>&1 && command -v matlab >/dev/null 2>&1; then
            echo "matlab     : $(command -v matlab)  (module $m)"
            return 0
        fi
    done
    echo "ERROR: no MATLAB module loaded. Tried: $MATLAB_MODULES" >&2
    echo "       'module avail 2>&1 | grep -i matlab' will show what exists here." >&2
    return 1
}

# Warn when the partition the job landed on does not match the detected cluster,
# which almost always means a forgotten "-p" on an Amarel submission.
tc_check_partition() {
    local got="${SLURM_JOB_PARTITION:-}"
    [ -n "$got" ] || return 0
    echo "cluster    : $TC_CLUSTER   partition: $got"
    case "$TC_CLUSTER:$got" in
        amarel:SOE_*|soe:main)
            echo "  ! partition '$got' looks wrong for cluster '$TC_CLUSTER'." >&2
            echo "    On Amarel submit with:  sbatch -p $PARTITION ..." >&2 ;;
    esac
}
