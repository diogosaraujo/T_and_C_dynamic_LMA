# Shared configuration for SOE cluster runs. Sourced by the other scripts in slurm/.
# Edit paths HERE, not in the individual scripts.

# Root for downloaded/derived model INPUT data, one subfolder per dataset
# (era5_land/ so far). Outside the git repo: the ERA5-Land download alone is ~1 GB,
# and GitHub rejects files over 100 MB.
TC_INPUT_DATA="${TC_INPUT_DATA:-/vol_efthymios/NFS07/dd1136/T_and_C/input_data}"

# Python environment built by setup_env.sh.
TC_VENV="${TC_VENV:-$HOME/envs/tc-preproc}"

# General-use CPU partition. SOE_legacy also works; SOE_nyg belongs to the GPU group.
PARTITION="${PARTITION:-SOE_main}"

# LMOD module providing python3.
PYTHON_MODULE="${PYTHON_MODULE:-Python/3.13.7}"

# Exported so download_era5_land.py picks up TC_INPUT_DATA as its default output root.
export TC_INPUT_DATA TC_VENV PARTITION PYTHON_MODULE

# AmeriFlux credentials are NOT set here -- never commit them. Export them in your shell
# before submitting; SLURM copies the submitting environment into the job:
#     export AMF_USER_ID=<your ameriflux username>
#     export AMF_USER_EMAIL=<your email>
# Put those two lines in ~/.bashrc on the cluster if you would rather not repeat them.
