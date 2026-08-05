# Running the preprocessing on the SOE cluster

SLURM wrappers for the Rutgers SOE HPC. Currently covers the ERA5-Land download; the
T&C model job-array wrapper comes later.

## Cluster quick reference

| | |
|---|---|
| Login | `ssh -p 222 <netid>@soemaster2.hpc.rutgers.edu` (port **222**; needs campus network or RU VPN) |
| File transfer | `scp <files> <netid>@soenfs1.hpc.rutgers.edu:` (port **22** — different host and port from login) |
| Home | `/volume/NFS/$USER` — backed up, slow |
| Scratch | `/mnt/beegfs/$USER` — fast, **not** backed up |
| Node-local | `/tmp/$USER/$SLURM_JOB_ID` — fastest, single-node runs |
| CPU partitions | `SOE_main` (new Epyc), `SOE_legacy` (older Xeon) |
| Walltime | 3 days default, up to 14 with `#SBATCH --time=` |
| Python | `ml Python/3.13.7` or `Python/3.14.6`, or conda at `/opt/apps/miniconda3` |

`SOE_nyg` and `--account=nyg` in the SOE docs are for the GPU-owning group — not us. The
general CPU partitions need no `--account`.

LMOD is only initialised for interactive shells that read `.bashrc`. If `module` is not
found, add the snippet from the SOE docs to your `.bashrc`; the scripts here also source
`/opt/apps/lmod/lmod/init/profile` defensively so they work either way.

## First-time setup

```bash
ssh -p 222 <netid>@soemaster2.hpc.rutgers.edu
git clone https://github.com/diogosaraujo/T_and_C_dynamic_LMA.git
cd T_and_C_dynamic_LMA

bash slurm/setup_env.sh          # venv at ~/envs/tc-preproc + cdsapi
```

Then create `~/.cdsapirc` with your key from <https://cds.climate.copernicus.eu/profile>:

```
url: https://cds.climate.copernicus.eu/api
key: <your-api-key>
```

`chmod 600 ~/.cdsapirc`, and accept the ERA5-Land licence on the dataset page.

## Check connectivity before submitting

```bash
bash slurm/check_cds_access.sh
```

**Already verified once (2026-08-04): compute nodes CAN reach the CDS**, confirmed by a
real ERA5-Land retrieval on `soeepyc16`. Rerun the check if downloads start failing, or
after any cluster networking change.

The script probes the login node, then a compute node via `srun`, then performs a real
CDS retrieval — and its verdict depends on that last step, not on the probes. Note
**`curl` is not installed on the SOE nodes**, so the probes use Python `urllib`; an
earlier curl-based version reported "NO route to the CDS" on a node that was downloading
successfully at that moment.

## Submit the download

```bash
sbatch slurm/submit_era5_download.sh
```

80 requests, 1985–2021, all 118 stations, ~1 GB. One node, 4 threads, 7-day walltime.

Extra flags pass straight through to the Python script:

```bash
sbatch slurm/submit_era5_download.sh --stations US-Ho1,US-MMS
sbatch slurm/submit_era5_download.sh --start-year 1985 --end-year 1990
JOBS=2 sbatch slurm/submit_era5_download.sh
```

Monitor:

```bash
squeue -u $USER
squeue -u $USER --start                 # if PENDING, when will it start
tail -f slurm/logs/era5_dl_<jobid>.out
sacct -j <jobid> --format=JobID,JobName,State,Elapsed,MaxRSS
```

If the job dies or hits the walltime, **just resubmit the same command** — completed
files are skipped, so it resumes where it stopped.

## Job array — gridded fallback only

```bash
sbatch slurm/submit_era5_download_array.sh --mode gridded
```

Only worth it for `--mode gridded` (~2,960 requests). For the default timeseries mode the
array is counterproductive: the CDS throttles per user, so extra tasks queue at ECMWF
anyway while holding SLURM allocations.

Array tasks take a round-robin slice of grid cells — never split mid-cell, so the
shared-cell copy always finds its source file. Each writes its own
`manifest.shard<NNN>.csv`. Merge afterwards:

```bash
cd preprocessing/data/era5_land
head -1 manifest.shard000.csv > manifest.csv
tail -q -n +2 manifest.shard*.csv | sort >> manifest.csv
rm manifest.shard*.csv
```

To change the shard count, edit `#SBATCH --array=0-N%M` — `NUM_SHARDS` is derived from
the array bounds, so the two stay in sync automatically.

## Why output goes to the shared filesystem, not `/tmp`

The SOE docs recommend staging to node-local `/tmp` and copying back at the end, which is
right for I/O-heavy compute. This job is different: it is network-bound and resumable.
Writing in place means a walltime kill loses only the file in flight, and a resubmit
continues. Staging to `/tmp` would discard days of downloads if the job were killed
before the copy-back step.

Output lands in `preprocessing/data/era5_land/`, which is gitignored — never push netCDF
to GitHub. Move it to `/mnt/beegfs/$USER` with `OUT_DIR=...` if home quota is tight, but
remember beegfs is not backed up.

## Expected timing

The bottleneck is the ECMWF queue, not the cluster. Observed on 2026-08-04, a trivial
probe request (1 variable, 2 days, 1 point) took **~2 minutes end to end** — ~80 s
queued in `accepted` before it even started running. That latency is per request and
largely independent of size, so the real 37-year × 7-variable requests will take
substantially longer each.

With 80 requests at `--jobs 4`, budget a day or more of mostly-idle walltime. That is why
the walltime is 7 days and `--jobs` stays low.
