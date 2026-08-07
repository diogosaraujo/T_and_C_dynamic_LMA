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
| Input data | `/vol_efthymios/NFS07/dd1136/T_and_C/input_data` — outside the repo, set in [config.sh](config.sh) |
| `curl` | **not installed** — use Python `urllib` for connectivity checks |

All shared paths live in [config.sh](config.sh), sourced by every script here. Change a
path there once rather than editing each script.

`SOE_nyg` and `--account=nyg` in the SOE docs are for the GPU-owning group — not us. The
general CPU partitions need no `--account`.

LMOD is only initialised for interactive shells that read `.bashrc`. If `module` is not
found, add the snippet from the SOE docs to your `.bashrc`; the scripts here also source
`/opt/apps/lmod/lmod/init/profile` defensively so they work either way.

## First-time setup

```bash
ssh -p 222 <netid>@soemaster2.hpc.rutgers.edu
git clone git@github.com:diogosaraujo/T_and_C_dynamic_LMA.git
cd T_and_C_dynamic_LMA

sbatch slurm/submit_setup_env.sh          # venv at ~/envs/tc-preproc
tail -f slurm/logs/setup_env_<jobid>.out
```

> **Everything runs through SLURM.** No work belongs in the login shell — not the
> downloads, not the verification, not even the `pip install`. Every step here has a
> submit script; `check_cds_access.sh` is the one script you invoke directly, and it only
> orchestrates `srun` calls that execute on compute nodes.

Re-run `submit_setup_env.sh` whenever `requirements.txt` changes; the venv is reused.

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

All work runs on compute nodes via `srun` — the script only orchestrates. It probes
reachability, then performs a real CDS retrieval, and its verdict depends on that last
step rather than the probe. Note **`curl` is not installed on the SOE nodes**, so the
probe uses Python `urllib`; an earlier curl-based version reported "NO route to the CDS"
on a node that was downloading successfully at that moment.

## Submit the download

```bash
sbatch slurm/submit_era5_download.sh
```

80 requests, 1985–2021, all 118 stations, ~2 GB, roughly 30 min. One node, 4 threads, 7-day walltime.

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

## Verify the download

```bash
sbatch slurm/submit_era5_verify.sh
```

Subset:

```bash
sbatch slurm/submit_era5_verify.sh --stations US-HBK,US-Ha2
```

Checks time axis, variable presence, coverage, physical ranges, and empirically confirms
the `tp`/`ssrd` accumulation convention. Writes a JSON report next to the data and exits
non-zero if any station fails. See [preprocessing/README.md](../preprocessing/README.md)
for what each check catches.

This needs `xarray` and `netCDF4`, added to `requirements.txt` after the download step
was built — **rerun `bash slurm/setup_env.sh`** to install them into the existing venv.

## AmeriFlux measurements + BADM

```bash
export AMF_USER_ID=<ameriflux username> AMF_USER_EMAIL=<email>
sbatch slurm/submit_ameriflux_download.sh --agree-policy --is-test --stations US-HBK,US-Ha2
sbatch slurm/submit_ameriflux_download.sh --agree-policy          # full run
```

Credentials travel via the submitting environment (SLURM copies it into the job) — never
put them in `config.sh` or any tracked file. Output goes to
`$TC_INPUT_DATA/ameriflux/`. The submit script runs the BADM inspector afterwards, so the
parameter-coverage table appears in the job log.

`--is-test` suppresses the emails AmeriFlux sends to site teams; use it while testing,
drop it for the real run. See [preprocessing/AMERIFLUX.md](../preprocessing/AMERIFLUX.md).

## Canopy height (GEDI/Landsat)

```bash
sbatch slurm/submit_canopy_height.sh                          # all stations
sbatch slurm/submit_canopy_height.sh --stations US-HBK,US-Ha2 # a subset
```

Samples Potapov et al. (2021) UMD GLAD 30 m canopy height — GEDI lidar calibrated onto
Landsat — at every station. The multi-GB continental mosaic is read over `/vsicurl/`, so
only the scanlines covering each station transfer; ~5 minutes for 118 stations.

Output goes to its own folder outside the repo, `$TC_INPUT_DATA/canopy_height/`, as
`canopy_height_gedi.csv` plus a JSON provenance sidecar.

This step **only downloads**. Where the AmeriFlux download is present it adds a validation
column comparing against measured BADM `HEIGHTC`, but it does not merge them — the
selection rule (BADM where present, GEDI as fallback) is applied later, when the `.mat`
files are built.

Needs `rasterio`, added to `requirements.txt` for this step: **re-run
`sbatch slurm/submit_setup_env.sh`** before the first use. The job checks for it and
exits with that instruction if it is missing.

## Rooting depth (Schenk & Jackson)

Needs a free [NASA Earthdata Login](https://urs.earthdata.nasa.gov), once:

```bash
printf 'machine urs.earthdata.nasa.gov login <user> password <pass>\n' >> ~/.netrc
chmod 600 ~/.netrc
sbatch slurm/submit_root_depth.sh
```

`EARTHDATA_USER`/`EARTHDATA_PASS` or `EARTHDATA_TOKEN` work too. The job checks before
doing any work and exits with instructions if none are visible.

Samples the ISLSCP II Ecosystem Rooting Depths grids (doi:10.3334/ORNLDAAC/929) and
writes `ZR95_H` and `ZR50_H` in mm to `$TC_INPUT_DATA/root_depth/`.

⚠️ **The grid is 1° (~110 km)** — only 54 distinct cells across 118 stations, so 15
CHEESEHEAD towers share one value, as do 8 Wisconsin and 7 Metolius sites. The output
carries `n_stations_sharing_cell` so this is visible rather than implied. Treat it as a
regional lookup, not a site measurement. **BADM has no rooting-depth variable at all**, so
unlike canopy height there is no measured alternative for any station.

⚠️ T&C aborts if `ZR95_H` exceeds the deepest `Zs` layer, so this has to be reconciled
with the soil-depth step when the `.mat` files are built.

## Time-varying LMA (PLSR temporal CV)

```bash
sbatch slurm/submit_lma_input.sh --reconstruct             # all stations, gap-free
sbatch slurm/submit_lma_input.sh --stations US-HBK,US-Ha2  # a subset
sbatch slurm/submit_lma_input.sh --dry-run                 # resolve inputs only
sbatch slurm/submit_lma_input.sh --audit                   # coverage, writes no series
```

Reads the PLSR temporal-CV output already on the cluster — nothing is downloaded:

```
$PLSR_ROOT/eco<ii>/time/PLSR_predictions_eco<ii>_<forest>_oofcv.mat
$PLSR_ROOT/eco<ii>/time/PLSR_fitting_coeff_eco<ii>_<forest>_oofcv_TEMPORAL.mat
$PREDICTOR_ROOT/LMA_ecoregion_no<ii>.csv
```

with the defaults

```
PLSR_ROOT      = $ECOREGION_ROOT/PLSR_temporal_cv_pixel_climatology_DOY/LMA
PREDICTOR_ROOT = $ECOREGION_ROOT/PLSR_inputs_pixel_climatology_DOY/LMA
```

`PREDICTOR_ROOT` is `opts.InputDir` from the MATLAB pipeline — the table the fit
itself read — and it supplies **both** the pixel coordinates and the predictors, so
the two can never come from different vintages. The older `ecoregion_no<ii>.csv` at
the ecoregions root carries no SSRD columns, which every fit selected as a champion;
it is used only when nothing else is present, and the log says so. Override any of
the three with `PLSR_ROOT=... ECOREGION_ROOT=... PREDICTOR_ROOT=... sbatch ...`.

Two series come out per station — `Y_plot_abs` (observed) and `yfit_plot_abs` (PLSR
modelled) — into `$TC_INPUT_DATA/lma/<station>/`, in **g/m²**.

`--reconstruct` rebuilds the modelled series by re-applying
`PLSR_fitting_coeff_*_TEMPORAL.mat` to every row of the predictor table instead of
reading the stored `yfit_plot_abs`. Use it: the pipeline drops rows whose response is
missing *before* predicting, so `yfit_plot_abs` exists only where an observation did
and the stored series is **not** gap-free. Only pixels in `uniq_pix_final` (those that
contributed to the fit) are used.

Each station takes the nearest ecoregion pixel **that appears in the predictions**;
matching against the ecoregion at large could land on a pixel the model never saw,
since `predict_with_fit` drops pixels absent from the training set. Where that pixel
has no usable data, or sits more than `--max-distance-km` away (default 50), the
ecoregion median is used instead and every row records which it was. `--fill-gaps`
extends the fallback to individual missing years — off by default, because the filler
comes from a different pixel with a different mean and introduces a step change.

Needs `h5py`: the PLSR files are MATLAB `-v7.3`, i.e. HDF5, which `scipy.io.loadmat`
cannot read. **Re-run `sbatch slurm/submit_setup_env.sh`** before first use; the job
checks and exits with that instruction if it is missing.

⚠️ The SLA conversion is deliberately **not** applied here. `Sl = 1/(LMA·f_C)` happens
at `.mat` build time so that revisiting `f_C` never means re-reading the PLSR output.
Note the shipped `LMA_US_xRM.mat` uses `SLA_H = 1/LMA` with **no** `f_C` — under the
project's `f_C = 0.5` every value there should double.

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

## Where the data lands

```
/vol_efthymios/NFS07/dd1136/T_and_C/input_data/era5_land/
    US-Ho1/
        US-Ho1_ERA5_Land_2m-temperature.nc          # t2m, d2m
        US-Ho1_ERA5_Land_pressure-precipitation.nc  # sp, tp
        US-Ho1_ERA5_Land_radiation-heat.nc          # ssrd
        US-Ho1_ERA5_Land_wind.nc                    # u10, v10
        US-Ho1_ERA5_Land.json                       # metadata sidecar
    ...
    manifest.csv
```

Four netCDFs per station, one per variable group — the CDS returns the time series split
that way rather than merged. All four share the same hourly time axis.

Outside the repo, so ~2 GB of netCDF can never be accidentally staged and pushed. The
path is set once in [config.sh](config.sh) as `TC_INPUT_DATA`; change it there if it ever
needs to move.

## Expected timing

Measured from job 35235 (2026-08-04), which completed all 320 CDS retrievals server-side:

| | |
|---|---|
| One full 37-year, 7-variable request | **~82 s** (queue + generate + transfer) |
| Payload per station | **~18.6 MB** compressed (4 netCDFs in one zip) |
| Full run: 80 requests at `--jobs 4` | **~30 min** |
| 2-station test | **~2 min** |

Add SLURM queue time before the job starts, which on `SOE_main` is usually short.

The 7-day walltime is deliberate headroom, not an expectation — CDS load varies and a
slow day should not kill a run mid-way. Latency is dominated by ECMWF queueing, so
raising `--jobs` past 4 buys little: the CDS throttles per user.
