# Preprocessing — ERA5-Land forcing download

Downloads hourly ERA5-Land data for every AmeriFlux station in the ecoregion pairing
lists, in native netCDF, as the first stage of the historical (1985–2014) T&C forcing
pipeline. Nothing is unit-converted here — that happens in the forcing-builder step.

## Setup

```bash
pip install -r requirements.txt
```

Then create a CDS API key file. Register at
<https://cds.climate.copernicus.eu>, accept the **ERA5-Land** licence on the dataset
page (downloads fail with a licence error otherwise), and copy your key from
<https://cds.climate.copernicus.eu/profile> into `~/.cdsapirc`:

```
url: https://cds.climate.copernicus.eu/api
key: <your-api-key>
```

On Windows that path is `C:\Users\<you>\.cdsapirc`.

> The legacy CDS was retired in September 2024. Old-style keys (`UID:key`) and the old
> `.../api/v2` URL no longer work, and `cdsapi < 0.7.2` cannot talk to the current store.

## Usage

```bash
# See exactly what would be requested, without contacting the CDS
python download_era5_land.py --dry-run

# Full run: 118 stations, 1985–2021, hourly
python download_era5_land.py

# A couple of sites, gently
python download_era5_land.py --stations US-Ho1,US-MMS --jobs 2
```

Re-running skips files that already exist, so an interrupted run resumes safely. Use
`--overwrite` to force a refetch.

| Flag | Purpose |
|---|---|
| `--site-list` | Station CSV (repeatable). Defaults to the deciduous + evergreen lists. |
| `--start-year` / `--end-year` | Period. Defaults 1985–2021. |
| `--mode` | `timeseries` (default) or `gridded` fallback. |
| `--stations` | Comma-separated StationIDs to restrict the run. |
| `--jobs` | Concurrent CDS requests (default 4). |
| `--no-dedup` | Download every station separately even when they share a grid cell. |
| `--dry-run` | Print the request plan and exit. |

## What gets downloaded

Seven variables, chosen to cover the T&C forcing fields that ERA5-Land has to supply:

| ERA5-Land variable | short | units | T&C field |
|---|---|---|---|
| `2m_temperature` | `t2m` | K | `Ta` |
| `surface_pressure` | `sp` | Pa | `Pre` |
| `2m_dewpoint_temperature` | `d2m` | K | `Tdew` (→ `ea`) |
| `total_precipitation` | `tp` | m (accumulated) | `Pr` |
| `10m_u_component_of_wind` | `u10` | m s⁻¹ | `Ws` |
| `10m_v_component_of_wind` | `v10` | m s⁻¹ | `Ws` |
| `surface_solar_radiation_downwards` | `ssrd` | J m⁻² (accumulated) | `Rsw` → radiation partition |

Timestamps come from the netCDF `time` coordinate (hourly, UTC).

Three things are deliberately *not* downloaded, and `era5_variables.py` records why so
they don't get re-added by mistake: **longwave** (T&C computes it internally from
`Ta, ea, N`), **cloud cover** (`N` is derived from the clearness index inside the
radiation partition), and **relative humidity** (`ea` comes from `Tdew` via Tetens).

## Output

```
<output dir>/
  US-Ho1/
    US-Ho1_ERA5_Land_2m-temperature.nc         # t2m, d2m
    US-Ho1_ERA5_Land_pressure-precipitation.nc # sp, tp
    US-Ho1_ERA5_Land_radiation-heat.nc         # ssrd
    US-Ho1_ERA5_Land_wind.nc                   # u10, v10
    US-Ho1_ERA5_Land.json                      # metadata sidecar
  ...
  manifest.csv                                 # one row per station
```

**Four files per station, not one.** The CDS time-series collection splits its response
by variable group and returns them zipped together; it does not merge them. The download
step unpacks the archive and routes each member to a predictable filename. All four share
the same hourly time axis, so the forcing builder can open them as one dataset
(`xarray.open_mfdataset`). The `timeseries_group` field in the sidecar maps each variable
to its file.

The output directory defaults to:

```
/vol_efthymios/NFS07/dd1136/T_and_C/input_data/era5_land/
```

Outside the repo, because ~1 GB of netCDF has no business in git. Override with
`--out <path>` for a one-off; to change it permanently, edit `TC_INPUT_DATA` in
`slurm/config.sh`.

The `.json` sidecar records station and grid coordinates, the variable table with units
and time conventions, the CDS collection and citation, the UTC/`DeltaGMT`/`t_bef`/`t_aft`
convention, and which stations shared a download.

## Two things to know before using the data

**Accumulated fields reset daily.** ERA5-Land `tp` and `ssrd` at hour H hold the
accumulation since **00 UTC that day**, not the flux over the preceding hour. They must
be de-accumulated with a within-day difference. Verify this against a known day the
first time you process a batch — getting it wrong silently inflates precipitation and
shortwave through the day.

**Stations are snapped to the 0.1° grid.** ERA5-Land resolves to ~11 km, so requests are
issued at the nearest grid point (offset recorded per station in the sidecar and
manifest). Stations closer together than one cell get identical data: the 18 CHEESEHEAD
towers and the 3 Howland towers each collapse onto a handful of cells. The script
downloads each cell once and copies, which turns 118 stations into 80 CDS requests.
Keep the **true** station coordinates and elevation for solar geometry and site
parameters — only the meteorology is grid-snapped.

## The two modes

`timeseries` (default) uses the CDS point collection
`reanalysis-era5-land-timeseries`: one request per station for the whole 1985–2021
period, so 80 requests total.

`gridded` uses the classic `reanalysis-era5-land` collection with a one-cell `area` box,
one request per station-year — 80 × 37 = 2,960 requests, and much slower. It exists
because the ECMWF product guide states the time-series entry "may be temporarily
disabled or completely deprecated at any point." If the default mode starts failing at
the dataset level, switch with `--mode gridded`; output lands in per-station
subdirectories as `<StationID>/<StationID>_ERA5_Land_<year>.nc`.

## Notes

- Downloads are **not** committed — `preprocessing/data/` is gitignored. Rerun the
  script on the cluster rather than pushing netCDFs through GitHub.
- CDS requests queue server-side; a full run takes hours to days depending on load.
  `--jobs 4` is deliberately modest, since the CDS throttles per user.
- ERA5-Land is land-only. A grid point that resolves over water returns missing values —
  the manifest's `grid_offset_km` column is the first place to look if a station's data
  comes back empty.
