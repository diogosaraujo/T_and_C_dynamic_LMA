# AmeriFlux measurements + BADM metadata

Second preprocessing stage. Two jobs:

1. **`download_ameriflux.py`** — pulls the BASE-BADM product (half-hourly/hourly flux and
   met measurements, plus that site's BADM metadata) for the stations in the ecoregion
   pairing CSVs, and saves public site-registry metadata for all of them.
2. **`inspect_ameriflux_badm.py`** — reports which T&C parameters the BADM can actually
   supply, per station, and what the fallback is where it cannot.

Output goes outside the repo, next to the ERA5-Land data:

```
/vol_efthymios/NFS07/dd1136/T_and_C/input_data/ameriflux/
    site_metadata.json / .csv     registry info for every requested station
    US-HBK/
        AMF_US-HBK_BASE_HH_<ver>.csv     measurements
        AMF_US-HBK_BIF_<policy>_<ver>.xlsx  BADM metadata
    _archives/                    the downloaded zips, kept so a re-unpack needs no refetch
    files_manifest.csv            every unpacked file, labelled base / badm / doc / other
    download_provenance.json      what was requested, by whom, under which policy
    badm_coverage.csv             station x parameter -> value found
    badm_values.csv               matched variables, long format
    badm_inventory.csv            EVERY variable reported, matched or not
```

## Account and policy

The public site registry needs no account:

```bash
sbatch slurm/submit_ameriflux_download.sh --metadata-only
```

Downloading measurements does. Register free at <https://ameriflux.lbl.gov>, then:

```bash
export AMF_USER_ID=<your ameriflux username>
export AMF_USER_EMAIL=<your email>
```

Downloading also means accepting the [AmeriFlux data use
policy](https://ameriflux.lbl.gov/data/data-policy/), so `--agree-policy` is required and
is recorded in `download_provenance.json`. Sites are shared under either **CC-BY-4.0** or
the **AmeriFlux Legacy** policy; the script splits the station list by policy
automatically (asking for the wrong one silently returns nothing for that site).

Under the legacy policy you are expected to notify site PIs and offer co-authorship where
appropriate. `download_provenance.json` records the per-site policy so this is traceable
when the paper is written.

## Test runs

```bash
# two stations, flagged as a test so site teams are not emailed
sbatch slurm/submit_ameriflux_download.sh --agree-policy --is-test --stations US-HBK,US-Ha2
```

The submit script runs the inspector afterwards, so the coverage table lands in the same
job log. All work runs on compute nodes -- nothing is executed in the login shell.

**Use `--is-test` while shaking out the pipeline.** The service emails site teams when
their data is downloaded; the flag exists so repeated test pulls don't spam them. Drop it
for the real run — that notification is the courtesy the policy expects.

`--dry-run` prints the exact request payload without contacting the service or needing
credentials.

## What the inspector looks for

Each target ties back to a parameter in CLAUDE.md §5–§7:

| Parameter | T&C use | Fallback when a site doesn't report it |
|---|---|---|
| `canopy_height` | `hc` | Potapov 2021 / Simard 2011 |
| `lai` | validates `LAI = Sl·B(1)`; checks spin-up | spin-up equilibrium only |
| `biomass` | seeds/validates `B_H(1:8)` | spin-up guess; NBCD/GEDI cross-check |
| `soil_texture` | `Psan`/`Pcla` for Saxton & Rawls | POLARIS / SSURGO / SoilGrids |
| `soil_chem` | `Porg` (OM = SOC × 1.72) | POLARIS / SoilGrids SOC |
| `soil_depth` | soil column depth | Pelletier 2016 / Shangguan 2017 |
| `root` | `ZR95` | Schenk & Jackson / Fan 2017 (rarely reported) |
| `species` | sanity-check the PFT and deciduous/evergreen split | IGBP alone |
| `disturbance` | flags harvest/fire sites | none — unflagged disturbance contaminates the LMA signal |
| `elevation` | `Zbas` | Copernicus GLO-30 |
| `igbp` | PFT and the `aSE` phenology switch | already in the pairing CSVs |
| `utc_offset` | aligns tower data with UTC ERA5-Land forcing | infer from longitude |

**It discovers rather than assumes.** BADM naming varies between sites and isn't fully
documented in one place, so the inspector matches on patterns and separately inventories
everything it *didn't* match, into `badm_inventory.csv`. A renamed or unexpected variable
shows up there as unmatched rather than being silently reported as absent — check that
file before concluding a parameter is unavailable.

Values of `-9999` are treated as missing, so a site that reports a variable with no data
counts as absent rather than supplying a nonsense number.

## Expected reality

Don't expect full coverage. BADM is contributed voluntarily and unevenly: canopy height
and IGBP are common, soil texture and biomass less so, rooting depth and depth-to-bedrock
rarely. The point of `badm_coverage.csv` is to make the split explicit — which sites are
parameterised from in-situ measurements, and which fall back to PFT lookups or gridded
products. That distinction belongs in the paper's methods, since it varies per site and
affects how much weight each tower's validation carries.

## On the cluster

```bash
export AMF_USER_ID=... AMF_USER_EMAIL=...
sbatch slurm/submit_ameriflux_download.sh --agree-policy --is-test --stations US-HBK,US-Ha2
```

The submit script runs the download and then the inspector, so the coverage table lands
in the job log. Credentials travel via the submitting environment — never put them in
`config.sh` or any tracked file.
