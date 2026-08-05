# LMA → Tethys-Chloris (T&C) Project

Context file for Claude Code. This project applies remote-sensing–derived **Leaf Mass per
Area (LMA)** estimates (from a companion PLSR project) into the **Tethys-Chloris (T&C)**
mechanistic ecohydrological model, to investigate how **fixed vs. dynamic LMA** changes
water, energy, and carbon fluxes across forested ecoregions of the CONUS.

---

## 1. Goal & experimental design

- **Core question:** does letting LMA vary in time (vs. holding it fixed) change simulated
  water/energy/carbon fluxes, and by how much?
- **Mechanistic channel:** LMA enters T&C only as **SLA** (`Sl = 1/LMA`), and SLA drives
  `LAI = Sl * B(1)` (leaf area from leaf-carbon pool). It does **not** feed photosynthetic
  capacity — `Vmax` is set directly (`Vmax_H=32` at US_xRM; the `Maximum_Rubisco_Capacity`
  call is commented out). So the LMA signal propagates through **leaf area → light
  interception, transpiring area**, not through Vmax. This is the effect the experiment
  isolates.
- **Scope:** deciduous and evergreen forest ecoregions (EPA Level III × forest type from
  the PLSR skill map), each paired to representative AmeriFlux tower(s) for validation.
- **Periods:** historical **1985–2014** (ERA5-Land forcing) and future **2015–2100**
  (GCM/RCP–SSP forcing). Multiple GCMs × SSP scenarios.
- **Validation:** AmeriFlux tower fluxes (GPP, ET, energy balance). Note the T&C papers
  were driven by *in-situ measured* meteorology; we substitute *reanalysis* (ERA5-Land) —
  a defensible but stated methodological difference, which is exactly why the tower
  comparison matters.

### Reference site (worked example)
US_xRM = high-elevation Colorado Front Range / Niwot area (40.28°, −105.55°, **2753 m**),
~960 mm/yr precip, ~1/3 snow. Modeled as a single **evergreen** conifer canopy (High
vegetation; Low/understory off).

---

## 2. Model overview (T&C)

- Two coupled timescales: **hourly** energy/water budget (`HYDROLOGIC_UNIT` → `SVAT_UNIT`)
  and **daily** vegetation carbon (`VEGGIE_UNIT` → `VEGETATION_DYNAMIC`, via `ode45`).
- Plot-scale = **single grid point**. Two vegetation layers: **High (`_H`)** / **Low
  (`_L`)**, one PFT each. 8 carbon pools `B(1:8)` (leaf, sapwood, fine root, reserve,
  fruit, heartwood, standing dead, idling).
- **SLA injection:** `MAIN_FRAME_SLA.m` overwrites the static `Sl_H` each year with the
  PLSR SLA (`Sl_H = 0.9*Sl_H + 0.1*SLA_ex.SLA_H(yr)`). ⚠️ **Known propagation concern:**
  `Restating_parameters.m` sets `VegH_Param_Dyn.Sl = Sl_H` **once before** the loop —
  verify the yearly SLA actually reaches the `LAI = Sl*B(1)` computation (plot `Sl_H_t` vs
  `LAI_H` to confirm dynamic propagation).
- **SLA units:** T&C needs SLA in **m²/gC**. LMA→SLA is `SLA = 1/(LMA · f_C)` with carbon
  fraction **f_C = 0.5** (decided; confirm with data provider whether LMA is dry-mass or
  carbon based). Paschalis's shared code applied only `SLA = 1/LMA` (missing f_C).

---

## 3. Repository layout & workflow

```
repo/
  preprocessing/   # Python: builds T&C-ready .mat forcing structs + site params
    download_era5_land.py   # DONE: hourly ERA5-Land per AmeriFlux station
    era5_variables.py       #       variable registry (units, T&C field mapping)
  slurm/           # SLURM submit scripts + job-array wrappers
    setup_env.sh, check_cds_access.sh, submit_era5_download*.sh
  T&C/             # MATLAB: T&C source + per-site run dirs (current layout)
  tc_model/        # planned: consolidated T&C source + GO_<site>.m / run_site.m
  config/          # planned: per-site config (coords, PFT, soil depth, scenario)
```

- **Preprocessing in Python** — developed and tested locally, but **executed on the SOE
  cluster** (SLURM) for the real runs: ERA5-Land / GCM pulls, humidity, radiation
  partition port, soil sampling, AmeriFlux metadata → writes the `.mat` structs T&C loads.
  Every step gets a `slurm/` wrapper, not just the MATLAB ones.
- **Modeling in MATLAB** (runs on SOE HPC): T&C proper.
- **git/GitHub is the bridge:** develop in VS Code (+ Claude Code), push; on the cluster
  `git clone`/`git pull` and run MATLAB via SLURM job arrays.

---

## 4. Meteorological forcing (input struct T&C reads)

The forcing `.mat` is **purely an input bundle** (nothing computed by T&C at load except
`Ds = max(esat−ea,0)`). Hourly. **20 of 26** stored variables are actually read.

### Required fields (per site, hourly unless noted)
`Date, Ta, Pre, Tdew, ea, esat, Pr, Ws, N, Ca, SAB1, SAB2, SAD1, SAD2, PARB, PARD` +
scalars `Lat, Lon, DeltaGMT, t_bef, t_aft`.

### Skip (recomputed or unused)
`Ds` (recomputed at load), `Rsw` (total SW — feeds the partition, not read by the run),
`U` (= relative humidity `ea/esat`; never referenced), `Zbas` (used only in radiation
preprocessing, not the run — run's reference height is `zatm=31` in MOD_PARAM),
`id_location` (set in launcher; set it **after** the `load`).

### How each field is produced — historical (ERA5-Land) vs future (GCM/RCP)
| Field | Historical (ERA5-Land) | Future (GCM/RCP) |
|---|---|---|
| Ta | `2m_temperature` (K→°C) | `tas` |
| Pre | `surface_pressure` | barometric from Zbas + mean Ta |
| Tdew | `2m_dewpoint` | from `hurs` (RH) via Magnus |
| ea | Tetens on Tdew | same |
| esat | Tetens on Ta | same |
| Pr | `total_precipitation` | `pr` (kg/m²/s→mm/day→hourly) |
| Ws | 10 m u/v → speed | `sfcWind` |
| N | derived (clearness index) | same routine |
| SAB/SAD/PAR | radiation partition of `ssrd` (hourly) | partition of `rsds` (daily→hourly first) |
| Ca | global annual-mean record (NOAA/Scripps) | SSP CO₂ pathway |

### Humidity (use ONE Tetens formula for both fields)
```
esat = 611.0 * exp(17.27 * Ta   / (Ta   + 237.3))   # Pa
ea   = 611.0 * exp(17.27 * Tdew / (Tdew + 237.3))   # Pa
Ds   = max(esat - ea, 0)
```
`esat←Ta`, `ea←Tdew` (correct: dewpoint is where current vapor saturates). Verified: `esat`
reproduces the stored values to 0.04 Pa; `U = ea/esat` (relative humidity) exactly. For the
GCM path, derive `Tdew` from `hurs` (or `huss`) then use the **same** Tetens set end-to-end
(the shared code mixed 17.27/237.3 for esat with 17.625/243.04 for the RH→Tdew step — a
~0.5°C inconsistency to avoid).

### Radiation partition (shared core, both pipelines)
`Automatic_Radiation_Partition.m` / `C_Automatic_Radiation_Partition.m` (use the vectorized
`C_` version): Gueymard clear-sky 2-band + Slingo(1989) cloudy → `SAB1/SAB2/SAD1/SAD2`,
`PARB/PARD`, and **N**. Inputs: `Date, Lat, Lon, Zbas, DeltaGMT, Pr, Tdew, Rsw`.
- **N is derived**, not downloaded: clearness index `A = Rsw/Rsw_clearsky`,
  `N = ((1/0.75)(1−A))^(1/3.4)`; `N=0` if clear, `N=1` if `Pr>0`, `0.15` where Rsw is NaN.
  No ERA5 cloud product needed. (N's external drivers: **shortwave, precip, dewpoint** +
  date + static site coords.)
- **ERA5-Land path skips the daily→hourly SW disaggregation** (ssrd is already hourly);
  feed hourly `ssrd` straight in. The GCM path (daily `rsds`) first builds a clear-sky-shaped
  hourly SW, then partitions.
- Bands sum exactly to total SW (`Ratio_Evaluator`). PAR ⊂ visible only (≈45% of total SW).
- Longwave is **not** an input — T&C computes it internally: `Incoming_Longwave(Ta,ea,N)`.
  So ERA5 `strd` / GCM `rlds` are unused.

### t_bef / t_aft (sub-hourly sun-averaging offsets) — a per-PRODUCT constant
Window `[H−t_bef, H+t_aft]` over which solar altitude is averaged; corrects for the
dataset's radiation timestamp convention. **ERA5-Land**: accumulations from 00 UTC,
deaccumulated → value at H = **preceding hour** → documented `t_bef=1, t_aft=0` (the
optimizer lands ~`0.75/0.25`, fine). Requires UTC (`DeltaGMT=0`) + correct Lon. **GCM**: no
native convention; impose `0/1` in the disaggregation. Set once per product, reuse across
sites (use the `C_` version's force option).

### CO₂ (Ca)
Globally well-mixed → **same across all sites** at a given time (not site-specific);
**varies in time**; **historical** = NOAA/Scripps global annual mean; **future** = SSP
pathway matching the scenario. Annual-step is fine.

### Elevation (Zbas) datum
Use **orthometric** (a.s.l./geoid, e.g. Copernicus GLO-30 or AmeriFlux metadata), not
ellipsoidal — but sensitivity is tiny (barometric scale height 8434.5 m). Bigger risk is
using the reanalysis grid-cell elevation vs the true site elevation in complex terrain.

---

## 5. Vegetation parameterization

Grounded in `MOD_PARAM_<site>.m`. **Most parameters are PFT-prescribed** (literature/trait
DB), not site-measured. Two layers, `_H`/`_L`, one PFT each; forest sites use High only,
`Ccrown=1`.

### Genuinely site-specific external inputs (short list)
1. **PFT / forest type** → selects the whole prescribed block + phenology switch `aSE`
   (0 evergreen, 1 deciduous, 2 grass) + optical class. **Source: AmeriFlux site description
   (deciduous vs evergreen)** — already shortlisted.
2. **Cover/story:** always forest, `Ccrown=1`, High active (as US_xRM).
3. **Canopy height `hc`:** AmeriFlux **BADM `HEIGHTC`** (partial coverage) else global
   canopy-height product (Potapov 2021 / Simard 2011).
4. **Rooting depth `ZR95`:** not in AmeriFlux → PFT/biome lookup (Schenk & Jackson 2002;
   Canadell 1996) or Fan et al. 2017; often just PFT-prescribed (US_xRM `ZR95_H=800 mm`).
   Keep **`ZR95` ≤ soil column depth**.
5. **SLA (`Sl_H`):** the dynamic PLSR input (see §2). Static MOD_PARAM value is a placeholder
   overwritten yearly.
6. **Initial carbon pools `B_H(1:8)`:** from **spin-up**, not measured (see below). Seed/
   validate with AmeriFlux BADM AGB/LAI or biomass products (NBCD/GEDI).

Everything else (photosynthesis FI/Vmax/a1/Do/CT/DSE/Ha/rjv/Nl/Knit, plant hydraulics
Psi/K/C/Axyl, phenology thresholds, allocation/turnover, optics) = **PFT-prescribed**.

### Spin-up
Seed the slow pools with any plausible guess, run over the forcing **repeatedly** until the
pools stop changing between cycles (state then set by climate+params, not the guess — a wrong
guess only slows convergence). Evergreen wood pools → decades–century of sim time. Spin up to
a **common baseline**, then branch fixed vs dynamic LMA so the treatment is the only
difference. Feed the final state back as the restart initial condition.

---

## 6. Soil / land-surface parameterization

Hydraulics are **derived**, not input: `Soil_parameters.m` (Saxton & Rawls, `SPAR=2`) takes
`Psan/Pcla/Porg` → `Osat, Ks, L, Pe, O33, alpVG, nVG`, thermal (`lan_*, cv_s`).

### Site-specific external pulls (only two big ones)
1. **Soil texture `Psan/Pcla/Porg`:** AmeriFlux BADM soil group where measured (best,
   in-situ), else **POLARIS** (CONUS 30 m, from SSURGO) / SSURGO / **SoilGrids** (global).
   Convert SOC/OM to `Porg` (OM ≈ SOC × 1.72).
2. **Soil depth to bedrock / column depth:** US_xRM uses **1 m free-draining** (`Kbot=NaN`).
   ⚠️ **1–2 m is often too shallow for deep-rooted FOREST sites** → biases toward drought
   stress, under-ET/GPP, contaminating the LMA flux signal. Set per site from depth-to-bedrock
   (Pelletier 2016 / Shangguan 2017 / SoilGrids `BDTICM`) or rooting depth; **deepen forests
   (often 2–5 m)**; keep `ZR95 ≤ column depth`. Watch shallow-water-table/riparian sites
   (free drainage wrong → need water-table lower boundary).

### Soil layering — **DECISION: use LAYERED (depth-resolved) soil**
Mesh `Zs = [0 10 20 50 100 150 200 300 400 500 600 700 800 1000] mm` → **`ms=13` layers
summing to 1 m** (graded: 10 mm at top → 200 mm at bottom, NOT 1 m each). Rule:
`length(Zs) == ms+1`.
Current US_xRM soil is **vertically uniform** (single texture replicated via `*ones(1,ms)`;
no Ks depth-decay). **We will instead build a layered profile:** sample SSURGO/POLARIS/
SoilGrids texture at **each `Zs` interval's depth**, run Saxton & Rawls **per layer**, and
fill `Osat/Ks/alpVG/nVG/...` arrays with **depth-varying** values. Optionally apply the
Ks-with-depth decay (`Ks(z)=Ks0·e^(−0.0011z)`, Lucky Hills style; off by default).

### Prescribed / not site data
Layer mesh, `SPAR`, `Kfc/Phy`, snow/ice physics (`TminS=−0.8, TmaxS=2.8, ros_max, Aice`,
freezing thresholds), interception. Initial soil moisture/SWE/temperatures = spin-up state.

---

## 7. Per-site external data sources (summary)

| Input | Source | Notes |
|---|---|---|
| Meteo (historical) | ERA5-Land (CDS API) | hourly; ssrd, tp, t2m, d2m, sp, 10u/10v |
| Meteo (future) | GCM/RCP–SSP daily | disaggregate to hourly |
| CO₂ | NOAA/Scripps (hist) + SSP pathway (future) | one series, all sites |
| PFT / forest type | AmeriFlux site description | deciduous vs evergreen |
| Canopy height | AmeriFlux BADM `HEIGHTC` → Potapov/Simard fallback | |
| Rooting depth | Schenk & Jackson / Fan 2017 / PFT lookup | ≤ column depth |
| Soil texture | AmeriFlux BADM → POLARIS/SSURGO/SoilGrids | layered by depth |
| Soil depth | Pelletier 2016 / Shangguan 2017 / SoilGrids BDTICM | deepen forests |
| Elevation | Copernicus GLO-30 / AmeriFlux metadata | orthometric |
| AGB / LAI (validation) | AmeriFlux BADM / NBCD / GEDI | seed/validate spin-up |

---

## 8. Compute environment — SOE HPC (Rutgers)

- **Access:** VPN or on-campus. Login: `ssh -p 222 dd1136@soemaster2.hpc.rutgers.edu`
  (port **222** = new cluster; 22 = old, decommissioned). Key auth set up (shared NFS home).
- **Filesystems:** home `/volume/NFS/$USER` (backed up, slow) · scratch `/mnt/beegfs/$USER`
  (fast, NOT backed up) · node-local `/tmp/$USER/$SLURM_JOB_ID`. **File transfer:**
  `soenfs1.hpc.rutgers.edu` **port 22** (FileZilla: soenfs1:22, key `id_rsa`).
- **Data lives OUTSIDE the repo:** model inputs under
  `/vol_efthymios/NFS07/dd1136/T_and_C/input_data/<dataset>/` (`era5_land/` so far), set
  once in `slurm/config.sh` as `$TC_INPUT_DATA` and honoured by the Python scripts. Keeps
  ~1 GB of netCDF out of git — repo checkout is `.../T_and_C/T_and_C_dynamic_LMA`, a
  sibling of `input_data`.
- **Modules (LMOD):** `ml Matlab/2025a` available → run T&C natively (no Compiler needed).
  ⚠️ verify Mapping Toolbox for the geospatial pairing (`matlab -batch "ver"`).
  LMOD is only initialised for shells that read `.bashrc` (snippet sourcing
  `/opt/apps/lmod/lmod/init/profile`) — batch scripts must source it defensively.
- **Python:** `ml Python/3.13.7` or `Python/3.14.6`, or shared conda at
  `/opt/apps/miniconda3` (`conda activate` / `conda activate py3146`). Preprocessing uses
  a venv at `~/envs/tc-preproc` built by `slurm/setup_env.sh`.
- **Partitions:** `SOE_main` (new Epyc) / `SOE_legacy` (older Xeon) for general CPU — no
  `--account` needed. `SOE_nyg` + `--account=nyg` in the SOE docs belongs to the
  GPU-owning group, not us. Interactive: `srun -p SOE_main --cpus-per-task=4 --mem=24G --pty bash`.
- ✅ **Compute nodes CAN reach the CDS** (verified 2026-08-04 on `soeepyc16` via `srun`:
  a real ERA5-Land retrieval succeeded). So download jobs can be submitted normally —
  no proxy, no login-node workaround needed.
- ⚠️ **`curl` is not installed** on the login or compute nodes. Use Python
  (`urllib.request`) for connectivity checks in scripts; a missing-`curl` error looks
  exactly like a firewall block and will send you chasing the wrong problem.
- **Run:** `matlab -nodisplay -nosplash -batch "GO_<site>"`. Single-node CPU on **SOE_main**.
  Stage inputs to `/tmp` scratch, copy `RES_*.mat` back. Default runtime 3 days (→14 with
  `#SBATCH --time=`).
- **Fleet:** SLURM **job array** — `#SBATCH --array=1-N`, `matlab -batch
  "run_site($SLURM_ARRAY_TASK_ID)"`; index each ecoregion/site/scenario.

---

## 9. Known issues & caveats

- **SLA propagation** (§2): confirm yearly SLA reaches `LAI=Sl*B(1)` (`VegH_Param_Dyn.Sl`
  set once before loop).
- **`Rd` spikes / `Rh` NaN** in results are **cosmetic** numerical-solver artifacts (mass
  still conserves; CK1 ≈ machine precision) — do not affect the valid carbon/energy/ET
  results. Optional fix: tighten the soil ODE solver.
- **Soil depth for forests** (§6): don't inherit 1 m blindly.
- **Reanalysis vs in-situ forcing** (§1): state it; validate against towers.
- **f_C = 0.5** (§2): confirm LMA basis with data provider.

## 10. Open TODOs
- [ ] Extend AmeriFlux fetch: `HEIGHTC`, `AG_BIOMASS`, `LAI`, BADM soil texture.
- [ ] POLARIS/SoilGrids depth-resolved texture + depth-to-bedrock sampling.
- [ ] Verify Mapping Toolbox on the cluster.
- [x] Verify compute-node outbound access to the CDS — confirmed working 2026-08-04.
- [x] ERA5-Land download (`preprocessing/download_era5_land.py` + `slurm/` wrappers).
- [ ] Build remaining Python preprocessing modules (GCM, radiation port, soil, params).
- [ ] Confirm the unit T&C expects for `Pre` (Pa vs mbar) against the US_xRM forcing.
- [ ] `run_site.m` + SLURM job-array wrapper.
- [ ] Spin-up protocol (common baseline, then fixed/dynamic branch).
