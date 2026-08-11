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
  fraction **f_C = 0.5**. ✅ **CONFIRMED 2026-08-11: the PLSR LMA is DRY MASS**, so the 0.5
  is required to reach gC. Not an open question any more — it is a factor of 2 on every
  station's leaf area, so do not revisit it casually. Paschalis's shared code applied only
  `SLA = 1/LMA` (missing f_C), which would double SLA and hence LAI.

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
  ~2 GB of netCDF out of git — repo checkout is `.../T_and_C/T_and_C_dynamic_LMA`, a
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
- ⚠️ **Nothing runs in the login shell** — every step (downloads, verification, even the
  `pip install`) goes through `sbatch`. Each preprocessing step has a wrapper in `slurm/`;
  `check_cds_access.sh` is invoked directly but only orchestrates `srun`.
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
- **6 stations ON HOLD (deferred 2026-08-10), 12 arms.** Everything else runs.
  - `fzero` in `Surface_Temperature_Snow` (`SVAT_UNIT` line 333), *"Initial function
    value must be finite and real"*: **US-NMj, US-Wi1, US-Wi2, US-Wi4, US-xSB**.
    Something reaching the snow surface-temperature solve is non-finite; check the
    forcing for NaNs at these sites before theorising.
  - **US-DPP** does not crash, it *stalls*: 8 h wall clock to reach `Iter: 2`,
    against 0.28–0.92 h for a whole run elsewhere. More wall time will not help.
    It also has 7 soil layers sitting at the no-silt boundary, which would push
    Saxton & Rawls to extreme conductivity and could be starving the soil ODE
    solver's step size — a hypothesis, not a diagnosis.
- **`B(6)` heartwood is INERT in this configuration — do not spend effort on it.**
  It is 0 in both IC vectors (Dr. Paschalis's US_xRM choice). That makes initial
  `TBio = 0.02*(B1+B2+B3+B4+B6)` ≈ **19–21 t DM/ha** against AmeriFlux BADM observed
  means of **203 (deciduous, n=8)** and **173 (evergreen, n=12)** t DM/ha
  (`badm_biomass.csv`, job 36168) — an 8–11x gap that looks alarming and is not.
  `TBio` has exactly one consumer, `Allocation_Coefficients`, which uses it only
  inside `if aSE ~= 2 && OPT_VCA >= 1`. **`OPT_VCA = 0` here**, so `so = 0.3` is
  constant and `TBio` changes nothing. `B(6)` feeds nothing else either
  (`OPT_SoilBiogeochemistry = 0` disables `BIOGEOCHEMISTRY_DYNAMIC3`). So the pool
  is a diagnostic accumulator only: it has no effect on any flux, pool or
  allocation coefficient. ⚠️ If `OPT_VCA` is ever switched on, this reverses and
  `B(6) = max(0, 50*TBio_target − (B1+B2+B3+B4))` becomes necessary.
  Related: heartwood never leaves (`Wm = 0` for all 8 PFTs), so `TBio` grows without
  bound — US-Ha2 passes the observed mean at year 19.6 and ends at 176% of it. Do
  not read modelled biomass as a validation target.
- **Soil depth for forests** (§6): don't inherit 1 m blindly.
- **Reanalysis vs in-situ forcing** (§1): state it; validate against towers.

## 10. Open TODOs
- [x] AmeriFlux BASE + BADM fetch (`preprocessing/download_ameriflux.py`) and a coverage
      report mapping BADM → T&C parameters (`inspect_ameriflux_badm.py`).
- [ ] Decide per-site parameter sources from `badm_coverage.csv` (which sites use in-situ
      values vs PFT/gridded fallbacks) — this split belongs in the methods.
- [ ] POLARIS/SoilGrids depth-resolved texture + depth-to-bedrock sampling.
- [ ] Verify Mapping Toolbox on the cluster.
- [x] Verify compute-node outbound access to the CDS — confirmed working 2026-08-04.
- [x] ERA5-Land download (`preprocessing/download_era5_land.py` + `slurm/` wrappers).
- [ ] Build remaining Python preprocessing modules (GCM, radiation port, soil, params).
- [ ] Confirm the unit T&C expects for `Pre` (Pa vs mbar) against the US_xRM forcing.
- [ ] `run_site.m` + SLURM job-array wrapper.
- [ ] Spin-up protocol (common baseline, then fixed/dynamic branch).
- [x] AmeriFlux BADM biomass coverage (`preprocessing/check_badm_biomass.py`,
      `slurm/submit_check_badm_biomass.sh`): 20/110 stations report usable standing
      tree AGB. Not wired into any input — biomass is a **validation** dataset here,
      not a parameter source. Richest unused field is **`BASAL_AREA` (315 values)`**,
      which pairs with T&C's `BA_H` output; also `AG_LIT_BIOMASS`, `WD_BIOMASS_*`,
      `SOIL_STOCK_C_ORG`, and `LAI` at 39 stations (the direct check on `LAI=Sl*B(1)`).
- [ ] **On hold:** the 5 `fzero` snow stations + US-DPP's stall (§9). 91/98 stations
      (182 arms) are complete without them, and both affected ecoregions survive the
      loss: *Northern Lakes and Forests* (US-Wi1, US-Wi2, US-Wi4, US-NMj) still has
      ~15 complete stations, and *Southern Coastal Plain* (US-DPP, US-xSB) keeps
      US-SP1/SP2/SP3/SP4. So this costs no ecoregion coverage.
- [x] **Ecoregion pairing VERIFIED (2026-08-10): 118/118 stations fall inside the EPA
      Level III polygon they are paired to.** `preprocessing/verify_ecoregion_pairing.py`
      (+ `slurm/submit_verify_pairing.sh`) does point-in-polygon against `us_eco_l3`,
      projecting with Albers parameters read from the shapefile's own `.prj` (matches
      pyproj to <1 mm) and aborting if control stations fail. Re-run it whenever the
      pairing file changes. Three apparent anomalies were checked and are all correct:
      US-NMj is *Northern* **Michigan** Jack Pine (not New Mexico); Southeastern Plains
      (65) genuinely reaches Maryland, so SERC/US-xSE belongs there; and US-HB5 (75) vs
      US-Sx2 (63) really are split by the boundary despite being 25 km apart.
      ⚠️ The old `gaftp.epa.gov` download path now 404s — EPA serves from S3.
