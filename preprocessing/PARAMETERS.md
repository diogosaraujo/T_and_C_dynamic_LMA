# T&C parameter lookup tables

Extracted from `T&C/Thanos_US_xRM/MOD_PARAM_US_xRM.m`, the working evergreen-conifer
parameterisation (US_xRM, Colorado Front Range, modelled as *Picea mariana*). The intent
is to reuse this set at the other **evergreen** stations, so a new site only needs its
forcing, soil and a few measured traits.

| File | Contents |
|---|---|
| `tc_evergreen_pft_parameters.csv` | 180 rows — the values to **reuse** at other evergreen sites |
| `tc_site_specific_parameters.csv` | 40 rows — what is deliberately **excluded**, and where each comes from instead |
| `tc_site_parameters_HBK_Ha2.csv` | 234 rows — **proposed** runnable sets for the two test stations (see below) |

## Proposed sets for US-HBK and US-Ha2

`tc_site_parameters_HBK_Ha2.csv` is long-format (`station_id, variable, proposed_value,
units, description, category, basis, source, confidence, notes`) so it extends to the
other 116 stations without changing shape. 117 rows per station.

| | **US-HBK** | **US-Ha2** |
|---|---|---|
| Site | Hubbard Brook, NH — 367 m, MAT 6 °C, MAP 1400 mm | Harvard Forest Hemlock, MA — 360 m, MAT 6.56 °C, MAP 1071 mm |
| Vegetation | northern hardwood (beech / sugar maple / yellow birch) | eastern hemlock (*Tsuga canadensis*) |
| Donor set | **`MOD_PARAM_ZURICH_SMA.m` High-layer deciduous block** | **`MOD_PARAM_US_xRM.m` evergreen conifer** |
| `aSE_H` | 1 (seasonal) | 0 (evergreen) |
| Optics | `Veg_Optical_Parameter(7)` BDT temperate | `Veg_Optical_Parameter(1)` NET temperate |

**41 vegetation parameters differ** between the two proposed sets — the whole phenology,
turnover, hydraulics and photosynthesis block. Non-vegetation parameters (snow density,
ice albedo, soil mesh, domain geometry) stay on one common baseline for both: those are
site/physics properties, and importing Zurich's just because it is the deciduous donor
would drag in an unrelated Swiss grassland configuration.

Range check against Table 3: **the deciduous set is clean except `rjv_H = 2.8`** (12%
above the 1.5–2.5 range — flagged, not changed, since it is the T&C authors' own value).
The evergreen set keeps its three known flags, and note `dc_C_H` stays harmless at US-Ha2
for the same reason as US_xRM: `Tcold_H = −50` never triggers cold shedding.

**Four values are mine, not sourced** — all marked `proposed_estimate` or
`needs_verification`:

- `Vmax_H = 45` at US-HBK. Zurich sets 0 (layer off) so a value *must* be supplied; 45 is
  a mid-range temperate-deciduous choice from Table 3's 10–150. **Calibrate against tower
  GPP** once the AmeriFlux data land.
- `PsiX50_H = −5.0` at US-Ha2, inherited and **probably too negative**: that is a boreal
  spruce value, *Tsuga canadensis* is notably drought-sensitive, and T&C's own temperate
  value is −3.5. Check Choat et al. (2012) before running — this controls simulated
  drought mortality.
- `ZR95_H = 800 mm` at both — PFT placeholder as agreed.
- `rjv_H = 2.8` at US-HBK, as above.

### What the downloaded BADM actually supplied

Coverage is **extremely uneven** between the two stations — 168 BADM rows at US-Ha2
against 68 at US-HBK — and that asymmetry is itself a result worth stating in the methods.

**US-Ha2 — nine parameters now come from site measurements:**

| | Measured | Was going to use | Impact |
|---|---|---|---|
| `Sl_H` | **0.01527** (LMA = 131 g/m²) | 0.017 (US_xRM) | in-situ ground truth for the study's central variable |
| `hc_H` | **21 m** | 8 m (US_xRM) | US_xRM's value would have been badly wrong |
| `LAI_H` | **4.4** | 4.03 | spin-up convergence target |
| `Nl_H` | **~35** (foliage N 1.18–1.67%) | 62 | boreal spruce value replaced; also brings it inside Table 3's 15–42 |
| `TBio_H` | **260 ton DM/ha** | 300 | seeds the carbon pools |
| `Slo_top` | **0** (TERRAIN = Flat) | 0 | confirms the flat-plot assumption |
| `Axyl_H` | SA_MAX = 230, **unit unstated** | 15.0 | check the unit — could be 15× |

`Sl_H = 1/(131 × 0.5) = 0.01527` lands within 10% of US_xRM's 0.017. That is a genuine
independent check on the whole LMA→SLA chain **including the f_C = 0.5 assumption**, and
the same number can validate the PLSR LMA product at this site.

Soil chemistry is also there — SOC 47–103 g/kg by horizon, bulk density, pH, horizon
depths to 35 cm, sourced to Compton & Boone (2000) — enough to build a layered `Porg`
profile without falling back to SoilGrids. Texture (`Psan`/`Pcla`) is still absent.

**US-HBK supplied almost nothing**: elevation, IGBP, UTC offset, disturbance, terrain. No
canopy height, LAI, biomass, soil, or species. It falls back to PFT defaults and gridded
products throughout — so the deciduous set stands or falls on the Zurich block plus
calibration against tower fluxes.

### Two site facts that change how these runs should be set up

⚠️ **US-HBK is not flat.** `TERRAIN = "Significant Slope (>5%)"`, `ASPECT = N`. US_xRM and
US-Ha2 are flat, and the plot-scale configuration assumes it (`Slo_top = 0`, `SvF = 1`). A
north-facing slope changes incoming shortwave and the sky-view factor, and enables lateral
subsurface flow. Either set `Slo_top`/`SvF`/`Asur` from a DEM, or state the flat-plot
assumption as a limitation. `SITE_SNOW_COVER_DAYS = 120` gives a free check on the snow
module.

⚠️ **US-Ha2's disturbance record reads `Pests and disease`** — the Harvard Forest hemlock
stand is affected by hemlock woolly adelgid. That matters more than a parameter: an
insect-driven canopy decline produces a multi-year LMA and LAI trajectory that has
**nothing to do with climate**. Using this site to test whether dynamic LMA changes fluxes
risks attributing an infestation signal to the LMA treatment. It is still a fine pipeline
test site; it is a questionable *inference* site, and the confound should be checked
against the flux record before it carries any weight in the results.

Still to come from elsewhere: `zatm`, soil texture, depth to bedrock, and all spin-up
state.

Together they account for every prescribed assignment in `MOD_PARAM_US_xRM.m` (281
assignments; 41 are computed inside the file and listed as derived). That completeness is
verified, not assumed — nothing in the source file is silently unaccounted for.

## Columns

`variable` (name exactly as in the model) · `value` · `units` · `description` ·
`category` · `layer` (H / L / both / n/a) · `active_at_US_xRM` · `source` ·
`source_confidence` · `notes`

### `source_confidence` — read this before citing anything

| Value | Meaning | Rows |
|---|---|---|
| `code_comment` | Stated verbatim in `MOD_PARAM_US_xRM.m` | 111 |
| `tc_source_code` | Read directly out of a T&C function (`Veg_Optical_Parameter.m`, `Root_Fraction_General.m`) | 4 |
| `tc_model_description` | Attributed to the T&C model description (Fatichi et al. 2012a) as the governing reference | 53 |
| `formulation_origin` | The underlying published formulation the parameter belongs to (e.g. Leuning 1995 for `a1`, `Do`, `go`) | 4 |
| `needs_verification` | **I could not attribute this confidently** — check the T&C manual or ask Fatichi/Paschalis before citing | 8 |

The `tc_model_description` rows name the right governing paper but I have **not** verified
each individual value against its table. Treat them as "correct source, unverified value"
until someone checks them against the T&C documentation. Do not paste these citations
into a manuscript unchecked.

## Two things found while extracting

**`DSE_H` and `Ha_H` are mislabelled in the source — RESOLVED.** The code reads:

```matlab
DSE_H = [0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H  = [89];     %% [kJ / mol K]  entropy factor - Plant Dependent
```

Table 3 of Fatichi et al. (2012a) settles it: `Ha` (activation energy) has range
40–95 kJ/mol and `ΔS` (entropy factor) 0.625–0.665 kJ/(mol K). Both US_xRM values sit
inside their correct ranges — 89 within 40–95, and 0.649 within 0.625–0.665. **The values
are right; only the comments are swapped.** No further checking needed before citing.

## Cross-check against the published parameter ranges

Table 3 of Fatichi et al. (2012a) tabulates typical ranges for every T&C parameter. Those
are now in the `typical_range_fatichi2012` column (68 of 180 rows — the table doesn't
cover plant hydraulics, which post-dates the 2012 papers). Every other value falls inside
its published range. Four do not:

| Parameter | Value | Published range | Verdict |
|---|---|---|---|
| `Nl_H` | 62 gC/gN | 15–42 | **~1.5× over.** Conifer needles are genuinely high-C:N, so the range may not span them — but this drives maintenance respiration (`r_H` is per unit N). Check before reuse. |
| `dc_C_H` | 0.214 d⁻¹°C⁻¹ | 0.0027–0.067 | **~3× over.** Moot at US_xRM since `Tcold_H = −50` never fires, but it *would* bite in a deciduous set where cold shedding is active. |
| `Tcold_H` | −50 °C | −12 to +10 | **Deliberate.** An off-switch, not a physical threshold. A deciduous set must bring it back in range. |
| `LAI_min_L`, `rjv_L` | 0.1, 2.6 | 0.001–0.05, 1.5–2.5 | Low layer, inactive. Marginal. |

This is a stronger form of validation than the citations alone: it confirms the US_xRM set
is internally consistent with the model's own documented parameter space, and isolates
exactly which entries are unusual.

**`Veg_Optical_Parameter(2)` is NET *Boreal*, not temperate.** For CONUS evergreen sites,
class 1 (NET Temperate) would be the natural choice — but rows 1, 2 and 3 of
`OPTICAL_PAR_VEG` are numerically identical, so it makes no difference to results. Noted
so nobody spends time "fixing" it.

## Parameters that need thought per site, not blind reuse

Everything in the PFT table is reusable in principle, but five entries deserve a check:

- **`ZR95_H` = 800 mm** — PFT placeholder for now, as agreed; overwrite if site data
  appear. Hard constraint: must be ≤ the deepest `Zs` layer, or `Root_Fraction_General`
  aborts the run.
- **`Zs`** ends at 1 m. CLAUDE.md §6 flags this as too shallow for deep-rooted forest;
  deepening it per site (often 2–5 m) also relaxes the `ZR95_H` constraint above.
- **`LDay_min_H` (14.1 h) and `LDay_cr_H` (11.8 h)** — day-length thresholds are latitude
  sensitive. The station set spans 28°N to 47°N, so these should be reviewed rather than
  copied verbatim to the far ends of that range.
- **`age_cr_H` = 1220 d** (~3.3 yr needle retention) is a defining evergreen trait and
  differs substantially between *Pinus*, *Picea*, *Tsuga* and *Abies*. BADM species
  records (the `species` column of `badm_coverage.csv`) can inform this.
- **`PsiX50_H` = −5.0 MPa** is a drought-resistant conifer value, plausibly too negative
  for mesic eastern conifers such as the *Tsuga canadensis* at US-Ha2.

**`Vmax_H` = 32 is set directly**, with the `Maximum_Rubisco_Capacity(Sl, PLNR, Nl)` call
commented out. That is exactly why dynamic LMA does not propagate to photosynthetic
capacity, and why the experiment isolates the leaf-area pathway (CLAUDE.md §1–§2). Leave
it hardcoded — restoring that call would change what the experiment measures.

## Low vegetation layer

The Low layer is off at US_xRM (`Ccrown` high-only, `Vmax_L = 0`, `ZR95_L = 0`, and
`Veg_Optical_Parameter(0)` returning NaN optics). Its parameters are still in the table,
flagged `layer = L` and `active_at_US_xRM = no`, so a complete `MOD_PARAM` can be
assembled without going back to the source file. Filter on
`active_at_US_xRM == yes` for the 114 rows that actually drive an evergreen forest run.

## Deciduous sites

This table is evergreen-specific, and **the repo contains no deciduous-forest
parameterisation to copy from.** Only two `MOD_PARAM` files exist:

| File | Active layer | Type |
|---|---|---|
| `Thanos_US_xRM/MOD_PARAM_US_xRM.m` | High | evergreen needleleaf conifer — the source of this table |
| `TeC_Source_Code-master/MOD_PARAM_ZURICH_SMA.m` | **Low** | managed C3 grassland |

✅ **Correction to an earlier reading of this file.** I first described Zurich's High layer
as unusable placeholders. That was wrong. Its *structural and state* variables are indeed
zeroed because the layer is switched off at that site (`Vmax_H = 0`, `ZR95_H = 0`,
`LAI_H = 0`, `B_H` zeros, `hc_H = 0`, optics class 0 → NaN). But its **physiological,
hydraulic and phenological block is fully and coherently configured for a deciduous
broadleaf tree**, and differs from US_xRM in exactly the ways a deciduous PFT should:

```
aSE_H = 1        age_cr_H = 150 d      Tcold_H = 7 degC     Tlo_H = 12.9 degC
d_leaf_H = 3.5   Klf_H = 1/15          LDay_cr_H = 12.30 h  Nl_H = 30
Do_H = 1000      a1_H = 7              PsiX50_H = -3.5      Wm_H = 1/16425
```

That is a usable deciduous starting point — see
`tc_site_parameters_HBK_Ha2.csv`. What it does **not** supply is `Vmax_H`, the optical
class, rooting depth, and the structural/state variables; those must be set explicitly.

### What the Fatichi et al. (2012) papers contain

Neither paper's *plot-scale* sites are deciduous forest — Lucky Hills (AZ) is desert
shrub and Reynolds Creek Mountain East (ID) is sagebrush, both simulated as a **single
Low-vegetation layer**. But the watershed-scale paper (Part 2) parameterises three PFTs
across the RME and Tollgate catchments, and one of them **is a deciduous broadleaf tree**:

| PFT | Type | Cover (RME / TOL) | Reported traits |
|---|---|---|---|
| Aspen (*Populus tremuloides*) | **deciduous broadleaf tree** | 26.1% / 8.7% | h = 9.5 m, LAI = 1.35, GPP ≈ 250–450 gC m⁻² yr⁻¹ |
| Douglas fir (*Pseudotsuga menziesii*) | evergreen needleleaf | 4.6% / 11.9% | LAI ≈ 2.0, GPP ≈ 350–650 |
| Low sagebrush (*Artemisia arbuscula*) | shrub | 69.3% / 79.4% | h ≈ 0.6 m, LAI = 0.77 |

Whitethorn acacia (*Acacia constricta*) at Lucky Hills is also deciduous, but a **shrub**
in the Low layer (`Ccrown` = 0.25), not a forest canopy.

⚠️ **The aspen parameter values are not printed in either paper.** Part 2 §3.2.1 says only
that traits were "inferred from literature… [White et al., 2000; Kattge and Knorr, 2007]",
and the numeric tables live in the auxiliary material (Text S1, Table 1) — which covers
the *plot-scale* runs, i.e. the shrub and sagebrush sites. So the papers confirm a
deciduous tree PFT exists in T&C and give plausibility targets (LAI, height, GPP), but not
a set to copy.

A deciduous set therefore has to be obtained or built. At minimum it must change:
`aSE_H` (0 → 1), the optical class (`Veg_Optical_Parameter(7)` = BDT temperate),
`age_cr_H` (needle longevity → one growing season), `Tcold_H`, `Tlo_H`, `Tls_H`,
`Bfac_ls_H`, `Klf_H`, `LDay_cr_H`, `d_leaf_H` (0.25 cm needle → broadleaf),
`Nl_H` (deciduous leaves are far less C-rich per N), `Vmax_H`, `PsiX50_H` and the other
hydraulic traits. That is most of the vegetation block — it is a new parameterisation,
not a tweak.

## References cited in the tables

- Fatichi, S., Ivanov, V.Y., Caporali, E. (2012a). A mechanistic ecohydrological model to
  investigate complex interactions in cold and warm water-controlled environments: 1.
  Theoretical framework and plot-scale analysis. *J. Adv. Model. Earth Syst.* 4, M05002.
- Arora, V.K., Boer, G.J. (2005). A parameterization of leaf phenology for the terrestrial
  ecosystem component of climate models. *Earth Interactions* 9, 1–17. — `CASE_ROOT = 1`
- Oleson, K.W. et al. (2010). Technical description of version 4.0 of the Community Land
  Model. NCAR/TN-478+STR; with Dorman & Sellers (1989), Asner et al. (1998), Ross (1975).
  — vegetation optics
- Leuning, R. (1995). A critical appraisal of a combined stomatal-photosynthesis model for
  C3 plants. *Plant Cell Environ.* 18, 339–355. — `a1`, `Do`, `go`
- Farquhar, G.D. et al. (1980); Collatz, G.J. et al. (1991). — C3 photosynthesis
- Kattge, J., Knorr, W. (2007). *Plant Cell Environ.* 30, 1176–1190. — `DSE`, `Ha`
  *(needs verification)*
- Bonan, G.B. (2003) — `r_H`, from the code comment
- Mahfouf, J.-F., Jacquemin, B. (1989) — `KcI`, from the code comment
- Saxton, K.E., Rawls, W.J. (2006). *SSSAJ* 70, 1569–1578. — `SPAR = 2`
- Schenk & Jackson (2002); Fan et al. (2017); Choat et al. (2012); Pelletier (2016);
  Shangguan (2017) — fallback sources named in the notes

---

## The T&C source papers (Fatichi, Ivanov & Caporali 2012, JAMES)

Both PDFs sit in the repo root but are **not committed** (13 MB combined; see
`.gitignore`). Both are open access, so a fresh clone can retrieve them from the DOIs
below.

- **Part 1** — *Theoretical framework and plot-scale analysis*, M05002, doi:10.1029/2011MS000086
- **Part 2** — *Spatiotemporal analyses*, M05003, doi:10.1029/2011MS000087

Part 1 **Table 3** is the full parameter list with the values used. It is not
transcribed here: `pdftotext` interleaves its columns and a mis-copied value would be
worse than no value. Look parameters up in the PDF directly.

### What the papers adjusted, and how

Part 2 §3.2 ¶[30] is the only description of parameter fitting anywhere in the pair,
and it is deliberately modest — fewer than a dozen test runs in total:

| Domain | Runs | What was changed |
|---|---|---|
| Lucky Hills (AZ) | 2 | soil **sealing** parameterization only |
| Reynolds Creek Mountain East (ID) | 5 | `Ks`, anisotropy ratio `aR`, **root depth**, `a1` (photosynthesis↔stomatal conductance), `Vmax` |
| Tollgate (ID) | 0 | not calibrated — reused the RME values |

¶[29]: *"no formal, traditional, or advanced calibration of the model can be carried
out… The choice of model parameter values has been made subjectively, based only on
available data or literature information."* ¶[31]: manual adjustment only, expert
judgement, parameters held inside physically realistic ranges, framed as *"a final
adjustment to refine the simulation skill, which is mainly dictated by the model
structure and boundary conditions."*

⚠️ **Root depth is in that list.** We prescribe `ZR95` from a 1° (~110 km) Schenk &
Jackson grid where 15 CHEESEHEAD towers share a single cell — a parameter the model's
own authors found needed site-level adjustment. State it in the methods; it is also
the strongest argument for the Fan et al. (2017) ~1 km upgrade.

### Soil hydraulics are DERIVED, never fetched

Both watersheds take hydraulic properties from the **Saxton & Rawls (2006)**
pedotransfer functions applied to sand/clay fractions — the same `SPAR = 2` path
`Soil_parameters.m` uses. `Ks`, `Osat`, `L`, `Pe`, `O33`, `alpVG`, `nVG` are all
outputs of texture, so the only soil quantities we fetch are `Psan`/`Pcla`/`Porg` and
the column depth.

- Lucky Hills: *"derived from the pedotransfer functions of Saxton and Rawls [2006]
  using a 0.75 fraction of sand and 0.10 fraction of clay"*
- RME: same functions, with *"spatially variable fractions of sand and clay derived
  from the soil map"*

What ¶[30] calls calibrating `Ks` is therefore adjusting the *pedotransfer output*,
not substituting a measured value. Useful for us: SSURGO ships `ksat_r` per horizon,
which `fetch_soil.py` records — an independent check on the derived `Ks`, not a
replacement for it.

### Site configurations worth comparing against ours

| | Lucky Hills | RME | Part 1 plot scale |
|---|---|---|---|
| Column depth | 2 m | 1 m | 1 m |
| Layers | 18, 10 mm → 400 mm | 10, 10 mm → 200 mm | — |
| Bottom BC | free drainage | **`Kbot` = 0.01 mm/h** (near-impermeable bedrock, from in-situ geology) | free drainage |
| Anisotropy `aR` | 1 | **140** (mimics topographic preferential flow) | n/a (flat plot) |
| `Ks` with depth | `Ks(z) = Ks(0)·e^(−0.0011z)` [Scott et al. 2000] | uniform | — |
| Texture | uniform 0.75 sand / 0.10 clay | spatially variable from soil map | — |

Three things follow for this project:

1. **`Kbot` is not always `NaN`.** US_xRM inherits free drainage, but RME used a
   near-impermeable 0.01 mm/h bedrock derived from in-situ geology. For the 18
   stations where SSURGO reports an actual bedrock contact, a near-impermeable base
   is arguably the better boundary condition than free drainage — worth deciding when
   those shallow sites are modelled.
2. **The 1 m column was flagged as a limitation by the authors themselves.** Part 2
   ¶[57] attributes a Tollgate discharge timing error partly to *"the uniform soil
   depth of 1 m assumed for the entire domain,"* noting variable depth *"can produce
   deep storages of water at the beginning of the melting season."* Direct support for
   the per-site depths in `fetch_soil.py`.
3. **`aR` barely matters for us.** Anisotropy drives lateral subsurface flow; our runs
   are single-point plot scale like Part 1, which is flat with no lateral exchange.
   The RME value of 140 is a distributed-domain device, not a transferable constant.

### Spin-up precedent

Part 1 ¶[100]: initial soil moisture and vegetation carbon pools *"are obtained after
spinning up the model with a simulation of the same duration as the analyzed period."*
One pass over the record, not the repeat-until-convergence protocol — worth knowing
when we set ours, since evergreen wood pools need far longer.
