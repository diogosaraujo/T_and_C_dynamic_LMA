#!/usr/bin/env python3
"""Build the model_run tree: per-station T&C run directories, fixed and dynamic LMA.

    model_run/
      Code/                       shared T&C source, copied ONCE
      GRAPH_MOD.m                 shared
      <STATION>/
        era5_land/
          fixed_lma/   GO_<ST>.m  MOD_PARAM_<ST>.m  LMA_<ST>.mat  Meteo_<ST>_*.mat
          dyn_lma/     same
        hist_gcm/  ssp126/  ssp245/  ssp370/  ssp585/     (created empty)

Everything is generated from T&C/Thanos_US_xRM as the template, with the
site-specific blocks substituted. Substitution is VERIFIED: every pattern must
match exactly once or the station is refused. A silently unsubstituted MOD_PARAM
would run happily with US_xRM's soil and canopy and produce plausible, wrong
numbers -- the failure mode this project has hit repeatedly.

WHAT IS SUBSTITUTED, and from where

    Zs, ms, Psan/Pcla/Porg   soil_profiles.csv   per-layer, Soil_parameters looped
    Kbot                     soil_sites.csv      0.01 mm/h where bedrock is reported,
                                                 NaN (free drainage) otherwise
    ZR95_H                   root_depth CSV      capped at the column depth
    zatm                     BASE_MeasurementHeight, max EC height, else hc + 12
    hc_H                     canopy_height CSV
    Sl_H                     LMA series          fixed arm: 1/(mean(LMA)*f_C)
                                                 dyn arm:   same value as the day-1
                                                            state; the yearly series
                                                            arrives via LMA_<ST>.mat
    aSE_H                    site list           0 evergreen, 1 deciduous

DECISIONS WORTH KNOWING

  * Station IDs contain a hyphen (US-HBK) which is not legal in a MATLAB
    identifier or file name. Directories keep the hyphen; MATLAB names use an
    underscore (MOD_PARAM_US_HBK.m), matching the template's own id_location
    convention ('US_xRM' for site US-xRM).
  * Code/ (145 files) and GRAPH_MOD.m are shared at the tree root rather than
    copied into every arm of every station. GO_<ST>.m addpaths up to them.
  * The fixed arm gets a CONSTANT Sl from the site's mean LMA, not the literature
    PFT value, so both arms share a mean and differ only in variability.
  * f_C = 0.5 is applied here: SLA = 1/(LMA*f_C). The shipped LMA_US_xRM.mat uses
    SLA = 1/LMA with no f_C, so its values are half what T&C should receive.
"""

from __future__ import annotations

import argparse
import csv
import filecmp
import json
import math
import os
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PREPROC = Path(__file__).resolve().parent
INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA",
                                 "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
DEFAULT_ROOT = INPUT_ROOT.parent / "model_run"
DEFAULT_TEMPLATE = REPO_ROOT / "T&C" / "Thanos_US_xRM"
DEFAULT_EXCLUDED = Path(__file__).resolve().parent / "excluded_stations.csv"
DEFAULT_SITE_LISTS = [
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv",
]
DEFAULT_HEIGHTS = PREPROC / "BASE_MeasurementHeight_20260715.csv"

F_C = 0.5                    # leaf carbon fraction; SLA = 1/(LMA*f_C)
KBOT_BEDROCK = 0.01          # [mm/h] near-impermeable, after Fatichi et al. at RME
ZATM_FALLBACK_OFFSET = 12.0  # [m] above canopy where no EC height is reported
SCENARIOS = ["hist_gcm", "ssp126", "ssp245", "ssp370", "ssp585"]
ARMS = ["fixed_lma", "dyn_lma"]
MISSING = {"", "NA", "N/A", "NaN", "nan", "-9999", None}


def mat_name(station: str) -> str:
    """US-HBK -> US_HBK. MATLAB identifiers cannot contain a hyphen."""
    return re.sub(r"[^0-9A-Za-z_]", "_", station)


def fnum(v, default=None):
    try:
        f = float(v)
        return f if math.isfinite(f) else default
    except (TypeError, ValueError):
        return default


# --------------------------------------------------------------------- inputs

def read_stations(paths, wanted):
    out = []
    seen = set()
    for p in paths:
        if not Path(p).is_file():
            print(f"  ! site list not found: {p}", file=sys.stderr)
            continue
        for r in csv.DictReader(open(p, newline="", encoding="utf-8-sig")):
            sid = (r.get("StationID") or "").strip()
            if not sid or sid in seen or (wanted and sid not in wanted):
                continue
            seen.add(sid)
            out.append({"station_id": sid,
                        "forest_type": (r.get("ForestType") or "").strip().lower(),
                        "lat": fnum(r.get("Lat")), "lon": fnum(r.get("Lon"))})
    return sorted(out, key=lambda s: s["station_id"])


def read_excluded(path):
    """Stations dropped from the study, from excluded_stations.csv.

    17 of the 118 are excluded -- 8 with no BADM archive, 8 young or disturbed
    stands that are not at equilibrium for spin-up, and US-KS1 with no rooting
    depth. Building forcing or run directories for them wastes work and, worse,
    puts them in run_list.txt where they would be simulated by accident.
    """
    out = {}
    p = Path(path) if path else None
    if not p or not p.is_file():
        return out
    for r in csv.DictReader(open(p, newline="", encoding="utf-8-sig")):
        sid = (r.get("station_id") or "").strip()
        if sid:
            out[sid] = (r.get("reason") or "excluded").strip()
    return out


def read_soil(soil_dir: Path):
    """(sites by station, profile layers by station) from fetch_soil.py output."""
    sites, layers = {}, defaultdict(list)
    sp = soil_dir / "soil_sites.csv"
    pp = soil_dir / "soil_profiles.csv"
    if sp.is_file():
        for r in csv.DictReader(open(sp, newline="", encoding="utf-8-sig")):
            sites[r["station_id"]] = r
    if pp.is_file():
        for r in csv.DictReader(open(pp, newline="", encoding="utf-8-sig")):
            layers[r["station_id"]].append(r)
    for v in layers.values():
        v.sort(key=lambda r: float(r["z_top_mm"]))
    return sites, layers


def read_lookup(path: Path, key_cols, val_cols):
    """station -> first parseable value among val_cols."""
    out = {}
    if not path or not Path(path).is_file():
        return out
    for r in csv.DictReader(open(path, newline="", encoding="utf-8-sig")):
        sid = next((r[k].strip() for k in key_cols if r.get(k)), None)
        if not sid:
            continue
        for c in val_cols:
            v = fnum(r.get(c))
            if v is not None:
                out[sid] = v
                break
    return out


def read_ec_heights(path: Path):
    """Max eddy-covariance measurement height per station, metres.

    The BASE measurement-height table lists one row per variable and version; the
    turbulent-flux variables carry the height that matters for zatm.
    """
    out = {}
    if not Path(path).is_file():
        return out
    fam = ("FC", "H", "LE", "USTAR", "WS", "TA")
    for r in csv.DictReader(open(path, newline="", encoding="utf-8-sig")):
        sid = (r.get("Site_ID") or "").strip()
        var = (r.get("Variable") or "").strip().upper()
        h = fnum(r.get("Height"))
        if not sid or h is None or h <= 0:
            continue
        base = re.sub(r"(_\d+)+$", "", re.sub(r"_PI(_F)?$", "", var))
        if base in fam:
            out[sid] = max(out.get(sid, 0.0), h)
    return out


def read_lma_series(lma_dir: Path, station: str):
    """[(year, LMA g/m2)] for a station, modelled series, from build_lma_input.py."""
    p = lma_dir / station / f"{station}_LMA_modelled.csv"
    if not p.is_file():
        return []
    out = []
    for r in csv.DictReader(open(p, newline="", encoding="utf-8-sig")):
        y, v = fnum(r.get("year")), fnum(r.get("LMA_g_m2"))
        if y is not None and v is not None:
            out.append((int(y), v))
    return sorted(out)


# ------------------------------------------------------------------ generation

CLAMPED = []          # (Psan, Pcla, Porg, Psil) of every layer that was adjusted


def texture_triple(layer):
    """Sand, clay and organic fractions, as emitted, leaving silt strictly positive.

    Soil_parameters (line 26) computes

        Psil = 1 - Psan - Pcla - Porg;
        if Psil < 0,  disp('SOIL PERCENTAGE INPUTS INCONSISTENT'),  return,  end

    ORGANIC MATTER IS IN THAT BUDGET. A layer of sand 0.95 / clay 0.04 / org 0.02
    has sand+clay well under 1 but Psil = -0.01, so guarding only sand+clay misses
    it -- which is what killed US-SP2 and US-SP4 in array job 35717. The failure
    presents as "Output argument Osat not assigned" rather than an error, because
    the check `return`s after a disp() instead of raising.

    Sand and clay are scaled down into the budget left by Porg; Porg itself is
    kept, since it comes from measured organic carbon and drives the Saxton &
    Rawls density adjustment. Values are rounded to the emitted precision FIRST
    and clamped after, because clamping at full precision and rounding afterwards
    can round the sum straight back over 1.

    The trigger is the EMITTED Psil actually being negative, not a safety margin
    around it. That distinction matters: Soil_parameters only rejects Psil<0, so
    every station that has already run is proof its own Psil >= 0, and a rule
    keyed on Psil<0 therefore cannot alter a profile behind a stored result. An
    earlier version clamped anything within 1e-4 of the limit and silently
    rewrote US-PFt, US-SP1 and US-SP3, which had run perfectly well on a genuine
    trace of silt.
    """
    r = lambda x: float(f"{float(x):.5g}")     # exactly what gets written out
    psan, pcla, porg = r(layer["Psan"]), r(layer["Pcla"]), r(layer["Porg"])

    if 1.0 - psan - pcla - porg >= 0.0:
        return psan, pcla, porg                # valid as-is; emitted unchanged

    # Shrink sand and clay proportionally until the ROUNDED values leave silt.
    # Re-rounding after scaling can drift by ~1e-5, so the result is checked
    # rather than assumed, widening the target if the first attempt misses.
    for target in (2e-5, 1e-4, 1e-3):
        budget = 1.0 - porg - target
        if budget <= 0.0:
            break                              # Porg alone fills the profile
        scale = budget / (psan + pcla)
        s, c = r(psan * scale), r(pcla * scale)
        if 1.0 - s - c - porg >= 0.0:
            CLAMPED.append((psan, pcla, porg, 1.0 - psan - pcla - porg))
            return s, c, porg
    # Unfixable by rescaling: leave it alone so Soil_parameters says so out loud
    # instead of running on a texture invented here.
    return psan, pcla, porg


def soil_block(layers, kbot):
    """MATLAB for the depth-resolved soil, replacing the single-triple template.

    Soil_parameters is scalar-only (it uses ^ and * on scalars, and `if Psil<0`),
    so the layered profile has to be a loop, not a vectorised call.
    """
    zs = [float(layers[0]["z_top_mm"])] + [float(l["z_bot_mm"]) for l in layers]
    ms = len(layers)
    tex = [texture_triple(l) for l in layers]
    psan = ", ".join(f"{s:.5g}" for s, _, _ in tex)
    pcla = ", ".join(f"{c:.5g}" for _, c, _ in tex)
    porg = ", ".join(f"{o:.5g}" for _, _, o in tex)
    zss = " ".join(f"{z:g}" for z in zs)
    return ms, zs, f"""%%%%%%%%%%% SOIL INPUT -- depth-resolved, generated by build_model_run.py
%%%% One triple per Zs layer from SSURGO/POLARIS. Saxton & Rawls runs PER LAYER:
%%%% Soil_parameters is scalar-only (it uses ^ and * on scalars and `if Psil<0`),
%%%% and so is the van Genuchten derivation below it -- alpVG and nVG come from
%%%% L and Pe, so they must be computed inside the loop too, not replicated from
%%%% layer 1.
Psan_Zs = [{psan}];
Pcla_Zs = [{pcla}];
Porg_Zs = [{porg}];
Psan = Psan_Zs(1); Pcla = Pcla_Zs(1); Porg = Porg_Zs(1); %%% surface, for reference
Color_Class = 0;
SPAR = 2; %%% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
%%%%%%%%%%%%%%%%%%%
Osat=zeros(1,ms); L=zeros(1,ms); Pe=zeros(1,ms); Ks=zeros(1,ms); O33=zeros(1,ms);
rsd=zeros(1,ms); lan_dry=zeros(1,ms); lan_s=zeros(1,ms); cv_s=zeros(1,ms);
K_usle=zeros(1,ms); nVG=zeros(1,ms); alpVG=zeros(1,ms);
for jk = 1:ms
    [Osat(jk),L(jk),Pe(jk),Ks(jk),O33(jk),rsd(jk),lan_dry(jk),lan_s(jk),cv_s(jk),K_usle(jk)] = ...
        Soil_parameters(Psan_Zs(jk),Pcla_Zs(jk),Porg_Zs(jk));
    p_jk = 3 + 2/L(jk);
    m_jk = 2/(p_jk-1);
    nVG(jk) = 1/(1-m_jk);
    alpVG(jk) = (((-101.9368*Pe(jk))*(2*p_jk*(p_jk-1))/(p_jk+3))* ...
        ((55.6+7.4*p_jk+p_jk^2)/(147.8+8.1*p_jk+0.092*p_jk^2)))^-1; %%% [1/mm]
end
K_usle = mean(K_usle);
Ks_Zs = Ks; %%% [mm/h] already per layer
%%%%%%%%%%%%%%%
"""


SUBS_REQUIRED = ("kbot", "zatm", "soil", "zs", "zr95", "sl", "ase", "hc")

# --------------------------------------------------------- deciduous PFT block
#
# The template is MOD_PARAM_US_xRM.m, an EVERGREEN subalpine conifer. Flipping
# aSE_H to 1 selects the deciduous code path in PHENOLOGY_STATE/VEGETATION_DYNAMIC
# but leaves every threshold that path reads at its conifer value, which killed
# all 38 deciduous stations: with Tcold_H = -50 C against a climatological minimum
# of -10 C, cold leaf shed can never fire, so leaf age accumulates past age_cr
# until the leaf pool is exhausted (US-HBK: 83.9% of days in senescence, LAI 4.5
# -> 0.01, reserves still full at ~340 gC).
#
# SOURCE: Code/PARAMETERS_ALL_CY_Cb.m, T&C's own 8-PFT table, column 8 -- the one
# column with aSE_H == 1. Preferred over MOD_PARAM_ZURICH_SMA.m (the only other
# deciduous parameterisation shipped) because the table is the model's own
# multi-PFT driver: internally consistent across all 8 PFTs and mutually
# consistent in its units and conventions, where ZURICH is one tuned Swiss site.
# Where they disagree the table is the more conservative choice, and ZURICH is
# noted per line below.
#
# NOT taken from the table: Sl_H (the PLSR LMA supplies it) and Wm_H (0 in the
# table for all 8 PFTs; heartwood B(6) is a pure accumulator here since
# OPT_SoilBiogeochemistry = 0, so it changes no flux).
DECIDUOUS_PFT = [
    # (name, LHS as written in the template, value, comment)
    ("Tcold_H",      "Tcold_H",      "5",       "[C] cold leaf shed (was -50: never fired). ZURICH 7"),
    ("age_cr_H",     "age_cr_H",     "110",     "[day] critical leaf age (was 1220). ZURICH 150"),
    # Tlo comes from the table's COLD-ADAPTED column 2, not the deciduous column 8.
    # It is compared against Tsmm, a 30-day running mean of the PROFILE-MEAN SOIL
    # temperature (VEGGIE_UNIT: Tsm = mean(Tdp)), which lags air by ~1 month. At
    # Hubbard Brook column 8's 13 C is not reached until DOY 178, putting leaf-out
    # in late June: the tower's own CO2 flux gives sustained uptake from DOY 149
    # (29 May, sd 6 d over 7 yr) and the model was +39 days late. Column 2 is an
    # evergreen column, but Tlo is a climate threshold rather than a leaf-form
    # trait, and column 2 is the only cold-adapted set in the table (Tlo 4,
    # LDay_min 9.5, low dc_C). Two independent checks corroborate ~4: calibrating
    # against the tower gives 4.8, and a 10 C air-temperature criterion gives 4.1.
    # One value for every deciduous station -- deliberately NOT calibrated per site.
    ("Tlo_H",        "Tlo_H",        "4.0",     "[C] mean T leaf onset; PFT 2 cold-adapted (col 8's 13 gives DOY 178 vs tower 149)"),
    ("dmg_H",        "dmg_H",        "30",      "[day] day of max growth"),
    ("Trr_H",        "Trr_H",        "3.0",     "[gC/m2 d] translocation rate (was 0.5). ZURICH 3.5"),
    ("LDay_min_H",   "LDay_min_H",   "11.5",    "[h] min day length for leaf onset (was 14.1). ZURICH 11.0"),
    ("LDay_cr_H",    "LDay_cr_H",    "12.0",    "[h] day length for senescence. ZURICH 12.30"),
    ("Klf_H",        "Klf_H",        "0.025",   "[1/d] dead leaf fall turnover. ZURICH 1/15"),
    ("eps_ac_H",     "eps_ac_H",     "0.3",     "[-] allocation to reserve. ZURICH 1"),
    # ZURICH, not the PFT table: see the go_H note below. The table's 0 is outside
    # Table 3's 1/600-1/30 and disables drought leaf mortality entirely.
    ("dd_max_H",     "dd_max_H",     "1/365",   "[1/d] max drought leaf mortality; ZURICH (table's 0 is outside Table 3)"),
    ("drn_H",        "drn_H",        "0.0020",  "[1/d] root turnover"),
    ("dsn_H",        "dsn_H",        "0.0027",  "[1/d] sapwood transfer"),
    ("LtR_H",        "LtR_H",        "1.5",     "[-] max leaf-to-root ratio"),
    ("Mf_H",         "Mf_H",         "0.0125",  "[1/d] fruit maturation turnover"),
    # Photosynthesis. Vmax stays PRESCRIBED -- Maximum_Rubisco_Capacity is left
    # commented out so LMA propagates only through leaf area (CLAUDE.md 1) -- but
    # at the deciduous value, not the conifer 32. ZURICH sets Vmax_H = 0, which
    # would re-enable the N-based route and change what the experiment isolates.
    ("Vmax_H",       "Vmax_H",       "50",      "[umol/m2 s] deciduous (evergreen keeps 32)"),
    ("Nl_H",         "Nl_H",         "35",      "[gC/gN] leaf C:N (was 62, a conifer value). ZURICH 30"),
    ("a1_H",         "a1_H",         "6",       "[-] stomatal slope (was 5). ZURICH 7"),
    ("Do_H",         "Do_H",         "600",     "[Pa] VPD sensitivity (was 700). ZURICH 1000"),
    ("rjv_H",        "rjv_H",        "2.1",     "[-] Jmax:Vmax scaling (was 1.8). ZURICH 2.8"),
    # Water stress. PsiL50 drives Bfac in BetaFactor.m (phenology triggers, drought
    # leaf mortality, growth limitation) -- it is NOT the xylem vulnerability, and
    # PsiX50 is unused here because OPT_PlantHydr = 0.
    ("Psi_sto_00_H", "Psi_sto_00_H", "-0.7",    "[MPa] stomatal, 2% loss (was -0.8)"),
    ("Psi_sto_50_H", "Psi_sto_50_H", "-2.5",    "[MPa] stomatal, 50% loss"),
    ("PsiL50_H",     "PsiL50_H",     "-6.5",    "[MPa] leaf, 50% loss -> Bfac (was -3.2). ZURICH -5.6"),
    ("KnitH",        "KnitH",        "0.4",     "[-] canopy nitrogen decay (was 0.35)"),
    # A deciduous canopy starts leafless, not with 1172-day-old needles that are
    # already at 96% of age_cr.
    ("AgeL_H",       r"AgeL_H\(1,:\)", "0",     "[day] initial leaf age (was 1172)"),
    # --- the rest of PFT 8, added after auditing the table column by column ---
    # A broadleaf is not a needle: d_leaf sets the leaf boundary-layer
    # conductance, so leaving it at the conifer 0.25 would give a 3 cm broadleaf
    # canopy the aerodynamic behaviour of spruce needles.
    ("d_leaf_H",     "d_leaf_H",     "7",       "[cm] leaf dimension (was 0.25, a needle)"),
    ("Bfac_lo_H",    "Bfac_lo_H",    "0.9",     "[-] leaf-onset water-stress threshold (was 0.99)"),
    ("mjDay_H",      "mjDay_H",      "230",     "[day] last day leaf onset may trigger (was 220)"),
    ("r_H",          "r_H",          "0.02",    "[-] respiration coefficient (was 0.055)"),
    ("gR_H",         "gR_H",         "0.22",    "[-] growth respiration (was 0.25)"),
    # THREE PARAMETERS COME FROM ZURICH, NOT THE PFT TABLE: go_H, dd_max_H, dc_C_H.
    #
    # Fatichi et al. (2012) Part 1 Table 3 publishes "expected realistic ranges"
    # spanning all vegetation types. On every parameter where the PFT table and
    # MOD_PARAM_ZURICH_SMA disagree AND Table 3 gives a bound, ZURICH is inside it
    # and the table is outside:
    #     go_H      ZURICH 0.01     table 0.001   range 0.005-0.04
    #     dd_max_H  ZURICH 1/365    table 0       range 1/600-1/30
    #     dc_C_H    ZURICH 2/365    table 78/365  range 1/365-1/15
    # dc_C_H is the clearest: both are written as N/365 and they differ by a factor
    # of 39, with only ZURICH's landing in range. So for these three the deciduous
    # site file is better vetted than the generic PFT column, and it is the right
    # source anyway -- ZURICH is a deciduous parameterisation, US_xRM is not.
    ("go_H",         "go_H",         "0.01",    "[mol/s m2] cuticular conductance; ZURICH (table's 0.001 is below Table 3)"),
    ("Ha_H",         "Ha_H",         "72",      "[kJ/mol] entropy factor (was 89)"),
    # PsiG50/PsiG99 set the growth-limitation curve in BetaFactor (Bfac_all);
    # PsiX50 belongs to the plant-hydraulics module, which is off here
    # (OPT_PlantHydr = 0), and is carried only to keep PFT 8 internally whole.
    ("PsiG50_H",     "PsiG50_H",     "-0.7",    "[MPa] growth limitation, 50% (was -0.8)"),
    ("PsiG99_H",     "PsiG99_H",     "-6.5",    "[MPa] growth limitation, 99% (was -2.5)"),
    ("PsiX50_H",     "PsiX50_H",     "-9",      "[MPa] xylem 50% (unused: OPT_PlantHydr=0)"),
    # PsiL00 pairs with PsiL50 in BetaFactor. The table writes it as
    # `PsiL00_H = -1-[-1.4 ...]`, so PFT 8 evaluates to +0.4 MPa. A POSITIVE water
    # potential at the 2%-loss anchor is unusual; the curve is still well formed
    # (p = log(49)/6.9 = 0.56, PLC = 0.025 at Psi_l = 0). Taken from the table on
    # instruction, flagged here because it is the one adopted value I cannot
    # independently justify.
    ("PsiL00_H",     "PsiL00_H",     "0.4",     "[MPa] leaf 2% loss = -1-(-1.4); POSITIVE, see note"),
    # Now that Tcold = 5 makes cold shedding live, dc_C sets how fast leaves drop.
    # 78/365 (PFT table AND US_xRM) is 3.2x above Table 3's 1/365-1/15; ZURICH's
    # 2/365 is inside. Same N/365 form, factor of 39 apart.
    ("dc_C_H",       "dc_C_H",       "2/365",   "[1/(d C)] cold foliage loss; ZURICH (78/365 is 3.2x above Table 3)"),
    # --- INITIAL CONDITIONS: not in the table, so chosen here -----------------
    # The table carries no state, and neither shipped IC set fits: US_xRM starts a
    # mature EVERGREEN canopy in full leaf on 1 January (LAI 4.03, 215 gC of
    # leaves, PHE_S=4 senescence), while ZURICH starts from BARE GROUND (all pools
    # zero), which needs a century of spin-up to become a forest. A deciduous
    # stand on 1 Jan 1985 is mature but LEAFLESS, so:
    #   - B(1) leaf = 0. That carbon is NOT carried into the reserve: at autumn
    #     senescence leaf carbon goes to litter, only nutrients are resorbed.
    #   - B(4) reserve = 350, sized from the model's own two constraints rather
    #     than from the evergreen number. VEGETATION_DYNAMIC holds the reserve at
    #     or above 0.67*B(2) = 219 [Friend et al. 1997, line 144], and the spring
    #     flush draws Tr = min(B(4),Trr) = 3 gC/m2/d only while PHE_S == 2, i.e.
    #     ~90 gC over dmg = 30 d. The reserve must therefore exceed 219 + 90 = 309
    #     to fund a flush without breaching the floor and triggering the model's
    #     reserve-refill diversion. US_xRM's 244 would breach it (244-90 = 154).
    #   - B(5) fruit, B(7) standing dead = 0: no reproductive allocation in
    #     January, and a deciduous canopy's leaves are on the ground by then.
    #   - B(6) heartwood = 0, as in Dr. Paschalis's US_xRM vector. That leaves
    #     TBio = 0.02*sum(B) = ~19 t DM/ha, so Allocation_Coefficients sees a
    #     young stand rather than a mature one -- a known bias, kept deliberately
    #     so both PFTs start in the same structural state. It applies equally to
    #     both arms of every station, so it largely cancels in the fixed-vs-dynamic
    #     difference, which is the reported quantity. Spin-up is the real fix.
    #   - sapwood, fine root, fruit, heartwood, standing dead kept from US_xRM: no
    #     deciduous source for them exists anywhere in the repo, and the one model
    #     constraint that touches them -- B(1) < LtR*B(3) = 1.5*262 = 393 gC --
    #     caps leaf carbon at LAI ~9.8, far above a hardwood canopy, so it never
    #     binds.
    #   - PHE_S = 1 (dormant, case 1 -> 2 on leaf onset), not 4 (senescence).
    # All of this is provisional until spin-up replaces it (CLAUDE.md 5).
    ("LAI_H",        r"LAI_H\(1,:\)", "0.0",    "leafless on 1 Jan (was 4.03, evergreen)"),
    ("B_H",          r"B_H\(1,:,:\)", "0 327 262 350 0 0 0 0",
     "leafless Jan; reserve 350 > 0.67*B(2)+Trr*dmg = 309; B(6)=0 as US_xRM"),
    ("PHE_S_H",      r"PHE_S_H\(1,:\)", "1",    "dormant (was 4 = senescence)"),
]
# ONE table entry deliberately NOT adopted:
#
# Cx_H -- the table has -10 in all 8 columns, but Plant_PV_Curve.m documents the
#   reference value in its own header as `%Cx = 150; %%% [kg/m^3 sapwood MPa]`,
#   which is exactly what US_xRM uses. A NEGATIVE capacitance would flip the sign
#   of a = -1000/(Cx*(1+(1-n_sap)/fwat)), so -10 reads as a placeholder rather
#   than a trait. It is unused here anyway (OPT_PlantHydr = 0), but where the code
#   and the table disagree the code wins.


def render_mod_param(template: str, st: dict) -> tuple[str, list[str]]:
    """Substitute the site blocks. Every pattern must fire exactly once."""
    txt = template
    fired = {}

    def sub(key, pattern, repl, count=1):
        nonlocal txt
        txt, n = re.subn(pattern, repl.replace("\\", "\\\\"), txt, count=count,
                         flags=re.M)
        fired[key] = fired.get(key, 0) + n

    ms, zs, block = st["_soil_block"]

    # Replace the ENTIRE soil block in one go: from the SOIL INPUT banner down to
    # and including Ks_Zs = Ks*ones(1,ms).
    #
    # Stopping at the Soil_parameters call, as this used to, leaves the van
    # Genuchten derivation behind -- SPAR, p = 3+2/L, nVG, alpVG. Those are scalar
    # expressions, and once the loop above makes L an array, p = 3+2/L becomes a
    # matrix right division and MATLAB refuses it (job 35704, "Matrix dimensions
    # must agree"). The generated block computes nVG and alpVG per layer itself,
    # so the original must go entirely, not partially.
    sub("soil",
        r"^%%%%%%%%%% SOIL INPUT.*\n(?:.*\n)*?^Ks_Zs\s*=\s*Ks\*ones\(1,ms\);.*\n",
        block)
    zss = " ".join(f"{z:g}" for z in zs)
    sub("zs", r"^Zs=\s*\[[^\]]*\];.*$",
        f"Zs = [{zss}]; %% ms+1 = {len(zs)}, from the SSURGO/POLARIS profile")
    # Literal, not a reference to something the soil block defines: Kbot is
    # assigned at template line 23 and the soil block lands at line 33, so an
    # indirection here is read before it exists (job 35696).
    kb = st["_kbot"]
    sub("kbot", r"^Kbot\s*=\s*[^;]*;.*$",
        f"Kbot = {'NaN' if kb is None else f'{kb:g}'}; %% [mm/h] "
        f"{'free drainage, no bedrock reported' if kb is None else 'near-impermeable bedrock (Fatichi et al., RME)'}")
    sub("zatm", r"^zatm\s*=\s*[^;]*;.*$",
        f"zatm = {st['zatm']:g}; %% [m] {st['zatm_src']}")
    sub("zr95", r"^ZR95_H\s*=\s*\[[^\]]*\];.*$",
        f"ZR95_H = [{st['zr95']:g}]; %% [mm] {st['zr95_src']}")
    sub("sl", r"^Sl_H\s*=\s*\[[^\]]*\];.*$",
        f"Sl_H = [{st['sl']:.6g}]; %% [m^2/gC] 1/(LMA*{F_C}), LMA={st['lma_mean']:.1f} g/m2")
    sub("ase", r"^aSE_H\s*=\s*\[[^\]]*\];.*$",
        f"aSE_H = [{st['ase']}]; %% {st['forest_type']}")
    sub("hc", r"^hc_H\(1,:\)\s*=\s*\[[^\]]*\];.*$",
        f"hc_H(1,:) = [{st['hc']:g}]; %% [m] {st['hc_src']}")

    # Evergreen stations keep the template block: the template IS an evergreen
    # conifer, so it is already self-consistent. Deciduous stations need the whole
    # PFT swapped, not just the aSE_H switch.
    if st["ase"] == 1:
        # Leaf optics. The template hardcodes Veg_Optical_Parameter(2) = NET Boreal
        # for BOTH PFTs, so a broadleaf deciduous canopy runs with needleleaf
        # optics: NIR reflectance 0.35 instead of 0.45 and a leaf-angle parameter
        # of 0.01 (near-vertical, needle-like) instead of 0.25 (horizontal). That
        # biases absorbed radiation and the light-interception profile.
        # Class 7 = BDT temperate (the file's own header lists the classes).
        # Evergreen is left alone deliberately: rows 1 (NET Temperate), 2 (NET
        # Boreal) and 3 are numerically IDENTICAL in OPTICAL_PAR_VEG, so the
        # template's class 2 is already correct for the needleleaf sites.
        sub("pft_optics", r"^\s*\[PFT_opt_H\(1\)\]\s*=\s*Veg_Optical_Parameter\([^)]*\);.*$",
            "[PFT_opt_H(1)]=Veg_Optical_Parameter(7); %% BDT temperate "
            "(was 2 = NET Boreal, a needleleaf class)")

        # Initial conditions carry no PFT-table provenance -- the table has no
        # state at all -- so they must not be stamped as if they did.
        IC_KEYS = {"LAI_H", "B_H", "PHE_S_H", "AgeL_H"}
        # These three come from the deciduous SITE file, not the PFT table: the
        # table puts all three outside Table 3's published ranges (see the note on
        # go_H above). Stamped accordingly so the provenance in the generated
        # MOD_PARAM is true rather than uniform.
        ZURICH_KEYS = {"go_H", "dd_max_H", "dc_C_H"}
        for key, lhs, val, note in DECIDUOUS_PFT:
            src = ("initial condition, chosen here: no state in the PFT table"
                   if key in IC_KEYS else
                   "MOD_PARAM_ZURICH_SMA.m deciduous" if key in ZURICH_KEYS else
                   "PARAMETERS_ALL_CY_Cb.m PFT 8")
            sub(f"pft_{key}", rf"^\s*{lhs}\s*=\s*[^;]*;.*$",
                f"{lhs.replace(chr(92), '')} = [{val}]; %% {note} [{src}]")

    required = list(SUBS_REQUIRED)
    if st["ase"] == 1:
        # Every deciduous parameter must land. A silently-missed one leaves a
        # conifer value in a deciduous run, which is exactly the failure that
        # produced 38 dying forests and passed every check we had.
        required += [f"pft_{k}" for k, _, _, _ in DECIDUOUS_PFT] + ["pft_optics"]
    problems = [f"{k}: matched {fired.get(k, 0)}x (expected 1)"
                for k in required if fired.get(k, 0) != 1]
    # Tdp(1,:) = Ta(1)*ones(1,ms) is the initial soil temperature and MUST survive.
    leftover = [ln for ln in txt.splitlines()
                if "*ones(1,ms)" in ln and not ln.strip().startswith("%")
                and "Tdp(1,:)" not in ln]
    if leftover:
        problems.append(f"{len(leftover)} un-layered *ones(1,ms) line(s) remain: "
                        + "; ".join(x.strip()[:40] for x in leftover[:3]))

    # Scalar expressions from the template that CANNOT survive once the soil loop
    # makes L, Pe and the rest per-layer arrays. p = 3+2/L is a matrix right
    # division on an array and MATLAB refuses it -- job 35704 -- but it parses
    # cleanly, so checkcode says nothing and the failure only appears at runtime.
    scalar_only = ("p=3+2/L", "m=2/(p-1)", "alpVG=(((", "nVG=L+1")
    survived = [s for s in scalar_only if s in txt]
    if survived:
        problems.append("scalar van Genuchten line(s) survived the soil block, "
                        "which the per-layer arrays break: " + ", ".join(survived))
    return txt, problems


GO_TEMPLATE = """clear
%% Generated by build_model_run.py -- do not edit by hand.
%% Station {station} ({forest_type}), {arm}
addpath('{code_rel}')
addpath('{root_rel}')   %%% shared GRAPH_MOD.m

dt = 3600; %%[s]
dth = 1;   %%[h]

ms = {ms};  %%% Soil layers (site-specific, from the SSURGO/POLARIS profile)
cc = 1;     %%% Crown area

SLA_ex = load('LMA_{mname}.mat');

id_location = '{mname}';
load('{meteo_name}')

x1 = 1;
x2 = length(Ta);
NN = x2-x1+1;

Date=Date(x1:x2); Pr=Pr(x1:x2); Ta=Ta(x1:x2);
Ws=Ws(x1:x2); ea=ea(x1:x2); SAD1=SAD1(x1:x2);
SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2); Pre=Pre(x1:x2);
SAB2=SAB2(x1:x2); N=N(x1:x2); Tdew=Tdew(x1:x2); esat=esat(x1:x2);
PARB=PARB(x1:x2); PARD=PARD(x1:x2);

Ds = esat-ea; Ds(Ds<0) = 0;   %% [Pa] vapour pressure deficit
Oa = 210000;                  %% [umolO2/mol]
Ws(Ws<=0) = 0.01;

[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1)=YE; Datam(:,2)=MO; Datam(:,3)=DA; Datam(:,4)=HO;
clear YE MO DA HO MI SE

PARAM_IC = strcat('MOD_PARAM_','{mname}');

{main_frame};

rmpath('{code_rel}')

save(['RES_', id_location], '-v7.3');

%%% NO FIGURES HERE -- deliberately. Plotting used to run at the end of this
%%% script, and GRAPH_MOD draws 315,576-point lines: jobs 36261/36262 finished the
%%% simulation in ~20 min, then spent 7.5 h rendering and were killed by the 8 h
%%% wall clock with the science already complete but the job marked failed.
%%% Figures are now a separate, cheap, independently retryable step:
%%%     sbatch slurm/submit_figures.sh <STATION> <fixed_lma|dyn_lma>
%%% which draws them from RES_*.mat. Keeping them out of here means a plotting
%%% problem can never cost a simulation, and the run time is the science alone.
fprintf('done. figures: sbatch slurm/submit_figures.sh %s <arm>\\n', id_location);
"""


def build_station(st, tmpl_txt, args, out_root):
    """Write both arms for one station. Returns (n_written, problems)."""
    sid, mname = st["station_id"], mat_name(st["station_id"])
    probs = []
    for arm in ARMS:
        d = out_root / sid / "era5_land" / arm
        d.mkdir(parents=True, exist_ok=True)
        # Depth from the arm directory to the shared Code/ at the tree root.
        code_rel = "../../../Code"

        txt, p = render_mod_param(tmpl_txt, st)
        probs += [f"{arm}: {x}" for x in p]
        (d / f"MOD_PARAM_{mname}.m").write_text(txt, encoding="utf-8")

        (d / f"GO_{mname}.m").write_text(GO_TEMPLATE.format(
            station=sid, forest_type=st["forest_type"], arm=arm, code_rel=code_rel,
            root_rel="../../..",
            ms=st["_soil_block"][0], mname=mname, meteo_name=st["meteo_name"],
            main_frame="MAIN_FRAME_SLA" if arm == "dyn_lma" else "MAIN_FRAME",
        ), encoding="utf-8")

        if st.get("meteo_src"):
            dst = d / st["meteo_name"]
            if not dst.exists():
                try:
                    os.symlink(st["meteo_src"], dst)
                except OSError:
                    shutil.copy2(st["meteo_src"], dst)
        write_lma_mat(d / f"LMA_{mname}.mat", st["lma"], st["sl"], arm)

    for s in SCENARIOS:
        (out_root / sid / s).mkdir(parents=True, exist_ok=True)
    return probs


def write_lma_mat(path: Path, series, sl_fixed, arm):
    """LMA_<ST>.mat with years, LMA, SLA_H, SLA_L -- the LMA_US_xRM.mat layout.

    f_C is applied here: SLA = 1/(LMA*f_C). The shipped LMA_US_xRM.mat stores
    SLA = 1/LMA with no carbon fraction, so its values are half what T&C needs.
    The fixed arm gets the site's mean SLA repeated, so both arms share a mean and
    differ only in year-to-year variability.
    """
    try:
        import numpy as np
        from scipy.io import savemat
    except ImportError:
        return
    if not series:
        return
    years = np.array([y for y, _ in series], dtype=np.uint16)
    lma = np.array([v for _, v in series], dtype=float)
    sla = 1.0 / (lma * F_C) if arm == "dyn_lma" else np.full(lma.size, sl_fixed)
    savemat(path, {"years": years, "LMA": lma,
                   "SLA_H": sla.reshape(-1, 1),
                   "SLA_L": np.zeros((lma.size, 1))})


# ------------------------------------------------------------------------ main

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    p.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    p.add_argument("--stations", default="US-HBK,US-Ha2",
                   help="comma-separated; 'all' for every station in the site lists")
    p.add_argument("--site-list", type=Path, action="append", default=None)
    p.add_argument("--exclude-file", type=Path, default=DEFAULT_EXCLUDED,
                   help="CSV of stations to drop (default: excluded_stations.csv)")
    p.add_argument("--soil", type=Path, default=INPUT_ROOT / "soil")
    p.add_argument("--lma", type=Path, default=INPUT_ROOT / "lma")
    p.add_argument("--meteo", type=Path, default=INPUT_ROOT / "meteo")
    p.add_argument("--canopy", type=Path,
                   default=INPUT_ROOT / "canopy_height" / "canopy_height_gedi.csv")
    p.add_argument("--root-depth", type=Path,
                   default=INPUT_ROOT / "root_depth" / "root_depth_schenk_jackson.csv")
    p.add_argument("--heights", type=Path, default=DEFAULT_HEIGHTS)
    p.add_argument("--meteo-years", default="1985_2020")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()

    tmpl = a.template / "MOD_PARAM_US_xRM.m"
    if not tmpl.is_file():
        print(f"ERROR: template not found: {tmpl}", file=sys.stderr)
        return 1
    tmpl_txt = tmpl.read_text(encoding="utf-8", errors="replace")

    wanted = None if a.stations.strip().lower() == "all" else {
        s.strip() for s in a.stations.split(",") if s.strip()}
    excluded = read_excluded(a.exclude_file)
    stations = [x for x in read_stations(a.site_list or DEFAULT_SITE_LISTS, wanted)
                if x["station_id"] not in excluded]
    if not stations:
        print("no stations selected", file=sys.stderr)
        return 1

    soil_sites, soil_layers = read_soil(a.soil)
    # AmeriFlux first, gridded fallback -- the rule agreed for every site
    # parameter. fetch_canopy_height.py writes the BADM value and the GEDI window
    # statistics side by side; the median is preferred over the mean because a
    # window containing a gap or a road drags the mean down.
    hc, hc_src = {}, {}
    if a.canopy.is_file():
        for r in csv.DictReader(open(a.canopy, newline="", encoding="utf-8-sig")):
            sid = (r.get("station_id") or "").strip()
            for col, lab in (("badm_heightc_m", "AmeriFlux BADM HEIGHTC"),
                             ("hc_gedi_median_m", "GEDI window median"),
                             ("hc_gedi_mean_m", "GEDI window mean")):
                v = fnum(r.get(col))
                if sid and v is not None and 0 < v < 120:
                    hc[sid], hc_src[sid] = v, lab
                    break
    zr = read_lookup(a.root_depth, ("station_id", "StationID"),
                     ("ZR95_H_mm", "ZR95_mm", "ZR95"))
    ec = read_ec_heights(a.heights)

    print(f"root      : {a.root}")
    print(f"template  : {a.template}")
    print(f"stations  : {len(stations)} ({len(excluded)} excluded)")
    print(f"soil      : {len(soil_sites)} sites, {len(soil_layers)} profiles")
    from collections import Counter as _C
    print(f"canopy    : {len(hc)} "
          f"({', '.join(f'{k} {v}' for k, v in _C(hc_src.values()).most_common())})"
          if hc else "canopy    : 0   <-- NONE FOUND")
    print(f"ZR95      : {len(zr)}   EC heights: {len(ec)}")
    print(f"meteo dir : {a.meteo}{'' if a.meteo.is_dir() else '   <-- NOT FOUND'}\n")

    ready, blocked, all_probs = [], [], {}
    clamped_stations = []
    for st in stations:
        sid = st["station_id"]
        miss = []

        layers = soil_layers.get(sid)
        srow = soil_sites.get(sid, {})
        if not layers:
            miss.append("soil profile")
        st["hc"] = hc.get(sid)
        st["hc_src"] = hc_src.get(sid, "")
        if st["hc"] is None:
            miss.append("canopy height")
        st["ase"] = 0 if st["forest_type"].startswith("ever") else 1

        # zatm: measured EC height, else canopy + offset.
        if sid in ec:
            st["zatm"], st["zatm_src"] = ec[sid], "max EC measurement height (AmeriFlux)"
        elif st["hc"] is not None:
            st["zatm"] = st["hc"] + ZATM_FALLBACK_OFFSET
            st["zatm_src"] = f"hc + {ZATM_FALLBACK_OFFSET:g} m (no EC height reported)"
        else:
            miss.append("zatm")

        st["lma"] = read_lma_series(a.lma, sid)
        if not st["lma"]:
            miss.append("LMA series")
        else:
            st["lma_mean"] = sum(v for _, v in st["lma"]) / len(st["lma"])
            st["sl"] = 1.0 / (st["lma_mean"] * F_C)

        # Column depth, bedrock boundary condition, and the ZR95 cap.
        depth = fnum(srow.get("column_depth_mm"))
        bed = fnum(srow.get("bedrock_cm"))
        kbot = KBOT_BEDROCK if bed is not None else None
        st["_kbot"] = kbot
        z = zr.get(sid)
        if z is None:
            miss.append("ZR95")
        elif depth is not None and z > depth:
            st["zr95"] = depth
            st["zr95_src"] = (f"Schenk & Jackson {z:g} mm CAPPED at the {depth:g} mm "
                              f"column ({srow.get('bedrock') or 'no bedrock'})")
        else:
            st["zr95"] = z
            st["zr95_src"] = "Schenk & Jackson D95"

        st["meteo_name"] = f"Meteo_{mat_name(sid)}_{a.meteo_years}.mat"
        cand = a.meteo / st["meteo_name"]
        st["meteo_src"] = cand if cand.is_file() else None
        if st["meteo_src"] is None:
            miss.append("meteo .mat")

        if layers:
            # Record which stations the silt clamp touched, so a rebuild says out
            # loud whose MOD_PARAM changed. Everyone else's is byte-identical to
            # the file that produced their stored RES, and needs no re-run.
            before = len(CLAMPED)
            st["_soil_block"] = soil_block(layers, kbot)
            if len(CLAMPED) > before:
                clamped_stations.append((sid, len(CLAMPED) - before))

        if miss:
            blocked.append((sid, miss))
            continue
        ready.append(st)

    print(f"{'=' * 66}\nREADY: {len(ready)}/{len(stations)} stations have every input\n{'=' * 66}")
    for sid, miss in blocked:
        print(f"  ! {sid:<8} missing: {', '.join(miss)}")
    if blocked:
        print()

    # Only these stations' MOD_PARAM texture differs from the previous build, so
    # only these need re-running on account of the silt clamp. Silence here means
    # every other MOD_PARAM is byte-identical to the one behind its stored RES.
    if clamped_stations:
        print("silt clamp applied (Psan+Pcla+Porg >= 1, so Soil_parameters' "
              "Psil=1-Psan-Pcla-Porg would be negative):")
        for sid, n in clamped_stations:
            print(f"  ~ {sid:<8} {n} layer(s); sand+clay scaled into the budget "
                  f"left by Porg")
        print()
    else:
        print("silt clamp: no station needed it -- textures unchanged\n")

    if a.dry_run:
        for st in ready:
            print(f"  - {st['station_id']:<8} ms={st['_soil_block'][0]:<3} "
                  f"Kbot={'0.01' if st['_kbot'] else 'NaN':<5} "
                  f"ZR95={st['zr95']:>6g} zatm={st['zatm']:>5g} hc={st['hc']:>5g} "
                  f"Sl={st['sl']:.5f}")
        return 0

    a.root.mkdir(parents=True, exist_ok=True)
    # Refresh Code/ every build rather than copying it only when absent. The T&C
    # source is edited as the project goes (Sl_min for the deciduous branch, for
    # one), and copy-once means the runs keep executing the stale tree with the
    # bug already fixed in the repo -- a failure that looks like the fix not
    # working. Report what changed so a rebuild says whether the fix landed.
    code_src, code_dst = a.template / "Code", a.root / "Code"
    src_files = sorted(p for p in code_src.rglob("*") if p.is_file())
    stale = [p.relative_to(code_src).as_posix() for p in src_files
             if not (code_dst / p.relative_to(code_src)).is_file()
             or not filecmp.cmp(p, code_dst / p.relative_to(code_src), shallow=False)]
    shutil.copytree(code_src, code_dst, dirs_exist_ok=True)
    if stale:
        shown = ", ".join(stale[:6]) + (" ..." if len(stale) > 6 else "")
        print(f"refreshed shared Code/ -- {len(stale)} of {len(src_files)} "
              f"file(s) updated: {shown}")
    else:
        print(f"shared Code/ already up to date ({len(src_files)} files)")
    g = a.template / "GRAPH_MOD.m"
    if g.is_file() and (not (a.root / "GRAPH_MOD.m").is_file()
                        or not filecmp.cmp(g, a.root / "GRAPH_MOD.m", shallow=False)):
        shutil.copy2(g, a.root / "GRAPH_MOD.m")

    for st in ready:
        probs = build_station(st, tmpl_txt, a, a.root)
        if probs:
            all_probs[st["station_id"]] = probs
        print(f"  {st['station_id']:<8} written  ms={st['_soil_block'][0]} "
              f"Kbot={'0.01' if st['_kbot'] else 'NaN'} ZR95={st['zr95']:g} "
              f"zatm={st['zatm']:g}"
              + ("   <-- SUBSTITUTION PROBLEM" if probs else ""))

    if all_probs:
        print("\nSUBSTITUTION PROBLEMS -- these MOD_PARAM files are NOT trustworthy:")
        for sid, ps in all_probs.items():
            for x in ps:
                print(f"   {sid:<8} {x}")

    (a.root / "build_manifest.json").write_text(json.dumps({
        "stations_ready": [s["station_id"] for s in ready],
        "stations_blocked": {s: m for s, m in blocked},
        "f_C": F_C, "Kbot_bedrock_mm_h": KBOT_BEDROCK,
        "zatm_fallback_offset_m": ZATM_FALLBACK_OFFSET,
        "arms": ARMS, "scenarios": SCENARIOS,
        "substitution_problems": all_probs,
        "note": "fixed_lma uses MAIN_FRAME with a constant Sl from the site mean LMA; "
                "dyn_lma uses MAIN_FRAME_SLA with the yearly series in LMA_<ST>.mat. "
                "SLA = 1/(LMA*f_C).",
    }, indent=2), encoding="utf-8")
    # The job array reads this: one "<station> <arm>" per line, in the order
    # submit_tc_run.sh indexes with SLURM_ARRAY_TASK_ID.
    (a.root / "run_list.txt").write_text(
        "".join(f"{s['station_id']} {arm}\n" for s in ready for arm in ARMS),
        encoding="utf-8")
    print(f"\nmanifest : {a.root / 'build_manifest.json'}")
    print(f"run list : {a.root / 'run_list.txt'}  ->  "
          f"sbatch --array=1-{len(ready) * len(ARMS)} slurm/submit_tc_run.sh")
    # Exit code, in dependency-chain terms: did this produce something a run can
    # trust? A station blocked for a missing input is the NORMAL outcome across a
    # whole network -- 3 of 101 in job 35712, two with no ERA5 and one with no
    # canopy height -- and it must not stop the chain. Only an untrustworthy
    # generated file, or nothing built at all, is a failure.
    if all_probs:
        return 1
    if not ready:
        print("nothing was built", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
