#!/usr/bin/env python3
"""Registry for the NEX-GDDP-CMIP6 forcing: models, scenarios, variables, paths.

One place for everything that is a fact about the dataset rather than a choice
about the science, so extract_gcm_stations.py and build_gcm_meteo.py cannot drift
apart on a variable name or a unit.

WHICH MODELS. The five here are the ISIMIP3b primary set, which CHELSA also
adopts. They were selected for structural independence of their ocean and
atmosphere components, fair-to-good process representation (CRESCENDO expert
survey), and availability of all required daily variables -- and the resulting
subset reproduces the CMIP6 ensemble mean ECS of 3.7 K and TCR of 2.0 K, with
three low-sensitivity models (GFDL-ESM4, MPI-ESM1-2-HR, MRI-ESM2-0) and two high
(IPSL-CM6A-LR, UKESM1-0-LL). Using a published subset means the choice is citable
rather than something we have to defend from scratch.

WHICH VARIABLES. T&C reads 7 of the 9 NEX-GDDP variables:

    tas tasmax tasmin   ->  Ta, and the diurnal amplitude
    hurs                ->  Tdew, ea      (see HUMIDITY below)
    pr                  ->  Pr
    rsds                ->  Rsw, then SAB/SAD/PAR/N in the MATLAB stage
    sfcWind             ->  Ws

rlds is NOT used: T&C computes incoming longwave internally from Ta, ea and N via
Incoming_Longwave(), so supplying it would be ignored. huss is not used either --
hurs is taken for all five models so the humidity route is identical across the
ensemble; see HUMIDITY below for what that costs and why it is the right trade.

Pre is absent from NEX-GDDP entirely and is computed barometrically from station
elevation (see build_gcm_meteo.py).

FILENAMES are globbed rather than constructed. The variant label is not constant
across models -- UKESM1-0-LL is r1i1p1f2 while the other four are r1i1p1f1 -- and
the grid label can be gn or gr. Building the name from parts would silently miss
files; globbing fails loudly instead.
"""
from __future__ import annotations

import os
from pathlib import Path

# Root of the NEX-GDDP store on the SOE cluster. Overridable so the same code can
# run against a local copy for testing.
NEXGDDP_ROOT = Path(os.environ.get(
    "NEXGDDP_ROOT", "/vol_efthymios/NFS07/Data/CMIP6/NEXGDDP"))

# ISIMIP3b / CHELSA primary five.
GCMS = ["GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"]

# Equilibrium climate sensitivity, for the methods table and for checking the
# subset still brackets the full ensemble. IPCC AR6 / CMIP6 reported values.
ECS_K = {"GFDL-ESM4": 2.6, "MPI-ESM1-2-HR": 3.0, "MRI-ESM2-0": 3.2,
         "IPSL-CM6A-LR": 4.6, "UKESM1-0-LL": 5.4}

# scenario -> (first year, last year). NEX-GDDP historical starts in 1950 but the
# drought layers this pairs with use 1980-2014, so we match that.
SCENARIOS = {"historical": (1980, 2014),
             "ssp126": (2015, 2100),
             "ssp585": (2015, 2100)}

# NEX-GDDP name -> what it becomes in the T&C forcing. Only these 7 are read.
VARIABLES = {
    "tas":     dict(units="K",        role="daily mean air temperature"),
    "tasmax":  dict(units="K",        role="daily maximum, sets diurnal amplitude"),
    "tasmin":  dict(units="K",        role="daily minimum, sets diurnal amplitude"),
    "hurs":    dict(units="%",        role="relative humidity -> Tdew, ea"),
    "pr":      dict(units="kg/m2/s",  role="precipitation -> Pr [mm/h]"),
    "rsds":    dict(units="W/m2",     role="downwelling shortwave -> Rsw"),
    "sfcWind": dict(units="m/s",      role="10 m wind speed -> Ws"),
}
UNUSED = {
    "rlds": "T&C computes incoming longwave internally via Incoming_Longwave(Ta,ea,N)",
    "huss": "hurs is used instead, for ensemble consistency -- see HUMIDITY below",
}

# HUMIDITY: hurs everywhere, huss only as a fallback.
#
# IPSL-CM6A-LR ships no huss in this archive, so a huss-first policy would treat
# one model of five differently from the other four. That is the wrong trade: a
# model-dependent humidity route puts part of the GCM-to-GCM spread down to
# method rather than climate, and cross-model spread is exactly what the
# five-model subset exists to measure. One route for all five is worth more than
# a marginally better route for four.
#
# What it costs, stated rather than waved away: from hurs the daily vapour
# pressure is e = (hurs/100) * esat(tas), and esat is convex in temperature, so
# esat evaluated at the daily MEAN sits a few percent below the mean of esat over
# the day. e therefore carries a small low bias that huss would not have. The
# STRUCTURE is unaffected -- both routes yield one daily vapour pressure held
# constant through the day, with esat following hourly Ta, so VPD still peaks in
# the afternoon either way. The earlier concern about hurs flattening VPD applies
# only to interpolating RH hourly, which is not what is done here.
HUMIDITY_PREFERENCE = ["hurs", "huss"]

# Physical constants, matching what the ERA5 path already uses so the two forcings
# are built with identical thermodynamics.
TETENS_A, TETENS_B, TETENS_C = 611.0, 17.27, 237.3   # esat = A*exp(B*T/(T+C)), Pa
P0_PA = 101325.0            # sea-level pressure
SCALE_HEIGHT_M = 8434.5     # barometric scale height, as used for Zbas elsewhere
EPS_RATIO = 0.622           # Rd/Rv, for q -> vapour pressure

# The GCM path has no native radiation timestamp convention -- the hourly series is
# something we construct -- so the convention is imposed rather than detected. The
# disaggregation puts the value at hour H over the interval (H, H+1], the opposite
# of de-accumulated ERA5-Land, hence 0/1 rather than 1/0.
T_BEF, T_AFT = 0.0, 1.0


def var_dir(gcm: str, scenario: str, var: str) -> Path:
    return NEXGDDP_ROOT / gcm / scenario / var


def find_year_files(gcm: str, scenario: str, var: str) -> dict[int, Path]:
    """{year: path} for one model/scenario/variable, by globbing.

    The variant (r1i1p1f1 vs r1i1p1f2) and grid label (gn vs gr) are not constant
    across models, so the filename is matched rather than built.
    """
    d = var_dir(gcm, scenario, var)
    if not d.is_dir():
        return {}
    out = {}
    for p in sorted(d.glob(f"{var}_day_{gcm}_{scenario}_*.nc")):
        stem = p.stem
        # ..._<grid>_<year>[_v...] -- the year is the last all-digit 4-char field
        for tok in reversed(stem.split("_")):
            if len(tok) == 4 and tok.isdigit():
                out[int(tok)] = p
                break
    return out


def expected_years(scenario: str) -> list[int]:
    y0, y1 = SCENARIOS[scenario]
    return list(range(y0, y1 + 1))


def tasks(gcms=None, scenarios=None, variables=None):
    """The (gcm, scenario, variable) work list, in a stable order for job arrays."""
    return [(g, s, v)
            for g in (gcms or GCMS)
            for s in (scenarios or SCENARIOS)
            for v in (variables or VARIABLES)]
