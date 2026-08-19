#!/usr/bin/env python3
"""Daily GCM station series -> hourly T&C forcing (stage 1 of 2).

    python build_gcm_meteo.py --index 3           # array form: one GCM x scenario
    python build_gcm_meteo.py --gcm GFDL-ESM4 --scenario ssp585
    python build_gcm_meteo.py --all

Writes Meteo_<ST>_<GCM>_<scen>_raw.mat, which finish_meteo.m then completes with
SAB1/SAB2/SAD1/SAD2, PARB/PARD and N via C_Automatic_Radiation_Partition -- the
same second stage the ERA5-Land path uses, so both forcings get identical solar
geometry and cloud physics.

=============================================================================
HUMIDITY: the GCMs give no dewpoint, so Tdew is DERIVED -- and with ONE formula
=============================================================================
ERA5-Land supplies d2m directly. NEX-GDDP does not, so Tdew comes from specific
humidity in two steps:

    e    = q * P / (eps + (1 - eps) * q)          eps = 0.622, P in Pa
    Tdew = C * ln(e/A) / (B - ln(e/A))            the EXACT inverse of Tetens

where A/B/C = 611.0 / 17.27 / 237.3 are the same constants used for esat(Ta).
Using the exact inverse rather than a second approximation is the point: the
shared code paired 17.27/237.3 for esat with 17.625/243.04 for the humidity ->
dewpoint step, which disagrees by ~0.5 K and biases Ds = esat - ea everywhere.
Here ea = tetens(Tdew) reproduces e to round-off by construction, and the build
asserts it.

hurs is used for ALL five models, via e = (hurs/100) * esat(tas), so the humidity
route is identical across the ensemble. IPSL-CM6A-LR ships no huss, and a
huss-first policy would have treated one model of five differently -- putting part
of the GCM-to-GCM spread down to method rather than climate, in exactly the
comparison the five-model subset exists to make.

Either route ends at ONE daily vapour pressure held constant through the day,
with esat following hourly Ta, so VPD peaks in the afternoon regardless. The cost
of hurs is that esat is convex in temperature, so esat(daily mean) sits a few
percent below the mean of esat over the day and e carries a small low bias. That
is a bias shared by all five models, which is the point.

The GCMs also do not share a CALENDAR: GFDL-ESM4 is 365-day, UKESM1-0-LL 360-day,
the rest standard. real_dates() places each on the real axis -- see its docstring
for what that costs.

=============================================================================
DISAGGREGATION: what is constructed rather than measured
=============================================================================
Ta      diurnal cycle from tasmin/tasmax placed on the solar day (minimum at
        sunrise, maximum 2/3 through daylight), then OFFSET so the 24-hour mean
        equals tas. Shape and mean cannot both be exact -- three daily numbers
        do not determine 24 hourly ones -- and the mean is what the energy
        balance integrates, so the mean wins and the amplitude is approximate.
Rsw     daily mean scaled onto the clear-sky cosine-of-zenith curve, so the
        hourly series has the right shape and reproduces the daily mean exactly.
Pr      see --precip-scheme. This is the weakest link and is a stated choice.
Ws      held constant through the day (no sub-daily information exists).
q       held constant through the day (see above).
Pre     NOT a GCM variable: barometric from station elevation, in MILLIBAR,
        which is what T&C expects -- the ERA5 path converts sp from Pa for the
        same reason.
Ca      annual SSP pathway from input4MIPs, interpolated to the hourly stamp.

Timestamps are UTC (DeltaGMT = 0) and hour H carries the interval (H, H+1], which
is why the GCM path forces t_bef/t_aft = 0/1 -- the opposite of de-accumulated
ERA5-Land. That is imposed, not detected: a constructed series has no convention
of its own to discover.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gcm_variables import (GCMS, SCENARIOS, VARIABLES, ECS_K,        # noqa: E402
                           TETENS_A, TETENS_B, TETENS_C, P0_PA,
                           SCALE_HEIGHT_M, EPS_RATIO, T_BEF, T_AFT)
from extract_gcm_stations import read_stations, EPOCH                # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA",
                                 "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
STATION_ROOT = INPUT_ROOT / "gcm_stations"
OUT_ROOT = INPUT_ROOT / "gcm_meteo"          # RAW staging only; see build_one
CO2_DIR = INPUT_ROOT / "co2"
# Where the finished forcing lands. Same default as config.sh: a sibling of
# input_data, so the run tree and the downloaded originals stay separate.
MODEL_RUN = Path(os.environ.get("MODEL_RUN", INPUT_ROOT.parent / "model_run"))
MATLAB_EPOCH_OFFSET = 719529.0            # datenum(1970,1,1)


def _count(d):
    import collections
    return collections.Counter(d.values())


def mat_name(s: str) -> str:
    return s.replace("-", "_")


def tetens(t_c):
    """Saturation vapour pressure [Pa] from temperature [C]."""
    return TETENS_A * np.exp(TETENS_B * t_c / (t_c + TETENS_C))


# The lowest vapour pressure the inverse is allowed to see. e -> 0 sends Tdew to
# -infinity, so a floor is unavoidable: the downscaled output reports hurs = 0
# exactly at some hours, and 0% relative humidity is itself unphysical.
#
# The floor is set at Tdew = -80 C, roughly the coldest dewpoint ever observed on
# Earth (Antarctic plateau), rather than at an arbitrary small number. The earlier
# 1e-3 Pa corresponds to Tdew = -103 C, which has no physical referent and looks
# alarming in a diagnostic.
#
# The choice barely matters downstream and that is the point: at Ta = -20 C, a
# floor of -80 C gives Ds = 99.92% of the maximum possible and -103 C gives
# 99.9991%. Anything below about -60 C means "no water vapour at all" to four
# decimal places. What matters is that the floor is finite, documented, and does
# not print a number that invites a reader to think it is a temperature.
#
# Tdew is not inert in T&C: C_Automatic_Radiation_Partition uses it for the
# clear-sky precipitable-water attenuation, so an absurd value would propagate
# into the Gueymard transmittance. At -80 C that term is already saturated.
TDEW_FLOOR_C = -80.0
E_FLOOR_PA = TETENS_A * np.exp(TETENS_B * TDEW_FLOOR_C / (TDEW_FLOOR_C + TETENS_C))


def tetens_inverse(e_pa):
    """Dewpoint [C] from vapour pressure [Pa] -- the exact inverse of tetens()."""
    e = np.clip(np.asarray(e_pa, dtype=float), E_FLOOR_PA, None)
    L = np.log(e / TETENS_A)
    return TETENS_C * L / (TETENS_B - L)


def vapour_pressure_from_q(q, p_pa):
    """e [Pa] from specific humidity [kg/kg] and pressure [Pa]."""
    q = np.clip(np.asarray(q, dtype=float), 1e-12, 0.05)
    return q * p_pa / (EPS_RATIO + (1.0 - EPS_RATIO) * q)


def barometric_pressure_pa(zbas_m, ta_mean_c=None):
    """Surface pressure [Pa] from elevation. Scale height as used elsewhere here."""
    return P0_PA * np.exp(-float(zbas_m) / SCALE_HEIGHT_M)


# ------------------------------------------------------------------ solar geometry
def solar_terms(doy, lat_deg):
    """(declination rad, sunset hour angle rad) -- FAO-56 Eqs 14 and 16."""
    dec = 0.409 * np.sin(2 * np.pi * doy / 365.0 - 1.39)
    phi = np.radians(lat_deg)
    x = np.clip(-np.tan(phi) * np.tan(dec), -1.0, 1.0)
    return dec, np.arccos(x)


def cos_zenith(doy, lat_deg, lon_deg, hours_utc):
    """cos(solar zenith) for each hour, clipped at 0. Shape (ndays, 24)."""
    dec, _ = solar_terms(doy, lat_deg)
    phi = np.radians(lat_deg)
    # solar time = UTC + lon/15; hour angle zero at local solar noon
    solar_h = hours_utc[None, :] + lon_deg / 15.0
    ha = np.radians(15.0 * (solar_h - 12.0))
    cz = (np.sin(phi) * np.sin(dec)[:, None]
          + np.cos(phi) * np.cos(dec)[:, None] * np.cos(ha))
    return np.clip(cz, 0.0, None)


def daylight_hours(doy, lat_deg):
    _, ws = solar_terms(doy, lat_deg)
    return 24.0 / np.pi * ws


# ------------------------------------------------------------------ disaggregation
def hourly_temperature(tmin, tmax, tmean, doy, lat, lon):
    """Hourly Ta [C]. Shape from tmin/tmax on the solar day; mean forced to tmean.

    Minimum at sunrise, maximum 2/3 of the way through daylight, half-cosine
    between, cosine relaxation overnight. The final offset makes the 24-hour mean
    equal tas exactly -- see the module docstring for why the mean is privileged.
    """
    n = len(doy)
    h = np.arange(24, dtype=float)
    dl = daylight_hours(doy, lat)
    noon = 12.0 - lon / 15.0                       # solar noon in UTC hours
    rise = noon[None] - dl[:, None] / 2.0 if np.ndim(noon) else noon - dl / 2.0
    rise = np.broadcast_to(np.atleast_1d(rise).reshape(n, 1), (n, 24))
    dl2 = np.broadcast_to(dl.reshape(n, 1), (n, 24))
    tpk = rise + 0.667 * dl2                       # time of daily maximum

    amp = (tmax - tmin).reshape(n, 1) / 2.0
    mid = (tmax + tmin).reshape(n, 1) / 2.0
    # phase: -pi at the minimum, 0 at the peak, then continuing round to the
    # next minimum. Using a single cosine of a piecewise-scaled phase keeps the
    # curve continuous at both joins.
    hh = np.broadcast_to(h.reshape(1, 24), (n, 24))
    dt = (hh - rise) % 24.0
    span_day = (tpk - rise) % 24.0
    span_day = np.where(span_day <= 0, 1e-6, span_day)
    span_night = 24.0 - span_day
    phase = np.where(dt <= span_day,
                     np.pi * (dt / span_day + 1.0),          # rise: pi -> 2pi
                     np.pi * ((dt - span_day) / span_night))  # fall: 0 -> pi
    # +cos, not -cos: at sunrise (phase = pi) this gives mid - amp = Tmin, and at
    # the peak (phase = 2pi) mid + amp = Tmax. The night branch then runs cos from
    # +1 back to -1, so the curve is continuous at both joins and returns to Tmin
    # at the next sunrise. Getting this sign backwards puts the daily maximum at
    # dawn, which still integrates to the right daily mean and so survives every
    # check except looking at the diurnal cycle itself.
    shape = mid + amp * np.cos(phase)
    shape += (tmean.reshape(n, 1) - shape.mean(axis=1, keepdims=True))
    return shape


def hourly_shortwave(rsds_daily, doy, lat, lon):
    """Hourly Rsw [W/m2]: the daily mean redistributed onto the clear-sky curve."""
    cz = cos_zenith(doy, lat, lon, np.arange(24, dtype=float))
    denom = cz.mean(axis=1, keepdims=True)
    out = np.where(denom > 0, rsds_daily.reshape(-1, 1) * cz / np.maximum(denom, 1e-9), 0.0)
    return np.clip(out, 0.0, None)


def hourly_precip(pr_daily_mm, scheme, doy, rng):
    """Hourly Pr [mm]. The daily total is preserved exactly by every scheme."""
    n = len(pr_daily_mm)
    if scheme == "uniform":
        return np.repeat(pr_daily_mm.reshape(n, 1) / 24.0, 24, axis=1)
    if scheme == "block":
        # Concentrate the day's total into a contiguous wet block, which restores
        # intensity. 6 h is the round-number compromise: long enough not to
        # manufacture extreme rates, short enough that canopy interception
        # saturates and throughfall happens, which uniform spreading prevents.
        nwet = 6
        out = np.zeros((n, 24))
        start = rng.integers(0, 24, size=n)
        idx = (start.reshape(n, 1) + np.arange(nwet).reshape(1, nwet)) % 24
        np.put_along_axis(out, idx,
                          np.repeat(pr_daily_mm.reshape(n, 1) / nwet, nwet, axis=1),
                          axis=1)
        return out
    raise ValueError(f"unknown precip scheme: {scheme}")


# ------------------------------------------------------------------------ inputs
def real_dates(ymd, doy, calendar):
    """Place a GCM calendar on the real one. Returns (datetime64[D], keep mask, note).

    The GCMs do not share a calendar. GFDL-ESM4 is 365-day, UKESM1-0-LL is
    360-day, the rest are standard, and the extraction stores date COMPONENTS
    precisely so this decision lives here rather than failing at read time.

      standard / noleap / all_leap  (year, month, day) IS a real date, so it is
            used directly. A 365-day model simply never emits 29 February, which
            leaves a one-day gap in leap years -- harmless, since T&C indexes the
            hourly loop rather than assuming a fixed year length.
      360_day  every month has 30 days, so 30 February exists and 31 January does
            not; no direct mapping is possible. The day-of-year is stretched onto
            the real year, which is the standard ISIMIP treatment.
    """
    y, m, d = ymd[:, 0], ymd[:, 1], ymd[:, 2]
    if str(calendar) in ("360_day", "360"):
        out = np.empty(len(y), dtype="datetime64[D]")
        for k, (yy, dd) in enumerate(zip(y, doy)):
            leap = (yy % 4 == 0 and (yy % 100 != 0 or yy % 400 == 0))
            ny = 366 if leap else 365
            rd = int(round((dd - 0.5) * ny / 360.0))
            out[k] = np.datetime64(f"{yy:04d}-01-01", "D") + min(max(rd, 0), ny - 1)
        return out, np.ones(len(y), bool), "360-day stretched onto the real year"
    out = np.empty(len(y), dtype="datetime64[D]")
    keep = np.ones(len(y), bool)
    for k in range(len(y)):
        try:
            out[k] = np.datetime64(f"{y[k]:04d}-{m[k]:02d}-{d[k]:02d}", "D")
        except ValueError:                    # e.g. 29 Feb from an all_leap model
            keep[k] = False
    note = None if keep.all() else f"{(~keep).sum()} date(s) absent from the real calendar"
    return out, keep, note


def load_station_series(gcm, scenario, station_root):
    """{var: (values, stations)} plus shared date axis, for one model/scenario."""
    out, meta = {}, {}
    for v in VARIABLES:
        p = station_root / gcm / scenario / f"{v}.npz"
        if not p.is_file():
            return None, f"missing {p}"
        d = np.load(p, allow_pickle=False)
        if "ymd" not in d:
            return None, (f"{p.name} predates the calendar fix -- re-extract with "
                          f"submit_gcm_extract.sh --force")
        out[v] = (d["ymd"], d["values"], [str(x) for x in d["stations"]])
        meta[v] = dict(calendar=str(d["calendar"]), source_var=str(d["source_var"]),
                       doy=d["doy"])
    # Every variable must share one date axis, or the column-wise alignment below
    # is silently wrong rather than merely misaligned.
    ref = out["tas"][0]
    for v, (yy, _, _) in out.items():
        if yy.shape != ref.shape or not np.array_equal(yy, ref):
            return None, (f"{v} has a different date axis than tas "
                          f"({len(yy)} vs {len(ref)} days)")
    out["_meta"] = meta
    return out, None


def read_co2(scenario, co2_dir):
    """(years, ppm) for one scenario, written by fetch_ssp_co2.py."""
    p = co2_dir / f"co2_{scenario}.csv"
    if not p.is_file():
        return None, None
    y, c = [], []
    for r in csv.DictReader(open(p, newline="", encoding="utf-8-sig")):
        try:
            y.append(int(r["year"])); c.append(float(r["co2_ppm"]))
        except (KeyError, TypeError, ValueError):
            continue
    if not y:
        return None, None
    o = np.argsort(y)
    return np.array(y)[o], np.array(c)[o]


def read_elevation(ameriflux_dir=None):
    """StationID -> Zbas [m], reusing the ERA5 path's reader.

    Not re-implemented here. Elevation lives in three different places -- the
    long-format badm_values.csv, the site registry in site_metadata.json, and
    badm_coverage.csv where it is stored as text like '2753 (site registry)' --
    and build_meteo_input.read_elevation already tries all three in order of
    directness. A second guess at the schema found zero of 101 stations (job
    37223) and skipped every one of them, which is exactly the drift that reusing
    the tested reader avoids.
    """
    from build_meteo_input import read_elevation as _read
    out, src = _read(ameriflux_dir or (INPUT_ROOT / "ameriflux"))
    return out, src


# ------------------------------------------------------------------------- build
def build_one(gcm, scenario, station, si, series, lat, lon, zbas, co2,
              scheme, seed, out_dir, dry, dest_root=MODEL_RUN):
    meta = series["_meta"]
    cal = meta["tas"]["calendar"]
    days, keep, cal_note = real_dates(series["tas"][0], meta["tas"]["doy"], cal)
    # Clip to the scenario's run window before anything else reads the arrays.
    # Folding it into `keep` rather than trimming afterwards means the day axis,
    # doy, and every variable through get() are cut by the same mask, so they
    # cannot drift apart. The extraction deliberately keeps a wider archive than
    # the run needs (historical NEX-GDDP goes back to 1950), so this is where the
    # window is actually imposed -- the year filter is on the stored year COLUMN,
    # which is exact for every calendar including the 360-day stretch.
    y0, y1 = SCENARIOS[scenario]
    yr = np.asarray(series["tas"][0][:, 0], dtype=int)
    keep = keep & (yr >= y0) & (yr <= y1)
    days = days[keep]
    doy = (days - days.astype("datetime64[Y]")).astype(int) + 1
    n = len(days)
    if n == 0:
        return None, f"no days inside {y0}-{y1}"

    get = lambda v: np.asarray(series[v][1][keep, si], dtype=float)
    tas, tmax, tmin = get("tas") - 273.15, get("tasmax") - 273.15, get("tasmin") - 273.15
    # tasmax < tasmin happens in a handful of downscaled cells; swap rather than
    # propagate a negative amplitude into the diurnal cycle.
    bad = tmax < tmin
    if bad.any():
        tmax[bad], tmin[bad] = tmin[bad], tmax[bad]

    Ta = hourly_temperature(tmin, tmax, tas, doy, lat, lon)
    Rsw = hourly_shortwave(get("rsds"), doy, lat, lon)
    Pr = hourly_precip(get("pr") * 86400.0, scheme, doy,
                       np.random.default_rng(seed + si))
    Ws = np.repeat(get("sfcWind").reshape(n, 1), 24, axis=1)

    p_pa = barometric_pressure_pa(zbas)
    # Whichever humidity variable the extraction found, the destination is the
    # same: a daily vapour pressure, then Tdew through the exact Tetens inverse.
    # huss is preferred; hurs is the fallback where a model ships no huss
    # (IPSL-CM6A-LR). Both give a constant Tdew through the day, so VPD still
    # follows hourly Ta -- the difference is that hurs was itself derived against
    # the model's own daily temperature, which adds a little noise but no bias.
    hum_src = meta["hurs"]["source_var"]
    if hum_src == "hurs":
        # NEX-GDDP occasionally reports RH slightly above 100 after downscaling;
        # the drought-layer methods clip it the same way.
        e_pa = np.clip(get("hurs"), 0.0, 100.0) / 100.0 * tetens(tas)
    else:
        e_pa = vapour_pressure_from_q(get("hurs"), p_pa)
    # Apply the floor HERE, not inside the inverse. hurs == 0 appears in the
    # downscaled output, which makes e exactly 0; flooring inside tetens_inverse
    # left the round-trip check comparing a floored value against an unfloored
    # one, and it failed by exactly the floor (1.00e-03 Pa) at 8-11 stations per
    # task in job 37232. Flooring once, before both uses, makes the identity hold
    # and keeps the check meaningful.
    n_floored = int(np.sum(e_pa < E_FLOOR_PA))   # days the GCM reported hurs ~ 0
    e_pa = np.maximum(e_pa, E_FLOOR_PA)
    Tdew_d = tetens_inverse(e_pa)
    # Dewpoint cannot exceed the air temperature; where the GCM's q implies it,
    # cap at saturation rather than emit a negative vapour-pressure deficit.
    Tdew = np.minimum(np.repeat(Tdew_d.reshape(n, 1), 24, axis=1), Ta)
    ea, esat = tetens(Tdew), tetens(Ta)

    # The exact-inverse claim in the docstring, asserted rather than trusted.
    chk = np.nanmax(np.abs(tetens(Tdew_d) - e_pa))
    if chk > 1e-6:
        return None, f"Tetens inverse disagrees by {chk:.2e} Pa"

    flat = lambda a: a.reshape(-1)
    hours = np.tile(np.arange(24), n)
    datenum = (MATLAB_EPOCH_OFFSET
               + np.repeat((days - np.datetime64("1970-01-01")).astype(int), 24)
               + hours / 24.0)

    out = dict(Date=datenum, Ta=flat(Ta), Tdew=flat(Tdew), Ws=flat(Ws),
               Pr=flat(Pr), Rsw=flat(Rsw), esat=flat(esat), ea=flat(ea))
    out["Ds"] = np.clip(out["esat"] - out["ea"], 0.0, None)
    out["U"] = out["ea"] / out["esat"]
    out["Pre"] = np.full(out["Date"].size, p_pa / 100.0)      # Pa -> MILLIBAR
    if co2[0] is not None:
        yr_frac = np.repeat(days.astype("datetime64[Y]").astype(int) + 1970, 24)
        out["Ca"] = np.interp(yr_frac, co2[0], co2[1])
    else:
        out["Ca"] = np.full(out["Date"].size, np.nan)
    out.update(Lat=float(lat), Lon=float(lon), Zbas=float(zbas), DeltaGMT=0.0,
               t_bef=T_BEF, t_aft=T_AFT, id_location=mat_name(station))

    dead = [k for k in ("Ta", "Rsw", "Pr", "Ws", "Tdew") if not np.isfinite(out[k]).any()]
    if dead:
        return None, f"all-NaN in {', '.join(dead)}"

    diag = dict(hours=int(out["Date"].size), calendar=cal, humidity=hum_src,
                cal_note=cal_note,
                zero_rh_days=n_floored,
                years=f"{int(days[0].astype('datetime64[Y]').astype(int))+1970}-"
                      f"{int(days[-1].astype('datetime64[Y]').astype(int))+1970}",
                Ta_C=[round(float(np.nanmin(out["Ta"])), 1),
                      round(float(np.nanmax(out["Ta"])), 1)],
                Ta_mean=round(float(np.nanmean(out["Ta"])), 2),
                Pre_mbar=round(float(out["Pre"][0]), 1),
                Pr_mm_yr=round(float(np.nansum(out["Pr"])) / max(1, n / 365.25), 1),
                Pr_max_h=round(float(np.nanmax(out["Pr"])), 2),
                Rsw_max=round(float(np.nanmax(out["Rsw"])), 1),
                Ds_mean=round(float(np.nanmean(out["Ds"])), 1),
                Ca=[round(float(np.nanmin(out["Ca"])), 1),
                    round(float(np.nanmax(out["Ca"])), 1)]
                if np.isfinite(out["Ca"]).any() else None)
    if dry:
        return diag, None

    from scipy.io import savemat
    # One directory per (scenario, GCM). Two reasons, and BOTH are needed:
    #   * finish_meteo.m globs every Meteo_*_raw.mat in a directory and stamps
    #     them all with one year tag, and historical (1985-2014) and the SSPs
    #     (2015-2100) do not share one.
    #   * the 15 array tasks run CONCURRENTLY, so a directory shared between
    #     GCMs means every task partitions every file that happens to be there
    #     and several write the same output at once. Job 37232 did exactly that:
    #     task 1 processed 283 files instead of its own 101, and MATLAB failed
    #     with "appears to be corrupt" on a file another task was writing.
    # Per-task ownership of the output directory is what makes the stage safe to
    # parallelise at all.
    #
    # This is the RAW staging area only. The finished file goes to dest_dir below,
    # inside model_run, and the raw file here is an intermediate that can be
    # deleted once finish_meteo.m has consumed it.
    out_dir = out_dir / scenario / mat_name(gcm)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Where finish_meteo.m must put the finished forcing: one copy per (station,
    # scenario, GCM), directly above the fixed_lma/dyn_lma pair that reads it. The
    # GCM directory is spelled as model_run spells it (GFDL-ESM4, hyphenated) --
    # mat_name's underscored form is for MATLAB identifiers and file names, not
    # for directories in the run tree.
    out["dest_dir"] = (dest_root / station / scenario / gcm).as_posix()

    savemat(out_dir / f"Meteo_{mat_name(station)}_{mat_name(gcm)}_{scenario}_raw.mat",
            {k: (np.asarray(v).reshape(-1, 1) if isinstance(v, np.ndarray) else v)
             for k, v in out.items()}, do_compression=True)
    return diag, None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", type=int, help="1-based index into (GCM x scenario)")
    ap.add_argument("--gcm", action="append")
    ap.add_argument("--scenario", action="append")
    ap.add_argument("--station", action="append")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--precip-scheme", choices=["uniform", "block"], default="block",
                    help="daily -> hourly precipitation. 'block' concentrates the "
                         "total into 6 contiguous hours and keeps intensity; "
                         "'uniform' spreads it over 24 h and systematically "
                         "under-produces throughfall and runoff. Default: block.")
    ap.add_argument("--seed", type=int, default=20260813,
                    help="fixes the block placement so builds are reproducible")
    ap.add_argument("--stations-root", type=Path, default=STATION_ROOT)
    ap.add_argument("--out", type=Path, default=OUT_ROOT,
                    help="RAW staging root; the finished forcing goes to "
                         "--model-run, not here")
    ap.add_argument("--model-run", type=Path, default=MODEL_RUN,
                    help="run tree the finished forcing is written into, at "
                         "<ST>/<scenario>/<GCM>/")
    ap.add_argument("--co2-dir", type=Path, default=CO2_DIR)
    ap.add_argument("--ameriflux", type=Path, default=None,
                    help="directory holding badm_values.csv / site_metadata.json "
                         "(default: $TC_INPUT_DATA/ameriflux)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    stations = read_stations()
    if a.station:
        stations = [s for s in stations if s["station"] in set(a.station)]
    elev, elev_src = read_elevation(a.ameriflux)

    work = [(g, s) for g in (a.gcm or GCMS) for s in (a.scenario or SCENARIOS)]
    if a.index is not None:
        if a.index < 1 or a.index > len(work):
            print(f"index {a.index} outside 1..{len(work)} -- nothing to do"); return 0
        work = [work[a.index - 1]]
    elif not (a.all or a.gcm or a.scenario):
        ap.error("give --index / --gcm / --scenario / --all")

    print(f"stations : {len(stations)}   elevation known for {len(elev)}"
          + (f" ({', '.join(f'{k}={v}' for k, v in sorted(_count(elev_src).items()))})"
             if elev else ""))
    if not elev:
        print("  ! no elevation for ANY station -- Pre and the radiation partition "
              "both need Zbas.\n"
              "    Check that $TC_INPUT_DATA/ameriflux holds badm_values.csv, "
              "site_metadata.json or badm_coverage.csv,\n"
              "    or pass --ameriflux <dir>.", file=sys.stderr)
    print(f"input    : {a.stations_root}")
    print(f"output   : {a.out}")
    print(f"precip   : {a.precip_scheme}   seed {a.seed}\n")

    rc = 0
    for g, sc in work:
        series, err = load_station_series(g, sc, a.stations_root)
        if series is None:
            print(f"  {g} {sc}: SKIPPED -- {err}"); rc = 1; continue
        co2 = read_co2(sc, a.co2_dir)
        if co2[0] is None:
            print(f"  ! no CO2 series for {sc} -- run fetch_ssp_co2.py; Ca will be NaN")
        names = series["tas"][2]
        print(f"  {g} {sc} (ECS {ECS_K.get(g,'?')} K): "
              f"{len(series['tas'][0])} days, {len(names)} stations")
        ok = 0
        for s in stations:
            if s["station"] not in names:
                print(f"    {s['station']}: not in the extraction"); continue
            z = elev.get(s["station"])
            if z is None:
                print(f"    {s['station']}: no elevation -- skipped"); rc = 1; continue
            diag, err = build_one(g, sc, s["station"], names.index(s["station"]),
                                  series, s["lat"], s["lon"], z, co2,
                                  a.precip_scheme, a.seed, a.out, a.dry_run,
                                  a.model_run)
            if diag is None:
                print(f"    {s['station']}: FAILED -- {err}"); rc = 1; continue
            ok += 1
            if a.dry_run or ok <= 2:
                print(f"    {s['station']}: {diag}")
        print(f"    -> {ok}/{len(stations)} built")
        # Machine-readable marker for the submit script. Stage 2 must partition
        # ONLY what this task wrote: finish_meteo.m globs a directory, and a shell
        # glob over the whole tree makes all 15 concurrent tasks partition all 15
        # directories and overwrite each other. Job 37293 did exactly that -- every
        # task reported 15 partition passes of 101 files, and MATLAB failed with
        # "appears to be corrupt" on a file another task was writing.
        if not a.dry_run and ok:
            print(f"BUILT_DIR: {(a.out / sc / mat_name(g)).resolve()}")
        print(flush=True)
    print("next: finish_meteo.m adds SAB/SAD/PAR/N (see slurm/submit_gcm_meteo.sh)")
    return rc


if __name__ == "__main__":
    sys.exit(main())
