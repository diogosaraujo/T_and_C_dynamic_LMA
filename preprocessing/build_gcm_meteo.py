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

huss is used rather than hurs deliberately. Specific humidity is conserved through
the day in the absence of mixing, so holding q constant gives a constant Tdew and
an esat that follows hourly Ta -- i.e. VPD peaks in the afternoon, as it must.
Interpolating hurs instead would hold RH constant and flatten VPD, and VPD drives
stomatal conductance directly.

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
OUT_ROOT = INPUT_ROOT / "gcm_meteo"
CO2_DIR = INPUT_ROOT / "co2"
MATLAB_EPOCH_OFFSET = 719529.0            # datenum(1970,1,1)


def mat_name(s: str) -> str:
    return s.replace("-", "_")


def tetens(t_c):
    """Saturation vapour pressure [Pa] from temperature [C]."""
    return TETENS_A * np.exp(TETENS_B * t_c / (t_c + TETENS_C))


def tetens_inverse(e_pa):
    """Dewpoint [C] from vapour pressure [Pa] -- the exact inverse of tetens()."""
    e = np.clip(np.asarray(e_pa, dtype=float), 1e-3, None)
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
def load_station_series(gcm, scenario, station_root):
    """{var: (dates, values, stations)} for one model/scenario."""
    out = {}
    for v in VARIABLES:
        p = station_root / gcm / scenario / f"{v}.npz"
        if not p.is_file():
            return None, f"missing {p}"
        d = np.load(p, allow_pickle=False)
        out[v] = (d["dates"], d["values"], [str(x) for x in d["stations"]])
    # every variable must share one calendar, or the alignment below is a lie
    ref = out["tas"][0]
    for v, (dd, _, _) in out.items():
        if dd.shape != ref.shape or not np.array_equal(dd, ref):
            return None, f"{v} has a different calendar than tas ({len(dd)} vs {len(ref)} days)"
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


def read_elevation():
    """StationID -> Zbas [m], from the same BADM table the ERA5 path uses."""
    p = INPUT_ROOT / "ameriflux" / "badm_values.csv"
    out = {}
    if not p.is_file():
        return out
    for r in csv.DictReader(open(p, newline="", encoding="utf-8-sig")):
        sid = (r.get("StationID") or r.get("SITE_ID") or "").strip()
        for k in ("Zbas", "elevation", "ELEV", "LOCATION_ELEV"):
            if sid and r.get(k):
                try:
                    out[sid] = float(r[k]); break
                except (TypeError, ValueError):
                    pass
    return out


# ------------------------------------------------------------------------- build
def build_one(gcm, scenario, station, si, series, lat, lon, zbas, co2,
              scheme, seed, out_dir, dry):
    dates = series["tas"][0]
    days = EPOCH + dates.astype("timedelta64[D]")
    doy = (days - days.astype("datetime64[Y]")).astype(int) + 1
    n = len(days)

    get = lambda v: np.asarray(series[v][1][:, si], dtype=float)
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
    e_pa = vapour_pressure_from_q(get("huss"), p_pa)
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

    diag = dict(hours=int(out["Date"].size),
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
    # Per-SCENARIO subdirectory, because finish_meteo.m globs every
    # Meteo_*_raw.mat in a directory and stamps them all with one year tag --
    # and historical (1980-2014) and the SSPs (2015-2100) do not share one.
    out_dir = out_dir / scenario
    out_dir.mkdir(parents=True, exist_ok=True)
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
    ap.add_argument("--out", type=Path, default=OUT_ROOT)
    ap.add_argument("--co2-dir", type=Path, default=CO2_DIR)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    stations = read_stations()
    if a.station:
        stations = [s for s in stations if s["station"] in set(a.station)]
    elev = read_elevation()

    work = [(g, s) for g in (a.gcm or GCMS) for s in (a.scenario or SCENARIOS)]
    if a.index is not None:
        if a.index < 1 or a.index > len(work):
            print(f"index {a.index} outside 1..{len(work)} -- nothing to do"); return 0
        work = [work[a.index - 1]]
    elif not (a.all or a.gcm or a.scenario):
        ap.error("give --index / --gcm / --scenario / --all")

    print(f"stations : {len(stations)}   elevation known for {len(elev)}")
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
                                  a.precip_scheme, a.seed, a.out, a.dry_run)
            if diag is None:
                print(f"    {s['station']}: FAILED -- {err}"); rc = 1; continue
            ok += 1
            if a.dry_run or ok <= 2:
                print(f"    {s['station']}: {diag}")
        print(f"    -> {ok}/{len(stations)} built\n", flush=True)
    print("next: finish_meteo.m adds SAB/SAD/PAR/N (see slurm/submit_gcm_meteo.sh)")
    return rc


if __name__ == "__main__":
    sys.exit(main())
