#!/usr/bin/env python3
"""Turn the downloaded ERA5-Land netCDFs into T&C forcing .mat files.

Two stages, because the radiation partition is not worth porting:

    build_meteo_input.py   ERA5-Land netCDF -> Meteo_<ST>_raw.mat
    finish_meteo.m         calls C_Automatic_Radiation_Partition, writes
                           Meteo_<ST>_<years>.mat with SAB/SAD/PAR/N filled in

C_Automatic_Radiation_Partition.m is several hundred lines of solar geometry,
Gueymard clear-sky and Slingo cloud physics, already written and tested. Porting
it to Python would repeat the era5_predictors exercise, where a careful port still
shipped two bugs that only a verification harness caught. MATLAB is on the cluster;
reuse the original.

UNITS -- read from the working Meteo_US_xRM_1985_2020.mat, not assumed:

    Ta, Tdew   degrees C          ERA5 t2m, d2m are K          -> -273.15
    Pre        MILLIBAR (hPa)     ERA5 sp is Pa                -> /100
    ea, esat   Pa                 Tetens, one formula for both
    Pr         mm/h               ERA5 tp is m accumulated     -> deaccumulate, *1000
    Ws         m/s                sqrt(u10^2 + v10^2)
    Rsw        W/m2               ERA5 ssrd is J/m2 accumulated-> deaccumulate, /3600
    Ca         ppm                Ca_Data.mat, hourly
    Date       MATLAB datenum

Pre in mbar is the one that would have been silently catastrophic: ERA5 delivers
~72000 Pa where T&C wants ~720, and nothing downstream would complain.

ACCUMULATION IS NOT UNIFORM ACROSS THE VARIABLES, which is what fixes t_bef/t_aft
downstream:

    ACCUMULATED from 00 UTC, reset daily   ssrd, tp
        the value at hour H is the total since midnight, not the hourly flux.
        De-accumulated here by within-day differencing (the 01 UTC value is
        already the hourly total), after which hour H means the mean over (H-1, H].
    INSTANTANEOUS at the timestamp         t2m, d2m, sp, u10, v10
        the value AT hour H. Nothing to de-accumulate.

finish_meteo.m therefore takes t_bef/t_aft from ssrd alone (1 and 0), since that
window exists to align solar geometry with the RADIATION timestamp. The residual
half-hour offset between instantaneous state variables and interval-mean fluxes is
left uncorrected and documented there.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA",
                                 "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
DEFAULT_ERA5 = INPUT_ROOT / "era5_land"
DEFAULT_OUT = INPUT_ROOT / "meteo"
# Ca_Data.mat lives with the T&C source inputs, not in T&C/Diogo. It is hourly
# 1975-2022 at 327-420 ppm; Ca_Data_Ann.mat sits beside it but stops in 2013 and
# is the wrong one for a 1985-2020 run.
DEFAULT_CA = REPO_ROOT / "T&C" / "TeC_Source_Code-master" / "Inputs" / "Ca_Data.mat"
DEFAULT_SITE_LISTS = [
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv",
]

# ERA5-Land short names -> what we need them for.
WANT = {"t2m": "Ta", "sp": "Pre", "d2m": "Tdew", "tp": "Pr",
        "u10": "u", "v10": "v", "ssrd": "Rsw"}
MATLAB_EPOCH_OFFSET = 719529.0    # datenum(1970,1,1)


def mat_name(s: str) -> str:
    return re.sub(r"[^0-9A-Za-z_]", "_", s)


def tetens(t_c):
    """Saturation vapour pressure [Pa] from temperature [degC]. One formula for
    both esat(Ta) and ea(Tdew) -- mixing 17.27/237.3 with 17.625/243.04 costs
    about 0.5 degC of consistency (CLAUDE.md section 4)."""
    return 611.0 * np.exp(17.27 * t_c / (t_c + 237.3))


def looks_accumulated(values, hours_utc, sample_days=60):
    """Is this field accumulated within the UTC day, or already an hourly flux?

    The gridded ERA5-Land product accumulates from 00 UTC; the newer time-series
    collection may deliver hourly fluxes directly. Getting it wrong is not loud:
    differencing an already-hourly field yields a plausible-looking series with
    roughly half the correct peak, which is what the first dry run produced
    (Rsw max ~550 W/m2 at 43 degrees N, where ~950 is expected).

    An accumulated field is non-decreasing between the 01 and 23 UTC steps of a
    day. An hourly flux is not: it falls every afternoon.
    """
    v = np.asarray(values, dtype=float)
    day = np.cumsum(hours_utc == 0)
    ok = drops = 0
    for d in np.unique(day)[1:sample_days + 1]:
        sel = (day == d) & (hours_utc >= 1)
        x = v[sel]
        if x.size < 6 or not np.isfinite(x).all():
            continue
        ok += 1
        if np.any(np.diff(x) < -1e-9 * max(1.0, np.nanmax(np.abs(x)))):
            drops += 1
    if ok == 0:
        return True, "undetermined, assumed accumulated"
    frac = drops / ok
    if frac < 0.1:
        return True, f"accumulated (within-day decrease on {drops}/{ok} days)"
    return False, f"already hourly (within-day decrease on {drops}/{ok} days)"


def deaccumulate(values, hours_utc):
    """ERA5-Land accumulations reset at 00 UTC: value(H) is the total since
    midnight. Difference within the day; the 01 UTC value already is the hourly
    total. Negatives from the reset boundary are clipped."""
    out = np.diff(values, prepend=np.nan)
    first = hours_utc <= 1
    out[first] = values[first]
    out[0] = values[0]
    return np.clip(np.nan_to_num(out, nan=0.0), 0.0, None)


def read_station_era5(station_dir: Path):
    """Merge a station's per-group netCDFs onto one hourly time axis."""
    try:
        import xarray as xr
    except ImportError:
        sys.exit("xarray is required: pip install -r requirements.txt")
    files = sorted(p for p in station_dir.glob("*.nc"))
    if not files:
        return None, "no .nc files"
    data, times = {}, None
    for f in files:
        with xr.open_dataset(f) as ds:
            tname = "valid_time" if "valid_time" in ds else (
                "time" if "time" in ds else None)
            if tname is None:
                continue
            t = ds[tname].values
            if times is None:
                times = t
            elif len(t) != len(times):
                return None, f"time axis mismatch in {f.name}: {len(t)} vs {len(times)}"
            for short in ds.data_vars:
                if short in WANT:
                    data[short] = np.asarray(ds[short].values, dtype=float).ravel()
    missing = [k for k in WANT if k not in data]
    if missing:
        return None, f"missing variable(s): {', '.join(missing)}"
    return (times, data), ""


def read_ca(path: Path):
    """Hourly CO2 [ppm] with its own datenum axis, from Ca_Data.mat."""
    if not Path(path).is_file():
        return None, None
    try:
        from scipy.io import loadmat
        d = loadmat(path, squeeze_me=True)
    except Exception:                                              # noqa: BLE001
        try:
            import h5py
            with h5py.File(path, "r") as f:
                d = {k: np.array(f[k]).ravel() for k in f if not k.startswith("#")}
        except Exception:                                          # noqa: BLE001
            return None, None
    ca = next((np.asarray(d[k], dtype=float).ravel()
               for k in ("Ca", "Ca_all", "ca") if k in d), None)
    # The axis is called Date_CO2 in the shipped file, not Date.
    dt = next((np.asarray(d[k], dtype=float).ravel()
               for k in ("Date_CO2", "Date", "Date_Ca", "t") if k in d), None)
    if ca is None or dt is None or dt.size != ca.size:
        return None, None
    return dt, ca


def read_elevation(ameriflux_dir: Path, override: Path | None = None):
    """Zbas [m] per station, from the AmeriFlux data already downloaded.

    inspect_ameriflux_badm.py tracks 'elevation' as a parameter and fills it from
    the site registry (GRP_LOCATION/LOCATION_ELEV) where BADM itself omits it, so
    three places can hold it. Tried in order of directness:

        1. badm_values.csv     long format, parameter == 'elevation'
        2. site_metadata.json  GRP_LOCATION.LOCATION_ELEV, the registry value
        3. badm_coverage.csv   which stores strings like '2753 (site registry)'

    Zbas should be ORTHOMETRIC (above the geoid). AmeriFlux reports site elevation
    that way; the sensitivity is small in any case, since the barometric scale
    height is 8434.5 m, but the ERA5 grid-cell elevation would be the wrong thing
    entirely in complex terrain.
    """
    out, src = {}, {}

    def take(sid, raw, where):
        if not sid or sid in out:
            return
        m = re.search(r"-?\d+(?:\.\d+)?", str(raw))
        if not m:
            return
        v = float(m.group(0))
        if -500 < v < 9000:            # a plausible land elevation
            out[sid], src[sid] = v, where

    vals = ameriflux_dir / "badm_values.csv"
    if vals.is_file():
        for r in csv.DictReader(open(vals, newline="", encoding="utf-8-sig")):
            if (r.get("parameter") or "").strip().lower() == "elevation" or                     "ELEV" in (r.get("variable") or "").upper():
                take((r.get("station_id") or "").strip(), r.get("value"), "badm_values")

    meta = ameriflux_dir / "site_metadata.json"
    if meta.is_file():
        try:
            j = json.loads(meta.read_text(encoding="utf-8"))
            recs = j.values() if isinstance(j, dict) else j
            for rec in recs:
                if not isinstance(rec, dict):
                    continue
                sid = (rec.get("SITE_ID") or rec.get("site_id") or "").strip()
                loc = rec.get("GRP_LOCATION") or {}
                if isinstance(loc, list):
                    loc = loc[0] if loc else {}
                take(sid, (loc or {}).get("LOCATION_ELEV"), "site_registry")
        except (json.JSONDecodeError, OSError, AttributeError):
            pass

    cov = ameriflux_dir / "badm_coverage.csv"
    if cov.is_file():
        for r in csv.DictReader(open(cov, newline="", encoding="utf-8-sig")):
            take((r.get("station_id") or "").strip(), r.get("elevation"), "badm_coverage")

    if override and Path(override).is_file():
        for r in csv.DictReader(open(override, newline="", encoding="utf-8-sig")):
            sid = (r.get("station_id") or r.get("StationID") or "").strip()
            for c in ("elevation_m", "elev_m", "Elevation", "Zbas"):
                if r.get(c) not in (None, ""):
                    out.pop(sid, None)
                    take(sid, r[c], "override")
                    break
    return out, src


def build(station, lat, lon, zbas, era5_dir, ca, years, out_dir, dry):
    sd = era5_dir / station
    if not sd.is_dir():
        return None, "no ERA5-Land directory"
    got, err = read_station_era5(sd)
    if got is None:
        return None, err
    times, raw = got

    # MATLAB datenum from numpy datetime64.
    secs = times.astype("datetime64[s]").astype(np.int64)
    datenum = MATLAB_EPOCH_OFFSET + secs / 86400.0
    dt64 = times.astype("datetime64[h]")
    hours = (dt64.astype(np.int64) % 24).astype(int)
    yrs = dt64.astype("datetime64[Y]").astype(int) + 1970

    keep = (yrs >= years[0]) & (yrs <= years[1])
    if not keep.any():
        return None, f"no hours inside {years[0]}-{years[1]}"

    Ta = raw["t2m"] - 273.15
    Tdew = raw["d2m"] - 273.15
    Pre = raw["sp"] / 100.0                       # Pa -> mbar
    Ws = np.hypot(raw["u10"], raw["v10"])
    # Detect rather than assume: the time-series collection may already deliver
    # hourly fluxes, and differencing those halves the peak without complaining.
    tp_acc, tp_why = looks_accumulated(raw["tp"], hours)
    sw_acc, sw_why = looks_accumulated(raw["ssrd"], hours)
    tp_h = deaccumulate(raw["tp"], hours) if tp_acc else np.clip(raw["tp"], 0.0, None)
    sw_h = deaccumulate(raw["ssrd"], hours) if sw_acc else np.clip(raw["ssrd"], 0.0, None)
    Pr = tp_h * 1000.0        # m -> mm over the hour
    Rsw = sw_h / 3600.0       # J/m2 over the hour -> W/m2
    esat, ea = tetens(Ta), tetens(Tdew)

    out = {k: v[keep] for k, v in
           dict(Date=datenum, Ta=Ta, Tdew=Tdew, Pre=Pre, Ws=Ws, Pr=Pr,
                Rsw=Rsw, esat=esat, ea=ea).items()}
    out["Ds"] = np.clip(out["esat"] - out["ea"], 0.0, None)
    out["U"] = out["ea"] / out["esat"]

    if ca[0] is not None:
        out["Ca"] = np.interp(out["Date"], ca[0], ca[1])
    else:
        out["Ca"] = np.full(out["Date"].size, np.nan)

    out.update(Lat=float(lat), Lon=float(lon), Zbas=float(zbas), DeltaGMT=0.0,
               id_location=mat_name(station))

    diag = {
        "hours": int(out["Date"].size),
        "years": f"{yrs[keep].min()}-{yrs[keep].max()}",
        "Ta_C": [round(float(np.nanmin(out['Ta'])), 1), round(float(np.nanmax(out['Ta'])), 1)],
        "Pre_mbar": [round(float(np.nanmin(out['Pre'])), 1), round(float(np.nanmax(out['Pre'])), 1)],
        "Pr_mm_h_max": round(float(np.nanmax(out["Pr"])), 2),
        "Rsw_W_m2_max": round(float(np.nanmax(out["Rsw"])), 1),
        "Ca_ppm": [round(float(np.nanmin(out['Ca'])), 1), round(float(np.nanmax(out['Ca'])), 1)]
        if np.isfinite(out["Ca"]).any() else None,
        "tp_why": tp_why, "ssrd_why": sw_why,
        "ssrd_raw_max": round(float(np.nanmax(raw["ssrd"])), 1),
        "Pr_total_mm_yr": round(float(np.nansum(out["Pr"])) /
                                max(1, (yrs[keep].max() - yrs[keep].min() + 1)), 1),
        "nan_fields": [k for k, v in out.items()
                       if isinstance(v, np.ndarray) and not np.isfinite(v).all()],
    }
    if dry:
        return diag, ""

    from scipy.io import savemat
    out_dir.mkdir(parents=True, exist_ok=True)
    savemat(out_dir / f"Meteo_{mat_name(station)}_raw.mat", out, do_compression=True)
    return diag, ""


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--era5", type=Path, default=DEFAULT_ERA5)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--ca", type=Path, default=DEFAULT_CA)
    p.add_argument("--site-list", type=Path, action="append", default=None)
    p.add_argument("--stations", default="US-HBK,US-Ha2",
                   help="comma-separated, or 'all'")
    p.add_argument("--ameriflux", type=Path, default=INPUT_ROOT / "ameriflux",
                   help="downloaded AmeriFlux dir; Zbas is read from badm_values.csv / "
                        "site_metadata.json / badm_coverage.csv found there")
    p.add_argument("--elevation", type=Path, default=None,
                   help="optional override CSV: station_id,elevation_m")
    p.add_argument("--start-year", type=int, default=1985)
    p.add_argument("--end-year", type=int, default=2020)
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()

    wanted = None if a.stations.strip().lower() == "all" else {
        s.strip() for s in a.stations.split(",") if s.strip()}
    stations = []
    for path in (a.site_list or DEFAULT_SITE_LISTS):
        if not Path(path).is_file():
            continue
        for r in csv.DictReader(open(path, newline="", encoding="utf-8-sig")):
            sid = (r.get("StationID") or "").strip()
            if not sid or (wanted and sid not in wanted):
                continue
            if any(s["station_id"] == sid for s in stations):
                continue
            try:
                stations.append({"station_id": sid, "lat": float(r["Lat"]),
                                 "lon": float(r["Lon"]),
                                 "elev": float(r.get("Elevation") or r.get("elev_m") or "nan")})
            except (KeyError, TypeError, ValueError):
                print(f"  ! {sid}: unusable Lat/Lon", file=sys.stderr)
    elev, elev_src = read_elevation(a.ameriflux, a.elevation)

    ca = read_ca(a.ca)
    print(f"era5     : {a.era5}{'' if a.era5.is_dir() else '   <-- NOT FOUND'}")
    print(f"out      : {a.out}")
    print(f"CO2      : {a.ca}{'' if ca[0] is not None else '   <-- NOT READ, Ca will be NaN'}")
    print(f"ameriflux: {a.ameriflux}{'' if a.ameriflux.is_dir() else '   <-- NOT FOUND'}")
    from collections import Counter as _C
    print(f"Zbas     : {len(elev)} station(s)"
          + (f" ({', '.join(f'{k} {v}' for k, v in _C(elev_src.values()).most_common())})"
             if elev else "   <-- NONE FOUND, no station can build"))
    print(f"stations : {len(stations)}   years {a.start_year}-{a.end_year}\n")

    ok, bad = 0, []
    for st in sorted(stations, key=lambda s: s["station_id"]):
        sid = st["station_id"]
        z = elev.get(sid, st.get("elev"))
        if z is None or not math.isfinite(z):
            bad.append((sid, "no elevation (Zbas) in badm_values.csv, site_metadata.json "
                             "or badm_coverage.csv -- needed for the radiation partition"))
            continue
        diag, err = build(sid, st["lat"], st["lon"], z, a.era5, ca,
                          (a.start_year, a.end_year), a.out, a.dry_run)
        if diag is None:
            bad.append((sid, err))
            continue
        ok += 1
        print(f"  {sid:<8} Zbas {z:>6.0f} m ({elev_src.get(sid, 'site list')})  "
              f"{diag['hours']:>7} h  {diag['years']}  "
              f"Ta {diag['Ta_C'][0]:>6.1f}..{diag['Ta_C'][1]:<5.1f} "
              f"Pre {diag['Pre_mbar'][0]:>6.1f}..{diag['Pre_mbar'][1]:<6.1f} mbar  "
              f"Rswmax {diag['Rsw_W_m2_max']:>6.1f}  P {diag['Pr_total_mm_yr']:>6.1f} mm/yr"
              + (f"   NaN in: {', '.join(diag['nan_fields'])}" if diag["nan_fields"] else ""))
        print(f"           ssrd {diag['ssrd_why']}, raw max {diag['ssrd_raw_max']:g}"
              f"; tp {diag['tp_why']}")

    print(f"\n{ok}/{len(stations)} stations written to {a.out}")
    for sid, why in bad:
        print(f"  ! {sid:<8} {why}")
    if ok:
        print("\nNext: MATLAB finish_meteo.m fills SAB/SAD/PAR/N via "
              "C_Automatic_Radiation_Partition and writes the final Meteo_*.mat")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
