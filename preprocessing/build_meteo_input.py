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

ACCUMULATED FIELDS. ERA5-Land tp and ssrd accumulate from 00 UTC and reset each
day, so the value at hour H is the total since midnight, not the flux over the
preceding hour. They are de-accumulated by within-day differencing, with the 01
UTC value used as-is.
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
DEFAULT_CA = REPO_ROOT / "T&C" / "Diogo" / "Ca_Data.mat"
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
    dt = next((np.asarray(d[k], dtype=float).ravel()
               for k in ("Date", "Date_Ca", "t") if k in d), None)
    return (dt, ca) if ca is not None else (None, None)


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
    Pr = deaccumulate(raw["tp"], hours) * 1000.0  # m -> mm per hour
    Rsw = deaccumulate(raw["ssrd"], hours) / 3600.0   # J/m2 -> W/m2
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
    p.add_argument("--elevation", type=Path,
                   default=INPUT_ROOT / "ameriflux" / "site_elevation.csv",
                   help="station_id,elevation_m; falls back to the site list")
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
    elev = {}
    if a.elevation.is_file():
        for r in csv.DictReader(open(a.elevation, newline="", encoding="utf-8-sig")):
            try:
                elev[(r.get("station_id") or r.get("StationID") or "").strip()] = \
                    float(r.get("elevation_m") or r.get("elev_m") or r.get("Elevation"))
            except (TypeError, ValueError):
                pass

    ca = read_ca(a.ca)
    print(f"era5     : {a.era5}{'' if a.era5.is_dir() else '   <-- NOT FOUND'}")
    print(f"out      : {a.out}")
    print(f"CO2      : {a.ca}{'' if ca[0] is not None else '   <-- NOT READ, Ca will be NaN'}")
    print(f"stations : {len(stations)}   years {a.start_year}-{a.end_year}\n")

    ok, bad = 0, []
    for st in sorted(stations, key=lambda s: s["station_id"]):
        sid = st["station_id"]
        z = elev.get(sid, st.get("elev"))
        if z is None or not math.isfinite(z):
            bad.append((sid, "no elevation (Zbas) -- needed for the radiation partition"))
            continue
        diag, err = build(sid, st["lat"], st["lon"], z, a.era5, ca,
                          (a.start_year, a.end_year), a.out, a.dry_run)
        if diag is None:
            bad.append((sid, err))
            continue
        ok += 1
        print(f"  {sid:<8} {diag['hours']:>7} h  {diag['years']}  "
              f"Ta {diag['Ta_C'][0]:>6.1f}..{diag['Ta_C'][1]:<5.1f} "
              f"Pre {diag['Pre_mbar'][0]:>6.1f}..{diag['Pre_mbar'][1]:<6.1f} mbar  "
              f"Rswmax {diag['Rsw_W_m2_max']:>6.1f}"
              + (f"   NaN in: {', '.join(diag['nan_fields'])}" if diag["nan_fields"] else ""))

    print(f"\n{ok}/{len(stations)} stations written to {a.out}")
    for sid, why in bad:
        print(f"  ! {sid:<8} {why}")
    if ok:
        print("\nNext: MATLAB finish_meteo.m fills SAB/SAD/PAR/N via "
              "C_Automatic_Radiation_Partition and writes the final Meteo_*.mat")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
