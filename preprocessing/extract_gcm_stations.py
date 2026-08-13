#!/usr/bin/env python3
"""Pull the station pixels out of the NEX-GDDP-CMIP6 store.

    python extract_gcm_stations.py --index 7          # array form, one task
    python extract_gcm_stations.py --gcm GFDL-ESM4    # one model, all scen/vars
    python extract_gcm_stations.py --all              # everything, serially
    python extract_gcm_stations.py --report           # what is on disk

One task = one (GCM, scenario, variable). It opens each yearly global file ONCE
and pulls all ~92 station pixels from it, then writes a single compact .npz of
daily values, shape (ndays, nstations).

That ordering is the whole point. The store is global daily 0.25-degree, so a
yearly file is O(100 MB) and we want O(100 bytes) of it. Looping stations on the
outside would reopen every file 92 times and turn a 7,000-file job into a
640,000-file one; on a shared filesystem that is the difference between hours and
days.

    5 GCMs x 3 scenarios x 7 variables            = 105 array tasks
    5 x 7 x (35 hist + 2 x 86 future) year-files  = 7,245 file reads total

GRID. NEX-GDDP is a regular lat/lon 0.25-degree grid with longitude on 0..360, so
CONUS station longitudes (negative) are wrapped before matching. Nearest-neighbour
selection, and the realised offset is recorded per station -- at 0.25 degrees the
worst case is ~14 km at this latitude, which is smaller than the ~9 km ERA5-Land
pixel only by a little and matters for the same reason: in complex terrain the
grid cell is not the tower.

Output, one file per task:
    $TC_INPUT_DATA/gcm_stations/<GCM>/<scenario>/<var>.npz
        values   float32 (ndays, nstation)   in the file's native units
        dates    int64   (ndays,)            days since 1850-01-01, proleptic
        stations <U8     (nstation,)
        lat lon  float32 (nstation,)         the GRID CELL centre actually used
        offset_km float32 (nstation,)        station-to-cell-centre distance
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gcm_variables import (GCMS, SCENARIOS, VARIABLES, NEXGDDP_ROOT,   # noqa: E402
                           find_year_files, expected_years, tasks)

REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA",
                                 "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
OUT_ROOT = INPUT_ROOT / "gcm_stations"
SITE_LISTS = [REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
              REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv"]
EXCLUDED = REPO_ROOT / "preprocessing" / "excluded_stations.csv"
EPOCH = np.datetime64("1850-01-01", "D")


def read_stations(exclude: Path | None = EXCLUDED) -> list[dict]:
    """Study stations with coordinates, minus the excluded list, sorted by ID."""
    drop = set()
    if exclude and exclude.is_file():
        for r in csv.DictReader(open(exclude, newline="", encoding="utf-8-sig")):
            sid = (r.get("station_id") or "").strip()
            if sid:
                drop.add(sid)
    out, seen = [], set()
    for p in SITE_LISTS:
        if not p.is_file():
            print(f"  ! site list missing: {p}", file=sys.stderr)
            continue
        for r in csv.DictReader(open(p, newline="", encoding="utf-8-sig")):
            sid = (r.get("StationID") or "").strip()
            if not sid or sid in seen or sid in drop:
                continue
            try:
                lat, lon = float(r["Lat"]), float(r["Lon"])
            except (KeyError, TypeError, ValueError):
                continue
            seen.add(sid)
            out.append(dict(station=sid, lat=lat, lon=lon,
                            forest_type=(r.get("ForestType") or "").strip().lower()))
    return sorted(out, key=lambda x: x["station"])


def locate(nc, stations):
    """Nearest grid indices for every station. Returns (iy, ix, cell_lat, cell_lon)."""
    lat = np.asarray(nc["lat"][:], dtype=float)
    lon = np.asarray(nc["lon"][:], dtype=float)
    wrap = lon.min() >= 0.0                      # NEX-GDDP uses 0..360
    iy, ix = [], []
    for s in stations:
        slon = s["lon"] % 360.0 if wrap else s["lon"]
        iy.append(int(np.argmin(np.abs(lat - s["lat"]))))
        ix.append(int(np.argmin(np.abs(lon - slon))))
    iy, ix = np.array(iy), np.array(ix)
    return iy, ix, lat[iy], lon[ix]


def haversine_km(lat1, lon1, lat2, lon2):
    lon1 = ((np.asarray(lon1) + 180) % 360) - 180
    lon2 = ((np.asarray(lon2) + 180) % 360) - 180
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def extract(gcm, scenario, var, stations, out_root, force=False):
    out = out_root / gcm / scenario / f"{var}.npz"
    if out.is_file() and not force:
        return "cached", None
    try:
        import netCDF4
    except ImportError:
        raise SystemExit("netCDF4 is required -- add it to requirements.txt")

    files = find_year_files(gcm, scenario, var)
    want = expected_years(scenario)
    have = [y for y in want if y in files]
    if not have:
        return "no files", f"{NEXGDDP_ROOT/gcm/scenario/var} empty or unreadable"
    missing = [y for y in want if y not in files]

    vals, dates = [], []
    iy = ix = clat = clon = None
    for y in have:
        with netCDF4.Dataset(files[y]) as nc:
            if iy is None:
                iy, ix, clat, clon = locate(nc, stations)
            v = nc.variables[var]
            # One read of the whole (time, lat, lon) block would pull the entire
            # global field into memory; slice per station index pair instead. The
            # file is open once either way, which is what actually costs.
            arr = np.empty((v.shape[0], len(stations)), dtype=np.float32)
            for k, (j, i) in enumerate(zip(iy, ix)):
                arr[:, k] = np.asarray(v[:, j, i], dtype=np.float32).ravel()
            vals.append(arr)
            t = nc.variables["time"]
            import netCDF4 as _n
            d = _n.num2date(t[:], t.units, getattr(t, "calendar", "standard"),
                            only_use_cftime_datetimes=False,
                            only_use_python_datetimes=True)
            dates.append(np.array([np.datetime64(x.date(), "D") for x in d]))

    values = np.concatenate(vals, axis=0)
    dd = np.concatenate(dates)
    order = np.argsort(dd)
    values, dd = values[order], dd[order]

    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out, values=values, dates=(dd - EPOCH).astype(np.int64),
        stations=np.array([s["station"] for s in stations]),
        lat=clat.astype(np.float32), lon=clon.astype(np.float32),
        offset_km=haversine_km([s["lat"] for s in stations],
                               [s["lon"] for s in stations],
                               clat, clon).astype(np.float32),
        units=np.array(VARIABLES[var]["units"]),
        years_missing=np.array(missing, dtype=np.int64))
    note = None
    if missing:
        note = (f"{len(missing)} year-file(s) missing: "
                f"{missing[0]}..{missing[-1]}" if len(missing) > 2 else str(missing))
    return f"{values.shape[0]} days x {values.shape[1]} stations", note


def report(out_root, stations):
    print(f"{'GCM':16s}{'scenario':12s}" + "".join(f"{v:>9s}" for v in VARIABLES))
    nmiss = 0
    for g in GCMS:
        for s in SCENARIOS:
            row = f"{g:16s}{s:12s}"
            for v in VARIABLES:
                p = out_root / g / s / f"{v}.npz"
                if p.is_file():
                    row += f"{'ok':>9s}"
                else:
                    row += f"{'-':>9s}"; nmiss += 1
            print(row)
    print(f"\n  {len(GCMS)*len(SCENARIOS)*len(VARIABLES)-nmiss} of "
          f"{len(GCMS)*len(SCENARIOS)*len(VARIABLES)} extracted, {nmiss} missing")
    # grid offsets, from whichever file exists
    for g in GCMS:
        for s in SCENARIOS:
            for v in VARIABLES:
                p = out_root / g / s / f"{v}.npz"
                if p.is_file():
                    d = np.load(p, allow_pickle=False)
                    o = d["offset_km"]
                    print(f"\n  grid offset ({g}): mean {o.mean():.1f} km, "
                          f"max {o.max():.1f} km ({d['stations'][o.argmax()]})")
                    return
    return


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", type=int, help="1-based index into the task list")
    ap.add_argument("--gcm", action="append")
    ap.add_argument("--scenario", action="append")
    ap.add_argument("--variable", action="append")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--out", type=Path, default=OUT_ROOT)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args(argv)

    stations = read_stations()
    print(f"nexgddp  : {NEXGDDP_ROOT}"
          f"{'' if NEXGDDP_ROOT.is_dir() else '   <-- NOT FOUND'}")
    print(f"output   : {a.out}")
    print(f"stations : {len(stations)}\n")

    if a.report:
        report(a.out, stations); return 0
    if not NEXGDDP_ROOT.is_dir():
        print(f"ERROR: {NEXGDDP_ROOT} not found", file=sys.stderr); return 1

    work = tasks(a.gcm, a.scenario, a.variable)
    if a.index is not None:
        if a.index < 1 or a.index > len(work):
            print(f"index {a.index} outside 1..{len(work)} -- nothing to do"); return 0
        work = [work[a.index - 1]]
    elif not (a.all or a.gcm or a.scenario or a.variable):
        ap.error("give --index / --gcm / --scenario / --variable / --all, or --report")

    rc = 0
    for i, (g, s, v) in enumerate(work, 1):
        try:
            st, note = extract(g, s, v, stations, a.out, a.force)
        except Exception as e:                                   # noqa: BLE001
            print(f"  [{i}/{len(work)}] {g} {s} {v}: FAILED -- "
                  f"{type(e).__name__}: {e}", flush=True)
            rc = 1
            continue
        print(f"  [{i}/{len(work)}] {g:14s} {s:11s} {v:8s} {st}"
              f"{'  ! ' + note if note else ''}", flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
