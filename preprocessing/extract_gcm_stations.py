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
        values     float32 (ndays, nstation) in the file's native units
        ymd        int32   (ndays, 3)         year, month, day IN THE MODEL'S OWN
                                              calendar -- not datetime64, because
                                              the models do not share one
        doy        int32   (ndays,)           day-of-year, needed to place a
                                              360-day calendar on the real one
        calendar   str                        as declared by the file
        source_var str                        which variable was actually read
                                              (huss, or hurs where a model has no
                                              huss -- IPSL-CM6A-LR)
        stations   <U8     (nstation,)
        lat lon    float32 (nstation,)        the GRID CELL centre actually used
        offset_km  float32 (nstation,)        station-to-cell-centre distance
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
                           HUMIDITY_PREFERENCE, find_year_files,
                           expected_years, tasks)

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

    # Humidity: hurs for every model, so the route is identical across the
    # ensemble. huss is kept only as a fallback for a model that ships no hurs.
    # Whichever was read is recorded in the npz and build_gcm_meteo.py converts
    # accordingly, so a substitution can never be silent.
    src_var = var
    alts = HUMIDITY_PREFERENCE if var in HUMIDITY_PREFERENCE else [var]
    files, want = {}, expected_years(scenario)
    for cand in alts:
        files = find_year_files(gcm, scenario, cand)
        if files:
            src_var = cand
            break
    have = [y for y in want if y in files]
    if not have:
        return "no files", (f"{NEXGDDP_ROOT/gcm/scenario/var} empty or unreadable"
                            + (f" (also tried {', '.join(alts[1:])})"
                               if len(alts) > 1 else ""))
    missing = [y for y in want if y not in files]

    vals, ymd, doys = [], [], []
    iy = ix = clat = clon = None
    calendar = "standard"
    for y in have:
        with netCDF4.Dataset(files[y]) as nc:
            if iy is None:
                iy, ix, clat, clon = locate(nc, stations)
            v = nc.variables[src_var]
            # ONE contiguous read of the bounding box holding every station, then
            # index inside it. The obvious version -- 101 strided v[:, j, i] point
            # reads per file -- walks the whole compressed chunk structure once
            # per station and measured 3.3-4 h per 35-year task (job 36896), which
            # put the 86-year SSP tasks over the wall clock. CONUS occupies about
            # 3% of a global 0.25-degree field, so the box costs little memory and
            # turns 101 reads into 1.
            wide = (ix.max() - ix.min()) > v.shape[2] // 2      # wraps the meridian
            if wide:
                arr = np.empty((v.shape[0], len(stations)), dtype=np.float32)
                for k, (j, i) in enumerate(zip(iy, ix)):
                    arr[:, k] = np.asarray(v[:, j, i], dtype=np.float32).ravel()
            else:
                y0, y1 = int(iy.min()), int(iy.max()) + 1
                x0, x1 = int(ix.min()), int(ix.max()) + 1
                box = np.asarray(v[:, y0:y1, x0:x1], dtype=np.float32)
                arr = box[:, iy - y0, ix - x0]
            vals.append(arr)

            t = nc.variables["time"]
            calendar = getattr(t, "calendar", "standard")
            import netCDF4 as _n
            # cftime objects ALWAYS work; Python datetimes do not. GFDL-ESM4 uses a
            # 365-day calendar and every one of its 21 tasks died on
            # "illegal calendar or reference date for python datetime". Store the
            # date COMPONENTS rather than a datetime64 so any calendar survives,
            # and let build_gcm_meteo.py decide how to place them on a real axis.
            d = _n.num2date(t[:], t.units, calendar,
                            only_use_cftime_datetimes=True)
            ymd.append(np.array([(x.year, x.month, x.day) for x in d], dtype=np.int32))
            doys.append(np.array([getattr(x, "dayofyr", 0) for x in d], dtype=np.int32))

    values = np.concatenate(vals, axis=0)
    YMD = np.concatenate(ymd, axis=0)
    DOY = np.concatenate(doys)
    order = np.lexsort((YMD[:, 2], YMD[:, 1], YMD[:, 0]))
    values, YMD, DOY = values[order], YMD[order], DOY[order]

    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out, values=values, ymd=YMD, doy=DOY,
        calendar=np.array(calendar), source_var=np.array(src_var),
        stations=np.array([s["station"] for s in stations]),
        lat=clat.astype(np.float32), lon=clon.astype(np.float32),
        offset_km=haversine_km([s["lat"] for s in stations],
                               [s["lon"] for s in stations],
                               clat, clon).astype(np.float32),
        units=np.array(VARIABLES[var]["units"]),
        years_missing=np.array(missing, dtype=np.int64))
    notes = []
    if src_var != var:
        notes.append(f"USED {src_var} (no {var} for this model)")
    if calendar not in ("standard", "gregorian", "proleptic_gregorian"):
        notes.append(f"calendar '{calendar}'")
    if missing:
        notes.append(f"{len(missing)} year-file(s) missing: {missing[0]}..{missing[-1]}"
                     if len(missing) > 2 else f"missing {missing}")
    return (f"{values.shape[0]} days x {values.shape[1]} stations",
            "; ".join(notes) or None)


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
