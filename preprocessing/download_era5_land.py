#!/usr/bin/env python3
"""Download hourly ERA5-Land forcing for the AmeriFlux stations of the LMA -> T&C study.

Reads the station lists produced by the ecoregion pairing step
(T&C/dynamic_lma_test/{deciduous,evergreen}_ameriflux.csv), and pulls hourly
ERA5-Land data at each station's coordinates for 1985-2021.

Two modes:

  timeseries (default)
      One request per station covering the whole period, using the CDS point
      time-series collection `reanalysis-era5-land-timeseries`. Output is a single
      <StationID>_ERA5_Land.nc per station.

  gridded
      Fallback using the classic gridded `reanalysis-era5-land` collection, one
      request per station-year over a small box around the site. Output is
      <StationID>/<StationID>_ERA5_Land_<year>.nc. Far more requests, but the
      time-series collection is documented as possibly being disabled or
      deprecated at any point, so keep this path working.

ERA5-Land snaps requests to its 0.1 degree grid, so stations closer together than
one grid cell resolve to identical data. Those are downloaded once and copied, which
for the current lists turns 118 stations into 80 requests.

Data is stored in native ERA5-Land netCDF, units unconverted. Each station also gets
a <StationID>_ERA5_Land.json sidecar describing the variables, their units, the time
convention, and the T&C forcing field each one feeds.

Setup:
    pip install -r requirements.txt
    # then put your CDS key in ~/.cdsapirc (see README.md)

Examples:
    python download_era5_land.py --dry-run
    python download_era5_land.py
    python download_era5_land.py --stations US-Ho1,US-MMS --jobs 2
    python download_era5_land.py --mode gridded --start-year 1985 --end-year 1986
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import sys
import threading
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

from era5_variables import (
    ACCUMULATION_NOTE,
    CDS_VARIABLE_NAMES,
    COORDINATE_NOTES,
    DELIBERATELY_OMITTED,
    GRIDDED_DATASET,
    TIMESERIES_DATASET,
    VARIABLES,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SITE_LISTS = [
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv",
]
DEFAULT_OUT = Path(__file__).resolve().parent / "data" / "era5_land"

GRID_STEP = 0.1  # ERA5-Land resolution, degrees
# A station sitting this close to a cell boundary is not safe to group with others,
# because our rounding and the CDS server's may disagree about which cell it lands in.
BOUNDARY_EPS = 1e-6

MONTHS = [f"{m:02d}" for m in range(1, 13)]
DAYS = [f"{d:02d}" for d in range(1, 32)]
HOURS = [f"{h:02d}:00" for h in range(24)]

_print_lock = threading.Lock()


def log(msg: str) -> None:
    with _print_lock:
        print(msg, flush=True)


# --------------------------------------------------------------------------------------
# Station list handling
# --------------------------------------------------------------------------------------


def snap_to_grid(value: float) -> float:
    """Round a coordinate to the nearest ERA5-Land 0.1 degree grid line."""
    return round(math.floor(value / GRID_STEP + 0.5) * GRID_STEP, 1)


def near_cell_boundary(value: float) -> bool:
    """True if the coordinate sits essentially exactly halfway between two grid lines."""
    frac = abs(value / GRID_STEP - math.floor(value / GRID_STEP))
    return abs(frac - 0.5) < BOUNDARY_EPS


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def read_stations(paths: list[Path], wanted: set[str] | None) -> list[dict]:
    """Read the pairing CSVs into one de-duplicated station list."""
    stations: dict[str, dict] = {}
    for path in paths:
        if not path.exists():
            raise SystemExit(f"site list not found: {path}")
        with path.open(newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                sid = (row.get("StationID") or "").strip()
                if not sid:
                    continue
                if wanted is not None and sid not in wanted:
                    continue
                try:
                    lat = float(row["Lat"])
                    lon = float(row["Lon"])
                except (KeyError, TypeError, ValueError):
                    log(f"  ! skipping {sid}: missing or unparseable Lat/Lon")
                    continue
                if sid in stations:
                    # The same tower can serve more than one ecoregion row.
                    stations[sid]["source_rows"] += 1
                    continue
                stations[sid] = {
                    "station_id": sid,
                    "station_name": (row.get("StationName") or "").strip(),
                    "lat": lat,
                    "lon": lon,
                    "forest_type": (row.get("ForestType") or "").strip(),
                    "igbp": (row.get("IGBP") or "").strip(),
                    "us_l3code": (row.get("US_L3CODE") or "").strip(),
                    "us_l3name": (row.get("US_L3NAME") or "").strip(),
                    "source_csv": path.name,
                    "source_rows": 1,
                }
    return sorted(stations.values(), key=lambda s: s["station_id"])


def group_by_grid_cell(stations: list[dict], dedup: bool) -> list[dict]:
    """Group stations that resolve to the same ERA5-Land grid cell.

    Returns one group per download. The first station in each group is the one whose
    file is actually fetched; the rest are copies.
    """
    groups: dict[str, dict] = {}
    for st in stations:
        glat = snap_to_grid(st["lat"])
        glon = snap_to_grid(st["lon"])
        st["grid_lat"] = glat
        st["grid_lon"] = glon
        st["grid_offset_km"] = round(haversine_km(st["lat"], st["lon"], glat, glon), 3)

        ambiguous = near_cell_boundary(st["lat"]) or near_cell_boundary(st["lon"])
        if not dedup or ambiguous:
            # Give it a private key so it is never merged with another station.
            key = f"solo:{st['station_id']}"
            if ambiguous:
                log(
                    f"  ! {st['station_id']} sits on a grid-cell boundary; "
                    "downloading it separately rather than sharing a cell"
                )
        else:
            key = f"{glat:.1f}_{glon:.1f}"

        groups.setdefault(key, {"grid_lat": glat, "grid_lon": glon, "members": []})
        groups[key]["members"].append(st)

    return [groups[k] for k in sorted(groups)]


# --------------------------------------------------------------------------------------
# Metadata sidecar
# --------------------------------------------------------------------------------------


def build_metadata(
    station: dict, group: dict, mode: str, start_year: int, end_year: int, files: list[str]
) -> dict:
    shared_with = [
        m["station_id"] for m in group["members"] if m["station_id"] != station["station_id"]
    ]
    dataset = TIMESERIES_DATASET if mode == "timeseries" else GRIDDED_DATASET
    return {
        "station_id": station["station_id"],
        "station_name": station["station_name"],
        "forest_type": station["forest_type"],
        "igbp": station["igbp"],
        "ecoregion": {
            "us_l3code": station["us_l3code"],
            "us_l3name": station["us_l3name"],
        },
        "coordinates": {
            "station_lat": station["lat"],
            "station_lon": station["lon"],
            "era5_land_grid_lat": station["grid_lat"],
            "era5_land_grid_lon": station["grid_lon"],
            "grid_offset_km": station["grid_offset_km"],
            "notes": (
                "The request is issued at the ERA5-Land grid point, not the tower "
                "coordinate; ERA5-Land resolves to a 0.1 degree grid (~11 km). Keep the "
                "true station coordinates for solar geometry and for site elevation."
            ),
        },
        "source": {
            "cds_dataset": dataset,
            "product": "ERA5-Land reanalysis",
            "provider": "Copernicus Climate Data Store (C3S/ECMWF)",
            "licence": "Licence to use Copernicus Products (CC-BY 4.0)",
            "citation": (
                "Munoz-Sabater, J. et al. (2021): ERA5-Land: a state-of-the-art global "
                "reanalysis dataset for land applications. Earth Syst. Sci. Data, 13, "
                "4349-4383."
            ),
            "format": "netCDF (native, unconverted)",
        },
        "period": {
            "start": f"{start_year}-01-01",
            "end": f"{end_year}-12-31",
            "temporal_resolution": "hourly",
        },
        "time_reference": COORDINATE_NOTES["time"],
        "accumulation_note": ACCUMULATION_NOTE,
        "variables": VARIABLES,
        "variables_deliberately_omitted": DELIBERATELY_OMITTED,
        "files": files,
        "provenance": {
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "script": Path(__file__).name,
            "mode": mode,
            "site_list": station["source_csv"],
            "downloaded_once_and_shared_with": shared_with,
        },
    }


# --------------------------------------------------------------------------------------
# Download
# --------------------------------------------------------------------------------------


def normalise_download(tmp: Path, final: Path) -> None:
    """Move a completed download into place, unwrapping a zip if CDS returned one."""
    with tmp.open("rb") as fh:
        magic = fh.read(4)

    if magic[:2] == b"PK":
        with zipfile.ZipFile(tmp) as zf:
            members = [n for n in zf.namelist() if n.lower().endswith(".nc")]
            if len(members) != 1:
                raise RuntimeError(
                    f"expected exactly one .nc inside the returned archive, got {members}"
                )
            with zf.open(members[0]) as src, final.open("wb") as dst:
                shutil.copyfileobj(src, dst)
        tmp.unlink()
        return

    if magic not in (b"CDF\x01", b"CDF\x02", b"\x89HDF"):
        raise RuntimeError(
            f"downloaded file is not netCDF (leading bytes {magic!r}); "
            "the CDS may have returned an error document"
        )
    tmp.replace(final)


def retrieve(client, dataset: str, request: dict, target: Path, retries: int) -> None:
    """Run one CDS retrieval with retries, writing atomically via a .part file."""
    tmp = target.with_suffix(target.suffix + ".part")
    if tmp.exists():
        tmp.unlink()

    last: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            client.retrieve(dataset, request, str(tmp))
            normalise_download(tmp, target)
            return
        except Exception as exc:  # cdsapi raises a variety of transport/server errors
            last = exc
            if tmp.exists():
                tmp.unlink()
            if attempt < retries:
                backoff = min(300, 20 * 2 ** (attempt - 1))
                log(f"  ! attempt {attempt}/{retries} failed ({exc}); retrying in {backoff}s")
                time.sleep(backoff)
    raise RuntimeError(f"all {retries} attempts failed; last error: {last}")


def timeseries_request(group: dict, start_year: int, end_year: int) -> dict:
    return {
        "variable": CDS_VARIABLE_NAMES,
        "location": {"latitude": group["grid_lat"], "longitude": group["grid_lon"]},
        "date": [f"{start_year}-01-01/{end_year}-12-31"],
        "data_format": "netcdf",
    }


def gridded_request(group: dict, year: int) -> dict:
    lat, lon = group["grid_lat"], group["grid_lon"]
    half = GRID_STEP / 2
    return {
        "variable": CDS_VARIABLE_NAMES,
        "year": str(year),
        "month": MONTHS,
        "day": DAYS,
        "time": HOURS,
        # [North, West, South, East] -- a single-cell box around the grid point.
        "area": [
            round(lat + half, 3),
            round(lon - half, 3),
            round(lat - half, 3),
            round(lon + half, 3),
        ],
        "data_format": "netcdf",
        "download_format": "unarchived",
    }


def station_files(station_id: str, out_dir: Path, args) -> list[Path]:
    """Output paths for one station, in the same order for every station in a cell."""
    if args.mode == "timeseries":
        return [out_dir / f"{station_id}_ERA5_Land.nc"]
    return [
        out_dir / station_id / f"{station_id}_ERA5_Land_{yr}.nc"
        for yr in range(args.start_year, args.end_year + 1)
    ]


def process_group(
    group: dict, out_dir: Path, args, client_factory
) -> list[dict]:
    """Download one grid cell and write files + sidecars for every station in it."""
    primary = group["members"][0]
    sid = primary["station_id"]
    results: list[dict] = []

    primary_files = station_files(sid, out_dir, args)
    if args.mode == "timeseries":
        requests = [timeseries_request(group, args.start_year, args.end_year)]
    else:
        requests = [gridded_request(group, yr) for yr in range(args.start_year, args.end_year + 1)]
        if not args.dry_run:
            primary_files[0].parent.mkdir(parents=True, exist_ok=True)
    targets = list(zip(primary_files, requests))

    pending = [(t, r) for t, r in targets if args.overwrite or not (t.exists() and t.stat().st_size > 0)]
    skipped = len(targets) - len(pending)
    if skipped:
        log(f"  = {sid}: {skipped} file(s) already present, skipping")

    if args.dry_run:
        for target, request in pending:
            log(f"  [dry-run] {target.name}  <-  {json.dumps(request, sort_keys=True)}")
    elif pending:
        client = client_factory()
        dataset = TIMESERIES_DATASET if args.mode == "timeseries" else GRIDDED_DATASET
        for target, request in pending:
            log(f"  > {sid}: requesting {target.name}")
            retrieve(client, dataset, request, target, args.retries)
            log(f"  + {sid}: wrote {target.name} ({target.stat().st_size / 1e6:.1f} MB)")

    # Every station in the cell gets its own file set, copied from the primary.
    for station in group["members"]:
        own_files = station_files(station["station_id"], out_dir, args)
        if station["station_id"] != sid and not args.dry_run:
            for src, dst in zip(primary_files, own_files):
                dst.parent.mkdir(parents=True, exist_ok=True)
                if src.exists() and (args.overwrite or not dst.exists()):
                    shutil.copy2(src, dst)

        if not args.dry_run:
            meta = build_metadata(
                station, group, args.mode, args.start_year, args.end_year,
                [f.name for f in own_files],
            )
            sidecar = own_files[0].parent / f"{station['station_id']}_ERA5_Land.json"
            sidecar.write_text(json.dumps(meta, indent=2), encoding="utf-8")

        results.append(
            {
                "station_id": station["station_id"],
                "station_name": station["station_name"],
                "forest_type": station["forest_type"],
                "lat": station["lat"],
                "lon": station["lon"],
                "grid_lat": station["grid_lat"],
                "grid_lon": station["grid_lon"],
                "grid_offset_km": station["grid_offset_km"],
                "n_files": len(own_files),
                "downloaded_from": sid,
                "shared_cell": station["station_id"] != sid,
            }
        )
    return results


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def parse_args(argv: list[str] | None = None):
    p = argparse.ArgumentParser(
        description="Download hourly ERA5-Land forcing for the study's AmeriFlux stations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--site-list", type=Path, action="append", default=None,
                   help="station CSV (repeatable); defaults to the deciduous + evergreen lists")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output directory")
    p.add_argument("--start-year", type=int, default=1985)
    p.add_argument("--end-year", type=int, default=2021)
    p.add_argument("--mode", choices=("timeseries", "gridded"), default="timeseries",
                   help="point time-series collection, or gridded fallback")
    p.add_argument("--stations", default=None,
                   help="comma-separated StationIDs to restrict the run to")
    p.add_argument("--jobs", type=int, default=4, help="concurrent CDS requests")
    p.add_argument("--shard", type=int, default=0,
                   help="0-based index of this shard, for SLURM job arrays")
    p.add_argument("--num-shards", type=int, default=1,
                   help="total number of shards; grid cells are split round-robin across them")
    p.add_argument("--retries", type=int, default=4, help="attempts per request")
    p.add_argument("--no-dedup", action="store_true",
                   help="download every station separately even when they share a grid cell")
    p.add_argument("--overwrite", action="store_true", help="re-download existing files")
    p.add_argument("--dry-run", action="store_true",
                   help="print the request plan without contacting the CDS")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if args.end_year < args.start_year:
        raise SystemExit("--end-year must be >= --start-year")

    site_lists = args.site_list or DEFAULT_SITE_LISTS
    wanted = {s.strip() for s in args.stations.split(",") if s.strip()} if args.stations else None

    stations = read_stations(site_lists, wanted)
    if not stations:
        raise SystemExit("no stations selected")
    if wanted:
        missing = wanted - {s["station_id"] for s in stations}
        if missing:
            log(f"  ! requested station(s) not found in the site lists: {sorted(missing)}")

    groups = group_by_grid_cell(stations, dedup=not args.no_dedup)
    all_cells = len(groups)

    if args.num_shards < 1 or not 0 <= args.shard < args.num_shards:
        raise SystemExit("--shard must satisfy 0 <= shard < --num-shards")
    if args.num_shards > 1:
        # Shard whole grid cells, never individual stations, so stations sharing a cell
        # stay in the same task and the copy step always finds its source file.
        groups = [g for i, g in enumerate(groups) if i % args.num_shards == args.shard]
        if not groups:
            log(f"shard {args.shard}/{args.num_shards}: no grid cells assigned, nothing to do")
            return 0

    per_group = 1 if args.mode == "timeseries" else (args.end_year - args.start_year + 1)
    log(f"stations           : {len(stations)}")
    if args.num_shards > 1:
        log(f"shard              : {args.shard} of {args.num_shards}")
        log(f"unique grid cells  : {len(groups)} of {all_cells} (this shard)")
    else:
        log(f"unique grid cells  : {len(groups)}")
    log(f"mode               : {args.mode} ({TIMESERIES_DATASET if args.mode == 'timeseries' else GRIDDED_DATASET})")
    log(f"period             : {args.start_year}-01-01 .. {args.end_year}-12-31 (hourly)")
    log(f"variables          : {', '.join(CDS_VARIABLE_NAMES)}")
    log(f"CDS requests       : {len(groups) * per_group}")
    log(f"output             : {args.out}")
    log("")

    if not args.dry_run:
        args.out.mkdir(parents=True, exist_ok=True)

    client_factory = None
    if not args.dry_run:
        try:
            import cdsapi
        except ImportError:
            raise SystemExit(
                "cdsapi is not installed. Run:  pip install -r requirements.txt"
            )
        # One client per worker thread; the client is not documented as thread-safe.
        local = threading.local()

        def client_factory():  # noqa: F811
            if not hasattr(local, "client"):
                local.client = cdsapi.Client()
            return local.client

    rows: list[dict] = []
    failures: list[tuple[str, str]] = []

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {
            pool.submit(process_group, g, args.out, args, client_factory): g for g in groups
        }
        for fut in as_completed(futures):
            group = futures[fut]
            sid = group["members"][0]["station_id"]
            try:
                rows.extend(fut.result())
            except Exception as exc:
                log(f"  X {sid}: FAILED -- {exc}")
                failures.append((sid, str(exc)))

    if not args.dry_run and rows:
        # Array tasks each write their own manifest, otherwise they would clobber
        # each other. Concatenate them afterwards (see slurm/README.md).
        name = "manifest.csv" if args.num_shards == 1 else f"manifest.shard{args.shard:03d}.csv"
        manifest = args.out / name
        with manifest.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(sorted(rows, key=lambda r: r["station_id"]))
        log(f"\nmanifest: {manifest}")

    expected = sum(len(g["members"]) for g in groups)
    log(f"\nstations completed: {len(rows)}/{expected}")
    if failures:
        log(f"failed grid cells : {len(failures)}")
        for sid, err in failures:
            log(f"  - {sid}: {err}")
        log("Re-run the same command to retry; completed files are skipped.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
