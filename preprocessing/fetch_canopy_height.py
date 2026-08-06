#!/usr/bin/env python3
"""Fetch canopy height (T&C `hc_H`) for the study's AmeriFlux stations.

Source: Potapov et al. (2021), "Mapping global forest canopy height through integration
of GEDI and Landsat data", Remote Sens. Environ. 253:112165 — the 2019 global 30 m
product from UMD GLAD, calibrated on GEDI lidar. Distributed as continental mosaics:

    https://glad.geog.umd.edu/Potapov/Forest_height_2019/Forest_height_2019_NAM.tif

Read remotely with GDAL's /vsicurl/, so the multi-GB mosaic is never downloaded — only
the scanlines covering each station. Use --download to cache it locally instead.

BADM reports HEIGHTC for only 44 of 110 stations, so this fills the other 66. Where BADM
*does* have a value the script reports both and their difference, which both validates
the product and tells you how much to trust it at the sites where BADM is silent.

Pixel encoding (from the GLAD product description):
    1-60   canopy height in metres
    0      no data
    101    water
    102    snow/ice
    103    no data

A single 30 m pixel is noisy relative to a flux-tower footprint, so the default samples a
window around the station and reports the mean, median and the fraction of valid pixels.

Examples:
    python fetch_canopy_height.py --dry-run
    python fetch_canopy_height.py
    python fetch_canopy_height.py --stations US-HBK,US-Ha2 --radius 150
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

try:
    import numpy as np
    import rasterio
    from rasterio.windows import Window
except ImportError:
    sys.exit("rasterio and numpy are required. Run:  pip install -r requirements.txt")

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SITE_LISTS = [
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv",
]
INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA", "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
DEFAULT_OUT = INPUT_ROOT / "canopy_height"

BASE_URL = "https://glad.geog.umd.edu/Potapov/Forest_height_2019"
CITATION = ("Potapov, P. et al. (2021). Mapping global forest canopy height through "
            "integration of GEDI and Landsat data. Remote Sens. Environ. 253, 112165.")

# Non-height sentinel values in the 8-bit raster.
NODATA = {0, 103}
WATER = 101
SNOW_ICE = 102
MAX_HEIGHT = 60


def log(msg: str = "") -> None:
    print(msg, flush=True)


# Stations dropped from the study for lack of data. Kept in one CSV so every step
# excludes the same set and the reason is recorded rather than remembered.
DEFAULT_EXCLUDED = Path(__file__).resolve().parent / "excluded_stations.csv"


def read_excluded(path: Path | None) -> dict[str, str]:
    """station_id -> reason, from the shared exclusion list. Missing file is not an error."""
    out: dict[str, str] = {}
    if not path or not Path(path).is_file():   # '' -> Path('.'), a directory
        return out
    with Path(path).open(newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            sid = (row.get("station_id") or "").strip()
            if sid:
                out[sid] = (row.get("reason") or "").strip()
    return out


def read_stations(paths: list[Path], wanted: set[str] | None) -> list[dict]:
    stations: dict[str, dict] = {}
    for path in paths:
        if not path.exists():
            raise SystemExit(f"site list not found: {path}")
        with path.open(newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                sid = (row.get("StationID") or "").strip()
                if not sid or sid in stations or (wanted is not None and sid not in wanted):
                    continue
                try:
                    lat, lon = float(row["Lat"]), float(row["Lon"])
                except (KeyError, TypeError, ValueError):
                    log(f"  ! {sid}: unusable Lat/Lon, skipping")
                    continue
                stations[sid] = {
                    "station_id": sid, "lat": lat, "lon": lon,
                    "forest_type": (row.get("ForestType") or "").strip(),
                    "igbp": (row.get("IGBP") or "").strip(),
                }
    return sorted(stations.values(), key=lambda s: s["station_id"])


def read_badm_heights(ameriflux_dir: Path) -> dict[str, float]:
    """Measured HEIGHTC per station, for validation. Takes the max reported value."""
    out: dict[str, float] = {}
    values = ameriflux_dir / "badm_values.csv"
    if not values.exists():
        return out
    with values.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row.get("variable") != "HEIGHTC":
                continue
            try:
                h = float(row["value"])
            except (KeyError, ValueError):
                continue
            sid = row["station_id"]
            out[sid] = max(out.get(sid, 0.0), h)
    return out


def sample_window(src, lon: float, lat: float, radius_m: float) -> dict:
    """Statistics over a square window of the given radius, centred on the station."""
    row, col = src.index(lon, lat)

    # Degrees per metre varies with latitude for longitude, not for latitude.
    deg_lat = radius_m / 111_320.0
    deg_lon = radius_m / (111_320.0 * max(math.cos(math.radians(lat)), 0.01))
    px_y = max(int(round(deg_lat / abs(src.transform.e))), 0)
    px_x = max(int(round(deg_lon / abs(src.transform.a))), 0)

    r0, c0 = max(row - px_y, 0), max(col - px_x, 0)
    r1 = min(row + px_y + 1, src.height)
    c1 = min(col + px_x + 1, src.width)
    if r1 <= r0 or c1 <= c0:
        return {"error": "window falls outside the mosaic"}

    data = src.read(1, window=Window(c0, r0, c1 - c0, r1 - r0)).astype(np.int16)
    total = data.size
    centre = int(data[min(row - r0, data.shape[0] - 1), min(col - c0, data.shape[1] - 1)])

    valid = data[(data >= 1) & (data <= MAX_HEIGHT)]
    return {
        "centre_px": centre,
        "n_px": total,
        "n_valid": int(valid.size),
        "pct_valid": round(100.0 * valid.size / total, 1) if total else 0.0,
        "pct_water": round(100.0 * int((data == WATER).sum()) / total, 1) if total else 0.0,
        "pct_nodata": round(100.0 * int(np.isin(data, list(NODATA)).sum()) / total, 1) if total else 0.0,
        "mean_m": round(float(valid.mean()), 2) if valid.size else "",
        "median_m": round(float(np.median(valid)), 2) if valid.size else "",
        "min_m": int(valid.min()) if valid.size else "",
        "max_m": int(valid.max()) if valid.size else "",
    }


def classify(stats: dict, badm: float | None) -> str:
    """A short verdict per station, so problems are visible without reading every column."""
    if "error" in stats:
        return "OUTSIDE_MOSAIC"
    if stats["n_valid"] == 0:
        if stats["pct_water"] > 50:
            return "WATER — no canopy height"
        return "NO_VALID_PIXELS"
    if stats["pct_valid"] < 50:
        return "SPARSE — under half the window is forest"
    if badm is not None:
        diff = abs(stats["mean_m"] - badm)
        if diff > 10:
            return f"DISAGREES with BADM by {diff:.0f} m"
        if diff > 5:
            return f"differs from BADM by {diff:.0f} m"
        return "ok (agrees with BADM)"
    return "ok"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Sample GEDI/Landsat canopy height at the study's AmeriFlux stations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--site-list", type=Path, action="append", default=None,
                   help="station CSV (repeatable); defaults to the deciduous + evergreen lists")
    p.add_argument("--stations", default=None,
                   help="comma-separated StationIDs to restrict the run to")
    p.add_argument("--exclude", default=None,
                   help="comma-separated StationIDs to skip, on top of --exclude-file")
    p.add_argument("--exclude-file", type=Path, default=DEFAULT_EXCLUDED,
                   help="CSV of stations dropped from the study (station_id,reason,...); "
                        "pass an empty string to ignore it")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output directory")
    p.add_argument("--region", default="NAM",
                   help="GLAD continental mosaic: NAM, SAM, AUS, NAFR, SAFR, NASIA, SASIA")
    p.add_argument("--radius", type=float, default=100.0,
                   help="half-width of the sampling window in metres (0 = single pixel)")
    p.add_argument("--ameriflux-dir", type=Path, default=None,
                   help="ameriflux output dir, to compare against measured BADM HEIGHTC")
    p.add_argument("--download", action="store_true",
                   help="download the mosaic first instead of streaming it (multi-GB)")
    p.add_argument("--dry-run", action="store_true",
                   help="list the stations and the source URL, then exit")
    args = p.parse_args(argv)

    wanted = {s.strip() for s in args.stations.split(",") if s.strip()} if args.stations else None
    dropped = read_excluded(args.exclude_file)
    excluded = set(dropped)
    if args.exclude:
        excluded |= {s.strip() for s in args.exclude.split(",") if s.strip()}
    stations = [s for s in read_stations(args.site_list or DEFAULT_SITE_LISTS, wanted)
                if s["station_id"] not in excluded]
    if not stations:
        raise SystemExit("no stations selected")

    url = f"{BASE_URL}/Forest_height_{args.region}.tif".replace(
        f"Forest_height_{args.region}", f"Forest_height_2019_{args.region}")

    log(f"stations   : {len(stations)}" + (f"  ({len(excluded)} excluded)" if excluded else ""))
    log(f"source     : {url}")
    log(f"window     : +/-{args.radius:g} m around each station")
    log(f"output     : {args.out}")
    log("")

    if args.dry_run:
        for s in stations:
            log(f"  [dry-run] {s['station_id']:8} {s['lat']:9.4f} {s['lon']:10.4f}")
        return 0

    badm = read_badm_heights(args.ameriflux_dir) if args.ameriflux_dir else {}
    if badm:
        log(f"comparing against BADM HEIGHTC for {len(badm)} station(s)\n")

    # Keep GDAL from listing the remote directory on every open.
    os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")

    args.out.mkdir(parents=True, exist_ok=True)
    target = url
    if args.download:
        local = args.out / Path(url).name
        if not local.exists():
            log(f"downloading {url} -> {local} (this is multi-GB) ...")
            import urllib.request
            urllib.request.urlretrieve(url, local)
        target = str(local)
    else:
        target = f"/vsicurl/{url}"

    rows = []
    with rasterio.open(target) as src:
        log(f"opened mosaic: {src.width} x {src.height}, {src.crs}\n")
        west, south, east, north = src.bounds
        for s in stations:
            sid, lat, lon = s["station_id"], s["lat"], s["lon"]
            if not (west <= lon <= east and south <= lat <= north):
                stats = {"error": "outside mosaic bounds"}
            else:
                try:
                    stats = sample_window(src, lon, lat, args.radius)
                except Exception as exc:
                    stats = {"error": f"{type(exc).__name__}: {exc}"}

            ref = badm.get(sid)
            verdict = classify(stats, ref)
            row = {
                "station_id": sid, "forest_type": s["forest_type"], "igbp": s["igbp"],
                "lat": lat, "lon": lon,
                "hc_gedi_mean_m": stats.get("mean_m", ""),
                "hc_gedi_median_m": stats.get("median_m", ""),
                "hc_gedi_centre_px_m": stats.get("centre_px", ""),
                "hc_gedi_min_m": stats.get("min_m", ""),
                "hc_gedi_max_m": stats.get("max_m", ""),
                "n_px": stats.get("n_px", ""), "n_valid_px": stats.get("n_valid", ""),
                "pct_valid": stats.get("pct_valid", ""),
                "pct_water": stats.get("pct_water", ""),
                "pct_nodata": stats.get("pct_nodata", ""),
                "badm_heightc_m": ref if ref is not None else "",
                "diff_gedi_minus_badm_m": (round(stats["mean_m"] - ref, 2)
                                           if ref is not None and stats.get("mean_m") != "" else ""),
                "flag": verdict,
                "error": stats.get("error", ""),
            }
            rows.append(row)
            mark = "!" if ("DISAGREE" in verdict or "NO_VALID" in verdict
                           or "WATER" in verdict or "OUTSIDE" in verdict) else " "
            log(f" {mark} {sid:8} hc={str(row['hc_gedi_mean_m']):>6} m"
                f"  valid={str(row['pct_valid']):>5}%"
                + (f"  BADM={ref:g} m  diff={row['diff_gedi_minus_badm_m']:+g}"
                   if ref is not None and row["diff_gedi_minus_badm_m"] != "" else "")
                + f"   {verdict}")

    out_csv = args.out / "canopy_height_gedi.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    meta = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "source_url": url,
        "citation": CITATION,
        "product": "Global Forest Canopy Height 2019, UMD GLAD (GEDI + Landsat), 30 m",
        "pixel_encoding": {"1-60": "canopy height (m)", "0": "no data",
                           "101": "water", "102": "snow/ice", "103": "no data"},
        "sampling_radius_m": args.radius,
        "n_stations": len(rows),
    }
    (args.out / "canopy_height_gedi.json").write_text(
        __import__("json").dumps(meta, indent=2), encoding="utf-8")

    got_value = [r for r in rows if r["hc_gedi_mean_m"] != ""]
    clean = sum(1 for r in rows if r["flag"].startswith("ok"))
    bad = [r for r in rows if r["flag"].startswith(("NO_VALID", "WATER", "OUTSIDE"))]
    disagree = [r for r in rows if "DISAGREE" in r["flag"]]
    compared = [r for r in rows if r["diff_gedi_minus_badm_m"] != ""]

    log(f"\n{len(got_value)}/{len(rows)} stations returned a height "
        f"({clean} of them flagged clean)")
    if compared:
        diffs = [abs(float(r["diff_gedi_minus_badm_m"])) for r in compared]
        log(f"validated against BADM at {len(compared)} station(s): "
            f"mean |difference| = {sum(diffs)/len(diffs):.1f} m, max = {max(diffs):.1f} m")
    if disagree:
        log(f"{len(disagree)} disagree with BADM by >10 m: "
            f"{', '.join(r['station_id'] for r in disagree)}")
        log("  Usually a DATE mismatch, not an error: this product is 2019, while BADM "
            "records a measurement date. Regenerating stands legitimately differ by "
            "10-20 m. Check HEIGHTC_DATE and the disturbance record before choosing.")
    if bad:
        log(f"{len(bad)} with no usable value: {', '.join(r['station_id'] for r in bad)}")
    log(f"\nwritten: {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
