#!/usr/bin/env python3
"""Fetch rooting depth (T&C `ZR95_H`) for the study's AmeriFlux stations.

Source: ISLSCP II Ecosystem Rooting Depths (Schenk & Jackson), ORNL DAAC,
doi:10.3334/ORNLDAAC/929. Two global 1-degree ASCII grids:

    95ecosys_rootdepth_1d.asc   depth containing 95% of roots (m)  -> ZR95_H
    50ecosys_rootdepth_1d.asc   depth containing 50% of roots (m)  -> ZR50_H

D95 maps directly onto T&C's ZR95_H; only a metres -> millimetres conversion is needed.

READ THIS BEFORE TRUSTING THE OUTPUT
------------------------------------
The grid is 1 degree, roughly 110 km. For this station network that resolves only
54 distinct cells across 118 stations: the 15 CHEESEHEAD towers all share one cell, as do
8 Wisconsin and 7 Metolius sites. So this is much closer to a regional lookup than to a
site measurement, and the output reports how many stations share each cell so that is
visible rather than implied. Fan et al. (2017) is ~1 km if a finer product is obtained
later.

Also: T&C requires ZR95_H <= the deepest soil layer, or Root_Fraction_General aborts the
run. The soil-depth step may end up constraining rooting depth more than this product
does, so the two have to be reconciled when the .mat files are built.

Credentials: the ORNL DAAC needs a free NASA Earthdata Login (https://urs.earthdata.nasa.gov).
Supply it by any of:
    ~/.netrc  with:  machine urs.earthdata.nasa.gov login <user> password <pass>
    EARTHDATA_USER / EARTHDATA_PASS environment variables
    EARTHDATA_TOKEN environment variable (bearer token)

Examples:
    python fetch_root_depth.py --dry-run
    python fetch_root_depth.py
    python fetch_root_depth.py --stations US-HBK,US-Ha2
"""

from __future__ import annotations

import argparse
import csv
import http.cookiejar
import json
import math
import netrc
import os
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SITE_LISTS = [
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv",
]
INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA", "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
DEFAULT_OUT = INPUT_ROOT / "root_depth"

EARTHDATA_HOST = "urs.earthdata.nasa.gov"
# The documentation names the files but not their directory; both layouts return 401
# (i.e. they exist behind auth), so try each.
BASE_CANDIDATES = [
    "https://daac.ornl.gov/daacdata/islscp_ii/vegetation/ecosystem_roots_1deg/comp",
    "https://daac.ornl.gov/daacdata/islscp_ii/vegetation/ecosystem_roots_1deg/data",
]
FILES = {"d95": "95ecosys_rootdepth_1d.asc", "d50": "50ecosys_rootdepth_1d.asc"}
CITATION = ("Schenk, H.J. and R.B. Jackson (2009). ISLSCP II Ecosystem Rooting Depths. "
            "ORNL DAAC, Oak Ridge, Tennessee, USA. doi:10.3334/ORNLDAAC/929. "
            "Underlying analysis: Schenk & Jackson (2002), Ecol. Monogr. 72:311-328.")


def log(msg: str = "") -> None:
    print(msg, flush=True)


# --------------------------------------------------------------------------------------
# Earthdata-authenticated download
# --------------------------------------------------------------------------------------


def build_opener() -> tuple[urllib.request.OpenerDirector, dict]:
    """Opener that can follow the Earthdata OAuth redirect chain."""
    headers = {"User-Agent": "T_and_C_dynamic_LMA/1.0 (research)"}
    token = os.environ.get("EARTHDATA_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
        return urllib.request.build_opener(
            urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar())), headers

    user = os.environ.get("EARTHDATA_USER")
    password = os.environ.get("EARTHDATA_PASS")
    if not (user and password):
        try:
            auth = netrc.netrc().authenticators(EARTHDATA_HOST)
            if auth:
                user, _, password = auth
        except (FileNotFoundError, netrc.NetrcParseError):
            pass
    if not (user and password):
        raise SystemExit(
            "No Earthdata credentials found. Register free at "
            f"https://{EARTHDATA_HOST} then either\n"
            f"  add to ~/.netrc:  machine {EARTHDATA_HOST} login <user> password <pass>\n"
            "  or export EARTHDATA_USER / EARTHDATA_PASS\n"
            "  or export EARTHDATA_TOKEN")

    mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    mgr.add_password(None, f"https://{EARTHDATA_HOST}", user, password)
    opener = urllib.request.build_opener(
        urllib.request.HTTPBasicAuthHandler(mgr),
        urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
    return opener, headers


def download(opener, headers: dict, name: str, dest: Path) -> Path:
    """Fetch one grid, trying each candidate directory."""
    if dest.exists() and dest.stat().st_size > 0:
        log(f"  = {name} already cached")
        return dest
    last = None
    for base in BASE_CANDIDATES:
        url = f"{base}/{name}"
        try:
            req = urllib.request.Request(url, headers=headers)
            with opener.open(req, timeout=180) as resp:
                data = resp.read()
            if data[:200].lstrip().lower().startswith(b"<!doctype") or b"<html" in data[:200].lower():
                last = f"{url} returned an HTML page (login likely failed)"
                continue
            dest.write_bytes(data)
            log(f"  + {name} ({len(data)/1024:.0f} KB) from {base}")
            return dest
        except urllib.error.HTTPError as exc:
            last = f"{url} -> HTTP {exc.code} {exc.reason}"
        except Exception as exc:
            last = f"{url} -> {type(exc).__name__}: {exc}"
    raise SystemExit(f"could not download {name}. Last attempt: {last}")


# --------------------------------------------------------------------------------------
# ASCII grid
# --------------------------------------------------------------------------------------


class Grid:
    """A global lat/lon ASCII grid, with or without an ArcInfo-style header."""

    def __init__(self, path: Path):
        tokens = path.read_text(encoding="utf-8", errors="replace").split()
        hdr = {}
        i = 0
        while i + 1 < len(tokens) and tokens[i].lower() in (
                "ncols", "nrows", "xllcorner", "yllcorner", "xllcenter", "yllcenter",
                "cellsize", "nodata_value"):
            hdr[tokens[i].lower()] = float(tokens[i + 1])
            i += 2

        values = [float(t) for t in tokens[i:]]
        if hdr:
            self.ncols, self.nrows = int(hdr["ncols"]), int(hdr["nrows"])
            self.cell = hdr.get("cellsize", 1.0)
            self.west = hdr.get("xllcorner", hdr.get("xllcenter", -180.0) - self.cell / 2)
            south = hdr.get("yllcorner", hdr.get("yllcenter", -90.0) - self.cell / 2)
            self.north = south + self.nrows * self.cell
            self.nodata = hdr.get("nodata_value", -99.0)
            self.header = True
        else:
            # ISLSCP II convention: headerless global grid, row 0 = northernmost.
            self.ncols, self.nrows, self.cell = 360, 180, 1.0
            self.west, self.north, self.nodata = -180.0, 90.0, -99.0
            self.header = False
            if len(values) != 360 * 180:
                raise SystemExit(
                    f"{path.name}: expected 64800 values for a global 1-degree grid, "
                    f"got {len(values)}. The layout differs from the assumed ISLSCP II "
                    "convention -- inspect the file before trusting any sample.")
        self.values = values

    def sample(self, lat: float, lon: float):
        col = int((lon - self.west) / self.cell)
        row = int((self.north - lat) / self.cell)
        if not (0 <= col < self.ncols and 0 <= row < self.nrows):
            return None
        v = self.values[row * self.ncols + col]
        return None if v <= self.nodata or v < 0 else v

    def describe(self) -> str:
        finite = [v for v in self.values if v > self.nodata and v >= 0]
        return (f"{self.ncols}x{self.nrows} @ {self.cell} deg, "
                f"header={self.header}, nodata={self.nodata:g}, "
                f"{len(finite)} valid cells, range {min(finite):.2f}-{max(finite):.2f} m")


# --------------------------------------------------------------------------------------


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
                    stations[sid] = {"station_id": sid, "lat": float(row["Lat"]),
                                     "lon": float(row["Lon"]),
                                     "forest_type": (row.get("ForestType") or "").strip(),
                                     "igbp": (row.get("IGBP") or "").strip()}
                except (KeyError, TypeError, ValueError):
                    log(f"  ! {sid}: unusable Lat/Lon, skipping")
    return sorted(stations.values(), key=lambda s: s["station_id"])


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Sample Schenk & Jackson rooting depth at the study's AmeriFlux stations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--site-list", type=Path, action="append", default=None)
    p.add_argument("--stations", default=None, help="comma-separated StationIDs")
    p.add_argument("--exclude", default=None,
                   help="comma-separated StationIDs to skip, on top of --exclude-file")
    p.add_argument("--exclude-file", type=Path, default=DEFAULT_EXCLUDED,
                   help="CSV of stations dropped from the study; '' to ignore")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--dry-run", action="store_true")
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

    log(f"stations : {len(stations)}" + (f"  ({len(excluded)} excluded)" if excluded else ""))
    log(f"source   : ISLSCP II Ecosystem Rooting Depths (Schenk & Jackson), 1 degree")
    log(f"output   : {args.out}")
    log("")
    if args.dry_run:
        for base in BASE_CANDIDATES:
            log(f"  [dry-run] would try {base}/{FILES['d95']}")
        for s in stations:
            log(f"  [dry-run] {s['station_id']:8} {s['lat']:9.4f} {s['lon']:10.4f}")
        return 0

    args.out.mkdir(parents=True, exist_ok=True)
    opener, headers = build_opener()
    grids = {}
    for key, name in FILES.items():
        path = download(opener, headers, name, args.out / name)
        grids[key] = Grid(path)
        log(f"    {name}: {grids[key].describe()}")
    log("")

    # How many stations land in each 1-degree cell -- the resolution caveat, made explicit.
    cell_of = {s["station_id"]: (math.floor(s["lat"]), math.floor(s["lon"])) for s in stations}
    per_cell: dict[tuple, list[str]] = {}
    for sid, cell in cell_of.items():
        per_cell.setdefault(cell, []).append(sid)

    rows = []
    for s in stations:
        sid = s["station_id"]
        d95 = grids["d95"].sample(s["lat"], s["lon"])
        d50 = grids["d50"].sample(s["lat"], s["lon"])
        share = sorted(x for x in per_cell[cell_of[sid]] if x != sid)
        rows.append({
            "station_id": sid, "forest_type": s["forest_type"], "igbp": s["igbp"],
            "lat": s["lat"], "lon": s["lon"],
            "ZR95_H_mm": round(d95 * 1000) if d95 is not None else "",
            "ZR50_H_mm": round(d50 * 1000) if d50 is not None else "",
            "d95_m": d95 if d95 is not None else "",
            "d50_m": d50 if d50 is not None else "",
            "grid_cell_1deg": f"{cell_of[sid][0]},{cell_of[sid][1]}",
            "n_stations_sharing_cell": len(per_cell[cell_of[sid]]),
            "stations_sharing_cell": ";".join(share),
            "flag": "ok" if d95 is not None else "NO DATA in the 1-degree cell",
        })
        log(f"  {sid:8} ZR95={str(rows[-1]['ZR95_H_mm']):>6} mm"
            f"  ZR50={str(rows[-1]['ZR50_H_mm']):>6} mm"
            f"  cell shared by {len(per_cell[cell_of[sid]])}"
            + ("" if d95 is not None else "   <- NO DATA"))

    out_csv = args.out / "root_depth_schenk_jackson.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    (args.out / "root_depth_schenk_jackson.json").write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "citation": CITATION,
        "product": "ISLSCP II Ecosystem Rooting Depths, 1 degree global",
        "variables": {"ZR95_H_mm": "depth containing 95% of roots, T&C ZR95_H",
                      "ZR50_H_mm": "depth containing 50% of roots, T&C ZR50_H"},
        "resolution_caveat": (
            f"1 degree (~110 km): {len(per_cell)} distinct cells for {len(rows)} stations. "
            "Closer to a regional lookup than a site measurement."),
        "constraint": "T&C requires ZR95_H <= the deepest Zs layer or the run aborts.",
        "n_stations": len(rows),
    }, indent=2), encoding="utf-8")

    got = [r for r in rows if r["ZR95_H_mm"] != ""]
    log(f"\n{len(got)}/{len(rows)} stations got a rooting depth")
    if got:
        vals = [r["ZR95_H_mm"] for r in got]
        log(f"ZR95 range {min(vals)}-{max(vals)} mm (median {sorted(vals)[len(vals)//2]} mm)")
    log(f"resolution: {len(per_cell)} distinct 1-degree cells for {len(rows)} stations")
    missing = [r["station_id"] for r in rows if r["ZR95_H_mm"] == ""]
    if missing:
        log(f"no value at: {', '.join(missing)}")
    log(f"\nwritten: {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
