#!/usr/bin/env python3
"""What the SPEI/SPI stacks actually contain at our stations, before trusting them.

    python inspect_drought.py --stations US-Ha2,US-NR1,US-Blo
    python inspect_drought.py --stations US-Ha2 --index SPEI12_ts --show-years

Reads through era5_predictors.Era5Monthly rather than opening the .mat files
directly. That class already encodes the parts that are easy to get wrong and
were derived from the MATLAB source: which stack is indexed [lat, lon] on the
shifted -180..180 grid and which is [lon, lat] on the native 0..360 one, the
1800-column longitude split, and the fact that h5py hands back MATLAB's
dimensions reversed. Re-deriving any of that here would just be a second chance
to get it wrong.

WHAT TO LOOK AT, AND WHY

  * time span      SPEI12 starts 1980-12, SPEI3 in 1980-03 -- each index has its
                   own start because an N-month index needs N months of history
                   before its first value. If the printed span does not match
                   DROUGHT_FILES, the registry and the files disagree.
  * NaN fraction   a station over water or outside the stack's domain comes back
                   all-NaN, and that is the failure mode most likely to be
                   mistaken for "no droughts here".
  * value range    SPEI is a standardised index, so |values| should mostly sit
                   under 3 and the mean should be near 0. A mean far from 0 means
                   the pixel lookup is landing somewhere unintended.
  * known droughts a US site should show 2012 (Midwest), 2002 and 2012 (Colorado
                   Front Range), 2002/2007/2012 (Southwest). Seeing them is the
                   real check that the pixel is the right one; the statistics
                   above only prove a number was read.

This writes nothing. It is a look before anything classifies years as drought.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from era5_predictors import (DEFAULT_ERA5_ROOT, DROUGHT_FILES,        # noqa: E402
                             SI_ORDER, Era5Monthly)

REPO_ROOT = Path(__file__).resolve().parents[1]
SITE_LISTS = [REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
              REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv"]


def read_sites(wanted=None):
    out = {}
    for p in SITE_LISTS:
        if not p.is_file():
            print(f"  ! site list not found: {p}", file=sys.stderr)
            continue
        for r in csv.DictReader(open(p, newline="", encoding="utf-8-sig")):
            sid = (r.get("StationID") or "").strip()
            if not sid or (wanted and sid not in wanted):
                continue
            try:
                out[sid] = (float(r["Lat"]), float(r["Lon"]),
                            (r.get("ForestType") or "").strip())
            except (KeyError, TypeError, ValueError):
                print(f"  ! {sid}: unusable Lat/Lon", file=sys.stderr)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--era5-root", type=Path, default=DEFAULT_ERA5_ROOT)
    ap.add_argument("--stations", default=None,
                    help="comma-separated; default is a spread of five sites")
    ap.add_argument("--index", default="SPEI12_ts",
                    help=f"which index to detail; one of {', '.join(SI_ORDER)}")
    ap.add_argument("--show-years", action="store_true",
                    help="print the annual mean of --index for every year")
    ap.add_argument("--threshold", type=float, default=-1.0,
                    help="annual mean below this counts as a drought year")
    a = ap.parse_args(argv)

    if a.index not in DROUGHT_FILES:
        print(f"ERROR: --index must be one of {', '.join(SI_ORDER)}", file=sys.stderr)
        return 1

    print("REGISTRY (era5_predictors.DROUGHT_FILES)")
    for k in SI_ORDER:
        fn, var, order, start = DROUGHT_FILES[k]
        path = a.era5_root / fn
        print(f"  {k:<11}{fn:<22}var={var:<6}{order:<8}start={start[0]}-{start[1]:02d}"
              f"   {'OK' if path.is_file() else 'MISSING'}")
    print()

    wanted = ([s.strip() for s in a.stations.split(",")] if a.stations
              else ["US-Ha2", "US-NR1", "US-Blo", "US-SRM", "US-Me2"])
    sites = read_sites(set(wanted))
    if not sites:
        print("ERROR: none of the requested stations are in the site lists",
              file=sys.stderr)
        return 1

    try:
        store = Era5Monthly(a.era5_root)
    except RuntimeError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    bad = 0
    for sid in sorted(sites):
        lat, lon, pft = sites[sid]
        print(f"{'=' * 70}\n{sid}  ({pft})  lat {lat:.4f}  lon {lon:.4f}")
        try:
            ser = store.pixel_series(lat, lon)
        except Exception as e:                                   # noqa: BLE001
            print(f"  FAILED: {type(e).__name__}: {e}"); bad += 1; continue

        clat, clon = store.cell_center(lat, lon)
        print(f"  nearest posting: lat {clat:.4f}  lon {clon:.4f}"
              f"   (indices lat={ser['lat_index']}, lon={ser['lon_index']})")

        print(f"  {'index':<11}{'n':>6}{'span':>18}{'NaN%':>7}"
              f"{'mean':>8}{'min':>8}{'max':>8}")
        for k in SI_ORDER:
            s = np.asarray(ser["si"][k], dtype=float)
            t = ser["si_time"][k]
            fin = np.isfinite(s)
            span = f"{t[0]:%Y-%m}..{t[-1]:%Y-%m}" if len(t) else "-"
            print(f"  {k:<11}{s.size:>6}{span:>18}"
                  f"{100*(1-fin.mean()):>6.1f}%"
                  f"{np.nanmean(s):>8.2f}{np.nanmin(s):>8.2f}{np.nanmax(s):>8.2f}")
            if not fin.any():
                print(f"    ! all-NaN -- this pixel carries no {k}")
                bad += 1

        # Annual means of the chosen index, and which years would be flagged.
        s = np.asarray(ser["si"][a.index], dtype=float)
        t = ser["si_time"][a.index]
        yrs = np.array([d.year for d in t])
        rows = []
        for y in range(int(yrs.min()), int(yrs.max()) + 1):
            v = s[yrs == y]
            if np.isfinite(v).any():
                rows.append((y, float(np.nanmean(v)), int(np.isfinite(v).sum())))
        dry = [r for r in rows if r[1] < a.threshold]
        print(f"\n  {a.index}: {len(rows)} years, "
              f"{len(dry)} below {a.threshold} "
              f"({100*len(dry)/max(1,len(rows)):.0f}%)")
        print("  driest: " + ", ".join(f"{y} ({v:+.2f})"
                                       for y, v, _ in sorted(rows, key=lambda r: r[1])[:6]))
        if a.show_years:
            for y, v, n in rows:
                mark = "  <-- drought" if v < a.threshold else ""
                print(f"    {y}  {v:+6.2f}  ({n:2d} months){mark}")
        print()

    print(f"{'=' * 70}")
    print("Sanity: 2012 should be dry across the Midwest and Front Range; 2002 and\n"
          "2012 in Colorado; 2002/2007/2012 in the Southwest. If the driest years\n"
          "look nothing like that, the pixel lookup is landing in the wrong place\n"
          "and no amount of thresholding will fix it.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
