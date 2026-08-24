#!/usr/bin/env python3
"""Label each station-year drought or normal, from SPEI.

    python classify_drought.py --out drought_years.csv
    python classify_drought.py --index SPEI6_ts --months 5,6,7,8,9
    python classify_drought.py --percentile 20 --out drought_years.csv

Writes the station,year,class CSV that analyze_daily_effect.py --drought consumes,
with the SPEI value kept alongside so a composite can be re-cut at a different
threshold without re-reading the stacks.

WHY SPEI AND NOT SPI. SPEI is precipitation minus potential evapotranspiration;
SPI is precipitation alone. A hot year with normal rainfall is a drought for a
vegetation model and invisible to SPI, and hot-dry years are exactly where the
LMA treatment is expected to act -- dynamic LMA sheds leaf area the fixed arm
cannot, so the arms should diverge most when evaporative demand is high.

WHY NOT DERIVE DROUGHT FROM THE MODEL'S OWN PRECIPITATION. The classification has
to be independent of the thing being tested. Both arms share one forcing, so a
model-derived index would be identical between them -- but it would still be
derived from the same series that produced the fluxes, and "the effect is larger
in years the model itself calls dry" is a weaker statement than "in years an
independent index calls dry".

TWO WAYS TO CUT IT

  --threshold   absolute, the conventional -1.0. Comparable across stations, but
                a wet station may contribute no drought years at all and an arid
                one may contribute fifteen, so the composite is unbalanced and
                some stations drop out of it entirely.
  --percentile  each station's own driest N% of years. Guarantees every station
                contributes, at the cost of "drought" meaning something different
                at each -- the driest 20% of Harvard Forest is not the driest 20%
                of a Southwest site.

Neither is right in general. Absolute is the honest default for a fleet-wide
claim; percentile is better for asking whether the effect scales with relative
dryness within a site. Both are written to the same schema, so a figure can be
redone either way.

ANNUAL MEANS ARE SMOOTHER THAN MONTHLY VALUES. P(monthly SPEI < -1) is about 16%
by construction, but averaging twelve months shrinks the variance, so -1.0 on an
annual mean is a rarer and more severe event -- around 9-12% of years at the
stations checked so far. Do not read the conventional monthly thresholds onto
these numbers without that adjustment.

NO SILENT FALLBACKS. A station whose pixel returns all-NaN, or which is missing
from the site lists, is named and skipped, and the exit code is non-zero.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from era5_predictors import (DEFAULT_ERA5_ROOT, SI_ORDER, Era5Monthly)  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
SITE_LISTS = [REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
              REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv"]


def read_sites(wanted=None) -> dict:
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
                out[sid] = (float(r["Lat"]), float(r["Lon"]))
            except (KeyError, TypeError, ValueError):
                print(f"  ! {sid}: unusable Lat/Lon", file=sys.stderr)
    return out


def annual(series: np.ndarray, keys: np.ndarray, months) -> dict:
    """{year: mean SPEI} over the requested months. keys are year*100+month."""
    k = np.asarray(keys, dtype=int)
    yr, mo = k // 100, k % 100
    sel = np.isin(mo, months) if months else np.ones(k.size, bool)
    out = {}
    for y in np.unique(yr[sel]):
        v = series[sel & (yr == y)]
        if np.isfinite(v).any():
            out[int(y)] = float(np.nanmean(v))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--era5-root", type=Path, default=DEFAULT_ERA5_ROOT)
    ap.add_argument("--out", type=Path, required=True, help="CSV to write")
    ap.add_argument("--index", default="SPEI12_ts",
                    help=f"one of {', '.join(SI_ORDER)} (default SPEI12_ts)")
    ap.add_argument("--months", default=None,
                    help="comma-separated months to average, e.g. 5,6,7,8,9 for "
                         "the growing season. Default: all twelve.")
    ap.add_argument("--threshold", type=float, default=-1.0,
                    help="annual mean below this is a drought year")
    ap.add_argument("--percentile", type=float, default=None,
                    help="instead of --threshold, take each station's own driest "
                         "N%% of years")
    ap.add_argument("--years", default=None,
                    help="restrict to YYYY-YYYY, e.g. 1985-2020 to match the "
                         "ERA5 runs")
    ap.add_argument("--stations", default=None, help="comma-separated subset")
    a = ap.parse_args(argv)

    if a.index not in SI_ORDER:
        print(f"ERROR: --index must be one of {', '.join(SI_ORDER)}", file=sys.stderr)
        return 1
    months = ([int(m) for m in a.months.split(",")] if a.months else None)
    y0, y1 = (None, None)
    if a.years:
        try:
            y0, y1 = (int(x) for x in a.years.split("-"))
        except ValueError:
            print("ERROR: --years must look like 1985-2020", file=sys.stderr)
            return 1

    sites = read_sites(set(s.strip() for s in a.stations.split(",")) if a.stations
                       else None)
    if not sites:
        print("ERROR: no stations", file=sys.stderr)
        return 1
    try:
        store = Era5Monthly(a.era5_root)
    except RuntimeError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    cut = (f"driest {a.percentile:.0f}% per station" if a.percentile is not None
           else f"annual mean < {a.threshold}")
    print(f"index     : {a.index}"
          f"{'   months ' + a.months if months else '   all months'}")
    print(f"criterion : {cut}")
    print(f"years     : {a.years or 'all'}\nstations  : {len(sites)}\n")

    rows, skipped = [], []
    for sid in sorted(sites):
        lat, lon = sites[sid]
        try:
            ser = store.pixel_series(lat, lon)
        except Exception as e:                                   # noqa: BLE001
            skipped.append((sid, f"{type(e).__name__}: {e}")); continue
        s = np.asarray(ser["si"][a.index], dtype=float)
        if not np.isfinite(s).any():
            skipped.append((sid, f"{a.index} is all-NaN at this pixel")); continue

        yv = annual(s, ser["si_time"][a.index], months)
        if y0 is not None:
            yv = {y: v for y, v in yv.items() if y0 <= y <= y1}
        if not yv:
            skipped.append((sid, "no years left after filtering")); continue

        if a.percentile is not None:
            lim = float(np.percentile(list(yv.values()), a.percentile))
        else:
            lim = a.threshold
        for y in sorted(yv):
            rows.append((sid, y, "drought" if yv[y] <= lim else "normal",
                         round(yv[y], 4), round(lim, 4)))

    if skipped:
        print(f"SKIPPED -- {len(skipped)} station(s):")
        for sid, why in skipped:
            print(f"  ! {sid:<10}{why}")
        print()
    if not rows:
        print("ERROR: nothing classified", file=sys.stderr)
        return 1

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["station", "year", "class", "spei", "cut"])
        w.writerows(rows)

    n_st = len({r[0] for r in rows})
    dry = [r for r in rows if r[2] == "drought"]
    per = {}
    for r in dry:
        per[r[0]] = per.get(r[0], 0) + 1
    none = n_st - len(per)
    print(f"{len(rows)} station-years over {n_st} stations; "
          f"{len(dry)} drought ({100*len(dry)/len(rows):.0f}%)")
    if per:
        c = sorted(per.values())
        print(f"drought years per station: min {c[0]}  median {c[len(c)//2]}  "
              f"max {c[-1]}"
              + (f"   -- {none} station(s) with NONE" if none else ""))
    if none and a.percentile is None:
        print("  A station contributing no drought years drops out of the "
              "composite entirely.\n  --percentile 20 guarantees every station "
              "contributes, at the cost of\n  'drought' meaning something "
              "different at each.")
    print(f"\n-> {a.out}")
    return 1 if skipped else 0


if __name__ == "__main__":
    sys.exit(main())
