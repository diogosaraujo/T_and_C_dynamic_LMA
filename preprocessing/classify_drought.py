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
from gcm_variables import GCMS, SCENARIOS, var_dir                     # noqa: E402
from results_dir import NoResultsDir, resolve_out                      # noqa: E402


def gcm_spei_file(gcm: str, scenario: str, months: int) -> Path | None:
    """The SPEI stack for one model/scenario, found by globbing.

    The variant (r1i1p1f1 vs r1i1p1f2) and grid label (gn vs gr) are not constant
    across models, so the filename is matched rather than built -- the same
    reasoning as gcm_variables.find_year_files. Unlike the daily variables there
    is ONE file per model/scenario, not one per year.
    """
    d = var_dir(gcm, scenario, "spei")
    if not d.is_dir():
        return None
    hits = sorted(d.glob(f"spei_{months}_{gcm}_{scenario}_*.nc"))
    return hits[0] if hits else None


def gcm_series(path: Path, sites: dict) -> dict:
    """{station: {year: annual mean SPEI}} for every site, from ONE open file.

    Opened once and looped over stations rather than the reverse: the array is
    (time, 600, 1440) float32, so re-opening per station would re-read chunks
    116 times over.

    Longitude is stored 0..360 here, against the -180..180 the site lists use.
    Latitude and longitude are 2-D but constant along the other axis, so the
    first row and column are the postings.
    """
    from netCDF4 import Dataset, num2date
    out = {}
    with Dataset(str(path)) as ds:
        t = ds.variables["time"]
        dates = num2date(t[:], t.units, t.calendar)
        yr = np.array([d.year for d in dates])
        mo = np.array([d.month for d in dates])
        la = np.asarray(ds.variables["latitude"][:, 0], dtype=float)
        lo = np.asarray(ds.variables["longitude"][0, :], dtype=float)
        v = ds.variables["spei"]
        for sid, (slat, slon) in sites.items():
            i = int(np.abs(la - slat).argmin())
            j = int(np.abs(lo - (slon % 360.0)).argmin())
            s = np.asarray(v[:, i, j], dtype=float)
            out[sid] = (s, yr * 100 + mo)
    return out

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
    """{year: (mean SPEI, n months used)}. keys are year*100+month.

    The COUNT is returned, not just the mean, because an N-month SPEI has no
    value for the first N-1 months and the published stacks are TRIMMED rather
    than padded: spei_12_..._1980-2014.nc holds 409 months, 420 - 11, and starts
    1980-12-16. So the first year of any stack is one month long.

    Historically that is harmless -- the 1985-2014 window drops 1980 anyway --
    but the SSP stacks run 2015-2100 and the SSP window IS 2015-2100, so 2015
    would otherwise be labelled from a single December while every other year
    used twelve. The caller enforces a minimum and reports what it dropped.
    """
    k = np.asarray(keys, dtype=int)
    yr, mo = k // 100, k % 100
    sel = np.isin(mo, months) if months else np.ones(k.size, bool)
    out = {}
    for y in np.unique(yr[sel]):
        v = series[sel & (yr == y)]
        n = int(np.isfinite(v).sum())
        if n:
            out[int(y)] = (float(np.nanmean(v)), n)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", choices=["era5", "gcm"], default="era5",
                    help="era5: SPEI from the ERA5-Land monthly stacks -- the "
                         "OBSERVED climate, valid only for the era5_land runs. "
                         "gcm: SPEI from each model's own NEX-GDDP output, which "
                         "is what the GCM runs need, because a GCM does not "
                         "reproduce actual weather and the observed dry years are "
                         "not its dry years.")
    ap.add_argument("--gcms", default=None, help="gcm source: comma-separated subset")
    ap.add_argument("--scenarios", default=None,
                    help="gcm source: comma-separated subset of "
                         + ", ".join(SCENARIOS))
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
    cut = (f"driest {a.percentile:.0f}% per station" if a.percentile is not None
           else f"annual mean < {a.threshold}")
    print(f"index     : {a.index}"
          f"{'   months ' + a.months if months else '   all months'}")
    print(f"criterion : {cut}")
    print(f"years     : {a.years or 'all'}\nstations  : {len(sites)}\n")

    # (station, gcm, scenario) -> {year: annual SPEI}. One entry per run family,
    # because a GCM's dry years are its own: applying the observed labels to a
    # GCM run would mislabel most years, since a GCM does not reproduce actual
    # weather. US-NR1's observed driest are 2012/2013/2002/2006; ACCESS-CM2's are
    # 2013/1984/2009.
    rows, skipped, dropped, series = [], [], [], {}

    if a.source == "era5":
        try:
            store = Era5Monthly(a.era5_root)
        except RuntimeError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 1
        for sid in sorted(sites):
            lat, lon = sites[sid]
            try:
                ser = store.pixel_series(lat, lon)
            except Exception as e:                               # noqa: BLE001
                skipped.append((f"{sid}", f"{type(e).__name__}: {e}")); continue
            series[(sid, "", "era5_land")] = (
                np.asarray(ser["si"][a.index], dtype=float),
                np.asarray(ser["si_time"][a.index], dtype=int))
    else:
        acc = int("".join(c for c in a.index if c.isdigit()) or 12)
        gcms = [g.strip() for g in a.gcms.split(",")] if a.gcms else list(GCMS)
        scens = ([x.strip() for x in a.scenarios.split(",")] if a.scenarios
                 else list(SCENARIOS))
        for g in gcms:
            for sc in scens:
                f = gcm_spei_file(g, sc, acc)
                if f is None:
                    skipped.append((f"{g}/{sc}", f"no spei_{acc}_{g}_{sc}_*.nc"))
                    continue
                print(f"  {g:<15}{sc:<12}{f.name}", flush=True)
                for sid, (sv, keys) in gcm_series(f, sites).items():
                    series[(sid, g, sc)] = (sv, keys)

    for (sid, g, sc), (sv, keys) in sorted(series.items()):
        who = f"{sid} {g} {sc}".strip()
        if not np.isfinite(sv).any():
            skipped.append((who, f"{a.index} is all-NaN at this pixel")); continue
        yc = annual(sv, keys, months)
        if y0 is not None:
            yc = {y: v for y, v in yc.items() if y0 <= y <= y1}
        # Require the full set of requested months. A year built from fewer is
        # not comparable to one built from twelve, and the threshold and the
        # percentile cut are both taken over these values -- a short year would
        # move the cut for every other year at the station.
        need = len(months) if months else 12
        short = {y: n for y, (_, n) in yc.items() if n < need}
        yv = {y: m for y, (m, n) in yc.items() if n >= need}
        for y, n in sorted(short.items()):
            dropped.append((who, y, n, need))
        if not yv:
            skipped.append((who, "no years left after filtering")); continue
        lim = (float(np.percentile(list(yv.values()), a.percentile))
               if a.percentile is not None else a.threshold)
        for y in sorted(yv):
            rows.append((sid, g, sc, y, "drought" if yv[y] <= lim else "normal",
                         round(yv[y], 4), round(lim, 4)))

    if dropped:
        # Not an error: the first year of every stack is short by construction,
        # because an N-month SPEI has no value for the first N-1 months. Said
        # out loud so a missing year is never a mystery later.
        yrs = sorted({y for _, y, _, _ in dropped})
        print(f"SHORT YEARS DROPPED -- {len(dropped)} station-year(s), "
              f"year(s) {', '.join(str(y) for y in yrs[:6])}"
              f"{' ...' if len(yrs) > 6 else ''}")
        for who, y, n, need in dropped[:6]:
            print(f"  - {who:<28}  {y}: {n} of {need} months")
        if len(dropped) > 6:
            print(f"  ... and {len(dropped) - 6} more")
        print()
    if skipped:
        print(f"SKIPPED -- {len(skipped)} station(s):")
        for who, why in skipped:
            # who is "<sid> <gcm> <scenario>", up to ~30 chars -- a :<10 field
            # left "US-HB2  era5_land" butted straight against the reason. The
            # trailing two spaces keep them apart however long the name runs.
            print(f"  ! {who:<28}  {why}")
        print()
    if not rows:
        print("ERROR: nothing classified", file=sys.stderr)
        return 1

    # Relative names go to $TC_RESULTS alongside the daily and annual tables.
    try:
        a.out = resolve_out(a.out)
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["station", "gcm", "scenario", "year", "class", "spei", "cut"])
        w.writerows(rows)

    n_st = len({r[0] for r in rows})
    dry = [r for r in rows if r[4] == "drought"]
    per = {}
    for r in dry:
        per[(r[0], r[1], r[2])] = per.get((r[0], r[1], r[2]), 0) + 1
    none = len({(r[0], r[1], r[2]) for r in rows}) - len(per)
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
