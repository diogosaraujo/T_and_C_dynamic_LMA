"""Fixed-vs-dynamic LMA per MONTH or SEASON, with the year kept.

WHY THIS EXISTS ALONGSIDE analyze_daily_effect.py. That script writes a day-of-
year CLIMATOLOGY: every row is a multi-year mean and there is no year column, so
it answers "when in the year does the treatment bite" and cannot answer "how did
2012 differ from 2011". Drought enters it only as separate climatologies per
class. This script keeps the year, so the output is a TIME SERIES and the
interannual questions -- trends, a specific drought year, a composite around an
event -- are all group-bys on it.

Daily resolution with the year kept was rejected on size: it would take the ERA5
table from 1.4M rows to ~50M, and roughly 15x that across the GCMs. A month is
short enough to resolve the seasonal cycle and the onset of a drought, and it
divides the file by ~30.

    station,key,year,period,variable,fixed,dyn,diff,rel_pct,n_days

'period' is 1-12 for --freq monthly, or DJF/MAM/JJA/SON for --freq seasonal.

NO DROUGHT ARGUMENT, DELIBERATELY. The year is in the output, so drought labels
join on (station, gcm, scenario, year) after the fact. Baking the classes in
would triple the rows to store a label that is one join away, and it would pin
the table to one threshold -- the fixed cut and the percentile cut could not
both be applied to the same file.

SEASONS SPAN THE YEAR BOUNDARY. DJF is assigned to the year of its January, so
December 2001 belongs to DJF 2002 -- the meteorological convention. The first
December and the last Jan-Feb of a record therefore fall in seasons that are one
month short; n_days shows it, and those seasons are the ones to drop before
plotting a trend.

RATIOS. Tfrac, Bowen and WUE are the mean of the DAILY ratios over the period,
with any day whose denominator flux has collapsed left out of that mean -- the
same sample-level rule analyze_daily_effect.py uses, and for the same reason: one
snowed-out day with QE near zero injects H/QE = 800 and drags the whole month
with it. n_days reports how many days actually contributed.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_daily_effect import (REPORT, RATIOS, Unusable,      # noqa: E402
                                  derive, read_run)
from check_treatment_effect import find_pairs                    # noqa: E402
from results_dir import NoResultsDir, resolve_out                # noqa: E402

# How a period value is formed from the daily series. Water and carbon fluxes
# accumulate; state and energy are averaged, because a monthly mean W/m2 is the
# quantity anyone plots. Matches analyze_daily_effect's daily aggregation.
PERIOD_SUM = {"GPP", "NPP", "ET", "T", "EG", "EIn", "Lk"}
PERIOD_MEAN = {"LAI_H", "QE", "H", "Rn"}

SEASON_OF = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
             6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}


def period_keys(year: np.ndarray, mo: np.ndarray, freq: str):
    """(period label per day, year per day) for the requested frequency."""
    if freq == "monthly":
        return mo.astype(int), year.astype(int)
    # DJF belongs to the year of its January: December rolls forward.
    lab = np.array([SEASON_OF[int(m)] for m in mo])
    yy = np.where(mo == 12, year + 1, year).astype(int)
    return lab, yy


def periods(fx: dict, dy: dict, freq: str):
    """Rows of (year, period, variable, fixed, dyn, diff, rel_pct, n)."""
    n = min(fx["LAI_H"].size, dy["LAI_H"].size)
    if n == 0:
        raise Unusable("no overlapping days")
    F = derive({k: v[:n] for k, v in fx.items() if k != "year"})
    D = derive({k: v[:n] for k, v in dy.items() if k != "year"})
    lab, yy = period_keys(fx["year"][:n], fx["mo"][:n], freq)

    rows = []
    # Sort periods in CALENDAR order, not lexicographic -- str() would file
    # month 10 between 1 and 2, which makes the CSV tedious to read.
    order = ({m: m for m in range(1, 13)} if freq == "monthly" else
             {"DJF": 1, "MAM": 2, "JJA": 3, "SON": 4})
    combos = sorted({(int(y), p if freq == "seasonal" else int(p))
                     for y, p in zip(yy, lab)},
                    key=lambda t: (t[0], order[t[1]]))
    for var in REPORT:
        if var in RATIOS:
            continue
        if var not in F or var not in D:
            continue
        fa, da_ = np.asarray(F[var], float), np.asarray(D[var], float)
        agg = np.nansum if var in PERIOD_SUM else np.nanmean
        for y, p in combos:
            m = (yy == y) & (lab == p)
            if not m.any() or not np.isfinite(fa[m]).any():
                continue
            f, v = float(agg(fa[m])), float(agg(da_[m]))
            if not (np.isfinite(f) and np.isfinite(v)):
                continue
            rel = 100.0 * (v - f) / abs(f) if abs(f) > 1e-12 else np.nan
            rows.append((y, p, var, f, v, v - f, rel, int(m.sum())))

    for var, (_, den) in RATIOS.items():
        if var not in F or var not in D or den not in F or den not in D:
            continue
        f_r, d_r = np.asarray(F[var], float), np.asarray(D[var], float)
        f_den, d_den = np.asarray(F[den], float), np.asarray(D[den], float)
        if not np.isfinite(f_den).any():
            continue
        floor = 0.01 * float(np.nanmean(np.abs(f_den)))
        ok = (np.abs(f_den) >= floor) & (np.abs(d_den) >= floor)
        for y, p in combos:
            m = (yy == y) & (lab == p) & ok
            if not m.any():
                continue
            if not (np.isfinite(f_r[m]).any() and np.isfinite(d_r[m]).any()):
                continue
            f, v = float(np.nanmean(f_r[m])), float(np.nanmean(d_r[m]))
            if not (np.isfinite(f) and np.isfinite(v)):
                continue
            rel = 100.0 * (v - f) / abs(f) if abs(f) > 1e-12 else np.nan
            rows.append((y, p, var, f, v, v - f, rel, int(m.sum())))
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--freq", choices=["monthly", "seasonal"], default="monthly")
    ap.add_argument("--pair", default=None,
                    help="glob over '<scenario>/<gcm>:<arm>', e.g. 'era5_land:*_ic'")
    ap.add_argument("--stations", default=None, help="comma-separated subset")
    ap.add_argument("--out", type=Path, default=None,
                    help="CSV; a bare name lands in $TC_RESULTS")
    a = ap.parse_args(argv)

    if not a.root.is_dir():
        print(f"ERROR: model_run root not found: {a.root}", file=sys.stderr)
        return 1
    try:
        out = resolve_out(a.out or f"period_effect_{a.freq}.csv")
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    stations = ([s.strip() for s in a.stations.split(",")] if a.stations else
                sorted(p.name for p in a.root.iterdir()
                       if p.is_dir() and next(p.glob("**/fixed_lma*"), None)))
    print(f"model_run : {a.root}\nfreq      : {a.freq}\n"
          f"pair      : {a.pair or '(all)'}\nstations  : {len(stations)}\n",
          flush=True)

    rows, skipped, npairs = [], [], 0
    for st in stations:
        for label, fxp, dyp in find_pairs(a.root, st, a.pair):
            try:
                fx, dy = read_run(fxp), read_run(dyp)
                got = periods(fx, dy, a.freq)
                if not got:
                    raise Unusable("no period produced a value")
                rows.extend((st, label) + r for r in got)
                npairs += 1
            except Unusable as e:
                skipped.append((st, label, str(e)))
            except Exception as e:                       # noqa: BLE001
                skipped.append((st, label, f"{type(e).__name__}: {e}"))

    if skipped:
        print(f"SKIPPED -- {len(skipped)} pair(s):")
        for st, label, why in skipped[:20]:
            print(f"  ! {st:<10}{label:<28}{why}")
        if len(skipped) > 20:
            print(f"  ... and {len(skipped) - 20} more")
        print()
    if not rows:
        print("ERROR: nothing was analysed.", file=sys.stderr)
        return 1

    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["station", "key", "year", "period", "variable",
                    "fixed", "dyn", "diff", "rel_pct", "n_days"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], r[3], r[4],
                        f"{r[5]:.6g}", f"{r[6]:.6g}", f"{r[7]:.6g}",
                        "" if not np.isfinite(r[8]) else f"{r[8]:.3f}", r[9]])

    # Summary: the treatment's size in each period, averaged over stations and
    # years. Read it as "when in the year", with the year-to-year detail left in
    # the file rather than collapsed here.
    print(f"{'variable':<9}" + "".join(f"{p:>9}" for p in
          (range(1, 13) if a.freq == "monthly" else ["DJF", "MAM", "JJA", "SON"])))
    per_list = list(range(1, 13)) if a.freq == "monthly" else ["DJF", "MAM", "JJA", "SON"]
    for var in REPORT:
        sel = [r for r in rows if r[4] == var and np.isfinite(r[8])]
        if not sel:
            continue
        by = {}
        for r in sel:
            by.setdefault(r[3], []).append(abs(r[8]))
        cells = "".join(f"{np.mean(by[p]):>8.2f}%" if p in by else f"{'-':>9}"
                        for p in per_list)
        print(f"{var:<9}{cells}")

    print(f"\n{npairs} pair(s), {len(rows)} row(s) -> {out}")
    return 1 if skipped else 0


if __name__ == "__main__":
    sys.exit(main())
