#!/usr/bin/env python3
"""Is the collapse a trajectory, or scattered single years?

The die-off run showed dynamic-only collapses outnumbering fixed-only by
4-8x under the SSPs, with nothing at all in ERA5 or GCM-historical. That is
either a dieback trajectory -- LMA rises, LAI falls, less carbon, LAI falls
further -- or the solver misbehaving at low LAI. The two look identical in a
count and completely different in time.

  DIEBACK          collapse starts, then persists to 2100. Long runs, last
                   year at or near the end of record, high duty cycle.
  NUMERICAL BLIPS  isolated years that recover. n_years small relative to
                   the span they are spread over, last year nowhere near
                   2100, duty cycle low.

  FORCED BY CLIMATE     most GCMs agree at that station
  ONE MODEL'S TRAJECTORY  n_gcm is 1

Reported per station-arm-variable, and summarised. Reads only the two CSVs
the die-off run already wrote; computes nothing new from the big tables.

longest_run is the key column: the most consecutive collapse years seen at
that station in any single GCM. A station with n_years=40 spread as forty
isolated years is a different object from one with a single 40-year run,
and only the run length tells them apart.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from results_dir import NoResultsDir, resolve_out                 # noqa: E402

END = 2100
DISP = {"T": "TR", "LAI_H": "LAI", "QE": "LE"}


def longest_run(years: np.ndarray) -> int:
    """Most consecutive years in a sorted unique array."""
    if years.size == 0:
        return 0
    y = np.unique(years.astype(int))
    best = run = 1
    for a, b in zip(y[:-1], y[1:]):
        run = run + 1 if b == a + 1 else 1
        best = max(best, run)
    return int(best)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--events", default="dieoff_events.csv")
    ap.add_argument("--out", default="dieoff_persistence.csv")
    a = ap.parse_args(argv)
    try:
        root = Path(a.results or resolve_out(".", create=False))
        out_p = resolve_out(a.out)
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    p = root / a.events
    if not p.is_file():
        print(f"ERROR: {p} not found -- run dieoff_summary.py first",
              file=sys.stderr)
        return 1
    E = pd.read_csv(p, low_memory=False)
    for c in ("dataset", "station", "variable", "arm", "year", "gcm"):
        if c not in E.columns:
            print(f"ERROR: {p.name} has no {c!r} column", file=sys.stderr)
            return 1
    E["year"] = pd.to_numeric(E["year"], errors="coerce")
    E = E[E["year"].notna()]
    E["gcm"] = E["gcm"].fillna("").astype(str)

    rows = []
    for keys, g in E.groupby(["dataset", "station", "pft", "variable", "arm"],
                             observed=True, dropna=False):
        yrs = g["year"].to_numpy(float)
        # Runs are per GCM: the same year appearing under three models is not
        # a three-year run, and concatenating them would invent persistence
        # that no single trajectory has.
        runs = [longest_run(sub["year"].to_numpy(float))
                for _, sub in g.groupby("gcm", observed=True)]
        first, last = int(np.nanmin(yrs)), int(np.nanmax(yrs))
        span = last - first + 1
        n_years = int(pd.Series(yrs).nunique())
        rows.append(dict(zip(["dataset", "station", "pft", "variable", "arm"],
                             keys))
                    | {"n_gcm": int(g["gcm"].nunique()),
                       "n_years": n_years, "first": first, "last": last,
                       "span": span,
                       "duty": round(n_years / span, 3) if span else np.nan,
                       "longest_run": int(max(runs)) if runs else 0,
                       "reaches_end": bool(last >= END - 5)})
    P = pd.DataFrame(rows)
    P.to_csv(out_p, index=False)
    print(f"-> {out_p}  ({len(P)} rows)\n")

    c = P[P["variable"].isin(["GPP", "LAI_H", "T"])]
    print("PERSISTENCE, canopy variables, by dataset and arm")
    print(c.groupby(["dataset", "arm"], observed=True)
           .agg(station_vars=("station", "size"),
                med_longest_run=("longest_run", "median"),
                max_longest_run=("longest_run", "max"),
                med_duty=("duty", "median"),
                pct_reaching_2100=("reaches_end",
                                   lambda s: round(100 * s.mean(), 1)),
                med_n_gcm=("n_gcm", "median")).to_string())

    print("\nHow the collapses are shaped (canopy variables):")
    for lo, hi, lab in ((1, 1, "single isolated years"),
                        (2, 5, "2-5 consecutive"),
                        (6, 20, "6-20 consecutive"),
                        (21, 10_000, "21+ consecutive")):
        s = c[(c["longest_run"] >= lo) & (c["longest_run"] <= hi)]
        if len(c):
            print(f"  {lab:<22} {len(s):>5}  ({100*len(s)/len(c):>5.1f}%)"
                  f"   reaching 2100: {100*s['reaches_end'].mean():>5.1f}%"
                  if len(s) else f"  {lab:<22} {0:>5}")

    print("\nTop 12 by longest run:")
    top = c.sort_values("longest_run", ascending=False).head(12)
    print(top[["dataset", "station", "pft", "variable", "arm", "n_gcm",
               "n_years", "first", "last", "longest_run",
               "reaches_end"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
