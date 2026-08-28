#!/usr/bin/env python3
"""Is a large percent change real, or just a small denominator?

rel_pct = 100*(dyn-fixed)/|fixed| explodes when the fixed arm is near zero,
and that is a division artefact rather than a large effect: a year with
fixed TR = 0.4 mm and dynamic TR = 0.8 mm reads +100% for a 0.4 mm
difference. This bins every station-year by how big its fixed-arm value is
RELATIVE TO ITS OWN SITE MEAN, and reports the spread of |rel_pct| in each
bin, per variable and dataset.

How to read it. If med_abs_rel is roughly flat across bins and the largest
|rel_pct| sits in the 75-125% bin, the range is real -- big changes are
happening at ordinary flux levels. If med_abs_rel climbs steeply as the
fraction falls and the maxima live in the low bins, those values are
denominator artefacts and a guard is doing real work.

The guard column is the point of the exercise: at each candidate threshold
it says how many station-years would be dropped and what the largest
surviving |rel_pct| becomes. A threshold that removes a handful of rows and
cuts the maximum from 1e8 to something physical is earning its place; one
that removes thousands is reshaping the result and should not be used.

Nothing here is filtered or written -- it only reports.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from results_dir import NoResultsDir, resolve_out                 # noqa: E402
from station_metrics import DATASETS, table_path                  # noqa: E402

VARS = ["GPP", "LAI_H", "T", "ET", "QE", "H"]
DISP = {"T": "TR", "LAI_H": "LAI", "QE": "LE"}
EDGES = [0, .01, .05, .25, .50, .75, 1.25, np.inf]
NAMES = ["<1%", "1-5%", "5-25%", "25-50%", "50-75%", "75-125%", ">125%"]
GUARDS = [0.001, 0.005, 0.01, 0.05]


def one(d: pd.DataFrame, ds: str, var: str) -> None:
    d = d[d["variable"] == var].copy()
    if d.empty:
        print(f"\n{ds} / {DISP.get(var, var)}: no rows")
        return
    f = pd.to_numeric(d["fixed"], errors="coerce")
    v = pd.to_numeric(d["dyn"], errors="coerce")
    with np.errstate(divide="ignore", invalid="ignore"):
        d["rel"] = 100 * (v - f) / f.abs()
    d["rel"] = d["rel"].replace([np.inf, -np.inf], np.nan)
    site_mean = f.groupby(d["station"]).transform("mean")
    d["frac"] = f / site_mean.replace(0, np.nan)
    d["bin"] = pd.cut(d["frac"], EDGES, labels=NAMES)

    a = d["rel"].abs()
    print(f"\n{ds} / {DISP.get(var, var)}   n={len(d)}   "
          f"median |rel|={np.nanmedian(a):.2f}%   max |rel|={np.nanmax(a):.4g}%")
    t = (d.groupby("bin", observed=True)
           .agg(n=("rel", "size"),
                med=("rel", lambda s: np.nanmedian(s.abs())),
                p95=("rel", lambda s: np.nanpercentile(s.abs(), 95)),
                mx=("rel", lambda s: np.nanmax(s.abs()))))
    if not t.empty:
        print(t.rename(columns={"med": "med_abs_rel", "mx": "max_abs"})
               .to_string(float_format=lambda x: f"{x:,.4g}"))
    for g in GUARDS:
        keep = d["frac"] >= g
        drop = int((~keep).sum())
        mx = np.nanmax(d.loc[keep, "rel"].abs()) if keep.any() else np.nan
        print(f"    guard {g*100:>5.1f}% of site mean -> drop {drop:>6} "
              f"({100*drop/len(d):>5.2f}%), max |rel| becomes {mx:,.4g}%")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--datasets", default=",".join(DATASETS))
    ap.add_argument("--variables", default=",".join(VARS))
    a = ap.parse_args(argv)
    try:
        root = Path(a.results or resolve_out(".", create=False))
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    missing = []
    for ds in [x.strip() for x in a.datasets.split(",")]:
        p = table_path(root, ds, "annual")
        if not Path(p).is_file():
            missing.append(f"{ds}: {Path(p).name} not found")
            continue
        d = pd.read_csv(p, usecols=lambda c: c in
                        {"station", "year", "variable", "fixed", "dyn"},
                        low_memory=False)
        print(f"\n{'=' * 66}\n{ds}   {len(d):,} rows   {Path(p).name}\n{'=' * 66}")
        for var in [x.strip() for x in a.variables.split(",")]:
            one(d, ds, var)
    if missing:
        print("\nNOT READ:", file=sys.stderr)
        for m in missing:
            print(f"  {m}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
