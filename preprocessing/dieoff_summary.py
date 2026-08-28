#!/usr/bin/env python3
"""Where the canopy collapses, and in which arm.

GPP, LAI and TR fall to near zero together in the SSP scenarios while ET,
LE and H never do -- soil evaporation and sensible heat continue after the
leaves are gone. So the huge percent changes are not a division artefact;
they are dieback, and the ratio cannot describe them.

THE QUESTION A RATIO CANNOT ANSWER. Three situations produce the same
enormous percent change and mean opposite things:

    fixed collapses, dynamic does not   dynamic LMA PREVENTS dieback
    both collapse                       the ratio is noise between two
                                        near-zero numbers, meaningless
    dynamic collapses, fixed does not   dynamic LMA CAUSES dieback

This reports both arms' absolute values and which arm crossed, so the
three are told apart instead of averaged together.

A station-year is "collapsed" for a variable when its value falls below
--thresh (default 1%) of that station's OWN mean for that variable, in the
same scenario. Relative to the site, because 1% of mean annual TR is a
different number of mm at a Colorado conifer than a Carolina pine.

Reports only -- filters nothing, writes nothing except the CSVs asked for.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from results_dir import NoResultsDir, resolve_out                 # noqa: E402
from station_metrics import DATASETS, table_path, read_sites      # noqa: E402

CANOPY = ["GPP", "LAI_H", "T"]
OTHER = ["ET", "QE", "H"]
DISP = {"T": "TR", "LAI_H": "LAI", "QE": "LE"}


def flag(d: pd.DataFrame, thresh: float) -> pd.DataFrame:
    """Per (station, gcm, variable, year): collapsed in fixed / dyn / both."""
    g = ["gcm", "station", "variable"]
    out = d.copy()
    for arm in ("fixed", "dyn"):
        v = pd.to_numeric(out[arm], errors="coerce")
        out[f"_m_{arm}"] = v.groupby([out[c] for c in g]).transform("mean")
        out[f"col_{arm}"] = v < thresh * out[f"_m_{arm}"]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--datasets", default=",".join(DATASETS))
    ap.add_argument("--thresh", type=float, default=0.01)
    ap.add_argument("--out-prefix", default="dieoff")
    a = ap.parse_args(argv)
    try:
        root = Path(a.results or resolve_out(".", create=False))
        out_ev = resolve_out(f"{a.out_prefix}_events.csv")
        out_st = resolve_out(f"{a.out_prefix}_by_station.csv")
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    sites = read_sites()
    events, missing = [], []
    for ds in [x.strip() for x in a.datasets.split(",")]:
        p = table_path(root, ds, "annual")
        if not Path(p).is_file():
            missing.append(f"{ds}: {Path(p).name} not found"); continue
        d = pd.read_csv(p, usecols=lambda c: c in
                        {"station", "key", "year", "variable", "fixed", "dyn"},
                        low_memory=False)
        d["gcm"] = d["key"].astype(str).str.extract(r"^[^/]+/([^:]+):",
                                                    expand=False).fillna("")
        d = d[d["variable"].isin(CANOPY + OTHER)]
        d = flag(d, a.thresh)
        d["dataset"] = ds
        hit = d[d["col_fixed"] | d["col_dyn"]]
        print(f"\n{'=' * 70}\n{ds}: {len(d):,} station-years, "
              f"{len(hit):,} with a collapse in either arm")
        if hit.empty:
            continue
        for var in CANOPY + OTHER:
            h = hit[hit["variable"] == var]
            if h.empty:
                continue
            both = int((h["col_fixed"] & h["col_dyn"]).sum())
            fo = int((h["col_fixed"] & ~h["col_dyn"]).sum())
            do = int((~h["col_fixed"] & h["col_dyn"]).sum())
            print(f"  {DISP.get(var, var):<4} both {both:>6}   "
                  f"fixed only {fo:>6}   dynamic only {do:>6}   "
                  f"stations {h['station'].nunique():>3}   "
                  f"years {int(h['year'].min())}-{int(h['year'].max())}")
        events.append(hit[["dataset", "gcm", "station", "year", "variable",
                           "fixed", "dyn", "col_fixed", "col_dyn"]])

    if not events:
        print("\nno collapse anywhere at this threshold", flush=True)
        if missing:
            for m in missing:
                print(f"  {m}", file=sys.stderr)
            return 1
        return 0

    E = pd.concat(events, ignore_index=True).merge(sites, on="station",
                                                   how="left")
    E["arm"] = np.where(E["col_fixed"] & E["col_dyn"], "both",
                np.where(E["col_fixed"], "fixed_only", "dynamic_only"))
    E.to_csv(out_ev, index=False)
    print(f"\n-> {out_ev}  ({len(E)} rows)")

    # Per station: is this one site failing repeatedly, or many sites once?
    # The ssp126 -> ssp585 near-doubling means one or the other, and they are
    # different findings.
    by = (E.groupby(["dataset", "station", "pft", "variable", "arm"],
                    observed=True)
           .agg(n_years=("year", "nunique"), first=("year", "min"),
                last=("year", "max"), n_gcm=("gcm", "nunique")).reset_index())
    by.to_csv(out_st, index=False)
    print(f"-> {out_st}  ({len(by)} rows)")

    print("\nstations affected per dataset (canopy variables):")
    c = E[E["variable"].isin(CANOPY)]
    print(c.groupby(["dataset", "arm"], observed=True)
           .agg(station_years=("year", "size"),
                stations=("station", "nunique")).to_string())
    if missing:
        print("\nNOT READ:", file=sys.stderr)
        for m in missing:
            print(f"  {m}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
