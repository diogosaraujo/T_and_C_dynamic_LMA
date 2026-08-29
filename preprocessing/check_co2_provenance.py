#!/usr/bin/env python3
"""Where each arm's CO2 actually came from, checked against the forcing files.

Read-only and quick -- a few small arrays, no compute -- so it is fine to run in
the login shell, unlike anything that touches the model or the archives.

    source slurm/config.sh
    source $TC_VENV/bin/activate
    python preprocessing/check_co2_provenance.py

THE TWO ARMS DO NOT SHARE A CO2 SOURCE, which is the thing worth confirming:

    ERA5-Land arm   Ca_Data.mat, shipped with the T&C source tree
                    (TeC_Source_Code-master/Inputs), hourly 1975-2022.
                    build_meteo_input.py interpolates it onto the run's stamps.
    GCM arms        $TC_INPUT_DATA/co2/co2_<scenario>.csv, written by
                    fetch_ssp_co2.py from the input4MIPs CMIP6 concentrations
                    (MAGICC7, Meinshausen et al. 2020). historical AND both SSPs.

So "historical" means two different CO2 records depending on which arm you are
looking at: the ERA5 historical arm uses the shipped observational file, the GCM
historical arm uses input4MIPs. They are close but not identical, and a methods
section that names only one of them is wrong for half the runs.

WHAT THIS CHECKS. Not just that the source files exist -- it opens the FORCING
.mat files the runs actually read and compares their Ca against each candidate
source. A file that matches Ca_Data.mat to a few tenths of a ppm came from
Ca_Data.mat; that is provenance, rather than a claim about which script was
supposed to have run.

Published anchors are printed for scale, not as a pass/fail: the concentration
pathways are ~397 ppm at 2014, ~446 at ssp126 2100 and ~1135 at ssp585 2100.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import sys
from pathlib import Path

import numpy as np

PREPROC = Path(__file__).resolve().parent
REPO_ROOT = PREPROC.parent
SHIPPED_CA = REPO_ROOT / "T&C" / "TeC_Source_Code-master" / "Inputs" / "Ca_Data.mat"

# Approximate published CMIP6 concentrations, for scale only.
ANCHORS = {("historical", 2014): 397.5, ("ssp126", 2100): 445.6,
           ("ssp585", 2100): 1135.2}


def matlab_to_date(x: float) -> dt.date:
    """MATLAB datenum -> date. Offset 366 for the year-0 origin."""
    return dt.date.fromordinal(int(x) - 366)


def read_mat_vars(path: Path, want) -> dict:
    """Named variables from a .mat, whatever version it is.

    h5py first: every Meteo_*.mat here is -v7.3, and scipy's failure mode on
    those differs between machines (NotImplementedError locally, OSError on the
    cluster), which is how an earlier version silently read nothing.
    """
    out = {}
    try:
        import h5py
        with h5py.File(path, "r") as f:
            for k in want:
                if k in f:
                    out[k] = np.array(f[k]).ravel()
        if out:
            return out
    except Exception:                                  # noqa: BLE001
        pass
    try:
        import scipy.io as sio
        m = sio.loadmat(path, squeeze_me=True)
        for k in want:
            if k in m:
                out[k] = np.asarray(m[k]).ravel()
    except Exception as e:                             # noqa: BLE001
        print(f"    ! cannot read {path.name}: {type(e).__name__}: {e}")
    return out


def shipped_series(path: Path):
    """(datenum, ppm) from Ca_Data.mat, or (None, None)."""
    if not path.is_file():
        print(f"  ! not found: {path}")
        return None, None
    d = read_mat_vars(path, ("Ca", "Ca_all", "ca", "Date_CO2", "Date"))
    ca = next((d[k] for k in ("Ca", "Ca_all", "ca") if k in d), None)
    ax = next((d[k] for k in ("Date_CO2", "Date") if k in d), None)
    if ca is None or ax is None:
        print(f"  ! {path.name} has no recognisable Ca/date pair: {list(d)}")
        return None, None
    print(f"  {path.name}: n={ca.size}, {np.nanmin(ca):.1f}-{np.nanmax(ca):.1f} ppm, "
          f"{matlab_to_date(ax.min())} to {matlab_to_date(ax.max())}")
    for yr in (1985, 2000, 2014, 2020):
        t = dt.date(yr, 7, 1).toordinal() + 366
        if ax.min() <= t <= ax.max():
            print(f"      {yr} mid-year: {ca[np.argmin(np.abs(ax - t))]:.1f} ppm")
    return ax, ca


def csv_series(path: Path):
    """(years, ppm) from co2_<scenario>.csv, or (None, None)."""
    if not path.is_file():
        return None, None
    y, c = [], []
    with path.open(newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            try:
                y.append(int(r["year"])); c.append(float(r["co2_ppm"]))
            except (KeyError, TypeError, ValueError):
                continue
    if not y:
        return None, None
    return np.array(y, float), np.array(c, float)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model-run", type=Path, default=None)
    ap.add_argument("--co2-dir", type=Path, default=None)
    ap.add_argument("--ca-file", type=Path, default=SHIPPED_CA)
    ap.add_argument("--stations", type=int, default=2,
                    help="how many forcing files to open per arm")
    a = ap.parse_args(argv)

    mr = a.model_run or Path(os.environ.get("MODEL_RUN", ""))
    co2_dir = a.co2_dir or Path(os.environ.get("TC_INPUT_DATA", "")) / "co2"

    print("=" * 72)
    print("1. ERA5-Land arm source -- Ca_Data.mat, shipped with the T&C source")
    print("=" * 72)
    ax, ca = shipped_series(Path(a.ca_file))

    print()
    print("=" * 72)
    print("2. GCM arm source -- co2_<scenario>.csv from input4MIPs")
    print("=" * 72)
    print(f"  {co2_dir}")
    series = {}
    for sc in ("historical", "ssp126", "ssp585"):
        y, c = csv_series(co2_dir / f"co2_{sc}.csv")
        if y is None:
            print(f"  {sc:12s} MISSING -- run fetch_ssp_co2.py")
            continue
        series[sc] = (y, c)
        print(f"  {sc:12s} {int(y.min())}-{int(y.max())}  "
              f"{c[0]:.1f} -> {c[-1]:.1f} ppm  (n={y.size})")
        for (s2, yr), ref in ANCHORS.items():
            if s2 == sc and yr in y:
                got = float(c[y == yr][0])
                flag = "ok" if abs(got - ref) < 5 else "DIFFERS"
                print(f"      {yr}: {got:.1f} ppm vs published ~{ref} ppm  [{flag}]")

    print()
    print("=" * 72)
    print("3. What the FORCING FILES actually contain (the decisive check)")
    print("=" * 72)
    if not mr or not Path(mr).is_dir():
        print(f"  ! MODEL_RUN='{mr}' is not a directory; skipping.")
        print("    Run 'source slurm/config.sh' first.")
        return 0

    # One ERA5 file and one per GCM scenario, then match each against both
    # candidate sources. Whichever it reproduces is where its CO2 came from.
    groups = {"era5_land": [], "historical": [], "ssp126": [], "ssp585": []}
    for f in Path(mr).glob("*/**/Meteo_*.mat"):
        s = str(f).lower().replace("\\", "/")
        for key in groups:
            if f"/{key}/" in s and len(groups[key]) < a.stations:
                groups[key].append(f)
                break

    for key, files in groups.items():
        if not files:
            print(f"  {key:12s} no forcing file found")
            continue
        for f in files:
            d = read_mat_vars(f, ("Ca", "Date"))
            if "Ca" not in d or "Date" not in d:
                print(f"  {key:12s} {f.name}: no Ca/Date")
                continue
            cav, dav = d["Ca"], d["Date"]
            ok = np.isfinite(cav)
            if not ok.any():
                print(f"  {key:12s} {f.name}: Ca is all NaN")
                continue
            print(f"  {key:12s} {f.name}")
            print(f"      Ca {np.nanmin(cav):.1f}-{np.nanmax(cav):.1f} ppm, "
                  f"{matlab_to_date(dav.min())} to {matlab_to_date(dav.max())}")
            # Compare against the shipped hourly record.
            if ax is not None:
                pred = np.interp(dav[ok], ax, ca)
                print(f"      vs Ca_Data.mat        : max |diff| "
                      f"{np.nanmax(np.abs(cav[ok] - pred)):.3f} ppm")
            # And against the input4MIPs annual series for this scenario.
            sc = key if key in series else None
            if sc:
                yrs = np.array([matlab_to_date(t).year +
                                (matlab_to_date(t).timetuple().tm_yday - 1) / 365.0
                                for t in dav[ok]])
                pred = np.interp(yrs, series[sc][0], series[sc][1])
                print(f"      vs co2_{sc}.csv{' ' * max(0, 9 - len(sc))}: max |diff| "
                      f"{np.nanmax(np.abs(cav[ok] - pred)):.3f} ppm")
    print()
    print("A max |diff| near zero identifies the source. The ERA5 files should")
    print("match Ca_Data.mat; the GCM files should match their co2_*.csv.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
