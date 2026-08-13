#!/usr/bin/env python3
"""Annual global-mean CO2 for the historical period and each SSP.

    python fetch_ssp_co2.py --from-netcdf /path/to/input4MIPs/*.nc
    python fetch_ssp_co2.py --report

Writes $TC_INPUT_DATA/co2/co2_<scenario>.csv with columns year,co2_ppm -- one
series per scenario, which build_gcm_meteo.py interpolates onto the hourly stamp.

WHY ONE SERIES PER SCENARIO AND NOT PER MODEL. CMIP6 ScenarioMIP runs are
CONCENTRATION-driven: every model is prescribed the same harmonised greenhouse-gas
concentrations produced with MAGICC7 (Meinshausen et al. 2020, GMD 13, 3571) and
distributed through input4MIPs. So all five GCMs share a scenario's CO2 exactly,
and there is nothing model-specific to fetch. NEX-GDDP-CMIP6 itself carries no CO2
variable at all -- it is nine surface weather fields and nothing else -- which is
why this has to come from outside.

T&C needs Ca globally well-mixed and annual resolution is sufficient (CLAUDE.md
section 4), so the latitudinal and monthly structure in the input4MIPs files is
averaged away here.

SOURCE FILES, from https://esgf-node.llnl.gov/search/input4mips/ or
http://greenhousegases.science.unimelb.edu.au :

    historical 1980-2014   input4MIPs.CMIP6.CMIP.UoM.UoM-CMIP-1-2-0
    ssp126     2015-2100   input4MIPs.CMIP6.ScenarioMIP.UoM.UoM-IMAGE-ssp126-1-2-1
    ssp585     2015-2100   input4MIPs.CMIP6.ScenarioMIP.UoM.UoM-REMIND-MAGPIE-ssp585-1-2-1

Download them once and point --from-netcdf at them; the scenario is taken from the
filename. This script deliberately does NOT ship hard-coded concentration values:
a wrong constant here would scale photosynthesis at every station and every year
without ever failing loudly, so the numbers have to come from the archive.

A two-column CSV (year, ppm) is accepted instead via --from-csv, for the case
where the values were pulled by hand.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path

import numpy as np

INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA",
                                 "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
OUT_DIR = INPUT_ROOT / "co2"
SCENARIOS = {"historical": (1980, 2014), "ssp126": (2015, 2100), "ssp585": (2015, 2100)}
CO2_VARS = ["mole_fraction_of_carbon_dioxide_in_air", "co2"]


def scenario_of(path: Path) -> str | None:
    n = path.name.lower()
    for s in ("ssp126", "ssp585", "ssp245", "ssp370"):
        if s in n:
            return s
    if "historical" in n or re.search(r"uom-cmip", n):
        return "historical"
    return None


def annual_from_netcdf(path: Path):
    """(years, ppm) -- global, annual mean of whatever CO2 field the file holds."""
    try:
        import netCDF4
    except ImportError:
        raise SystemExit("netCDF4 is required -- add it to requirements.txt")
    with netCDF4.Dataset(path) as nc:
        name = next((v for v in CO2_VARS if v in nc.variables), None)
        if name is None:
            cand = [v for v in nc.variables if "co2" in v.lower()]
            if not cand:
                return None, None, f"no CO2 variable in {path.name}"
            name = cand[0]
        var = nc.variables[name]
        vals = np.asarray(var[:], dtype=float)
        units = (getattr(var, "units", "") or "").lower()
        t = nc.variables["time"]
        import netCDF4 as _n
        dates = _n.num2date(t[:], t.units, getattr(t, "calendar", "standard"),
                            only_use_cftime_datetimes=False,
                            only_use_python_datetimes=True)
        years = np.array([d.year for d in dates])
    # collapse every non-time dimension (latitude, sector) to a global mean
    while vals.ndim > 1:
        vals = np.nanmean(vals, axis=-1)
    if "1e-6" in units or "ppm" not in units and vals.max() < 1e-3:
        vals = vals * 1e6                        # mole fraction -> ppm
    out_y, out_c = [], []
    for y in range(int(years.min()), int(years.max()) + 1):
        m = years == y
        if m.any():
            out_y.append(y); out_c.append(float(np.nanmean(vals[m])))
    return np.array(out_y), np.array(out_c), None


def write(scenario, years, ppm, out_dir):
    y0, y1 = SCENARIOS[scenario]
    m = (years >= y0) & (years <= y1)
    if not m.any():
        return 0, f"file covers {years.min()}-{years.max()}, need {y0}-{y1}"
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"co2_{scenario}.csv"
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["year", "co2_ppm"])
        for y, c in zip(years[m], ppm[m]):
            w.writerow([int(y), f"{c:.3f}"])
    return int(m.sum()), None


def report(out_dir):
    print(f"{'scenario':12s}{'years':>14s}{'first ppm':>11s}{'last ppm':>10s}{'change':>9s}")
    for s, (y0, y1) in SCENARIOS.items():
        p = out_dir / f"co2_{s}.csv"
        if not p.is_file():
            print(f"{s:12s}{'MISSING':>14s}"); continue
        rows = list(csv.DictReader(open(p, newline="", encoding="utf-8-sig")))
        y = [int(r["year"]) for r in rows]; c = [float(r["co2_ppm"]) for r in rows]
        print(f"{s:12s}{f'{min(y)}-{max(y)}':>14s}{c[0]:>11.1f}{c[-1]:>10.1f}"
              f"{c[-1]-c[0]:>+9.1f}")
    miss = [s for s in SCENARIOS if not (out_dir / f"co2_{s}.csv").is_file()]
    if miss:
        print(f"\n  missing: {', '.join(miss)} -- see the module docstring for the "
              f"input4MIPs dataset names")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from-netcdf", type=Path, nargs="*", default=None)
    ap.add_argument("--from-csv", type=Path, nargs="*", default=None,
                    help="two-column year,ppm CSV; scenario taken from the filename")
    ap.add_argument("--out", type=Path, default=OUT_DIR)
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args(argv)

    if a.report or not (a.from_netcdf or a.from_csv):
        report(a.out)
        return 0 if a.report else 1

    rc = 0
    for p in list(a.from_netcdf or []) + list(a.from_csv or []):
        s = scenario_of(p)
        if s is None or s not in SCENARIOS:
            print(f"  {p.name}: cannot tell which scenario -- skipped"); rc = 1; continue
        if p.suffix.lower() == ".nc":
            y, c, err = annual_from_netcdf(p)
        else:
            rows = list(csv.reader(open(p, newline="", encoding="utf-8-sig")))
            rows = [r for r in rows if len(r) >= 2 and r[0].strip().isdigit()]
            y = np.array([int(r[0]) for r in rows], float)
            c = np.array([float(r[1]) for r in rows]); err = None
        if err or y is None:
            print(f"  {p.name}: {err}"); rc = 1; continue
        n, werr = write(s, y, c, a.out)
        if werr:
            print(f"  {p.name} -> {s}: {werr}"); rc = 1
        else:
            print(f"  {p.name} -> co2_{s}.csv  ({n} years, "
                  f"{c.min():.1f}-{c.max():.1f} ppm)")
    print()
    report(a.out)
    return rc


if __name__ == "__main__":
    sys.exit(main())
