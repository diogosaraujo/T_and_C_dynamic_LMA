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
CO2_VARS = ["mole_fraction_of_carbon_dioxide_in_air",
            "mole-fraction-of-carbon-dioxide-in-air", "co2"]

# ------------------------------------------------------------------ ESGF download
# The input4MIPs source_id for each scenario. These are the harmonised MAGICC7
# concentrations every CMIP6 model was prescribed, so one per scenario is enough.
SOURCE_IDS = {"historical": "UoM-CMIP-1-2-0",
              "ssp126": "UoM-IMAGE-ssp126-1-2-1",
              "ssp585": "UoM-REMIND-MAGPIE-ssp585-1-2-1"}

# ESGF is mid-migration to ESGF-NG, so more than one index is tried in order.
# esgf-node.llnl.gov now 302s to the ORNL 1.5 bridge, which is queried directly.
ESGF_INDEXES = ["https://esgf-node.ornl.gov/esgf-1-5-bridge",
                "https://esgf-node.llnl.gov/esg-search/search",
                "https://esgf.ceda.ac.uk/esg-search/search",
                "https://esgf-data.dkrz.de/esg-search/search"]

# The facet spellings that actually work, established by querying the live index:
#   * variable_id is HYPHENATED here, not underscored. Using the CF standard-name
#     spelling (mole_fraction_...) returns HTTP 422, not an empty result, so the
#     failure is loud but the cause is not obvious.
#   * grid_label gr1-GMNHSH with frequency yr is the annual global-mean/NH/SH
#     product -- a few kB, and exactly the quantity T&C wants. The alternative
#     gn-15x360deg is monthly on 15 latitude bands and much larger for no gain.
# BOTH spellings are needed. The archive is not self-consistent: the CMIP
# historical collection indexes the variable with HYPHENS, while the ScenarioMIP
# collection uses UNDERSCORES -- verified against the live index, where the wrong
# spelling returns zero documents rather than an error. The filenames are
# hyphenated in both cases, so only the search facet differs.
CO2_VARIABLE_IDS = ["mole-fraction-of-carbon-dioxide-in-air",
                    "mole_fraction_of_carbon_dioxide_in_air"]
PREFERRED = dict(frequency="yr", grid_label="gr1-GMNHSH")


def _esgf_query(index, params, timeout=60):
    import json
    import urllib.parse
    import urllib.request
    url = index + "?" + urllib.parse.urlencode(params)
    # urllib, not curl: curl is not installed on the SOE nodes and its absence
    # looks exactly like a firewall block.
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


def _http_urls(docs):
    """HTTPServer download links from ESGF file docs, newest version first."""
    out = []
    for d in docs:
        for u in d.get("url", []) or []:
            parts = u.split("|")
            if len(parts) >= 3 and parts[2].strip().lower() == "httpserver":
                out.append((d.get("title", ""), parts[0]))
    return out


def discover(scenario, verbose=True):
    """[(filename, url)] for one scenario, trying each index until one answers."""
    base = dict(project="input4MIPs", source_id=SOURCE_IDS[scenario],
                type="File", limit=50, format="application/solr+json")
    for index in ESGF_INDEXES:
        for vid in CO2_VARIABLE_IDS:
            for extra in (PREFERRED, {}):      # narrow first, then widen
                try:
                    js = _esgf_query(index, {**base, "variable_id": vid, **extra})
                except Exception as e:                           # noqa: BLE001
                    if verbose:
                        print(f"    {index.split('/')[2]}: {type(e).__name__}: {e}")
                    break                                        # try next spelling
                docs = js.get("response", {}).get("docs", [])
                urls = _http_urls(docs)
                if urls:
                    if verbose:
                        print(f"    {index.split('/')[2]}: {len(urls)} file(s), "
                              f"variable_id with "
                              f"{'hyphens' if '-' in vid else 'underscores'}"
                              f"{', narrowed' if extra else ', unfiltered'}")
                    return urls
    return []


def download(scenario, raw_dir, force=False, verbose=True):
    """Fetch the input4MIPs CO2 file for one scenario. Returns the path or None."""
    import urllib.request
    raw_dir.mkdir(parents=True, exist_ok=True)
    urls = discover(scenario, verbose)
    if not urls:
        return None, "no HTTPServer URL found on any index"
    name, url = urls[0]
    dest = raw_dir / (name or f"co2_{scenario}.nc")
    if dest.is_file() and not force:
        return dest, "cached"
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as fh:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
    except Exception as e:                                       # noqa: BLE001
        tmp.unlink(missing_ok=True)
        return None, f"{type(e).__name__}: {e}"
    tmp.replace(dest)
    return dest, f"{dest.stat().st_size/1024:.0f} kB"


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
        nc_dims = tuple(var.dimensions)
        t = nc.variables["time"]
        import netCDF4 as _n
        cal = getattr(t, "calendar", "standard")
        # The historical GHG file runs from year 0 (…_gr1-GMNHSH_0000-2014.nc) and
        # its time units reference year 0, which cftime rejects unless year zero is
        # declared legal. Ask for it, and fall back to cftime objects: dates before
        # 1678 cannot be represented as Python datetimes at all, and we only need
        # the integer year.
        try:
            dates = _n.num2date(t[:], t.units, cal, has_year_zero=True,
                                only_use_cftime_datetimes=False,
                                only_use_python_datetimes=True)
        except (ValueError, TypeError):
            dates = _n.num2date(t[:], t.units, cal, has_year_zero=True)
        years = np.array([d.year for d in dates])
    # The gr1-GMNHSH product carries three SECTORS -- global mean, NH, SH -- along
    # a second axis. Averaging them is not the global mean, it is the mean of a
    # mean and its two halves; it happens to land close because GM is roughly
    # (NH+SH)/2, which is exactly why the mistake would never show up in the
    # numbers. Select the global sector explicitly instead.
    if vals.ndim > 1:
        sector = None
        for cand in ("sector", "sector_name", "region"):
            if cand in nc_dims:
                sector = cand
                break
        if sector is not None and vals.shape[-1] == 3:
            vals = vals[..., 0]                  # UoM order: 0 global, 1 NH, 2 SH
            how = "sector 0 (global mean) of 3"
        else:
            how = f"mean over trailing axis of size {vals.shape[-1]}"
            vals = np.nanmean(vals, axis=-1)     # latitude bands: a mean is fair
    else:
        how = "1-D, no reduction"
    while vals.ndim > 1:
        vals = np.nanmean(vals, axis=-1)
    print(f"    {path.name}: dims {nc_dims} -> {how}")
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
    ap.add_argument("--download", action="store_true",
                    help="fetch the input4MIPs files from ESGF, then convert")
    ap.add_argument("--scenario", action="append",
                    help="limit --download to these scenarios (default: all three)")
    ap.add_argument("--raw-dir", type=Path, default=None,
                    help="where downloads land (default: <out>/raw)")
    ap.add_argument("--from-netcdf", type=Path, nargs="*", default=None)
    ap.add_argument("--from-csv", type=Path, nargs="*", default=None,
                    help="two-column year,ppm CSV; scenario taken from the filename")
    ap.add_argument("--out", type=Path, default=OUT_DIR)
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args(argv)

    raw_dir = a.raw_dir or (a.out / "raw")
    fetched = []
    if a.download:
        want = a.scenario or list(SCENARIOS)
        print(f"downloading input4MIPs CO2 for: {', '.join(want)}")
        print(f"  raw dir: {raw_dir}\n")
        for s in want:
            if s not in SOURCE_IDS:
                print(f"  {s}: no input4MIPs source_id known -- skipped"); continue
            print(f"  {s}  ({SOURCE_IDS[s]})")
            p, how = download(s, raw_dir, force=a.force)
            if p is None:
                print(f"    FAILED: {how}")
            else:
                print(f"    -> {p.name}  [{how}]")
                fetched.append(p)
        print()
        if not fetched:
            print("Nothing downloaded. ESGF is mid-migration to ESGF-NG and its\n"
                  "index endpoints move; fetch the three files by hand instead --\n"
                  "dataset names are in the module docstring -- and re-run with\n"
                  "--from-netcdf. The conversion below is unaffected.", file=sys.stderr)
            return 1

    if a.report or not (a.from_netcdf or a.from_csv or fetched):
        report(a.out)
        return 0 if a.report else 1

    rc = 0
    for p in fetched + list(a.from_netcdf or []) + list(a.from_csv or []):
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
