#!/usr/bin/env python3
"""Verify downloaded ERA5-Land station files before they feed the T&C forcing builder.

The download step only checks that a file is netCDF and non-empty. That is not enough:
a file can be well-formed and still be wrong for our purposes -- truncated period,
missing hours, an all-NaN land-mask miss, or an accumulation convention that differs
from what the pipeline assumes. Each of those produces plausible-looking forcing and
surfaces much later as bad ET or GPP, so check them here instead.

Per station it checks:

  structure     all expected group files + sidecar present and readable
  variables     every requested short name present, in the group the registry claims
  time axis     hourly, strictly increasing, no gaps or duplicates, covers the period,
                and identical across all four group files
  coverage      fraction of missing values (an all-NaN variable means the grid point
                resolved over water -- ERA5-Land is land-only)
  ranges        physically plausible values, which catches unit surprises
  accumulation  EMPIRICAL confirmation that tp and ssrd reset at 00 UTC rather than
                being per-hour fluxes (see check_accumulation for the logic)

Exit status is 0 only if every station passes. Warnings do not fail the run.

    python verify_era5_land.py                       # all stations found on disk
    python verify_era5_land.py --stations US-HBK,US-Ha2
    python verify_era5_land.py --report report.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

try:
    import xarray as xr
except ImportError:
    sys.exit("xarray is required. Run:  pip install -r requirements.txt")

from era5_variables import TIMESERIES_GROUPS, VARIABLES

INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA", "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
DEFAULT_DIR = INPUT_ROOT / "era5_land"

# Names the CDS has used for the hourly time coordinate across API generations.
TIME_NAMES = ("valid_time", "time", "forecast_reference_time")

# Plausible physical ranges in the NATIVE ERA5-Land units. Deliberately wide: these
# catch unit changes and corruption, not subtle bias.
RANGES = {
    "t2m": (180.0, 340.0, "K"),
    "d2m": (180.0, 340.0, "K"),
    "sp": (40000.0, 110000.0, "Pa"),
    "tp": (0.0, 2.0, "m accumulated/day"),
    "u10": (-120.0, 120.0, "m/s"),
    "v10": (-120.0, 120.0, "m/s"),
    "ssrd": (0.0, 5.0e7, "J/m2 accumulated/day"),
}

ACCUMULATED = {v["short_name"] for v in VARIABLES if v["time_convention"].startswith("accumulated")}
GROUP_OF = {v["short_name"]: v["timeseries_group"] for v in VARIABLES}


class Result:
    def __init__(self, station_id: str):
        self.station_id = station_id
        self.errors: list[str] = []
        self.warnings: list[str] = []
        self.info: dict = {}

    def error(self, msg: str) -> None:
        self.errors.append(msg)

    def warn(self, msg: str) -> None:
        self.warnings.append(msg)

    @property
    def status(self) -> str:
        if self.errors:
            return "FAIL"
        return "WARN" if self.warnings else "PASS"


def find_time(ds) -> str | None:
    for name in TIME_NAMES:
        if name in ds.coords or name in ds.dims:
            return name
    for name in ds.coords:
        if np.issubdtype(ds[name].dtype, np.datetime64):
            return str(name)
    return None


def check_time_axis(res: Result, group: str, times: np.ndarray, args) -> None:
    if times.size == 0:
        res.error(f"{group}: time axis is empty")
        return

    diffs = np.diff(times).astype("timedelta64[s]").astype(np.int64)
    if diffs.size:
        if (diffs <= 0).any():
            res.error(f"{group}: time axis is not strictly increasing "
                      f"({int((diffs <= 0).sum())} non-positive steps)")
        off = diffs[diffs != 3600]
        if off.size:
            uniq = sorted({int(v) for v in off})[:5]
            res.error(f"{group}: {off.size} non-hourly step(s), e.g. {uniq} s "
                      "-- gaps or duplicated timestamps")

    start = times[0].astype("datetime64[s]").item()
    end = times[-1].astype("datetime64[s]").item()
    res.info.setdefault("period", {})[group] = f"{start:%Y-%m-%d %H:%M} .. {end:%Y-%m-%d %H:%M}"

    # Compare against the exact requested bounds, not just the years: a file truncated
    # a few days into the final year still has the right end year and would slip past
    # a year-granular check.
    want_start = datetime(args.start_year, 1, 1, 0)
    want_end = datetime(args.end_year, 12, 31, 23)
    if start > want_start:
        res.error(f"{group}: starts {start:%Y-%m-%d %H:%M}, after requested "
                  f"{want_start:%Y-%m-%d %H:%M}")
    if end < want_end:
        res.error(f"{group}: ends {end:%Y-%m-%d %H:%M}, before requested "
                  f"{want_end:%Y-%m-%d %H:%M} -- the series is truncated")

    expected = int((want_end - want_start).total_seconds() // 3600) + 1
    if times.size != expected and not res.errors:
        res.warn(f"{group}: {times.size} timesteps, expected {expected} for an hourly "
                 f"{args.start_year}-{args.end_year} series")


def check_accumulation(res: Result, name: str, values: np.ndarray, times: np.ndarray) -> None:
    """Empirically confirm the daily-reset accumulation convention.

    ERA5-Land accumulations restart at 00 UTC, so within a day the series rises
    monotonically and then drops once at the boundary. A per-hour flux series would
    instead fall about as often as it rises. Counting the fraction of negative steps
    separates the two cleanly: ~1/24 = 4% for daily accumulation, near 50% for a flux.
    """
    finite = np.isfinite(values)
    if finite.sum() < 48:
        res.warn(f"{name}: too few finite values to test the accumulation convention")
        return

    v = values.astype("float64")
    diffs = np.diff(v)
    ok = np.isfinite(diffs)
    if ok.sum() == 0:
        return
    neg_frac = float((diffs[ok] < -1e-12).sum()) / float(ok.sum())
    res.info.setdefault("negative_step_fraction", {})[name] = round(neg_frac, 4)

    if neg_frac > 0.20:
        res.error(
            f"{name}: {neg_frac:.1%} of hourly steps decrease -- this does NOT look like "
            "the documented accumulate-from-00-UTC convention (expected ~4%). The "
            "de-accumulation logic in the forcing builder would be wrong; inspect a day "
            "of raw values before proceeding."
        )
        return
    if neg_frac < 0.01:
        res.warn(
            f"{name}: only {neg_frac:.2%} of steps decrease (expected ~4% for a daily "
            "reset). Possibly a longer accumulation period -- verify before de-accumulating."
        )
        return

    # The drops should land on the 00 UTC boundary, not scattered through the day.
    drop_idx = np.nonzero(diffs < -1e-12)[0]
    if drop_idx.size:
        hours = times[drop_idx + 1].astype("datetime64[h]").astype(np.int64) % 24
        at_midnight = float((hours == 0).sum()) / float(hours.size)
        res.info.setdefault("drops_at_00utc", {})[name] = round(at_midnight, 4)
        if at_midnight < 0.9:
            res.warn(
                f"{name}: only {at_midnight:.1%} of decreases occur at 00 UTC; the reset "
                "may not be on the day boundary this pipeline assumes."
            )


def verify_station(station_dir: Path, args) -> Result:
    sid = station_dir.name
    res = Result(sid)

    sidecar = station_dir / f"{sid}_ERA5_Land.json"
    if not sidecar.exists():
        res.warn("metadata sidecar is missing")

    axes: dict[str, np.ndarray] = {}
    seen: set[str] = set()

    for group in TIMESERIES_GROUPS:
        path = station_dir / f"{sid}_ERA5_Land_{group}.nc"
        if not path.exists():
            res.error(f"{group}: file missing ({path.name})")
            continue
        if path.stat().st_size == 0:
            res.error(f"{group}: file is empty")
            continue

        try:
            ds = xr.open_dataset(path)
        except Exception as exc:
            res.error(f"{group}: cannot open ({exc})")
            continue

        with ds:
            tname = find_time(ds)
            if tname is None:
                res.error(f"{group}: no recognisable time coordinate "
                          f"(coords: {list(ds.coords)})")
                continue
            times = ds[tname].values
            axes[group] = times
            check_time_axis(res, group, times, args)

            for var in ds.data_vars:
                short = str(var)
                if short not in GROUP_OF:
                    continue
                seen.add(short)
                if GROUP_OF[short] != group:
                    res.warn(f"{short} found in group '{group}' but registry says "
                             f"'{GROUP_OF[short]}'")

                values = np.asarray(ds[var].values).squeeze()
                if values.ndim > 1:
                    values = values.reshape(values.shape[0], -1)[:, 0]

                finite = np.isfinite(values)
                missing = 1.0 - (float(finite.sum()) / float(values.size or 1))
                res.info.setdefault("missing_fraction", {})[short] = round(missing, 4)
                if missing >= 1.0:
                    res.error(f"{short}: all values missing -- the grid point most likely "
                              "resolved over water (ERA5-Land is land-only)")
                    continue
                if missing > 0.01:
                    res.warn(f"{short}: {missing:.1%} of values missing")

                lo, hi, unit = RANGES.get(short, (None, None, ""))
                if lo is not None:
                    vals = values[finite]
                    if vals.size and (vals.min() < lo or vals.max() > hi):
                        res.error(
                            f"{short}: values outside the plausible range "
                            f"[{lo}, {hi}] {unit} (got {vals.min():.4g} .. {vals.max():.4g}) "
                            "-- check whether the units are what the registry claims"
                        )

                if short in ACCUMULATED and not args.quick:
                    check_accumulation(res, short, values, times)

    expected = {v["short_name"] for v in VARIABLES}
    for missing_var in sorted(expected - seen):
        res.error(f"variable '{missing_var}' not found in any group file")

    # All group files must share one time axis, or they cannot be merged downstream.
    if len(axes) > 1:
        ref_group, ref = next(iter(axes.items()))
        for group, times in axes.items():
            if group == ref_group:
                continue
            if times.shape != ref.shape or not np.array_equal(times, ref):
                res.error(f"{group}: time axis differs from '{ref_group}' -- the group "
                          "files cannot be combined into one dataset")
                break
    if axes:
        res.info["n_timesteps"] = int(next(iter(axes.values())).size)

    return res


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Verify downloaded ERA5-Land station files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--dir", type=Path, default=DEFAULT_DIR, help="era5_land output directory")
    p.add_argument("--stations", default=None, help="comma-separated StationIDs to check")
    p.add_argument("--start-year", type=int, default=1985)
    p.add_argument("--end-year", type=int, default=2021)
    p.add_argument("--quick", action="store_true",
                   help="skip the accumulation-convention check (faster)")
    p.add_argument("--report", type=Path, default=None, help="write a JSON report here")
    args = p.parse_args(argv)

    if not args.dir.is_dir():
        raise SystemExit(f"not a directory: {args.dir}")

    wanted = {s.strip() for s in args.stations.split(",") if s.strip()} if args.stations else None
    dirs = sorted(d for d in args.dir.iterdir() if d.is_dir())
    if wanted is not None:
        dirs = [d for d in dirs if d.name in wanted]
        for missing in sorted(wanted - {d.name for d in dirs}):
            print(f"  ! {missing}: no directory in {args.dir}")
    if not dirs:
        raise SystemExit(f"no station directories found in {args.dir}")

    print(f"verifying {len(dirs)} station(s) in {args.dir}")
    print(f"expected period: {args.start_year}-01-01 .. {args.end_year}-12-31, hourly\n")

    results = [verify_station(d, args) for d in dirs]

    for res in sorted(results, key=lambda r: (r.status == "PASS", r.station_id)):
        if res.status == "PASS":
            print(f"  PASS  {res.station_id}  ({res.info.get('n_timesteps', '?')} timesteps)")
            continue
        print(f"  {res.status}  {res.station_id}")
        for msg in res.errors:
            print(f"          ERROR: {msg}")
        for msg in res.warnings:
            print(f"          warn : {msg}")

    n_fail = sum(r.status == "FAIL" for r in results)
    n_warn = sum(r.status == "WARN" for r in results)
    n_pass = sum(r.status == "PASS" for r in results)
    print(f"\n{n_pass} passed, {n_warn} with warnings, {n_fail} failed "
          f"(of {len(results)})")

    if args.report:
        args.report.write_text(json.dumps(
            [{"station_id": r.station_id, "status": r.status, "errors": r.errors,
              "warnings": r.warnings, "info": r.info} for r in results],
            indent=2), encoding="utf-8")
        print(f"report: {args.report}")

    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
