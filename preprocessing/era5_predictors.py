#!/usr/bin/env python3
"""Recompute PLSR climate predictors straight from the ERA5-Land monthly stacks.

A port of build_predictor_names and build_climate_predictor_row from
PLSR_PREPROCESS_PIXEL_CLIM_DOY_CORE.m. It exists because the preprocessed tables
(LMA_ecoregion_no<ii>.csv) only emit a row where the pixel's dominant NLCD class
that year was forest AND a DOY-climatology key existed. Every other pixel-year is
absent, so a station's series has holes -- or its own ERA5-Land cell is missing
from the table altogether -- even though the ERA5-Land forcing behind it exists.
This module goes back to that forcing.

    store = Era5Monthly(Path('/vol_efthymios/NFS07/Data/ERA5_Land/monthly'))
    series = store.pixel_series(45.20, -68.74)
    vals   = store.predictor_row(series, date(2001, 7, 15), ['P - 3mo', 'SSRD - 7mo'])

Everything is derived from the MATLAB source, so the conventions that are easy to
get wrong are kept deliberately visible:

  * PET predictors are NEGATED (`row = -1 * pet`), so a larger value is less
    evaporative demand. Copying the raw variable would flip the sign of beta.
  * Two georeferences are in play. SPEI12/SPI3/SPI6/SPI12 are indexed [lat, lon]
    on a -180..180 grid; SPEI3/SPEI6/tp/tas/pet/ssrd are indexed [lon, lat] on the
    native 0..360 grid. Mixing them silently returns a pixel on the wrong continent.
  * The drought indices start on different months (SPEI12/SPI12 in 1980-12,
    SPEI3/SPI3 in 1980-03, SPEI6/SPI6 in 1980-06), so each has its own time axis.
  * Seasonal means drop the previous year's occurrence of the sampling season
    before averaging.

Verify before trusting: `python era5_predictors.py verify <table.csv>` recomputes
predictors for rows that ARE in a preprocessed table and reports the mismatch.
Run that on the cluster before any fallback value is used.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from datetime import date, timedelta
from pathlib import Path

import numpy as np

DEFAULT_ERA5_ROOT = Path("/vol_efthymios/NFS07/Data/ERA5_Land/monthly")

# opts.DroughtPath in the MATLAB source. Each entry is (filename, variable,
# axis order). 'latlon' means the array is [lat, lon, time] on the -180..180
# reference R; 'lonlat' means [lon, lat, time] on the native 0..360 reference
# R_tp_tas. This split is in the MATLAB indexing and is not negotiable.
DROUGHT_FILES = {
    "SPEI12_ts": ("Global_SPEI12.mat", "R_SI", "latlon", (1980, 12)),
    "SPEI3_ts": ("Global_SPEI_3.mat", "spei", "lonlat", (1980, 3)),
    "SPEI6_ts": ("Global_SPEI_6.mat", "spei", "lonlat", (1980, 6)),
    "SPI12_ts": ("Global_SPI12.mat", "R_SI", "latlon", (1980, 12)),
    "SPI3_ts": ("Global_SPI3.mat", "R_SI", "latlon", (1980, 3)),
    "SPI6_ts": ("Global_SPI6.mat", "R_SI", "latlon", (1980, 6)),
}
CLIMATE_FILES = {
    "tp": ("tp_1980_2022_monthly.mat", ("data_all",), "lonlat"),
    "tas": ("t2m_1980_2022_monthly.mat", ("data_all",), "lonlat"),
    "pet": ("e_1980_2022_monthly.mat", ("data_all",), "lonlat"),
    "ssrd": ("ssrd_1980_2022_monthly.mat", ("data_all", "ssrd", "SSRD"), "lonlat"),
}
CLIMATE_START = (1980, 1)          # t_climate = datetime(1980,1,30):calmonths(1):...
GRID_NC = "1980.nc"
SI_ORDER = ["SPEI12_ts", "SPEI3_ts", "SPEI6_ts", "SPI12_ts", "SPI3_ts", "SPI6_ts"]
SEASON_BY_MONTH = ["DJF", "DJF", "MAM", "MAM", "MAM", "JJA",
                   "JJA", "JJA", "SON", "SON", "SON", "DJF"]


def predictor_names() -> list[str]:
    """The 146 predictor names, in table-column order (build_predictor_names)."""
    names: list[str] = [""] * 146
    for i, si in enumerate(SI_ORDER):
        base = 13 * i
        names[base] = si
        for lag in range(1, 12):
            names[base + lag] = f"{si}-{lag}mo"
        names[base + 12] = f"{si}-sev"
    names[78], names[79], names[80] = "MAR", "MAT", "PET"
    for i in range(1, 13):
        names[80 + i] = f"P - {i - 1}mo"
        names[92 + i] = f"T - {i - 1}mo"
        names[104 + i] = f"PET - {i - 1}mo"
    for off, var in ((117, "P"), (118, "T"), (119, "PET")):
        for k, season in enumerate(("Summer", "Spring", "Fall", "Winter")):
            names[off + 3 * k] = f"{var} - {season}"
    names[129] = "SSRD"
    for i in range(1, 13):
        names[129 + i] = f"SSRD - {i - 1}mo"
    for k, season in enumerate(("Summer", "Spring", "Fall", "Winter")):
        names[142 + k] = f"SSRD - {season}"
    return names


PREDICTOR_NAMES = predictor_names()
NAME_INDEX = {n: i for i, n in enumerate(PREDICTOR_NAMES)}


def ym_key(y: int, m: int) -> int:
    return y * 100 + m


def month_axis(start: tuple[int, int], n: int) -> np.ndarray:
    """n consecutive monthly year*100+month keys from `start`."""
    y, m = start
    out = np.empty(n, dtype=int)
    for i in range(n):
        out[i] = ym_key(y, m)
        m += 1
        if m == 13:
            y, m = y + 1, 1
    return out


def si_to_si_droughts(si: np.ndarray, threshold: float = -1.0) -> np.ndarray:
    """Drought-severity series behind the '-sev' predictors -- NOT VERIFIED.

    SI_to_SIdroughts lives in /vol_efthymios/NFS07/dd1136/functions/ and is not
    reproduced here. The obvious reading -- keep the index where it sits at or
    below the threshold, zero elsewhere -- was tested against the real eco1 table
    (job 35425) and is WRONG: it returns 0 for windows the table gives severities
    of 0.2 to 5.6, because those windows never touch -1 at all. The real function
    is evidently run-based (a drought event spanning the months either side of the
    crossing, in the theory-of-runs sense) rather than a per-month threshold.

    Kept only so the shape of the calculation is documented. era5_fill refuses to
    fill a group whose fit selected a '-sev' predictor, so this is not reached in
    normal use; SEV_VERIFIED flips once the real function is ported.
    """
    out = np.zeros_like(si, dtype=float)
    mask = np.isfinite(si) & (si <= threshold)
    out[mask] = si[mask]
    return out


# Job 35425: 97 of 150 '-sev' values disagreed with the table, while all 3,500
# ported values matched. Until SI_to_SIdroughts is ported, a fit that selected a
# '-sev' predictor cannot be filled from ERA5-Land.
SEV_VERIFIED = False


def unverifiable(predictors) -> list[str]:
    """Selected predictors this module cannot yet reproduce."""
    if SEV_VERIFIED:
        return []
    return [p for p in predictors if p.endswith("-sev")]


class Era5Monthly:
    """Lazy reader over the monthly ERA5-Land stacks used by the PLSR preprocessing."""

    def __init__(self, root: Path = DEFAULT_ERA5_ROOT):
        self.root = Path(root)
        if not self.root.is_dir():
            raise RuntimeError(f"ERA5-Land monthly root not found: {self.root}")
        self._handles: dict[Path, object] = {}
        self._grid_ready = False

    # ------------------------------------------------------------------ grid
    def _load_grid(self) -> None:
        """Latitude/longitude postings, exactly as load_climate_inputs builds them."""
        if self._grid_ready:
            return
        nc = self.root / GRID_NC
        if not nc.is_file():
            raise RuntimeError(f"grid file not found: {nc}")
        try:
            from netCDF4 import Dataset
            with Dataset(str(nc)) as ds:
                lat = np.asarray(ds.variables["latitude"][:], dtype=float)
                lon_native = np.asarray(ds.variables["longitude"][:], dtype=float)
        except ImportError:
            import xarray as xr
            with xr.open_dataset(nc) as ds:
                lat = np.asarray(ds["latitude"].values, dtype=float)
                lon_native = np.asarray(ds["longitude"].values, dtype=float)

        # longitude = [lon(1801:end)-360, lon(1:1800)] -- MATLAB is 1-based, so the
        # split is at index 1800 here. Guard the assumption rather than trust it.
        split = 1800
        if lon_native.size <= split:
            raise RuntimeError(
                f"{GRID_NC}: expected >{split} longitudes for the -180..180 shift, "
                f"found {lon_native.size}")
        lon_shift = np.concatenate([lon_native[split:] - 360.0, lon_native[:split]])

        self.lat = lat
        self.lon_native = lon_native
        self.lon_shift = lon_shift
        # R (shifted) and R_tp_tas (native): postings spanning min..max.
        self._ref_latlon = _Postings(lat, lon_shift)
        self._ref_lonlat = _Postings(lat, lon_native)
        self._grid_ready = True

    def cell_center(self, lat: float, lon: float) -> tuple[float, float]:
        """Centre of the ERA5-Land posting nearest (lat, lon), in -180..180."""
        self._load_grid()
        i, j = self._ref_latlon.index(lat, lon)
        return float(self.lat[self._ref_latlon.lat_pos(i)]), \
            float(self.lon_shift[self._ref_latlon.lon_pos(j)])

    # ------------------------------------------------------------------- io
    def _dataset(self, fname: str, candidates: tuple[str, ...]):
        """Open a .mat stack lazily. v7.3 is HDF5, which h5py can slice in place."""
        path = self.root / fname
        if path not in self._handles:
            if not path.is_file():
                raise RuntimeError(f"missing ERA5-Land stack: {path}")
            import h5py
            if not h5py.is_hdf5(str(path)):
                raise RuntimeError(
                    f"{path.name} is not MATLAB -v7.3; a whole-array read would be "
                    f"tens of GB. Re-save it with '-v7.3' to allow lazy slicing.")
            self._handles[path] = h5py.File(path, "r")
        fh = self._handles[path]
        for name in candidates:
            if name in fh:
                return fh[name]
        raise RuntimeError(f"{path.name}: none of {candidates} present; has "
                           f"{list(fh.keys())[:8]}")

    @staticmethod
    def _series_at(dset, i: int, j: int, order: str) -> np.ndarray:
        """One pixel's time series. h5py reverses MATLAB's dimension order."""
        # MATLAB (a, b, time) is stored by h5py as (time, b, a).
        a, b = (i, j) if order == "latlon" else (j, i)
        if dset.ndim != 3:
            raise RuntimeError(f"expected a 3-D stack, got shape {dset.shape}")
        return np.asarray(dset[:, b - 1, a - 1], dtype=float).ravel()

    # -------------------------------------------------------------- extract
    def pixel_series(self, lat: float, lon: float) -> dict:
        """Monthly series at the posting nearest (lat, lon).

        Mirrors extract_pixel_climate_series: SPEI12/SPI* on the shifted grid
        indexed [lat, lon]; SPEI3/SPEI6 and tp/tas/pet/ssrd on the native grid
        indexed [lon, lat].
        """
        self._load_grid()
        i_s, j_s = self._ref_latlon.index(lat, lon)
        i_n, j_n = self._ref_lonlat.index(lat, lon)

        out: dict = {"si": {}, "si_time": {}}
        for name in SI_ORDER:
            fname, var, order, start = DROUGHT_FILES[name]
            dset = self._dataset(fname, (var,))
            i, j = (i_s, j_s) if order == "latlon" else (i_n, j_n)
            s = self._series_at(dset, i, j, order)
            out["si"][name] = s
            out["si_time"][name] = month_axis(start, s.size)

        for key, (fname, cands, order) in CLIMATE_FILES.items():
            dset = self._dataset(fname, cands)
            s = self._series_at(dset, i_n, j_n, order)
            out[key] = s
        n = min(out[k].size for k in CLIMATE_FILES)
        out["climate_time"] = month_axis(CLIMATE_START, n)
        out["lat_index"], out["lon_index"] = i_s, j_s
        return out

    # ------------------------------------------------------------ predictors
    def predictor_row(self, series: dict, sample: date,
                      wanted: list[str] | None = None) -> dict[str, float]:
        """Predictors for one sampling date (build_climate_predictor_row).

        Returns {name: value}; `wanted` restricts the work to the handful of
        predictors a fit actually selected.
        """
        want = set(wanted) if wanted is not None else set(PREDICTOR_NAMES)
        unknown = want - set(NAME_INDEX)
        if unknown:
            raise RuntimeError(f"not ERA5-derivable predictor(s): {sorted(unknown)}")
        row: dict[str, float] = {n: math.nan for n in want}
        key = ym_key(sample.year, sample.month)

        for name in SI_ORDER:
            fam = [n for n in want if n == name or n.startswith(name + "-")]
            if not fam:
                continue
            ts, axis = series["si"][name], series["si_time"][name]
            hit = np.nonzero(axis == key)[0]
            if not hit.size:
                continue
            jj = int(hit[0])
            for n in fam:
                if n == name:
                    row[n] = ts[jj]
                elif n.endswith("-sev"):
                    lo = max(0, jj - 11)
                    row[n] = -1.0 * np.nansum(si_to_si_droughts(ts[lo:jj + 1]))
                else:
                    lag = int(n.rsplit("-", 1)[1].removesuffix("mo"))
                    src = jj - lag
                    if 0 <= src < ts.size:
                        row[n] = ts[src]

        axis = series["climate_time"]
        hit = np.nonzero(axis == key)[0]
        if not hit.size or int(hit[0]) - 11 < 0:
            return row                      # MATLAB returns the partial row here too
        jj = int(hit[0])
        win = slice(jj - 11, jj + 1)
        tp, tas = series["tp"], series["tas"]
        pet, ssrd = -1.0 * series["pet"], series["ssrd"]

        for n, v in (("MAR", np.nanmean(tp[win])), ("MAT", np.nanmean(tas[win])),
                     ("PET", np.nanmean(pet[win])), ("SSRD", np.nanmean(ssrd[win]))):
            if n in want:
                row[n] = float(v)
        for lag in range(12):
            src = jj - lag
            for pre, arr in (("P", tp), ("T", tas), ("PET", pet), ("SSRD", ssrd)):
                n = f"{pre} - {lag}mo"
                if n in want:
                    row[n] = float(arr[src])

        # Seasonal means over the same 12-month window, minus the previous year's
        # occurrence of the sampling season (drop_prev_peak_season).
        idx = np.arange(jj, jj - 12, -1)
        months = [(axis[k] % 100) for k in idx]
        yrs = [(axis[k] // 100) for k in idx]
        labels = [SEASON_BY_MONTH[m - 1] for m in months]
        peak = SEASON_BY_MONTH[sample.month - 1]
        keep = [t for t, (lab, yr) in enumerate(zip(labels, yrs))
                if not (lab == peak and yr == sample.year - 1)]
        idx, labels = idx[keep], [labels[t] for t in keep]

        season_of = {"Summer": "JJA", "Spring": "MAM", "Fall": "SON", "Winter": "DJF"}
        for pre, arr in (("P", tp), ("T", tas), ("PET", pet), ("SSRD", ssrd)):
            for season, code in season_of.items():
                n = f"{pre} - {season}"
                if n not in want:
                    continue
                sel = [k for k, lab in zip(idx, labels) if lab == code]
                if sel:
                    row[n] = float(np.nanmean(arr[sel]))
        return row

    def close(self) -> None:
        for fh in self._handles.values():
            try:
                fh.close()
            except Exception:                                      # noqa: BLE001
                pass
        self._handles.clear()


class _Postings:
    """MATLAB georefpostings index lookup.

    georefpostings spans [min, max] with N postings, so the step is
    (max-min)/(N-1) and the posting is the nearest one. ColumnsStartFrom decides
    whether row 1 is the northernmost latitude, RowsStartFrom whether column 1 is
    the westernmost longitude -- both taken from the order of the source vectors,
    the same test load_climate_inputs makes.
    """

    def __init__(self, lat: np.ndarray, lon: np.ndarray):
        self.nlat, self.nlon = lat.size, lon.size
        self.latmin, self.latmax = float(lat.min()), float(lat.max())
        self.lonmin, self.lonmax = float(lon.min()), float(lon.max())
        self.from_north = lat[0] >= lat[-1]      # colsFrom = 'north' when descending
        self.from_west = lon[0] <= lon[-1]       # rowsFrom = 'west' when ascending
        self.dlat = (self.latmax - self.latmin) / (self.nlat - 1)
        self.dlon = (self.lonmax - self.lonmin) / (self.nlon - 1)

    def index(self, lat: float, lon: float) -> tuple[int, int]:
        # A -180..180 station longitude on a 0..360 reference has to be wrapped,
        # which is what geographicToDiscrete does internally.
        if lon < self.lonmin - 0.5 * self.dlon:
            lon += 360.0
        elif lon > self.lonmax + 0.5 * self.dlon:
            lon -= 360.0
        i = (self.latmax - lat) / self.dlat if self.from_north else (lat - self.latmin) / self.dlat
        j = (lon - self.lonmin) / self.dlon if self.from_west else (self.lonmax - lon) / self.dlon
        i, j = int(round(i)) + 1, int(round(j)) + 1
        if not (1 <= i <= self.nlat and 1 <= j <= self.nlon):
            raise RuntimeError(f"({lat}, {lon}) falls outside the ERA5-Land grid")
        return i, j

    def lat_pos(self, i: int) -> int:
        return i - 1 if self.from_north else self.nlat - i

    def lon_pos(self, j: int) -> int:
        return j - 1 if self.from_west else self.nlon - j


def doy_to_date(year: int, doy: float) -> date:
    """doy_to_year_datetime: clamp into the year, then offset from 1 January."""
    if not math.isfinite(doy):
        raise ValueError("non-finite DOY")
    leap = (year % 4 == 0 and year % 100 != 0) or year % 400 == 0
    d = max(1, min(365 + int(leap), int(round(doy))))
    return date(year, 1, 1) + timedelta(days=d - 1)


# --------------------------------------------------------------- verification

def verify(table: Path, root: Path, n_rows: int, tol: float) -> int:
    """Recompute predictors for rows that ARE in a preprocessed table and compare.

    This is the only end-to-end check that the port -- grids, axis order, time
    axes, sign conventions -- reproduces what the MATLAB wrote. Run it on the
    cluster before any fallback value is used in anger.
    """
    with open(table, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        print(f"{table}: no rows", file=sys.stderr)
        return 1
    have = [n for n in PREDICTOR_NAMES if n in rows[0]]
    print(f"table   : {table}")
    print(f"rows    : {len(rows)}, {len(have)}/146 predictor columns present")
    missing = [n for n in PREDICTOR_NAMES if n not in rows[0]]
    if missing:
        print(f"  ! absent from the table: {', '.join(missing[:6])}"
              f"{f' (+{len(missing) - 6} more)' if len(missing) > 6 else ''}")

    step = max(1, len(rows) // n_rows)
    sample = rows[::step][:n_rows]
    store = Era5Monthly(root)
    worst: list[tuple[float, str, float, float, str]] = []
    checked = nan_both = 0
    for r in sample:
        lat, lon = float(r["lat"]), float(r["lon"])
        d = doy_to_date(int(r["Year"]), float(r["DOY"]))
        got = store.predictor_row(store.pixel_series(lat, lon), d, have)
        for n in have:
            try:
                ref = float(r[n])
            except (TypeError, ValueError):
                ref = math.nan
            mine = got[n]
            if not math.isfinite(ref) and not math.isfinite(mine):
                nan_both += 1
                continue
            checked += 1
            scale = max(abs(ref), abs(mine), 1e-12)
            rel = abs(ref - mine) / scale
            worst.append((rel, n, ref, mine, f"{r['Pixel ID']}/{r['Year']}"))
    store.close()

    worst.sort(reverse=True)
    bad = [w for w in worst if w[0] > tol or not math.isfinite(w[0])]
    print(f"compared: {checked} value(s) over {len(sample)} row(s) "
          f"({nan_both} NaN in both)")
    print(f"mismatch: {len(bad)} above rel tol {tol:g}")

    # Split the verdict by predictor, and separate the ported predictors from the
    # one inferred function. '-sev' failing says nothing about the grids, time axes
    # or sign conventions -- it is a different, isolated problem -- so reporting a
    # single pass/fail hides the result that actually matters.
    sev_bad = [w for w in bad if w[1].endswith("-sev")]
    other_bad = [w for w in bad if not w[1].endswith("-sev")]
    n_sev = sum(1 for n in have if n.endswith("-sev")) * len(sample)
    print(f"  ported predictors : {len(other_bad)} mismatch(es) of "
          f"{checked - n_sev} compared")
    print(f"  '-sev' (inferred) : {len(sev_bad)} mismatch(es) of {n_sev} compared")

    per = {}
    for rel, n, _, _, _ in bad:
        p = per.setdefault(n, [0, 0.0])
        p[0] += 1
        p[1] = max(p[1], rel)
    if per:
        print("\n  mismatches by predictor:")
        for n, (cnt, mx) in sorted(per.items(), key=lambda kv: -kv[1][0]):
            print(f"    {n:<18} {cnt:>4} row(s)   worst rel {mx:.2e}")

    show = other_bad[:10] or worst[:6]
    print("\n  worst values:")
    for rel, n, ref, mine, where in show:
        flag = "  <-- MISMATCH" if rel > tol else ""
        print(f"    {n:<18} table {ref:>16.6g}  recomputed {mine:>16.6g}  "
              f"rel {rel:.2e}  [{where}]{flag}")

    if other_bad:
        fams = sorted({n.split(" - ")[0].split("-")[0] for _, n, _, _, _ in other_bad})
        print(f"\nFAILED on ported predictors. Families affected: {', '.join(fams)}")
        return 1
    if sev_bad:
        print("\nPORT VERIFIED for every predictor except '-sev'. The grids, time axes,\n"
              "lag indexing, seasonal windows and the PET sign all reproduce the table.\n"
              "'-sev' needs the real SI_to_SIdroughts from\n"
              "/vol_efthymios/NFS07/dd1136/functions/ -- the fallback refuses to fill a\n"
              "group whose fit selected one rather than substituting a guess.")
        return 2
    print("\nOK -- the port reproduces the table.")
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("verify", help="recompute predictors for existing table rows")
    v.add_argument("table", type=Path)
    v.add_argument("--era5-root", type=Path, default=DEFAULT_ERA5_ROOT)
    v.add_argument("--rows", type=int, default=25)
    v.add_argument("--tol", type=float, default=1e-6)
    n = sub.add_parser("names", help="print the 146 predictor names in column order")
    n.add_argument("--check", type=Path, default=None,
                   help="compare against a preprocessed table's header")
    a = p.parse_args()

    if a.cmd == "names":
        if a.check:
            with open(a.check, newline="", encoding="utf-8-sig") as fh:
                header = [c.strip() for c in next(csv.reader(fh))]
            ok = header[:146] == PREDICTOR_NAMES
            print(f"table header first 146 columns match the port: {ok}")
            if not ok:
                for i, (h, m) in enumerate(zip(header[:146], PREDICTOR_NAMES)):
                    if h != m:
                        print(f"  col {i + 1}: table {h!r} != port {m!r}")
                return 1
            return 0
        for i, n_ in enumerate(PREDICTOR_NAMES, 1):
            print(f"{i:3d}  {n_}")
        return 0
    return verify(a.table, a.era5_root, a.rows, a.tol)


if __name__ == "__main__":
    sys.exit(main())
