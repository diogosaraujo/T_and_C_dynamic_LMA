#!/usr/bin/env python3
"""Build per-station time-varying LMA series from the PLSR temporal-CV output.

Two series are produced for every AmeriFlux station, both in g/m2:

    observed  <- Y_plot_abs      the iLMA observations that entered the fit
    modelled  <- yfit_plot_abs   the PLSR reconstruction (pixel climatology + anomaly)

Each station takes its nearest pixel -- ERA5-Land 0.1 deg cells, so the nearest is
the cell the tower sits in -- and that pixel's annual series is written out.

The absolute values come back from the fit's own baselines, none of which are
re-estimated here: predictors are demeaned by the PIXEL mean, normalised by the
ECOREGION mean and standard deviation, run through beta, then rescaled by sigma_Y
and re-centred on the pixel mean. Where the nearest pixel was never mapped as the
station's forest type the fit holds no mean for it, so the ecoregion mean stands in
for the pixel mean; the predictors still come from that pixel, so only the level is
substituted, and every row records which baseline applied.

Inputs (on the cluster, outside this repo):

    <plsr-root>/eco<ii>/time/PLSR_predictions_eco<ii>_<forest>_oofcv.mat
    <plsr-root>/eco<ii>/time/PLSR_fitting_coeff_eco<ii>_<forest>_oofcv_TEMPORAL.mat
    <predictor-root>/LMA_ecoregion_no<ii>.csv    -- Pixel ID, lat, lon, LU, time, predictors

The predictor table is the one the MATLAB fit itself read (resolve_input_table:
<OutRoot>/PLSR_inputs_pixel_climatology_DOY/LMA/). It is also the coordinate
source, so pixel IDs, coordinates and predictors always come from a single
vintage. The older ecoregion_no<ii>.csv at the ecoregions root lacks the SSRD
columns every fit selected and may index pixels differently; it is used only if
nothing else is present, and then it is flagged in the log.

Output:

    $TC_INPUT_DATA/lma/<station>/<station>_LMA_observed.csv
    $TC_INPUT_DATA/lma/<station>/<station>_LMA_modelled.csv
    $TC_INPUT_DATA/lma/lma_manifest.csv

The SLA conversion is deliberately NOT applied here. LMA -> SLA depends on the
carbon fraction f_C, which is a modelling decision, and it belongs with the .mat
build so a change to f_C does not mean re-reading the PLSR output.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

import era5_predictors

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SITE_LISTS = [
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv",
]
DEFAULT_EXCLUDED = Path(__file__).resolve().parent / "excluded_stations.csv"
INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA", "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
DEFAULT_OUT = INPUT_ROOT / "lma"

ECOREGION_ROOT = Path("/vol_efthymios/NFS07/dd1136/ecoregions")
DEFAULT_PLSR_ROOT = ECOREGION_ROOT / "PLSR_temporal_cv_pixel_climatology_DOY" / "LMA"
# opts.InputDir in the MATLAB pipeline: fullfile(OutRoot, this, target_name).
PLSR_INPUT_SUBDIR = Path("PLSR_inputs_pixel_climatology_DOY")
DEFAULT_PREDICTOR_ROOT = ECOREGION_ROOT / PLSR_INPUT_SUBDIR / "LMA"

# PLSR_TEMPORAL_CV_OOF_METRICS_PIXEL_DOY_CORE.m spells the first class 'deciduos'
# (missing a 'u'). Try the typo first, then the correct spelling, so the script
# keeps working if the upstream name is ever fixed.
FOREST_LABELS = {
    "deciduous": ["deciduos", "deciduous"],
    "evergreen": ["evergreen"],
    "mixed": ["mixed"],
}
# lu_value = 40 + lu_id, with lu_id 1/2/3 = deciduous/evergreen/mixed. The ecoregion
# table holds several land-use classes, so reconstruction must filter to the right one.
LU_BASE = 40
FOREST_LU = {"deciduous": 1, "evergreen": 2, "mixed": 3}

# Physically implausible LMA is a sign of a unit or column mix-up, not a real leaf.
LMA_MIN, LMA_MAX = 10.0, 600.0


# --------------------------------------------------------------------------- io

def read_excluded(path: Path | None) -> dict[str, str]:
    """Station -> reason. An unset or absent file simply excludes nothing."""
    out: dict[str, str] = {}
    if not path or not Path(path).is_file():   # '' becomes Path('.'), a directory
        return out
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            sid = (row.get("station_id") or "").strip()
            if sid:
                out[sid] = (row.get("reason") or "excluded").strip()
    return out


def read_stations(paths, wanted: set[str] | None) -> list[dict]:
    stations, seen = [], set()
    for path in paths:
        if not Path(path).is_file():
            print(f"  ! site list not found: {path}", file=sys.stderr)
            continue
        with open(path, newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                sid = (row.get("StationID") or "").strip()
                if not sid or sid in seen:
                    continue
                if wanted and sid not in wanted:
                    continue
                try:
                    lat, lon = float(row["Lat"]), float(row["Lon"])
                    eco = int(row["ECO_IDX"])
                except (KeyError, TypeError, ValueError):
                    print(f"  ! {sid}: unusable Lat/Lon/ECO_IDX, skipped", file=sys.stderr)
                    continue
                seen.add(sid)
                stations.append({
                    "station_id": sid,
                    "lat": lat,
                    "lon": lon,
                    "eco_idx": eco,
                    "forest_type": (row.get("ForestType") or "").strip().lower(),
                    "eco_name": (row.get("US_L3NAME") or "").strip(),
                    "plsr_q2": (row.get("PLSR_TemporalQ2") or "").strip(),
                })
    return stations


HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"


def _is_hdf5(path: Path) -> bool:
    """True for a MATLAB -v7.3 file.

    A v7.3 .mat is HDF5 behind a 512-byte userblock holding the ASCII
    'MATLAB 7.3 MAT-file...' header, so the signature is NOT at offset 0 --
    checking only there misreports every real file as legacy v7 and hands it to
    scipy.io.loadmat, which refuses it. h5py.is_hdf5 understands userblocks
    (HDF5 probes offset 0 then successive powers of two from 512); the manual
    check is the fallback for when h5py is not installed.
    """
    try:
        import h5py
        return bool(h5py.is_hdf5(str(path)))
    except ImportError:
        pass
    with open(path, "rb") as fh:
        if fh.read(len(HDF5_SIGNATURE)) == HDF5_SIGNATURE:
            return True
        offset = 512
        while offset <= 8192:
            fh.seek(offset)
            if fh.read(len(HDF5_SIGNATURE)) == HDF5_SIGNATURE:
                return True
            offset *= 2
    return False


def load_prediction_mat(path: Path) -> dict[str, np.ndarray]:
    """Read the numeric fields of PLSR_predictions_*.mat.

    The pipeline saves with '-v7.3', i.e. HDF5, which scipy.io.loadmat cannot
    read; older MATLAB versions would produce a v7 file, so both are handled.
    'time' is skipped -- it is a MATLAB datetime object and awkward to decode,
    while 'time_yrs' carries the same information as a plain number.
    """
    wanted = ("Y_plot_abs", "yfit_plot_abs", "time_yrs", "pixel_id", "plot_has_observed")
    out: dict[str, np.ndarray] = {}
    if _is_hdf5(path):
        import h5py  # imported lazily so the rest of the script runs without it
        with h5py.File(path, "r") as fh:
            for key in wanted:
                if key in fh:
                    out[key] = np.asarray(fh[key][()]).squeeze().ravel()
    else:
        from scipy.io import loadmat
        raw = loadmat(path, squeeze_me=True)
        for key in wanted:
            if key in raw:
                out[key] = np.atleast_1d(raw[key]).ravel()
    missing = [k for k in ("Y_plot_abs", "yfit_plot_abs", "time_yrs", "pixel_id") if k not in out]
    if missing:
        raise RuntimeError(f"{path.name}: missing field(s) {', '.join(missing)}")
    n = len(out["time_yrs"])
    bad = {k: len(v) for k, v in out.items() if len(v) != n}
    if bad:
        raise RuntimeError(f"{path.name}: field length mismatch against time_yrs={n}: {bad}")
    return out


def _mat_array(fh, key, hdf5: bool) -> np.ndarray:
    """Read one numeric field, undoing HDF5's dimension reversal.

    MATLAB stores column-major, so an (m, n) matrix comes back from h5py as
    (n, m). Vectors are unaffected once ravelled; genuine matrices are not.
    """
    a = np.asarray(fh[key][()] if hdf5 else fh[key], dtype=float)
    if hdf5 and a.ndim == 2 and min(a.shape) > 1:
        a = a.T
    return a


def _mat_strings(fh, key, hdf5: bool) -> list[str]:
    """Read a MATLAB cellstr. Under HDF5 each cell is a ref to uint16 char codes."""
    if not hdf5:
        return [str(s).strip() for s in np.atleast_1d(fh[key]).ravel()]
    out = []
    for ref in np.asarray(fh[key][()]).ravel():
        out.append("".join(chr(int(c)) for c in np.asarray(fh[ref][()]).ravel()))
    return out


def read_fit_mat(path: Path) -> dict:
    """Read PLSR_fitting_coeff_*_TEMPORAL.mat -- everything predict_with_fit needs."""
    num = ("beta", "mu_X", "sigma_X", "mu_Y", "sigma_Y",
           "mu_X_pix_final", "mu_Y_pix_final", "uniq_pix_final", "r2", "ncomp")
    out: dict = {}
    if _is_hdf5(path):
        import h5py
        with h5py.File(path, "r") as fh:
            for k in num:
                if k in fh:
                    out[k] = _mat_array(fh, k, True)
            for k in ("predictor_name_VIP", "predictor_name_NS"):
                if k in fh:
                    out["predictors"] = _mat_strings(fh, k, True)
                    break
    else:
        from scipy.io import loadmat
        raw = loadmat(path, squeeze_me=True)
        for k in num:
            if k in raw:
                out[k] = np.atleast_1d(np.asarray(raw[k], dtype=float))
        for k in ("predictor_name_VIP", "predictor_name_NS"):
            if k in raw:
                out["predictors"] = _mat_strings(raw, k, False)
                break
    missing = [k for k in ("beta", "mu_X", "sigma_X", "mu_Y", "sigma_Y", "mu_X_pix_final",
                           "mu_Y_pix_final", "uniq_pix_final") if k not in out]
    if missing or "predictors" not in out:
        raise RuntimeError(f"{path.name}: missing {', '.join(missing + ['predictors'] * ('predictors' not in out))}")

    out["beta"] = out["beta"].ravel()
    out["uniq_pix_final"] = out["uniq_pix_final"].ravel().astype(int)
    out["mu_Y_pix_final"] = out["mu_Y_pix_final"].ravel()
    out["mu_X"] = out["mu_X"].ravel()
    out["sigma_X"] = out["sigma_X"].ravel()
    out["mu_Y"] = float(np.ravel(out["mu_Y"])[0])
    out["sigma_Y"] = float(np.ravel(out["sigma_Y"])[0])
    npred, npix = len(out["predictors"]), len(out["uniq_pix_final"])
    mx = np.atleast_2d(out["mu_X_pix_final"])
    if mx.shape != (npix, npred):
        if mx.shape == (npred, npix):
            mx = mx.T
        else:
            raise RuntimeError(f"{path.name}: mu_X_pix_final is {mx.shape}, expected "
                               f"({npix}, {npred}) = (pixels, predictors)")
    out["mu_X_pix_final"] = mx
    if len(out["beta"]) != npred + 1:
        raise RuntimeError(f"{path.name}: beta has {len(out['beta'])} terms, expected "
                           f"{npred + 1} (intercept + {npred} predictors)")

    # Ecoregion means in RAW units, for pixels the fit never saw. mu_X/mu_Y are NOT
    # these: the fit demeans per pixel before standardising, so mu_X and mu_Y are
    # means of the anomalies and come out at 0 (eco1 evergreen: mu_X = [0 -0 -0 0 -0],
    # mu_Y = 0.000000). Substituting those would put a station at 0 g/m2. The
    # ecoregion mean of the variable itself is the mean of the per-pixel means the
    # fit stored -- 88.27 g/m2 for eco1 evergreen.
    out["mu_X_eco"] = out["mu_X_pix_final"].mean(axis=0)
    out["mu_Y_eco"] = float(np.mean(out["mu_Y_pix_final"]))
    return out


def predict_from_raw(fit: dict, X_raw: np.ndarray, mu_X_row: np.ndarray,
                     mu_Y_row: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """predict_with_fit, given raw predictors and the baselines to demean by.

    Returns (absolute values, mask of non-finite predictors). Non-finite
    predictors are zeroed after standardising, i.e. replaced by the mean, exactly
    as the MATLAB does -- which quietly pulls a prediction toward its baseline, so
    the mask is reported rather than discarded.
    """
    X = (X_raw - mu_X_row - fit["mu_X"]) / fit["sigma_X"]
    bad = ~np.isfinite(X)
    X = np.where(bad, 0.0, X)
    yfit = fit["beta"][0] + X @ fit["beta"][1:]
    return yfit * fit["sigma_Y"] + fit["mu_Y"] + mu_Y_row, bad


def reconstruct_modelled(fit: dict, eco_csv: Path, lu_value: int) -> tuple[dict, dict]:
    """Re-apply the fitted model to EVERY row of the predictor table.

    The pipeline only predicted rows that survived rmmissing on the response, so
    yfit_plot_abs inherits the observations' gaps -- 26 of 37 years for eco1
    evergreen. The predictors are ~99% complete in every year, so evaluating the
    same coefficients over the full table fills the rest. This reproduces
    predict_with_fit, demeaning by the pixel mean and then standardising by the
    ecoregion mean and standard deviation:

        X_anom   = X_raw - mu_X_pix_final(pix,:)
        X        = (X_anom - mu_X) ./ sigma_X ;  non-finite -> 0
        yfit_abs = ([1 X]*beta) .* sigma_Y + mu_Y + mu_Y_pix_final(pix)

    Every baseline comes from the fit file untouched; nothing is re-estimated here.

    TWO BASELINES. A pixel in uniq_pix_final has its own mean, and that mean was
    computed by the fit over the years the pixel carried lu_value ONLY -- the fit
    stores LU and was run on those rows -- so a pixel mapped evergreen 1985-2002 and
    mixed after contributes evergreen-period statistics, which is what a pixel that
    changed class should contribute. Its series still spans every year, since the
    predictors are the same climate record however the map labelled the pixel.

    A pixel the fit never saw (never lu_value) has no mean of its own, so the
    ECOREGION mean stands in for the pixel mean -- mu_X_eco / mu_Y_eco, the means of
    the per-pixel means the fit stored. Predictors still come from that pixel: only
    the baseline is substituted, so the station keeps its local climate signal and
    loses only the local level.

    Returns (series, diagnostics). series has 'years', 'pixels', 'values' and
    'pixel_mean' (True where the pixel's own mean was used).
    """
    preds = fit["predictors"]
    pix_row = {int(p): i for i, p in enumerate(fit["uniq_pix_final"])}

    years, pixels, rows = [], [], []
    class_years: dict[int, set[int]] = {}
    seen_class: dict[int, set[str]] = {}
    doy_by_pixel: dict[int, float] = {}
    with open(eco_csv, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        cols = {c.strip(): c for c in (reader.fieldnames or [])}
        missing = [p for p in preds if p not in cols]
        if missing:
            raise RuntimeError(f"{eco_csv.name}: predictor column(s) absent: {', '.join(missing)}")
        pid_col = cols.get("Pixel ID") or cols.get("Pixel_ID")
        lu_col = cols.get("LU")
        # 'time' is a date string ('15-Jul-1985'); 'Year' is the plain year. Prefer it.
        t_col = cols.get("Year") or cols.get("year") or cols.get("time")
        doy_col = cols.get("DOY")
        if not (pid_col and lu_col and t_col):
            raise RuntimeError(f"{eco_csv.name}: need Pixel ID, LU and Year/time columns")
        for r in reader:
            try:
                pid = int(float(r[pid_col]))
                yr = int(str(r[t_col]).strip()[-4:])
                lu = int(float(r[lu_col]))
            except (TypeError, ValueError):
                continue
            seen_class.setdefault(pid, set()).add(str(r[lu_col]).strip())
            if lu == lu_value:
                class_years.setdefault(pid, set()).add(yr)
                # DOY is a per-(pixel, LU) constant -- it comes from a climatology
                # table, not from the year -- so one value per pixel is enough to
                # re-sample a year the table never emitted.
                if doy_col and pid not in doy_by_pixel:
                    try:
                        doy_by_pixel[pid] = float(r[doy_col])
                    except (TypeError, ValueError):
                        pass
            vals = []
            for p in preds:
                try:
                    vals.append(float(r[cols[p]]))
                except (TypeError, ValueError):
                    vals.append(np.nan)
            years.append(yr)
            pixels.append(pid)
            rows.append(vals)

    diag = {"n_rows": len(rows),
            "n_pixels": len(seen_class),
            "n_pixels_in_fit": len(pix_row),
            "n_pixels_class_switched": sum(1 for p, v in seen_class.items()
                                           if len(v) > 1 and p in pix_row),
            "class_years": {p: len(v) for p, v in class_years.items()},
            "doy_by_pixel": doy_by_pixel,
            "eco_doy": (float(np.median(list(doy_by_pixel.values())))
                        if doy_by_pixel else float("nan"))}
    if not rows:
        return {"years": np.array([]), "pixels": np.array([]),
                "values": np.array([]), "pixel_mean": np.array([], bool)}, diag

    X_raw = np.asarray(rows, dtype=float)
    years, pixels = np.asarray(years, int), np.asarray(pixels, int)
    in_fit = np.array([p in pix_row for p in pixels])
    idx = np.array([pix_row.get(p, 0) for p in pixels])

    mu_X_row = np.where(in_fit[:, None], fit["mu_X_pix_final"][idx, :], fit["mu_X_eco"])
    mu_Y_row = np.where(in_fit, fit["mu_Y_pix_final"][idx], fit["mu_Y_eco"])

    vals, bad = predict_from_raw(fit, X_raw, mu_X_row, mu_Y_row)
    diag["n_rows_with_imputed_predictor"] = int(np.any(bad, axis=1).sum())
    diag["n_rows_all_predictors_imputed"] = int(np.all(bad, axis=1).sum())

    # A pixel whose predictors are non-finite in EVERY year gets every predictor
    # zeroed, so its "prediction" is a constant equal to its baseline -- a flat series
    # that would enter the experiment as a dynamic input while carrying no dynamics at
    # all. eco1 has three (166, 182, 513). Drop them so no station can match one; the
    # next-nearest pixel is then used, and the run log names them.
    all_bad = np.all(bad, axis=1)
    dead = sorted({int(p) for p in np.unique(pixels[all_bad])
                   if np.all(all_bad[pixels == p])})
    diag["dead_pixels"] = dead
    if dead:
        keep = ~np.isin(pixels, dead)
        years, pixels, vals, in_fit = years[keep], pixels[keep], vals[keep], in_fit[keep]
        diag["n_rows"] = int(keep.sum())

    diag["n_pixels_kept"] = int(len(np.unique(pixels)))
    diag["n_rows_ecoregion_baseline"] = int((~in_fit).sum())
    return {"years": years, "pixels": pixels, "values": vals, "pixel_mean": in_fit}, diag


def read_pixel_coords(path: Path) -> dict[int, tuple[float, float]]:
    """Pixel ID -> (lat, lon). The table repeats each pixel once per year."""
    coords: dict[int, tuple[float, float]] = {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        cols = {c.strip().lower(): c for c in (reader.fieldnames or [])}
        pid_col = cols.get("pixel id") or cols.get("pixel_id") or cols.get("pixelid")
        lat_col = cols.get("lat") or cols.get("latitude")
        lon_col = cols.get("lon") or cols.get("longitude")
        if not (pid_col and lat_col and lon_col):
            raise RuntimeError(
                f"{path.name}: need Pixel ID/lat/lon columns, found {reader.fieldnames}")
        for row in reader:
            try:
                pid = int(float(row[pid_col]))
                lat, lon = float(row[lat_col]), float(row[lon_col])
            except (TypeError, ValueError):
                continue
            if math.isfinite(lat) and math.isfinite(lon):
                coords.setdefault(pid, (lat, lon))
    return coords


# ---------------------------------------------------------------------- helpers

def haversine_km(lat0: float, lon0: float, lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    r = 6371.0088
    p0, p1 = math.radians(lat0), np.radians(lat)
    dp, dl = p1 - p0, np.radians(lon - lon0)
    a = np.sin(dp / 2) ** 2 + math.cos(p0) * np.cos(p1) * np.sin(dl / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def annual_median(years: np.ndarray, values: np.ndarray) -> dict[int, tuple[float, int]]:
    """year -> (median, n contributing rows), ignoring NaN.

    The source tables are annual, so this is normally a no-op; it exists so a
    pixel appearing twice in one year cannot silently pick an arbitrary value.
    The PLSR pipeline aggregates the same way (add_oof_year_medians).
    """
    out: dict[int, tuple[float, int]] = {}
    finite = np.isfinite(values)
    for yr in np.unique(years[finite]):
        v = values[(years == yr) & finite]
        if v.size:
            out[int(yr)] = (float(np.median(v)), int(v.size))
    return out


def resolve_prediction_file(plsr_root: Path, eco: int, forest: str) -> Path | None:
    for label in FOREST_LABELS.get(forest, [forest]):
        cand = plsr_root / f"eco{eco}" / "time" / f"PLSR_predictions_eco{eco}_{label}_oofcv.mat"
        if cand.is_file():
            return cand
    return None


def resolve_tables(predictor_root: Path, ecoregion_root: Path, eco: int,
                   target: str = "LMA") -> list[Path]:
    """Candidate pixel/predictor tables for one ecoregion, best first.

    Mirrors resolve_input_table in the MATLAB pipeline, which reads from
    OutRoot/PLSR_inputs_pixel_climatology_DOY/<TARGET>/ and tries
    '<TARGET>_ecoregion_no<ii>.csv' then 'ecoregion_no<ii>.csv'. Those two are the
    real inputs; the copies at the ecoregions root are an older vintage kept only
    as a last resort (see is_legacy_table).
    """
    cands = [predictor_root / f"{target}_ecoregion_no{eco}.csv",
             predictor_root / f"ecoregion_no{eco}.csv",
             ecoregion_root / f"{target}_ecoregion_no{eco}.csv",
             ecoregion_root / f"ecoregion_no{eco}.csv"]
    seen, out = set(), []
    for c in cands:
        if c not in seen and c.is_file():
            seen.add(c)
            out.append(c)
    return out


def is_legacy_table(path: Path, predictor_root: Path) -> bool:
    """True for the older tables outside the directory the fit read.

    Those lack the SSRD columns every fit selected, so they cannot support
    --reconstruct, and there is no guarantee they index pixels the same way --
    which matters because the pixel IDs in the .mat are resolved against them.
    """
    return path.parent.resolve() != predictor_root.resolve()


def table_columns(path: Path) -> list[str]:
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return [c.strip() for c in next(csv.reader(fh), [])]


def missing_columns(path: Path, cols) -> list[str]:
    have = set(table_columns(path))
    return [c for c in cols if c not in have]


def has_coord_columns(path: Path) -> bool:
    cols = {c.lower() for c in table_columns(path)}
    return (bool(cols & {"pixel id", "pixel_id", "pixelid"})
            and bool(cols & {"lat", "latitude"})
            and bool(cols & {"lon", "longitude"}))


def pick_predictor_table(cands: list[Path], predictors) -> tuple[Path, dict[str, list[str]]]:
    """First candidate carrying every predictor. Raises with per-file detail if none."""
    misses: dict[str, list[str]] = {}
    for c in cands:
        gap = missing_columns(c, predictors)
        if not gap:
            return c, misses
        misses[str(c)] = gap
    detail = "; ".join(f"{Path(k).name} lacks {', '.join(v[:3])}"
                       + (f" (+{len(v) - 3} more)" if len(v) > 3 else "")
                       for k, v in misses.items())
    raise RuntimeError(
        f"no table carries all {len(predictors)} predictor(s): {detail}. "
        f"Point --predictor-root at the directory the fit actually used "
        f"(PLSR_inputs_pixel_climatology_DOY/LMA).")


def resolve_fit_file(plsr_root: Path, eco: int, forest: str) -> Path | None:
    """The '_TEMPORAL' and '_TEMPORAL_ALLYEARS' files are the same struct saved twice."""
    for label in FOREST_LABELS.get(forest, [forest]):
        for suffix in ("TEMPORAL", "TEMPORAL_ALLYEARS"):
            cand = (plsr_root / f"eco{eco}" / "time" /
                    f"PLSR_fitting_coeff_eco{eco}_{label}_oofcv_{suffix}.mat")
            if cand.is_file():
                return cand
    return None


def era5_fill(store, fit: dict, rows: list[dict], lat: float, lon: float, doy: float,
              pix_index: int | None, tag: str) -> tuple[int, str]:
    """Fill a station's missing years from the ERA5-Land monthly stacks.

    The preprocessed table only emits a pixel-year whose dominant NLCD class was
    forest and which had a DOY-climatology key, so a year can be absent even
    though the ERA5-Land forcing behind it exists. This recomputes the predictors
    that the fit selected, at the same sampling DOY, and predicts with the same
    baseline the rest of the station's series uses -- the pixel mean when the fit
    has one (pix_index), otherwise the ecoregion mean.

    Returns (n filled, note).
    """
    missing = [r for r in rows if r["LMA_g_m2"] == ""]
    if not missing or not math.isfinite(doy):
        return 0, ""
    preds = fit["predictors"]
    sev = [p for p in preds if p.endswith("-sev")]
    if pix_index is None:
        mu_X_row, mu_Y_row = fit["mu_X_eco"], fit["mu_Y_eco"]
        label = "ecoregion_mean"
    else:
        mu_X_row, mu_Y_row = fit["mu_X_pix_final"][pix_index], fit["mu_Y_pix_final"][pix_index]
        label = "pixel_mean"

    series = store.pixel_series(lat, lon)
    n = 0
    for r in missing:
        try:
            when = era5_predictors.doy_to_date(int(r["year"]), doy)
            got = store.predictor_row(series, when, preds)
        except Exception as exc:                                   # noqa: BLE001
            print(f"      ! {tag} {r['year']}: {type(exc).__name__}: {exc}")
            continue
        X_raw = np.array([[got[p] for p in preds]], dtype=float)
        if not np.any(np.isfinite(X_raw)):
            continue
        val, _ = predict_from_raw(fit, X_raw, mu_X_row, mu_Y_row)
        r.update(LMA_g_m2=round(float(val[0]), 4), source=f"{label}_era5", n_values=1)
        n += 1
    note = ""
    if n:
        note = f"{n} year(s) filled from ERA5-Land"
        if sev:
            note += (f" using inferred SI_to_SIdroughts for {', '.join(sev)} "
                     f"-- verify against MATLAB")
    return n, note


def write_series(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["year", "LMA_g_m2", "source", "pixel_id", "n_values"])
        w.writeheader()
        w.writerows(rows)


# ------------------------------------------------------------------------- core

def build_station(st: dict, series: dict[str, dict], coords: dict[int, tuple[float, float]],
                  years_wanted: list[int], max_dist_km: float, fill_gaps: bool,
                  class_years: dict[int, int] | None = None) -> dict:
    """Return {'observed': rows, 'modelled': rows, **manifest fields} for one station.

    The station takes its NEAREST pixel outright -- these are ERA5-Land 0.1 deg cells,
    so the nearest is the one the tower sits in. No search for a better-classified
    pixel further away: a cell 20 km off is a different stand under a different
    climate, and pretending otherwise would be worse than losing the local level.

    What varies is the baseline. If that pixel is in the fit, its own mean is used.
    If it was never mapped as the station's forest type, the fit has no mean for it
    and the ecoregion mean stands in -- reconstruct_modelled has already made that
    substitution, and 'pixel_mean' records which applied.
    """
    m = series["modelled"]
    have = sorted({int(p) for p in np.unique(m["pixels"])} & set(coords)) if len(m["pixels"]) else []
    info = {"n_pixels": len(have), "pixel_id": "", "pixel_lat": "", "pixel_lon": "",
            "distance_km": "", "pixel_class_years": "", "baseline": "", "note": ""}

    if not have:
        info["note"] = "no pixel in the predictor table has coordinates"
        best = None
    else:
        plat = np.array([coords[p][0] for p in have])
        plon = np.array([coords[p][1] for p in have])
        d = haversine_km(st["lat"], st["lon"], plat, plon)
        k = int(np.argmin(d))
        best, dist = have[k], float(d[k])
        n_class = class_years.get(best, 0) if class_years is not None else ""
        info.update(pixel_id=best, pixel_lat=round(plat[k], 5), pixel_lon=round(plon[k], 5),
                    distance_km=round(dist, 3), pixel_class_years=n_class)
        if dist > max_dist_km:
            info["note"] = f"nearest pixel is {dist:.1f} km away (> {max_dist_km} km)"

    # Was the pixel's own mean available, or did the ecoregion mean stand in?
    if best is not None and len(m["pixels"]):
        sel = m["pixels"] == best
        used_pixel_mean = bool(sel.any() and m["pixel_mean"][sel][0])
        info["baseline"] = "pixel_mean" if used_pixel_mean else "ecoregion_mean"
        if not used_pixel_mean:
            info["note"] = ((info["note"] + "; ") if info["note"] else "") + \
                f"nearest pixel is never mapped {st['forest_type']}; ecoregion mean used"

    out = {**info}
    for kind in ("observed", "modelled"):
        k = series[kind]
        k_yrs = np.asarray(k["years"], dtype=int)
        k_pix = np.asarray(k["pixels"], dtype=int)
        vals = np.asarray(k["values"], dtype=float)
        base = k.get("pixel_mean")
        sel = k_pix == best if best is not None else np.zeros(len(k_yrs), bool)
        pix_series = annual_median(k_yrs[sel], vals[sel]) if best is not None else {}

        label = info["baseline"] or "pixel_mean"
        if base is not None and best is not None and sel.any():
            label = "pixel_mean" if base[sel][0] else "ecoregion_mean"
        source = label if pix_series else "none"
        rows = []
        for yr in years_wanted:
            if yr in pix_series:
                v, n = pix_series[yr]
                rows.append({"year": yr, "LMA_g_m2": round(v, 4), "source": label,
                             "pixel_id": best, "n_values": n})
            else:
                rows.append({"year": yr, "LMA_g_m2": "", "source": "missing",
                             "pixel_id": "", "n_values": 0})

        present = [r["LMA_g_m2"] for r in rows if r["LMA_g_m2"] != ""]
        out[kind] = rows
        out[f"{kind}_source"] = source if present else "none"
        out[f"{kind}_n_years"] = len(present)
        out[f"{kind}_mean"] = round(float(np.mean(present)), 3) if present else ""
        out[f"{kind}_min"] = round(float(np.min(present)), 3) if present else ""
        out[f"{kind}_max"] = round(float(np.max(present)), 3) if present else ""
        oor = [v for v in present if not (LMA_MIN <= v <= LMA_MAX)]
        if oor:
            out["note"] = (out["note"] + "; " if out["note"] else "") + \
                f"{len(oor)} {kind} value(s) outside {LMA_MIN}-{LMA_MAX} g/m2"
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--site-list", type=Path, action="append", default=None)
    p.add_argument("--stations", default=None, help="comma-separated StationIDs")
    p.add_argument("--exclude", default=None, help="comma-separated StationIDs to drop")
    p.add_argument("--exclude-file", type=Path, default=DEFAULT_EXCLUDED,
                   help="CSV of stations to drop (default: excluded_stations.csv; '' to disable)")
    p.add_argument("--plsr-root", type=Path, default=DEFAULT_PLSR_ROOT,
                   help="directory holding eco<ii>/time/PLSR_predictions_*.mat")
    p.add_argument("--ecoregion-root", type=Path, default=ECOREGION_ROOT,
                   help="PLSR output root; only used to derive --predictor-root and as a "
                        "last-resort source of the older ecoregion_no<ii>.csv")
    p.add_argument("--predictor-root", type=Path, default=None,
                   help="directory holding LMA_ecoregion_no<ii>.csv -- the tables the fit "
                        "read, and this script's source of pixel coordinates AND predictors. "
                        f"Default: <ecoregion-root>/{PLSR_INPUT_SUBDIR.as_posix()}/LMA")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--start-year", type=int, default=1985)
    p.add_argument("--end-year", type=int, default=2021)
    p.add_argument("--max-distance-km", type=float, default=50.0,
                   help="beyond this, fall back to the ecoregion median (default: 50)")
    p.add_argument("--fill-gaps", action="store_true",
                   help="(retained for compatibility; no effect since the fallback is now "
                        "a baseline substitution, not a separate median series)")
    p.add_argument("--reconstruct", action="store_true",
                   help="rebuild the modelled series by re-applying "
                        "PLSR_fitting_coeff_*_TEMPORAL.mat to the full predictor table, "
                        "instead of reading the gapped yfit_plot_abs. Uses only pixels in "
                        "uniq_pix_final, i.e. those that contributed to the fit.")
    p.add_argument("--era5-root", type=Path, default=era5_predictors.DEFAULT_ERA5_ROOT,
                   help="ERA5-Land monthly stacks, used to fill years the preprocessed "
                        "table never emitted (it only keeps forest pixel-years)")
    p.add_argument("--no-era5-fallback", action="store_true",
                   help="leave those years blank instead of recomputing the predictors")
    p.add_argument("--dry-run", action="store_true", help="resolve inputs, write nothing")
    p.add_argument("--audit", action="store_true",
                   help="read everything and report year coverage; write only lma_audit.csv. "
                        "Use this to check completeness before committing to a series.")
    args = p.parse_args()

    if args.predictor_root is None:
        args.predictor_root = args.ecoregion_root / PLSR_INPUT_SUBDIR / "LMA"

    wanted = {s.strip() for s in args.stations.split(",") if s.strip()} if args.stations else None
    dropped = read_excluded(args.exclude_file)
    excluded = set(dropped)
    if args.exclude:
        excluded |= {s.strip() for s in args.exclude.split(",") if s.strip()}

    stations = [s for s in read_stations(args.site_list or DEFAULT_SITE_LISTS, wanted)
                if s["station_id"] not in excluded]
    if not stations:
        print("no stations to process", file=sys.stderr)
        return 1
    years_wanted = list(range(args.start_year, args.end_year + 1))

    print(f"stations   : {len(stations)} ({len(excluded)} excluded)")
    print(f"years      : {args.start_year}-{args.end_year}")
    print(f"plsr root  : {args.plsr_root}")
    print(f"eco root   : {args.ecoregion_root}")
    print(f"pred root  : {args.predictor_root}/LMA_ecoregion_no<ii>.csv"
          f"{'' if args.predictor_root.is_dir() else '   <-- NOT FOUND'}")
    print(f"output     : {args.out}")
    print(f"modelled   : {'RECONSTRUCTED from fit coefficients' if args.reconstruct else 'yfit_plot_abs as stored'}")
    print("matching   : nearest pixel outright; its own mean if the fit has one, "
          "else the ecoregion mean\n")

    # One (ecoregion, forest) pair serves many stations -- read each file once.
    groups: dict[tuple[int, str], list[dict]] = {}
    for st in stations:
        groups.setdefault((st["eco_idx"], st["forest_type"]), []).append(st)

    manifest, failures, filled_total = [], 0, 0
    year_hits: dict[int, int] = {}
    store = None
    if args.reconstruct and not args.no_era5_fallback and not (args.dry_run or args.audit):
        try:
            store = era5_predictors.Era5Monthly(args.era5_root)
            print(f"era5 fill  : {args.era5_root}")
        except Exception as exc:                                   # noqa: BLE001
            print(f"era5 fill  : DISABLED -- {type(exc).__name__}: {exc}")
    pix_lookup: dict[int, int] = {}
    for (eco, forest), members in sorted(groups.items()):
        tag = f"eco{eco} {forest} ({len(members)} station{'s' if len(members) > 1 else ''})"
        mat = resolve_prediction_file(args.plsr_root, eco, forest)
        tables = resolve_tables(args.predictor_root, args.ecoregion_root, eco)
        if mat is None or not tables:
            why = ("no prediction .mat" if mat is None else
                   f"no LMA_ecoregion_no{eco}.csv under {args.predictor_root}")
            print(f"  ! {tag}: {why} -- skipped")
            failures += 1
            for st in members:
                manifest.append({**st, "status": "missing_input", "note": why})
            continue
        if args.dry_run:
            cols = table_columns(tables[0])
            print(f"  - {tag}: would read {mat.name} + {tables[0]}")
            print(f"      {len(cols)} columns, coords {'OK' if has_coord_columns(tables[0]) else 'MISSING'}"
                  f", {sum(c.upper().startswith('SSRD') for c in cols)} SSRD column(s)"
                  f"{'   <-- LEGACY TABLE' if is_legacy_table(tables[0], args.predictor_root) else ''}")
            continue

        try:
            pred = load_prediction_mat(mat)
            obs = {"years": pred["time_yrs"].astype(int),
                   "pixels": pred["pixel_id"].astype(int),
                   "values": np.asarray(pred["Y_plot_abs"], dtype=float)}
            class_years = None
            if args.reconstruct:
                fit_path = resolve_fit_file(args.plsr_root, eco, forest)
                if fit_path is None:
                    raise RuntimeError("no PLSR_fitting_coeff_*_TEMPORAL.mat to reconstruct from")
                fit = read_fit_mat(fit_path)
                ptab, _ = pick_predictor_table(tables, fit["predictors"])
                # Coordinates come from the predictor table itself, so pixel IDs,
                # coordinates and predictors are guaranteed to be one vintage.
                eco_csv = ptab
                modelled, diag = reconstruct_modelled(fit, ptab, LU_BASE + FOREST_LU[forest])
                if not len(modelled["values"]):
                    raise RuntimeError(f"reconstruction produced no rows ({diag})")
                class_years = diag.pop("class_years")
                doy_by_pixel, eco_doy = diag.pop("doy_by_pixel"), diag.pop("eco_doy")
                pix_lookup = {int(p): i for i, p in enumerate(fit["uniq_pix_final"])}
                src = (f"reconstructed from {fit_path.name} over {ptab.name}: "
                       f"{diag['n_rows']} pixel-years over {diag['n_pixels_kept']} pixel(s) "
                       f"({diag['n_pixels_in_fit']} with a pixel mean, the rest on the "
                       f"ecoregion mean {fit['mu_Y_eco']:.1f} g/m2), "
                       f"{len(fit['predictors'])} predictor(s), "
                       f"{diag['n_rows_with_imputed_predictor']} row(s) with an imputed "
                       f"predictor, {diag['n_pixels_class_switched']} fit pixel(s) changed class"
                       + (f", dropped {len(diag['dead_pixels'])} pixel(s) with no usable "
                          f"predictor in any year" if diag["dead_pixels"] else ""))
            else:
                modelled = {"years": pred["time_yrs"].astype(int),
                            "pixels": pred["pixel_id"].astype(int),
                            "values": np.asarray(pred["yfit_plot_abs"], dtype=float)}
                src = f"{mat.name}"
                fit, doy_by_pixel, eco_doy = None, {}, float("nan")
                eco_csv = next((t for t in tables if has_coord_columns(t)), None)
                if eco_csv is None:
                    raise RuntimeError(
                        "no candidate table carries Pixel ID/lat/lon: "
                        + ", ".join(str(t) for t in tables))
            coords = read_pixel_coords(eco_csv)
            if is_legacy_table(eco_csv, args.predictor_root):
                print(f"    ! coordinates from the legacy {eco_csv} -- pixel IDs may "
                      f"not match the fit; check {args.predictor_root}")
        except Exception as exc:                                   # noqa: BLE001
            print(f"  ! {tag}: {type(exc).__name__}: {exc}")
            failures += 1
            for st in members:
                manifest.append({**st, "status": "read_error", "note": str(exc)[:200]})
            continue

        print(f"  {tag}: {len(coords)} pixels, {src}")
        series = {"observed": obs, "modelled": modelled}
        for st in members:
            res = build_station(st, series, coords, years_wanted, args.max_distance_km,
                                args.fill_gaps, class_years)
            obs_rows, mod_rows = res.pop("observed"), res.pop("modelled")

            # Years the preprocessed table never emitted -- the pixel was not forest
            # that year, so no row exists even though the ERA5-Land forcing does.
            if store is not None and fit is not None and res.get("pixel_id") != "":
                pid = int(res["pixel_id"])
                pix_index = pix_lookup.get(pid) if res["baseline"] == "pixel_mean" else None
                n_fill, note = era5_fill(
                    store, fit, mod_rows, coords[pid][0], coords[pid][1],
                    doy_by_pixel.get(pid, eco_doy), pix_index, f"{st['station_id']}")
                if n_fill:
                    res["note"] = ((res["note"] + "; ") if res["note"] else "") + note
                    present = [r["LMA_g_m2"] for r in mod_rows if r["LMA_g_m2"] != ""]
                    res.update(modelled_n_years=len(present),
                               modelled_mean=round(float(np.mean(present)), 3),
                               modelled_min=round(float(np.min(present)), 3),
                               modelled_max=round(float(np.max(present)), 3))
                    filled_total += n_fill
            if args.audit:
                # Which years does the modelled series actually carry? This is the
                # question -- the pipeline drops rows with a missing response before
                # predicting, so 'modelled' is not necessarily gap-free.
                res["missing_years"] = " ".join(
                    str(r["year"]) for r in mod_rows if r["LMA_g_m2"] == "")
                for r in mod_rows:
                    year_hits.setdefault(r["year"], 0)
                    year_hits[r["year"]] += r["LMA_g_m2"] != ""
            else:
                sdir = args.out / st["station_id"]
                write_series(sdir / f"{st['station_id']}_LMA_observed.csv", obs_rows)
                write_series(sdir / f"{st['station_id']}_LMA_modelled.csv", mod_rows)
            status = "ok" if res["modelled_n_years"] else "no_data"
            if status != "ok":
                failures += 1
            manifest.append({**st, **res, "status": status})
            flag = "" if status == "ok" else "  <-- " + status
            print(f"      {st['station_id']:9} pixel {str(res['pixel_id']):>6} "
                  f"{str(res['distance_km']):>8} km  obs {res['observed_n_years']:>2}y "
                  f"mod {res['modelled_n_years']:>2}y  [{res['modelled_source']}]{flag}")

    if args.dry_run:
        print(f"\ndry run: {len(groups)} ecoregion x forest group(s), {len(stations)} stations")
        return 0

    args.out.mkdir(parents=True, exist_ok=True)
    cols = ["station_id", "forest_type", "eco_idx", "eco_name", "plsr_q2", "lat", "lon",
            "pixel_id", "pixel_lat", "pixel_lon", "distance_km", "n_pixels",
            "pixel_class_years", "baseline",
            "observed_source", "observed_n_years", "observed_mean", "observed_min", "observed_max",
            "modelled_source", "modelled_n_years", "modelled_mean", "modelled_min", "modelled_max",
            "status", "note"]

    if args.audit:
        n_st = len(manifest)
        want = len(years_wanted)
        complete = [m for m in manifest if m.get("modelled_n_years") == want]
        print(f"\n{'=' * 62}\nAUDIT: modelled series coverage, {args.start_year}-{args.end_year}"
              f" ({want} years)\n{'=' * 62}")
        print(f"stations with a COMPLETE modelled series : {len(complete)}/{n_st}")
        if year_hits:
            worst = sorted(year_hits.items(), key=lambda kv: kv[1])
            print(f"\nyears no station has  : "
                  f"{' '.join(str(y) for y, n in worst if n == 0) or '(none)'}")
            print("thinnest years        : " +
                  ", ".join(f"{y}:{n}/{n_st}" for y, n in worst[:8] if n))
            full = [y for y, n in year_hits.items() if n == n_st]
            print(f"years all stations have: {len(full)}/{want}")
        with open(args.out / "lma_audit.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols + ["missing_years"], extrasaction="ignore")
            w.writeheader()
            w.writerows(manifest)
        print(f"\nper-station detail: {args.out / 'lma_audit.csv'}")
        if len(complete) < n_st:
            print("\nThe modelled series is NOT gap-free. The pipeline drops rows whose\n"
                  "response (iLMA) is missing before predicting, so yfit_plot_abs exists\n"
                  "only where an observation existed. Regenerating a complete series means\n"
                  "re-applying PLSR_fitting_coeff_*_TEMPORAL.mat to the full predictor\n"
                  "table, which is ~98% complete in every year.")
        return 0 if len(complete) == n_st else 1

    with open(args.out / "lma_manifest.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(manifest)

    with open(args.out / "lma_provenance.json", "w", encoding="utf-8") as fh:
        json.dump({
            "plsr_root": str(args.plsr_root),
            "ecoregion_root": str(args.ecoregion_root),
            "predictor_root": str(args.predictor_root),
            "predictor_table": f"{args.predictor_root}/LMA_ecoregion_no<ii>.csv",
            "years": [args.start_year, args.end_year],
            "observed_field": "Y_plot_abs",
            "modelled_field": ("reconstructed: beta/mu_X/sigma_X/mu_Y/sigma_Y and the "
                               "per-pixel means from PLSR_fitting_coeff_*_TEMPORAL.mat, "
                               "re-applied to every row of LMA_ecoregion_no<ii>.csv. "
                               "Predictors are demeaned by the pixel mean, then normalised "
                               "by the ecoregion mean and standard deviation.")
                              if args.reconstruct else "yfit_plot_abs",
            "units": "g/m2 (leaf dry mass per unit leaf area)",
            "baseline_rule": ("pixel mean (mu_X_pix_final/mu_Y_pix_final, computed by the "
                              "fit over the years that pixel carried the station's forest "
                              "type) where the fit covers the nearest pixel; otherwise the "
                              "ecoregion mean, i.e. the mean of those per-pixel means"),
            "max_distance_km": args.max_distance_km,
            "note": "iLMA is the dataset name, not an inverse. SLA conversion is applied "
                    "at .mat build time, not here.",
        }, fh, indent=2)

    ok = sum(1 for m in manifest if m.get("status") == "ok")
    if store is not None:
        store.close()
    fell_back = sum(1 for m in manifest if m.get("baseline") == "ecoregion_mean")
    print(f"\nwrote {ok}/{len(manifest)} stations to {args.out}")
    if filled_total:
        print(f"  {filled_total} station-year(s) recomputed from ERA5-Land because the "
              f"preprocessed table had no row")
    if fell_back:
        print(f"  {fell_back} station(s) sit on a pixel the fit never saw as their forest "
              f"type -- ecoregion mean used as the baseline")
    print(f"  manifest: {args.out / 'lma_manifest.csv'}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
