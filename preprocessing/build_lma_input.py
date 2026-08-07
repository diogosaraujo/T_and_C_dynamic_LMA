#!/usr/bin/env python3
"""Build per-station time-varying LMA series from the PLSR temporal-CV output.

Two series are produced for every AmeriFlux station, both in g/m2:

    observed  <- Y_plot_abs      the iLMA observations that entered the fit
    modelled  <- yfit_plot_abs   the PLSR reconstruction (pixel climatology + anomaly)

Each station is matched to the nearest pixel of its ecoregion, and the pixel's
annual series is written out. Where that pixel carries no usable data the
ecoregion median for the same year/forest type is used instead; every row records
which of the two it came from, so the fallback is never silent.

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
    return out


def reconstruct_modelled(fit: dict, eco_csv: Path, lu_value: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """Re-apply the fitted model to the FULL predictor table.

    The pipeline only predicted rows that survived rmmissing on the response, so
    yfit_plot_abs inherits the observations' gaps. The predictors are ~98% complete
    in every year, so evaluating the same coefficients over every row of
    LMA_ecoregion_no<ii>.csv fills those years. This reproduces predict_with_fit exactly:

        X_anom   = X_raw - mu_X_pix_final(pix,:)
        X        = (X_anom - mu_X) ./ sigma_X ;  non-finite -> 0
        yfit_abs = ([1 X]*beta) .* sigma_Y + mu_Y + mu_Y_pix_final(pix)

    Only pixels in uniq_pix_final are used -- those are the pixels that contributed
    to the fit, i.e. the ones that had at least one LMA observation.

    Returns (years, pixels, values, diagnostics).
    """
    preds = fit["predictors"]
    pix_row = {int(p): i for i, p in enumerate(fit["uniq_pix_final"])}

    years, pixels, rows = [], [], []
    n_lu_skip = n_pix_skip = 0
    with open(eco_csv, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        cols = {c.strip(): c for c in (reader.fieldnames or [])}
        missing = [p for p in preds if p not in cols]
        if missing:
            raise RuntimeError(f"{eco_csv.name}: predictor column(s) absent: {', '.join(missing)}")
        pid_col = cols.get("Pixel ID") or cols.get("Pixel_ID")
        lu_col, t_col = cols.get("LU"), cols.get("time")
        if not (pid_col and lu_col and t_col):
            raise RuntimeError(f"{eco_csv.name}: need Pixel ID, LU and time columns")
        for r in reader:
            try:
                if int(float(r[lu_col])) != lu_value:
                    n_lu_skip += 1
                    continue
                pid = int(float(r[pid_col]))
                yr = int(str(r[t_col])[-4:])
            except (TypeError, ValueError):
                continue
            if pid not in pix_row:          # pixel contributed nothing to the fit
                n_pix_skip += 1
                continue
            vals = []
            for p in preds:
                try:
                    vals.append(float(r[cols[p]]))
                except (TypeError, ValueError):
                    vals.append(np.nan)
            years.append(yr)
            pixels.append(pid)
            rows.append(vals)

    diag = {"n_rows": len(rows), "n_rows_wrong_lu": n_lu_skip, "n_rows_pixel_not_in_fit": n_pix_skip}
    if not rows:
        return np.array([]), np.array([]), np.array([]), diag

    X_raw = np.asarray(rows, dtype=float)
    idx = np.array([pix_row[p] for p in pixels])
    X = (X_raw - fit["mu_X_pix_final"][idx, :] - fit["mu_X"]) / fit["sigma_X"]
    # predict_with_fit zeroes non-finite predictors, i.e. substitutes the mean. That
    # silently degrades a prediction toward the pixel climatology, so it is counted.
    bad = ~np.isfinite(X)
    X[bad] = 0.0
    diag["n_rows_with_imputed_predictor"] = int(np.any(bad, axis=1).sum())
    diag["n_rows_all_predictors_imputed"] = int(np.all(bad, axis=1).sum())

    yfit = fit["beta"][0] + X @ fit["beta"][1:]
    vals = yfit * fit["sigma_Y"] + fit["mu_Y"] + fit["mu_Y_pix_final"][idx]
    return np.asarray(years), np.asarray(pixels), vals, diag


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


def write_series(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["year", "LMA_g_m2", "source", "pixel_id", "n_values"])
        w.writeheader()
        w.writerows(rows)


# ------------------------------------------------------------------------- core

def build_station(st: dict, series: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
                  coords: dict[int, tuple[float, float]],
                  years_wanted: list[int], max_dist_km: float, fill_gaps: bool) -> dict:
    """Return {'observed': rows, 'modelled': rows, **manifest fields} for one station.

    `series` maps each kind to (years, pixels, values). The two need not share a row
    index: the observed series is sparse, while a reconstructed modelled series covers
    every pixel-year.
    """
    # Match on the modelled series -- that is what gets forced, and its pixel set is
    # the one that contributed to the fit (pixels with no LMA at all never appear in
    # uniq_pix_final). Matching against the ecoregion at large could otherwise land
    # on a pixel the model never saw.
    pix = series["modelled"][1].astype(int) if len(series["modelled"][1]) else np.array([], int)
    have = sorted({int(p) for p in np.unique(pix)} & set(coords))
    info = {"n_pixels": len(have), "pixel_id": "", "pixel_lat": "", "pixel_lon": "",
            "distance_km": "", "note": ""}
    if not have:
        info["note"] = "no prediction pixel has coordinates in the ecoregion table"
        best = None
    else:
        plat = np.array([coords[p][0] for p in have])
        plon = np.array([coords[p][1] for p in have])
        d = haversine_km(st["lat"], st["lon"], plat, plon)
        k = int(np.argmin(d))
        best, dist = have[k], float(d[k])
        info.update(pixel_id=best, pixel_lat=round(plat[k], 5), pixel_lon=round(plon[k], 5),
                    distance_km=round(dist, 3))
        if dist > max_dist_km:
            info["note"] = f"nearest pixel {dist:.1f} km away (> {max_dist_km} km); using ecoregion median"
            best = None

    out = {**info}
    for kind in ("observed", "modelled"):
        k_yrs, k_pix, k_vals = series[kind]
        k_yrs = np.asarray(k_yrs, dtype=int)
        k_pix = np.asarray(k_pix, dtype=int)
        vals = np.asarray(k_vals, dtype=float)
        # Fallback median is over the valid pixels only, since those are the only ones
        # present in either series.
        eco_series = annual_median(k_yrs, vals)
        sel = k_pix == best if best is not None else np.zeros(len(k_yrs), bool)
        pix_series = annual_median(k_yrs[sel], vals[sel]) if best is not None else {}

        source = "pixel" if pix_series else "ecoregion_median"
        rows = []
        for yr in years_wanted:
            if yr in pix_series:
                v, n = pix_series[yr]
                rows.append({"year": yr, "LMA_g_m2": round(v, 4), "source": "pixel",
                             "pixel_id": best, "n_values": n})
            elif (fill_gaps or not pix_series) and yr in eco_series:
                v, n = eco_series[yr]
                rows.append({"year": yr, "LMA_g_m2": round(v, 4), "source": "ecoregion_median",
                             "pixel_id": "", "n_values": n})
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
                   help="fill individual missing years from the ecoregion median "
                        "(default: fall back only when the pixel has no data at all)")
    p.add_argument("--reconstruct", action="store_true",
                   help="rebuild the modelled series by re-applying "
                        "PLSR_fitting_coeff_*_TEMPORAL.mat to the full predictor table, "
                        "instead of reading the gapped yfit_plot_abs. Uses only pixels in "
                        "uniq_pix_final, i.e. those that contributed to the fit.")
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
    print(f"fallback   : {'per-year' if args.fill_gaps else 'whole-series'} ecoregion median\n")

    # One (ecoregion, forest) pair serves many stations -- read each file once.
    groups: dict[tuple[int, str], list[dict]] = {}
    for st in stations:
        groups.setdefault((st["eco_idx"], st["forest_type"]), []).append(st)

    manifest, failures = [], 0
    year_hits: dict[int, int] = {}
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
            obs = (pred["time_yrs"].astype(int), pred["pixel_id"].astype(int),
                   np.asarray(pred["Y_plot_abs"], dtype=float))
            if args.reconstruct:
                fit_path = resolve_fit_file(args.plsr_root, eco, forest)
                if fit_path is None:
                    raise RuntimeError("no PLSR_fitting_coeff_*_TEMPORAL.mat to reconstruct from")
                fit = read_fit_mat(fit_path)
                ptab, _ = pick_predictor_table(tables, fit["predictors"])
                # Coordinates come from the predictor table itself, so pixel IDs,
                # coordinates and predictors are guaranteed to be one vintage.
                eco_csv = ptab
                ry, rp, rv, diag = reconstruct_modelled(fit, ptab, LU_BASE + FOREST_LU[forest])
                if not len(rv):
                    raise RuntimeError(f"reconstruction produced no rows (LU filter: {diag})")
                modelled = (ry, rp, rv)
                src = (f"reconstructed from {fit_path.name} over {ptab.name}: "
                       f"{diag['n_rows']} pixel-years, {len(fit['predictors'])} predictor(s), "
                       f"{diag['n_rows_with_imputed_predictor']} row(s) with an imputed predictor")
            else:
                modelled = (pred["time_yrs"].astype(int), pred["pixel_id"].astype(int),
                            np.asarray(pred["yfit_plot_abs"], dtype=float))
                src = f"{mat.name}"
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
                                args.fill_gaps)
            obs_rows, mod_rows = res.pop("observed"), res.pop("modelled")
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
            "modelled_field": ("reconstructed: beta/mu/sigma from "
                               "PLSR_fitting_coeff_*_TEMPORAL.mat re-applied to the full "
                               "LMA_ecoregion_no<ii>.csv predictor table, restricted to "
                               "uniq_pix_final and LU = 40 + forest id")
                              if args.reconstruct else "yfit_plot_abs",
            "units": "g/m2 (leaf dry mass per unit leaf area)",
            "fallback": "per-year ecoregion median" if args.fill_gaps
                        else "whole-series ecoregion median",
            "max_distance_km": args.max_distance_km,
            "note": "iLMA is the dataset name, not an inverse. SLA conversion is applied "
                    "at .mat build time, not here.",
        }, fh, indent=2)

    ok = sum(1 for m in manifest if m.get("status") == "ok")
    fell_back = sum(1 for m in manifest if m.get("modelled_source") == "ecoregion_median")
    print(f"\nwrote {ok}/{len(manifest)} stations to {args.out}")
    if fell_back:
        print(f"  {fell_back} used the ecoregion median instead of a pixel series")
    print(f"  manifest: {args.out / 'lma_manifest.csv'}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
