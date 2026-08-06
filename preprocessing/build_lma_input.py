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
    <ecoregion-root>/ecoregion_no<ii>.csv        -- Pixel ID, lat, lon

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

# PLSR_TEMPORAL_CV_OOF_METRICS_PIXEL_DOY_CORE.m spells the first class 'deciduos'
# (missing a 'u'). Try the typo first, then the correct spelling, so the script
# keeps working if the upstream name is ever fixed.
FOREST_LABELS = {
    "deciduous": ["deciduos", "deciduous"],
    "evergreen": ["evergreen"],
    "mixed": ["mixed"],
}

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


def _is_hdf5(path: Path) -> bool:
    with open(path, "rb") as fh:
        return fh.read(8) == b"\x89HDF\r\n\x1a\n"


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


def write_series(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["year", "LMA_g_m2", "source", "pixel_id", "n_values"])
        w.writeheader()
        w.writerows(rows)


# ------------------------------------------------------------------------- core

def build_station(st: dict, pred: dict[str, np.ndarray], coords: dict[int, tuple[float, float]],
                  years_wanted: list[int], max_dist_km: float, fill_gaps: bool) -> dict:
    """Return {'observed': rows, 'modelled': rows, **manifest fields} for one station."""
    yrs = pred["time_yrs"].astype(int)
    pix = pred["pixel_id"].astype(int)

    # Only pixels that actually appear in the predictions can be matched -- picking
    # the nearest pixel of the ecoregion at large could land on one the model never
    # saw (predict_with_fit drops pixels absent from the training set).
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
    for kind, field in (("observed", "Y_plot_abs"), ("modelled", "yfit_plot_abs")):
        vals = np.asarray(pred[field], dtype=float)
        eco_series = annual_median(yrs, vals)                      # fallback, all pixels
        pix_series = annual_median(yrs[pix == best], vals[pix == best]) if best is not None else {}

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
                   help="directory holding ecoregion_no<ii>.csv")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--start-year", type=int, default=1985)
    p.add_argument("--end-year", type=int, default=2021)
    p.add_argument("--max-distance-km", type=float, default=50.0,
                   help="beyond this, fall back to the ecoregion median (default: 50)")
    p.add_argument("--fill-gaps", action="store_true",
                   help="fill individual missing years from the ecoregion median "
                        "(default: fall back only when the pixel has no data at all)")
    p.add_argument("--dry-run", action="store_true", help="resolve inputs, write nothing")
    args = p.parse_args()

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
    print(f"output     : {args.out}")
    print(f"fallback   : {'per-year' if args.fill_gaps else 'whole-series'} ecoregion median\n")

    # One (ecoregion, forest) pair serves many stations -- read each file once.
    groups: dict[tuple[int, str], list[dict]] = {}
    for st in stations:
        groups.setdefault((st["eco_idx"], st["forest_type"]), []).append(st)

    manifest, failures = [], 0
    for (eco, forest), members in sorted(groups.items()):
        tag = f"eco{eco} {forest} ({len(members)} station{'s' if len(members) > 1 else ''})"
        mat = resolve_prediction_file(args.plsr_root, eco, forest)
        eco_csv = args.ecoregion_root / f"ecoregion_no{eco}.csv"
        if mat is None or not eco_csv.is_file():
            why = "no prediction .mat" if mat is None else f"no {eco_csv.name}"
            print(f"  ! {tag}: {why} -- skipped")
            failures += 1
            for st in members:
                manifest.append({**st, "status": "missing_input", "note": why})
            continue
        if args.dry_run:
            print(f"  - {tag}: would read {mat.name} + {eco_csv.name}")
            continue

        try:
            pred = load_prediction_mat(mat)
            coords = read_pixel_coords(eco_csv)
        except Exception as exc:                                   # noqa: BLE001
            print(f"  ! {tag}: {type(exc).__name__}: {exc}")
            failures += 1
            for st in members:
                manifest.append({**st, "status": "read_error", "note": str(exc)[:200]})
            continue

        print(f"  {tag}: {mat.name}, {len(coords)} pixels")
        for st in members:
            res = build_station(st, pred, coords, years_wanted, args.max_distance_km,
                                args.fill_gaps)
            sdir = args.out / st["station_id"]
            write_series(sdir / f"{st['station_id']}_LMA_observed.csv", res.pop("observed"))
            write_series(sdir / f"{st['station_id']}_LMA_modelled.csv", res.pop("modelled"))
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
    with open(args.out / "lma_manifest.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(manifest)

    with open(args.out / "lma_provenance.json", "w", encoding="utf-8") as fh:
        json.dump({
            "plsr_root": str(args.plsr_root),
            "ecoregion_root": str(args.ecoregion_root),
            "years": [args.start_year, args.end_year],
            "observed_field": "Y_plot_abs", "modelled_field": "yfit_plot_abs",
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
