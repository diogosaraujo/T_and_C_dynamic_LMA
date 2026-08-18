#!/usr/bin/env python3
"""Build the GCM arms of the model_run tree: 5 GCMs x 3 scenarios x fixed/dynamic LMA.

    python build_gcm_model_run.py --dry-run          # what is ready, writes nothing
    python build_gcm_model_run.py                    # build everything available
    python build_gcm_model_run.py --gcm GFDL-ESM4 --scenario ssp585

    model_run/<STATION>/<scenario>/<GCM>/
        fixed_lma/   GO_<ST>.m  MOD_PARAM_<ST>.m  LMA_<ST>.mat  Meteo_<ST>_<GCM>_<scen>_*.mat
        dyn_lma/     same

WHY THIS IS SHORT, AND WHY IT DOES NOT DUPLICATE build_model_run.py

MOD_PARAM is a property of the SITE, not of the forcing: soil layers, canopy
height, rooting depth, zatm and the PFT block do not change because a different
meteorology drives the site. The one field that does change is Sl_H, which is
derived from the LMA series. So this reads the already-built, already-verified
era5_land MOD_PARAM and patches exactly that one line, rather than re-running the
soil/root/canopy substitution machinery. Re-implementing it would mean two copies
of the logic that decides Kbot, caps ZR95 at the column depth and clamps the silt
budget -- and the failure mode of those drifting apart is a plausible, wrong
MOD_PARAM, which is the failure this project has hit most often.

That makes build_model_run.py a PREREQUISITE. Run it first; this refuses any
station whose era5_land arm is missing.

THE LMA SERIES

The PLSR projections arrive per ecoregion x GCM x scenario, already computed:

    <plsr-root>/LMA_ecoregion_no<eco>_<GCM>_<scenario>_projection.csv
        one row per (pixel, year), 1985-2100, with
        LMA_Future = LMA_Baseline + LMA_Anomaly   (exact; verified on read)

Each station takes THE CELL IT SITS IN. Land-cover class (LU 41 deciduous, 42
evergreen) breaks ties within that cell but does not select across cells: a
0.25 degree cell is labelled by its DOMINANT cover, and a tower can sit in a cell
whose majority class differs from the stand it measures. Selecting on class first
sent US-Ha2 -- evergreen hemlock inside a deciduous landscape -- to a pixel
187.9 km away, about seven cells from where it stands. The LMA series represents
the cell's climate, and the station shares the climate of its own cell whatever
the label says.

Note the files span 1985-2100, i.e. they already include the historical period,
and the two SSPs are identical over 1985-2014 because the scenarios only diverge
after 2015. The historical arm therefore reads its series from the ssp126 file
and cross-checks it against ssp585 where both are present; a disagreement means
the projections were not produced from a common historical baseline and is
reported rather than averaged away.

THE FIXED ARM

Sl_H = 1/(mean(LMA_Future over 1985-2014) * f_C), computed PER GCM, and held at
that value for the historical run AND both futures.

Per GCM matters: it keeps the fixed/dynamic pairing inside a single model, so a
GCM's own LMA offset cancels in the difference instead of contaminating it.
Historical-mean matters more: using each scenario's own mean would absorb the
future LMA trend into both arms and null out precisely the signal the experiment
exists to measure.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_model_run import GO_TEMPLATE, write_lma_mat, mat_name, F_C  # noqa: E402
from gcm_variables import GCMS, SCENARIOS                              # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA",
                                 "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
MODEL_RUN = Path(os.environ.get("MODEL_RUN", INPUT_ROOT.parent / "model_run"))
ECOREGION_ROOT = Path(os.environ.get("ECOREGION_ROOT",
                                     "/vol_efthymios/NFS07/dd1136/ecoregions"))
PLSR_ROOT = (ECOREGION_ROOT / "GCMs"
             / "PLSR_future_trait_predictions_GCM_clim_DOY_python" / "LMA")
GCM_METEO = INPUT_ROOT / "gcm_meteo"
SITE_LISTS = [REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
              REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv"]
EXCLUDED = Path(__file__).resolve().parent / "excluded_stations.csv"
ARMS = ["fixed_lma", "dyn_lma"]

HIST_YEARS = (1985, 2014)      # the projections start in 1985, not 1980
FUT_YEARS = (2015, 2100)
YEAR_TAG = {"historical": "1980_2014", "ssp126": "2015_2100", "ssp585": "2015_2100"}
# LU code -> forest type, the 40 + lu_id convention (1 deciduous, 2 evergreen).
LU_OF = {"deciduous": 41, "evergreen": 42}
SL_LINE = re.compile(r"^Sl_H\s*=\s*\[[^\]]*\];.*$", re.M)


def haversine_km(lat0, lon0, lat, lon):
    lon0 = ((lon0 + 180) % 360) - 180
    lon = ((np.asarray(lon) + 180) % 360) - 180
    p0, p = np.radians(lat0), np.radians(np.asarray(lat))
    dp, dl = p - p0, np.radians(lon - lon0)
    a = np.sin(dp / 2) ** 2 + np.cos(p0) * np.cos(p) * np.sin(dl / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def read_stations(wanted=None):
    drop = set()
    if EXCLUDED.is_file():
        for r in csv.DictReader(open(EXCLUDED, newline="", encoding="utf-8-sig")):
            if (r.get("station_id") or "").strip():
                drop.add(r["station_id"].strip())
    out, seen = [], set()
    for p in SITE_LISTS:
        if not p.is_file():
            print(f"  ! site list not found: {p}", file=sys.stderr)
            continue
        for r in csv.DictReader(open(p, newline="", encoding="utf-8-sig")):
            sid = (r.get("StationID") or "").strip()
            if not sid or sid in seen or sid in drop or (wanted and sid not in wanted):
                continue
            try:
                eco = int(r["ECO_IDX"]); lat = float(r["Lat"]); lon = float(r["Lon"])
            except (KeyError, TypeError, ValueError):
                print(f"  ! {sid}: unusable ECO_IDX/Lat/Lon, skipped", file=sys.stderr)
                continue
            seen.add(sid)
            out.append(dict(station_id=sid, eco=eco, lat=lat, lon=lon,
                            forest_type=(r.get("ForestType") or "").strip().lower()))
    return sorted(out, key=lambda s: s["station_id"])


# --------------------------------------------------------------------- LMA series
_CACHE: dict[tuple, dict] = {}


def load_projection(eco, gcm, scenario, plsr_root):
    """{(lat, lon, lu): {year: LMA_Future}} for one ecoregion/GCM/scenario file."""
    key = (eco, gcm, scenario, str(plsr_root))
    if key in _CACHE:
        return _CACHE[key], None
    p = plsr_root / f"LMA_ecoregion_no{eco}_{gcm}_{scenario}_projection.csv"
    if not p.is_file():
        return None, f"no projection file {p.name}"
    px: dict = {}
    bad = 0
    with open(p, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            if (r.get("ValidPrediction") or "").strip().lower() not in ("true", "1", ""):
                continue
            try:
                lat, lon = float(r["lat"]), float(r["lon"])
                yr = int(r["year"])
                base, anom = float(r["LMA_Baseline"]), float(r["LMA_Anomaly"])
                fut = float(r["LMA_Future"])
                lu = int(float(r["LU"]))
            except (KeyError, TypeError, ValueError):
                bad += 1
                continue
            # The stated identity, checked rather than assumed -- if a future file
            # ever stores an anomaly in different units this is where it shows.
            if abs(base + anom - fut) > 1e-6:
                bad += 1
                continue
            px.setdefault((lat, lon, lu), {})[yr] = fut
    _CACHE[key] = px
    return px, (f"{bad} unusable row(s)" if bad else None)


def grid_step(keys):
    """Grid spacing in degrees, inferred from the pixel coordinates themselves."""
    for i in (0, 1):
        v = np.unique(np.round([k[i] for k in keys], 6))
        if len(v) > 1:
            d = np.diff(v)
            d = d[d > 1e-6]
            if len(d):
                return float(np.min(d))
    return 0.25                                   # NEX-GDDP, if it cannot be told


def station_series(st, gcm, scenario, plsr_root):
    """[(year, LMA)] for one station, plus a diagnostic dict.

    THE PIXEL IS THE ONE THE STATION SITS IN, not the nearest one of a matching
    land-cover class.

    The first version took the nearest pixel whose LU equalled the station's
    AmeriFlux forest type. On the ERA5-Land 0.1 degree grid the historical path
    uses, those are almost always the same cell. On the 0.25 degree GCM grid they
    are not: a ~25 km cell is classified by its DOMINANT cover, and a tower can
    easily sit in a cell whose majority class differs from the stand it measures.
    US-Ha2 -- an evergreen hemlock site inside a deciduous-dominated landscape --
    was matched 187.9 km away, roughly seven cells from where it stands.

    So: take the containing cell when one exists, preferring a matching LU only to
    break ties within it, and fall back to nearest-matching-LU only when the
    station falls outside every cell in the table. The land-cover class of a 25 km
    cell is a property of the landscape; the LMA series is a property of the
    climate, and it is the climate the station shares with its own cell.
    """
    px, note = load_projection(st["eco"], gcm, scenario, plsr_root)
    if px is None:
        return None, {"error": note}
    if not px:
        return None, {"error": f"projection file for ecoregion {st['eco']} "
                               f"({gcm} {scenario}) holds no usable pixel"}
    lu = LU_OF.get(st["forest_type"])
    keys = list(px)
    half = grid_step(keys) / 2.0 + 1e-6

    slon360 = st["lon"] % 360.0
    inside = [k for k in keys
              if abs(k[0] - st["lat"]) <= half
              and min(abs(k[1] - slon360), 360 - abs(k[1] - slon360)) <= half]

    if inside:
        match = [k for k in inside if k[2] == lu]
        chosen = (match or inside)[0]
        how = "containing cell" if match else f"containing cell, LU {chosen[2]} not {lu}"
        pool = len(inside)
    else:
        cand = [k for k in keys if lu is None or k[2] == lu] or keys
        lats = np.array([c[0] for c in cand]); lons = np.array([c[1] for c in cand])
        d = haversine_km(st["lat"], st["lon"], lats, lons)
        chosen = cand[int(np.argmin(d))]
        how = "NEAREST -- station outside every cell in the table"
        pool = len(cand)

    d0 = float(haversine_km(st["lat"], st["lon"],
                            np.array([chosen[0]]), np.array([chosen[1]]))[0])
    return (sorted(px[chosen].items()),
            {"pixel_km": d0, "n_pixels": pool, "note": note, "how": how,
             "lu_used": chosen[2], "lu_wanted": lu, "grid_deg": round(half * 2, 4),
             "pixel_lat": chosen[0], "pixel_lon": ((chosen[1] + 180) % 360) - 180})


def clip(series, lo, hi):
    return [(y, v) for y, v in series if lo <= y <= hi]


# ------------------------------------------------------------------------- build
def build_one(st, gcm, scenario, series_fixed_mean, series, meteo_src, out_root,
              era5_mod_param, dry):
    sid, mname = st["station_id"], mat_name(st["station_id"])
    sl_fixed = 1.0 / (series_fixed_mean * F_C)
    tmpl = era5_mod_param.read_text(encoding="utf-8")
    patched, n = SL_LINE.subn(
        f"Sl_H = [{sl_fixed:.6g}]; %% [m^2/gC] 1/(LMA*{F_C}), "
        f"LMA={series_fixed_mean:.1f} g/m2 ({gcm} 1985-2014 mean)", tmpl)
    if n != 1:
        return 0, f"Sl_H matched {n} times in {era5_mod_param.name}, expected once"
    # ms (soil layer count) lives in GO_<ST>.m, not MOD_PARAM. Take it from the
    # era5_land GO rather than recomputing it from Zs: that file is what the
    # completed ERA5 runs actually used, so reading it cannot disagree with them.
    era5_go = era5_mod_param.parent / f"GO_{mname}.m"
    ms = None
    if era5_go.is_file():
        m = re.search(r"^ms\s*=\s*(\d+)", era5_go.read_text(encoding="utf-8"), re.M)
        ms = m.group(1) if m else None
    if ms is None:                       # fall back to the mesh: ms = len(Zs) - 1
        m = re.search(r"^Zs\s*=\s*\[([^\]]*)\]", patched, re.M)
        if m:
            ms = str(len(m.group(1).split()) - 1)
    if ms is None:
        return 0, f"cannot determine ms from {era5_go.name} or Zs in {era5_mod_param.name}"
    if dry:
        return 0, None

    for arm in ARMS:
        d = out_root / sid / scenario / gcm / arm
        d.mkdir(parents=True, exist_ok=True)
        (d / f"MOD_PARAM_{mname}.m").write_text(patched, encoding="utf-8")
        (d / f"GO_{mname}.m").write_text(GO_TEMPLATE.format(
            station=sid, forest_type=st["forest_type"], arm=f"{gcm} {scenario} {arm}",
            code_rel="../../../../Code", root_rel="../../../..",
            ms=ms, mname=mname, meteo_name=meteo_src.name,
            main_frame="MAIN_FRAME_SLA" if arm == "dyn_lma" else "MAIN_FRAME",
        ), encoding="utf-8")
        dst = d / meteo_src.name
        if not dst.exists():
            try:
                os.symlink(meteo_src, dst)
            except OSError:
                shutil.copy2(meteo_src, dst)
        write_lma_mat(d / f"LMA_{mname}.mat", series, sl_fixed, arm)
    return len(ARMS), None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gcm", action="append")
    ap.add_argument("--scenario", action="append")
    ap.add_argument("--station", action="append")
    ap.add_argument("--root", type=Path, default=MODEL_RUN)
    ap.add_argument("--plsr-root", type=Path, default=PLSR_ROOT)
    ap.add_argument("--meteo", type=Path, default=GCM_METEO)
    ap.add_argument("--max-pixel-km", type=float, default=50.0,
                    help="report stations whose nearest same-forest-type PLSR pixel "
                         "is further than this (default 50 km). Reporting only; "
                         "use --exclude-far to actually block them.")
    ap.add_argument("--exclude-far", action="store_true",
                    help="block combinations beyond --max-pixel-km instead of "
                         "reporting them")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    stations = read_stations(set(a.station) if a.station else None)
    gcms = a.gcm or GCMS
    scens = a.scenario or list(SCENARIOS)

    print(f"model_run : {a.root}")
    print(f"plsr      : {a.plsr_root}"
          f"{'' if a.plsr_root.is_dir() else '   <-- NOT FOUND'}")
    print(f"gcm meteo : {a.meteo}{'' if a.meteo.is_dir() else '   <-- NOT FOUND'}")
    print(f"stations  : {len(stations)}   gcms {len(gcms)}   scenarios {len(scens)}")
    print(f"target    : {len(stations)*len(gcms)*len(scens)*len(ARMS)} run directories\n")

    import collections
    how = collections.Counter()
    lu_mismatch = []
    written, ok_arms, blocked, runs = 0, 0, [], []
    offsets = []          # (station, km) so far matches can be named, not just counted
    for st in stations:
        sid, mname = st["station_id"], mat_name(st["station_id"])
        era5_mp = a.root / sid / "era5_land" / "fixed_lma" / f"MOD_PARAM_{mname}.m"
        if not era5_mp.is_file():
            blocked.append((sid, "-", "-", "no era5_land MOD_PARAM; run "
                                           "build_model_run.py first"))
            continue
        for gcm in gcms:
            # The fixed value is the 1985-2014 mean of THIS GCM, taken once and
            # reused for the historical run and both futures.
            hist_src = "ssp126" if "ssp126" in SCENARIOS else scens[0]
            hs, hd = station_series(st, gcm, hist_src, a.plsr_root)
            if hs is None:
                blocked.append((sid, gcm, "-", hd["error"]))
                continue
            hist = clip(hs, *HIST_YEARS)
            if len(hist) < 20:
                blocked.append((sid, gcm, "-",
                                f"only {len(hist)} historical years in the projection"))
                continue
            if a.exclude_far and hd["pixel_km"] > a.max_pixel_km:
                blocked.append((sid, gcm, "-", f"nearest {st['forest_type']} pixel is "
                                               f"{hd['pixel_km']:.0f} km away"))
                continue
            fixed_mean = float(np.mean([v for _, v in hist]))
            offsets.append((sid, hd["pixel_km"]))
            how[hd.get("how", "?")] += 1
            if hd.get("lu_used") != hd.get("lu_wanted"):
                lu_mismatch.append((sid, st["forest_type"], hd.get("lu_used"),
                                    hd.get("pixel_km")))

            for scen in scens:
                if scen == "historical":
                    series, diag = hist, hd
                else:
                    s, diag = station_series(st, gcm, scen, a.plsr_root)
                    if s is None:
                        blocked.append((sid, gcm, scen, diag["error"]))
                        continue
                    series = clip(s, *FUT_YEARS)
                if not series:
                    blocked.append((sid, gcm, scen, "empty series after clipping"))
                    continue
                mfile = (a.meteo / scen / mat_name(gcm) /
                         f"Meteo_{mname}_{mat_name(gcm)}_{scen}_{YEAR_TAG[scen]}.mat")
                if not mfile.is_file():
                    blocked.append((sid, gcm, scen, f"no forcing {mfile.name}"))
                    continue
                n, err = build_one(st, gcm, scen, fixed_mean, series, mfile,
                                   a.root, era5_mp, a.dry_run)
                if err:
                    blocked.append((sid, gcm, scen, err))
                    continue
                written += n
                ok_arms += len(ARMS)
                for arm in ARMS:
                    runs.append(f"{sid} {scen} {gcm} {arm}")

    print(f"{'=' * 72}")
    verb = "would be built" if a.dry_run else "written"
    print(f"{'DRY RUN -- ' if a.dry_run else ''}{ok_arms} run directories {verb}, "
          f"{len(blocked)} combination(s) blocked")
    print(f"{'=' * 72}")

    if offsets:
        o = np.array([d for _, d in offsets])
        print("\nSTATION-TO-PLSR-PIXEL DISTANCE")
        print("  " + "  ".join(f"p{p}={np.percentile(o, p):.1f}" for p in (50, 75, 90, 95, 99))
              + f"  max={o.max():.1f} km  (n={len(o)})")
        far = sorted({s for s, d in offsets if d > a.max_pixel_km},
                     key=lambda s: -max(d for ss, d in offsets if ss == s))
        if far:
            print(f"  beyond {a.max_pixel_km:.0f} km -- the LMA series comes from a pixel that far "
                  f"from the tower:")
            for s in far[:15]:
                d = max(d for ss, d in offsets if ss == s)
                print(f"    {s:<9} {d:6.1f} km")
            print(f"  These are NOT blocked. A distant pixel still has the right forest type and\n"
                  f"  ecoregion, but it is not the tower's own climate -- decide per station\n"
                  f"  whether that is acceptable, or re-run with --max-pixel-km to exclude them.")

    if blocked:
        # Group by reason, because one cause usually explains many rows and the
        # per-row list buries that. Every row is still printed underneath.
        import collections
        cat = collections.Counter()
        for sid, gcm, scen, why in blocked:
            key = re.sub(r"\d+", "N", why)
            key = re.sub(r"(for |in ecoregion )\S+", r"\1...", key)
            cat[key] += 1
        print(f"\nBLOCKED, by cause ({len(blocked)} total)")
        for why, n in cat.most_common():
            print(f"  {n:>4d} x  {why}")
        print(f"\nBLOCKED, in full")
        for sid, gcm, scen, why in sorted(blocked):
            print(f"  ! {sid:<9} {gcm:<14} {scen:<11} {why}")
    if runs and not a.dry_run:
        lst = a.root / "run_list_gcm.txt"
        lst.write_text("".join(r + "\n" for r in runs), encoding="utf-8")
        print(f"\nrun list : {lst}  ({len(runs)} arms)")
        print(f"           sbatch --array=1-{len(runs)}%NN slurm/submit_tc_run.sh "
              f"(with RUN_LIST=run_list_gcm.txt)")
    return 0 if ok_arms else 1


if __name__ == "__main__":
    sys.exit(main())
