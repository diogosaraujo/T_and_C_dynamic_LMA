#!/usr/bin/env python3
"""Build the GCM arms of the model_run tree: 5 GCMs x 3 scenarios x fixed/dynamic LMA.

    python build_gcm_model_run.py --dry-run          # what is ready, writes nothing
    python build_gcm_model_run.py                    # build everything available
    python build_gcm_model_run.py --gcm GFDL-ESM4 --scenario ssp585

    model_run/<STATION>/<scenario>/<GCM>/
        Meteo_<ST>_<GCM>_<scen>_<years>.mat     <- the forcing, ONE copy
        fixed_lma/   GO_<ST>.m  MOD_PARAM_<ST>.m  LMA_<ST>.mat
        dyn_lma/     same

The forcing sits above the arms, not inside them, and GO loads '../<name>'. Both
arms read the same file, model_run carries everything a run needs, and the tree
moves between clusters with a plain rsync -- none of which was true while each
arm held its own absolute symlink into input_data.

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

Station -> series, in two steps, the same rule build_lma_input.py applies to the
ERA5-Land runs:

  1. the cell the station SITS IN, whole, if that cell is mapped as the station's
     own forest type (LU 41 deciduous, 42 evergreen);
  2. if the cell is mapped as another type: its DYNAMICS, with the LEVEL replaced
     by the ecoregion mean baseline of the station's own type;
  3. if the station falls outside every cell in the table: the ecoregion mean of
     its own type supplies both.

Step 2 is build_lma_input.py's substitution transposed. There:

    mu_Y_row = np.where(in_fit, mu_Y_pix_final[idx], mu_Y_eco)

only the pixel MEAN is replaced, while the predictors still come from the
station's own pixel, so the level moves and the variability does not. Here
LMA_Baseline IS that per-cell mean and LMA_Anomaly is what the predictors
produced, so replacing the baseline alone reproduces the split exactly.

Selecting instead on nearest-same-type -- the first version of this -- sent
US-Ha2, evergreen hemlock inside a deciduous-dominated landscape, to a cell
187.9 km away, about seven cells from where it stands.

ONE DEPARTURE, stated rather than hidden: in the ERA5 path the station's own
forest-type beta is applied to the pixel's climate. These projections ship as
finished values, so under step 2 the anomaly carries the model of the cell's own
type. The level is corrected; the sensitivity is not. Re-fitting would be the
only way to close that, and it is not something this script can do.

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
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_model_run import (GO_TEMPLATE, write_lma_mat, mat_name, F_C,  # noqa: E402
                             is_dynamic, read_ic_table, apply_ic)
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

# DERIVED from gcm_variables.SCENARIOS, never restated. The LMA series written
# here and the forcing built by build_gcm_meteo.py must agree year for year,
# because MAIN_FRAME_SLA looks the simulated year up in the LMA file and stops if
# it is not there. Job 37524 is what a restated constant costs: the forcing said
# 1980-2014, this file said 1985-2014, and the historical dyn_lma arm died 27 s in
# on 1 January 1980. The window now has exactly one home.
HIST_YEARS = SCENARIOS["historical"]
FUT_YEARS = SCENARIOS["ssp585"]
YEAR_TAG = {s: f"{y0}_{y1}" for s, (y0, y1) in SCENARIOS.items()}
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
            # Baseline and anomaly are kept SEPARATE, because the fallback
            # recombines them: LMA_Future = LMA_Baseline + LMA_Anomaly, and the
            # ERA5-Land rule substitutes only the baseline. LMA_Baseline is
            # constant per cell, so the dict holds it once.
            e = px.setdefault((lat, lon, lu), {"base": base, "anom": {}})
            e["anom"][yr] = anom
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

    match = [k for k in inside if k[2] == lu]
    same = [k for k in keys if k[2] == lu]
    if not same:
        return None, {"error": f"no LU {lu} ({st['forest_type']}) cell anywhere in "
                               f"ecoregion {st['eco']} for {gcm} {scenario} -- the "
                               f"{st['forest_type']} PLSR was not fitted here"}

    # ---- 1. own cell, mapped as the station's own forest type: use it whole.
    if match:
        c = match[0]
        d0 = float(haversine_km(st["lat"], st["lon"],
                                np.array([c[0]]), np.array([c[1]]))[0])
        series = sorted((y, px[c]["base"] + a) for y, a in px[c]["anom"].items())
        return series, {"pixel_km": d0, "n_pixels": len(inside), "note": note,
                        "how": "own cell", "lu_used": c[2], "lu_wanted": lu,
                        "grid_deg": round(half * 2, 4), "pixel_lat": c[0],
                        "pixel_lon": ((c[1] + 180) % 360) - 180}

    # ---- 2. own cell exists but is mapped as another type: keep its DYNAMICS,
    #         take the LEVEL from the ecoregion mean of the station's own type.
    #         This is build_lma_input.py's substitution transposed:
    #             mu_Y_row = np.where(in_fit, mu_Y_pix_final[idx], mu_Y_eco)
    #         There, only the pixel MEAN is replaced and the predictors still come
    #         from the station's own pixel. Here, LMA_Baseline is that mean and
    #         LMA_Anomaly is what the predictors produced, so replacing the
    #         baseline alone reproduces the same split exactly.
    eco_base = float(np.mean([px[k]["base"] for k in same]))
    if inside:
        c = inside[0]
        d0 = float(haversine_km(st["lat"], st["lon"],
                                np.array([c[0]]), np.array([c[1]]))[0])
        series = sorted((y, eco_base + a) for y, a in px[c]["anom"].items())
        return series, {
            "pixel_km": d0, "n_pixels": len(same), "note": note,
            "how": "own cell dynamics + ecoregion baseline",
            "lu_used": c[2], "lu_wanted": lu, "eco_base": eco_base,
            "cell_base": px[c]["base"], "grid_deg": round(half * 2, 4),
            "pixel_lat": c[0], "pixel_lon": ((c[1] + 180) % 360) - 180}

    # ---- 3. the station falls outside every cell in the table: nothing local to
    #         borrow dynamics from, so the ecoregion mean supplies both.
    #
    #         How far the nearest cell is decides whether this is real. A station
    #         genuinely between forest cells sits tens of km from the nearest one;
    #         a station a few km away means the containing test is wrong -- a
    #         coordinate convention (cell corner vs centre) or a grid step
    #         mis-inferred from sparse coordinates.
    all_lat = np.array([k[0] for k in keys]); all_lon = np.array([k[1] for k in keys])
    d_any = float(haversine_km(st["lat"], st["lon"], all_lat, all_lon).min())
    s_lat = np.array([k[0] for k in same]); s_lon = np.array([k[1] for k in same])
    d_same = float(haversine_km(st["lat"], st["lon"], s_lat, s_lon).min())
    years = sorted({y for k in same for y in px[k]["anom"]})
    series = []
    for y in years:
        a = [px[k]["anom"][y] for k in same if y in px[k]["anom"]]
        if a:
            series.append((y, eco_base + float(np.mean(a))))
    return series, {
        "pixel_km": float("nan"), "n_pixels": len(same), "note": note,
        "how": "ecoregion mean (station outside every cell)",
        "d_nearest_any": d_any, "d_nearest_same": d_same,
        "lu_used": lu, "lu_wanted": lu, "eco_base": eco_base,
        "grid_deg": round(half * 2, 4),
        "pixel_lat": float(np.mean([k[0] for k in same])),
        "pixel_lon": float(np.mean([((k[1] + 180) % 360) - 180 for k in same]))}


def clip(series, lo, hi):
    return [(y, v) for y, v in series if lo <= y <= hi]


# ------------------------------------------------------------------------- build
def build_one(st, gcm, scenario, series_fixed_mean, series, meteo_src, out_root,
              era5_mod_param, dry, arms=None, ic_table=None, ic_key=""):
    sid, mname = st["station_id"], mat_name(st["station_id"])
    arms = arms or ARMS
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

    for arm in arms:
        d = out_root / sid / scenario / gcm / arm
        d.mkdir(parents=True, exist_ok=True)
        txt = patched
        if ic_table is not None:
            # No fallback to the template pools: a combination with no harvested
            # state is refused in main(), where it can be named.
            key = ic_key.format(station=sid, scenario=scenario, gcm=gcm)
            txt = apply_ic(txt, ic_table[(sid, key)], f"{sid}/{scenario}/{gcm}/{arm}")
        (d / f"MOD_PARAM_{mname}.m").write_text(txt, encoding="utf-8")
        (d / f"GO_{mname}.m").write_text(GO_TEMPLATE.format(
            station=sid, forest_type=st["forest_type"], arm=f"{gcm} {scenario} {arm}",
            code_rel="../../../../Code", root_rel="../../../..",
            ms=ms, mname=mname, meteo_path=f"../{meteo_src.name}",
            main_frame="MAIN_FRAME_SLA" if is_dynamic(arm) else "MAIN_FRAME",
        ), encoding="utf-8")
        # No symlink and no copy. meteo_src already IS d.parent/<name> -- the
        # forcing sits at <ST>/<scenario>/<GCM>/, one level above both arms, put
        # there by finish_meteo.m. Symlinking it into each arm is what made the
        # tree depend on absolute paths into input_data and stopped it moving
        # between clusters; it also duplicated a 40-115 MB file per arm.
        write_lma_mat(d / f"LMA_{mname}.mat", series, sl_fixed, arm)
    return len(arms), None


def inspect(plsr_root, stations, gcm, scenario):
    """What each ecoregion's projection file actually contains.

    The selection rule depends on whether a mixed ecoregion's file carries BOTH
    forest types or only the dominant one, and that cannot be inferred from the
    filename -- there is one file per (ecoregion, GCM, scenario) and no forest
    type in it. This answers the question directly, per ecoregion, naming the
    stations that would be affected.
    """
    import collections
    want = collections.defaultdict(set)
    for st in stations:
        want[st["eco"]].add(st["forest_type"])
    print(f"inspecting {gcm} {scenario}\n")
    print(f"  {'eco':>4s}  {'stations':28s}{'ModelForest in file':26s}"
          f"{'LU':10s}{'pixels':>7s}  verdict")
    trouble = []
    for eco in sorted(want):
        px, err = load_projection(eco, gcm, scenario, plsr_root)
        sts = [s for s in stations if s["eco"] == eco]
        lbl = ",".join(sorted({s["forest_type"][:4] for s in sts})) + f" ({len(sts)} st)"
        if px is None:
            print(f"  {eco:>4d}  {lbl:28s}{'-- ' + str(err):26s}")
            trouble.append((eco, str(err)))
            continue
        lus = sorted({k[2] for k in px})
        have = {41: "deciduous", 42: "evergreen"}
        types_present = {have.get(l, str(l)) for l in lus}
        missing = want[eco] - types_present
        v = "ok" if not missing else f"MISSING {','.join(sorted(missing))}"
        if missing:
            trouble.append((eco, v))
        print(f"  {eco:>4d}  {lbl:28s}{','.join(sorted(types_present)):26s}"
              f"{','.join(str(l) for l in lus):10s}{len(px):>7d}  {v}")
    print()
    if trouble:
        print(f"{len(trouble)} ecoregion(s) cannot serve one of their stations' forest")
        print("types from this product. Those stations need the PLSR to project their")
        print("own type over the containing cell, which the fit can do but this file")
        print("does not carry.")
    else:
        print("Every ecoregion carries every forest type its stations need, so the")
        print("station's own type can be selected on ModelForest directly.")
    return 0


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
    ap.add_argument("--inspect", action="store_true",
                    help="report what each ecoregion's projection file contains "
                         "(forest types, LU classes, pixel counts) and stop")
    ap.add_argument("--arms", type=lambda x: [y for y in x.split(",") if y],
                    default=None,
                    help="arm directory names (default fixed_lma,dyn_lma). Use "
                         "--arms spinup for the disposable spin-up pass. An arm "
                         "named dyn_lma* gets the yearly SLA series.")
    ap.add_argument("--ic", type=Path, default=None,
                    help="initial_state.csv from harvest_state.py. With this, "
                         "MOD_PARAM's LAI_H/B_H/PHE_S_H/AgeL_H come from the "
                         "harvested state and combinations without one are REFUSED.")
    ap.add_argument("--ic-key", default="era5_land/fixed_lma",
                    help="which harvested state to use, as it appears in the 'key' "
                         "column. {station}/{scenario}/{gcm} are substituted, so "
                         "'historical/{gcm}/spinup' picks each GCM's own spin-up.")
    ap.add_argument("--require-era5-state", action="store_true",
                    help="drop any station with no 'era5_land/fixed_lma' row in "
                         "--ic, so the GCM fleet covers exactly the stations the "
                         "ERA5 stage completed")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    stations = read_stations(set(a.station) if a.station else None)
    gcms = a.gcm or GCMS
    if a.inspect:
        return inspect(a.plsr_root, stations, gcms[0],
                       (a.scenario or ["ssp126"])[0])
    scens = a.scenario or list(SCENARIOS)

    print(f"model_run : {a.root}")
    print(f"plsr      : {a.plsr_root}"
          f"{'' if a.plsr_root.is_dir() else '   <-- NOT FOUND'}")
    print(f"gcm meteo : {a.meteo}{'' if a.meteo.is_dir() else '   <-- NOT FOUND'}")
    print(f"stations  : {len(stations)}   gcms {len(gcms)}   scenarios {len(scens)}")
    arms = a.arms or ARMS
    print(f"arms      : {', '.join(arms)}")
    print(f"target    : {len(stations)*len(gcms)*len(scens)*len(arms)} run directories\n")

    import collections
    how = collections.Counter()
    lu_mismatch = []
    written, ok_arms, blocked, runs = 0, 0, [], []
    ic_table = read_ic_table(a.ic) if a.ic else None
    if ic_table is not None:
        print(f"initial state: {a.ic}  ({len(ic_table)} row(s)), key '{a.ic_key}'")
    # Consistency over coverage: a GCM station with no completed ERA5 run has no
    # state to seed the spin-up from, so it would silently be the only station in
    # the fleet whose pools were never equilibrated. Drop it, by name.
    if a.require_era5_state:
        if ic_table is None:
            print("ERROR: --require-era5-state needs --ic", file=sys.stderr)
            return 1
        have = {st for (st, k) in ic_table if k == "era5_land/fixed_lma"}
        dropped = [s["station_id"] for s in stations if s["station_id"] not in have]
        stations = [s for s in stations if s["station_id"] in have]
        print(f"era5 states  : {len(have)}   dropped {len(dropped)} station(s) "
              f"with no completed ERA5 run")
        for d in dropped:
            print(f"  - {d}")
        print()
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
            if hd.get("how") == "own cell":
                offsets.append((sid, hd["pixel_km"]))
            how[hd.get("how", "?").split(" (")[0]] += 1
            if hd.get("how") != "own cell":
                lu_mismatch.append((sid, st["forest_type"], hd.get("how"),
                                    hd.get("n_pixels"), hd.get("cell_base"),
                                    hd.get("eco_base"), hd.get("lu_used"),
                                    hd.get("d_nearest_any"), hd.get("d_nearest_same"),
                                    hd.get("grid_deg")))

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
                # The forcing lives IN the run tree, at <ST>/<scenario>/<GCM>/,
                # one level above the two arms that share it. --meteo is the
                # legacy staging path, still consulted so a tree that has not
                # been through migrate_forcing.py says so instead of silently
                # blocking every combination.
                fname = f"Meteo_{mname}_{mat_name(gcm)}_{scen}_{YEAR_TAG[scen]}.mat"
                mfile = a.root / sid / scen / gcm / fname
                if not mfile.is_file():
                    legacy = a.meteo / scen / mat_name(gcm) / fname
                    blocked.append((sid, gcm, scen,
                                    f"forcing still in {legacy.parent} -- run "
                                    f"migrate_forcing.py" if legacy.is_file()
                                    else f"no forcing {fname}"))
                    continue
                if ic_table is not None:
                    key = a.ic_key.format(station=sid, scenario=scen, gcm=gcm)
                    if (sid, key) not in ic_table:
                        blocked.append((sid, gcm, scen,
                                        f"no harvested state '{key}' in {a.ic}"))
                        continue
                n, err = build_one(st, gcm, scen, fixed_mean, series, mfile,
                                   a.root, era5_mp, a.dry_run,
                                   a.arms, ic_table, a.ic_key)
                if err:
                    blocked.append((sid, gcm, scen, err))
                    continue
                written += n
                ok_arms += len(arms)
                for arm in arms:
                    runs.append(f"{sid} {scen} {gcm} {arm}")

    print(f"{'=' * 72}")
    verb = "would be built" if a.dry_run else "written"
    print(f"{'DRY RUN -- ' if a.dry_run else ''}{ok_arms} run directories {verb}, "
          f"{len(blocked)} combination(s) blocked")
    print(f"{'=' * 72}")

    if how:
        print("\nHOW EACH (station, GCM) PAIR GOT ITS SERIES")
        tot = sum(how.values())
        for k, n in how.most_common():
            print(f"  {n:>5d}  ({100*n/tot:4.1f}%)  {k}")

    if lu_mismatch:
        seen_lu = {}
        for sid, ft, howstr, npx, cb, eb, luu, da, ds, gd in lu_mismatch:
            seen_lu[sid] = (ft, howstr, npx, cb, eb, luu, da, ds, gd)
        print(f"\nBASELINE SUBSTITUTED FROM THE ECOREGION ({len(seen_lu)} station(s))")
        print("  The cell the station sits in is not mapped as its forest type. Its")
        print("  DYNAMICS are still used; only the LEVEL comes from the ecoregion mean")
        print("  of the station's own type -- the substitution build_lma_input.py makes")
        print("  for the ERA5-Land runs:")
        print("      mu_Y_row = np.where(in_fit, mu_Y_pix_final[idx], mu_Y_eco)")
        print(f"\n    {'station':10s}{'type':11s}{'cell LU':>8s}{'cell base':>11s}"
              f"{'eco base':>10s}{'shift':>8s}  n_cells")
        for sid, (ft, howstr, npx, cb, eb, luu, da, ds, gd) in sorted(seen_lu.items()):
            if cb is None:
                print(f"    {sid:10s}{ft:11s}{'-':>8s}{'-':>11s}"
                      f"{(eb if eb is not None else float('nan')):>10.1f}{'-':>8s}"
                      f"  {npx:>4d}   OUTSIDE: nearest cell {da:.1f} km "
                      f"(same type {ds:.1f} km), grid {gd} deg")
            else:
                print(f"    {sid:10s}{ft:11s}{luu:>8d}{cb:>11.1f}{eb:>10.1f}"
                      f"{eb - cb:>+8.1f}  {npx:>4d}")
        print("\n  One departure from the ERA5 path: there the station's OWN forest-type"
              "\n  beta is applied to the pixel's climate. These projections ship as"
              "\n  finished values, so the anomaly carries the model of the cell's own"
              "\n  type. The level is corrected; the sensitivity is not.")

    if offsets:
        o = np.array([d for _, d in offsets])
        print("\nSTATION-TO-CELL-CENTRE DISTANCE (own-cell matches only)")
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
        name = ("run_list_gcm.txt" if arms == ARMS
                else f"run_list_gcm_{'_'.join(arms)}.txt")
        lst = a.root / name
        lst.write_text("".join(r + "\n" for r in runs), encoding="utf-8")
        print(f"\nrun list : {lst}  ({len(runs)} arms)")
        print(f"           sbatch --array=1-{len(runs)}%NN slurm/submit_tc_run.sh "
              f"(with RUN_LIST={name})")
    return 0 if ok_arms else 1


if __name__ == "__main__":
    sys.exit(main())
