#!/usr/bin/env python3
"""Where in CONUS each LMA metric is large, one map figure per metric.

Five figures -- lma, flux, sd, spei, ta -- each a 6x4 grid of CONUS maps:

    rows    GPP  LAI  TR  ET  LE  H          the six fluxes
    columns ERA5-Land | GCM historical | ssp126 | ssp585

THIS IS THE SCATTER FIGURES' DATA ON A MAP. The per-station value is produced
by figure_dryness.y_values, imported rather than reimplemented, so a station's
colour here and its height in the corresponding scatter panel are the same
number by construction. Reimplementing the aggregation would let the two drift
apart silently, which is exactly the failure the median-across-GCMs rule exists
to prevent. Median across GCMs at the METRIC level, as everywhere else.

COLOUR IS PuOr, SYMMETRIC ABOUT ZERO, ONE SCALE PER ROW. Zero is white because
zero is the null result -- dynamic LMA changed nothing -- and a diverging scale
is the only kind that shows sign at a glance. The four columns of a row share
one scale so they can be read against each other, which is the entire reason
they sit side by side; rows do not, because a GPP slope and an LAI slope are
different units. Limits are the 98th percentile of |value| pooled over the row's
four datasets, with the bar extended at both ends, so a single extreme station
cannot flatten the other 81.

MARKER SIZE IS THE STATION'S FLUX MAGNITUDE, IN GLOBAL QUARTILES. The quartile
edges are computed per row over every station in all four datasets at once --
"global", not per panel -- so a large marker means the same flux in the ssp585
column as in the ERA5 one. Per-panel quartiles would put a quarter of the
stations in the top bin of every panel by construction and the size channel
would carry no information across columns. Magnitude is |mean annual flux| from
the FIXED arm, the control; the lma metric uses the DYNAMIC arm instead, since
dFlux/dLMA is a dynamic-arm quantity and has no fixed-arm counterpart.

A station with a metric but no flux magnitude is drawn at the MEDIAN size and
counted in the log -- the median is the bin that claims least about it. Dropping
it would remove a real result from the map over a missing decoration.

ONLY THE 82 STATIONS COMPLETE IN ALL FOUR DATASETS ARE DRAWN (--fleet common,
the default). The datasets do not cover the same stations -- 92 / 85 / 82 / 82 --
and on a map an unequal sample is not a footnote: the ERA5 column would carry ten
dots the ssp columns lack, in particular US-CPk, US-CZ2, US-MtB and US-SHC, and
a reader comparing columns would see that absence as a geographic pattern in the
result. Comparing columns is what the layout is for, so the columns must hold the
same stations. --fleet all keeps every station and prints the per-dataset counts.
NOTE that the scatter figures do not apply this filter, so a map and its scatter
can differ by those ten stations in the ERA5 and historical panels.

Writes map_stations.csv beside the figures: every point in every panel, which is
the table_map_stations product. It is emitted here rather than from
figure_tables.py so that the table and the figures cannot disagree.

Prerequisites: station_metrics.csv and station_sensitivity.csv from
station_metrics.py, and the annual effect tables for the flux metric.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_dryness import (Missing, PFTS, VARS, YAXIS,           # noqa: E402
                            FLUX_UNITS, y_values)
from results_dir import NoResultsDir, resolve_figure, resolve_out  # noqa: E402
from station_metrics import SITE_LISTS                            # noqa: E402

PANELS = [("era5", "ERA5-Land"), ("historical", "GCM historical"),
          ("ssp126", "GCM ssp126"), ("ssp585", "GCM ssp585")]

# ColorBrewer PuOr anchors, white at the centre. Orange and purple stay
# separable under deuteranopia and protanopia where a red-green pair does not.
PUOR = ["#7f3b08", "#b35806", "#e08214", "#fee0b6", "#ffffff",
        "#d8daeb", "#b2abd2", "#8073ac", "#542788", "#2d004b"]

# Marker areas in points^2 for flux-magnitude quartiles Q1..Q4.
SIZES = [8.0, 19.0, 36.0, 60.0]

# A FIXED CONUS FRAME, in EPSG:5070 metres, on every panel. Letting matplotlib
# autoscale to the points makes each panel's extent depend on which stations it
# happens to hold, so the same station sits at a different place in the ERA5 and
# ssp585 columns and the panels cannot be read against each other -- and with a
# station set that spans the continent east-west but not north-south, the
# autoscaled panels come out as flat strips. The 82 stations span
# x -2.34e6..2.11e6, y 6.7e5..2.81e6, so this frame contains all of them with
# room for the coastline.
XLIM = (-2.45e6, 2.30e6)
YLIM = (4.0e5, 3.15e6)


def read_site_coords() -> pd.DataFrame:
    """station -> lat, lon, pft, from the same two lists the runs were built from."""
    rows = []
    for path in SITE_LISTS:
        if not path.exists():
            raise Missing(f"site list not found: {path}")
        with path.open(newline="", encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                sid = (r.get("StationID") or "").strip()
                try:
                    lat, lon = float(r["Lat"]), float(r["Lon"])
                except (KeyError, TypeError, ValueError):
                    continue
                if sid and np.isfinite(lat) and np.isfinite(lon):
                    rows.append({"station": sid, "lat": lat, "lon": lon,
                                 "site_pft": (r.get("ForestType") or "").strip().lower()})
    if not rows:
        raise Missing("no station coordinates could be read from the site lists")
    return pd.DataFrame(rows).drop_duplicates("station")


def flux_magnitude(root: Path, arm: str) -> pd.DataFrame:
    """|mean annual flux| per (dataset, station, variable), median over GCMs.

    arm is "fixed" or "dyn". Annual, subset "all" -- the same cell the metrics
    themselves come from, so the size and the colour describe one population of
    station-years rather than two.
    """
    p = Path(root) / "station_metrics.csv"
    if not p.is_file():
        raise Missing("station_metrics.csv not found -- run station_metrics.py first")
    col = f"mean_{arm}"
    m = pd.read_csv(p, low_memory=False)
    for need in ("dataset", "station", "variable", "freq", "subset", col):
        if need not in m.columns:
            raise Missing(f"station_metrics.csv has no {need} column")
    m = m[(m["freq"] == "annual") & (m["subset"] == "all")].copy()
    if m.empty:
        raise Missing("station_metrics.csv has no annual/all rows")
    m["mag"] = pd.to_numeric(m[col], errors="coerce").abs()
    return (m.groupby(["dataset", "station", "variable"], as_index=False)["mag"]
             .median())


def quartile_bins(mag: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Bin index 0..3 per value, and the three global edges used.

    Edges come from the pooled distribution, so they are identical in all four
    panels of a row. NaN magnitudes return bin 1 -- the median bin, chosen
    because it asserts least -- and are counted by the caller.
    """
    good = mag[np.isfinite(mag)]
    if good.size < 4 or np.unique(good).size < 2:
        return np.full(mag.shape, 1, int), np.array([np.nan] * 3)
    edges = np.nanpercentile(good, [25, 50, 75])
    idx = np.digitize(mag, edges, right=False)
    idx = np.where(np.isfinite(mag), idx, 1)
    return np.clip(idx, 0, 3).astype(int), edges


def read_outline(path: Path, max_pts: int = 220) -> list:
    """CONUS land polygons in EPSG:5070, as a list of (N,2) rings.

    pyshp AND pyproj, NOT geopandas. requirements.txt excludes geopandas on
    purpose -- "optional and only used for the --basemap outline, so it is
    deliberately NOT listed" -- so the geopandas call this replaced could never
    run on the cluster: passing a valid --basemap raised ImportError and killed
    the job, and passing none printed a note. Either way there was no outline,
    and no combination of arguments could produce one. pyshp and pyproj are both
    listed, so this draws with the venv as it already is.

    NO DISSOLVE. geopandas got the silhouette by dissolving ~1250 Level III
    polygons into one; without shapely that is expensive to do properly. Filling
    every polygon in the same flat grey gives the identical picture, because the
    union of the fills IS the landmass and the seams are invisible when the
    edges are not drawn.

    Rings are decimated to max_pts. At 220 dpi a panel is ~440 px across 4750 km,
    so one pixel is ~11 km and vertices finer than that cannot be seen; keeping
    all ~500k of them would render 24 panels of invisible detail.
    """
    import shapefile                                            # pyshp
    from pyproj import CRS, Transformer

    p = Path(path)
    if not p.is_file():
        raise Missing(f"--basemap {p} does not exist. The EPA shapefile is "
                      f"downloaded by slurm/submit_verify_pairing.sh to "
                      f"$TC_INPUT_DATA/ecoregions/us_eco_l3.shp")
    prj = p.with_suffix(".prj")
    if not prj.is_file():
        raise Missing(f"no .prj beside {p.name}; the shapefile's projection "
                      f"cannot be read and the outline would be misplaced")
    # The EPA distributes us_eco_l3 in USA Contiguous Albers (GRS80), which is
    # NOT the same datum as EPSG:5070 (NAD83). Reproject from the file's own
    # .prj rather than assuming they match.
    tr = Transformer.from_crs(CRS.from_wkt(prj.read_text()), "EPSG:5070",
                              always_xy=True)
    rings = []
    for shp in shapefile.Reader(str(p)).iterShapes():
        pts = np.asarray(shp.points, float)
        if pts.size == 0:
            continue
        bounds = list(shp.parts) + [len(pts)]
        for a, b in zip(bounds[:-1], bounds[1:]):
            ring = pts[a:b]
            if len(ring) < 3:
                continue
            if len(ring) > max_pts:
                idx = np.linspace(0, len(ring) - 1, max_pts).astype(int)
                ring = ring[idx]
            x, y = tr.transform(ring[:, 0], ring[:, 1])
            rings.append(np.column_stack([x, y]))
    if not rings:
        raise Missing(f"{p.name} yielded no polygons")
    print(f"  basemap: {len(rings)} rings from {p.name}")
    return rings


def albers(lat, lon):
    """CONUS Albers Equal Area (EPSG:5070), one transformer for the whole run."""
    global _TR
    try:
        tr = _TR
    except NameError:
        from pyproj import Transformer
        tr = _TR = Transformer.from_crs("EPSG:4326", "EPSG:5070", always_xy=True)
    x, y = tr.transform(np.asarray(lon, float), np.asarray(lat, float))
    return np.asarray(x, float), np.asarray(y, float)


def build(d: pd.DataFrame, yk: dict, out_png: Path, figsize, outline: list | None,
          vmax_pct: float) -> pd.DataFrame:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.lines import Line2D

    cmap = LinearSegmentedColormap.from_list("puor", PUOR)

    fig, axes = plt.subplots(len(VARS), len(PANELS), figsize=figsize,
                             constrained_layout=True, squeeze=False)
    drawn, unsized, records = 0, 0, []
    for i, (disp, col) in enumerate(VARS):
        row = d[d["variable"] == col]
        v = row["value"].to_numpy(float)
        v = v[np.isfinite(v)]
        # Symmetric about zero: a diverging scale whose white is not at zero
        # would put the null result somewhere in the orange.
        vmax = float(np.nanpercentile(np.abs(v), vmax_pct)) if v.size else 1.0
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = float(np.nanmax(np.abs(v))) if v.size and np.isfinite(v).any() else 1.0
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = 1.0
        norm = Normalize(vmin=-vmax, vmax=vmax)

        # ONE set of quartile edges for the whole row, all four datasets pooled.
        bins, edges = quartile_bins(row["mag"].to_numpy(float))
        row = row.assign(size_bin=bins)
        unsized += int((~np.isfinite(row["mag"].to_numpy(float))).sum())

        for j, (ds, label) in enumerate(PANELS):
            ax = axes[i, j]
            if outline is not None:
                # A new collection per axes -- matplotlib will not share one --
                # but the projected vertex arrays behind them are read once.
                ax.add_collection(PolyCollection(
                    outline, facecolors="#eaeef2", edgecolors="none", zorder=1))
            sub = row[row["dataset"] == ds]
            for pft, _colr, mk in PFTS:
                w = sub[sub["pft"] == pft]
                if w.empty:
                    continue
                vals = w["value"].to_numpy(float)
                ok = np.isfinite(vals)
                if not ok.any():
                    continue
                x, y = albers(w["lat"].to_numpy(float)[ok],
                              w["lon"].to_numpy(float)[ok])
                s = np.array(SIZES)[w["size_bin"].to_numpy(int)[ok]]
                ax.scatter(x, y, c=vals[ok], cmap=cmap, norm=norm, marker=mk,
                           s=s, edgecolor="0.2", linewidth=0.3, zorder=3)
                drawn += int(ok.sum())
                for st, val, sb in zip(w["station"].to_numpy()[ok], vals[ok],
                                       w["size_bin"].to_numpy(int)[ok]):
                    records.append({"metric": yk["prefix"], "dataset": ds,
                                    "variable": disp, "station": st, "pft": pft,
                                    "value": float(val), "size_bin": int(sb) + 1})
            # Equal aspect, or Albers CONUS stretches to fill the panel and the
            # geography is quietly wrong. Fixed limits, so every panel is the
            # same map.
            ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
            ax.set_aspect("equal")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_linewidth(0.4); sp.set_color("0.7")
            if i == 0:
                ax.set_title(label, fontsize=8.5, pad=3)
            if j == 0:
                ax.set_ylabel(disp, fontsize=9, labelpad=4)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cb = fig.colorbar(sm, ax=list(axes[i, :]), fraction=0.02, pad=0.008,
                          shrink=0.88, aspect=22, extend="both")
        cb.ax.tick_params(labelsize=5.6)
        unit = ("" if yk["per_unit"].startswith(" (")
                else FLUX_UNITS.get(disp, "")) + yk["per_unit"]
        cb.set_label(unit.strip(), fontsize=5.6)
        if np.isfinite(edges).all():
            print(f"  {disp:<4} vmax {vmax:>10.4g}   size edges "
                  f"{edges[0]:.4g} / {edges[1]:.4g} / {edges[2]:.4g}")
        else:
            print(f"  {disp:<4} vmax {vmax:>10.4g}   size edges unavailable")

    if drawn == 0:
        raise Missing(f"{yk['prefix']}: every panel empty")

    pft_h = [Line2D([], [], marker=mk, ls="", color="0.3", markersize=5, label=p)
             for p, _c, mk in PFTS]
    size_h = [Line2D([], [], marker="o", ls="", color="0.55",
                     markersize=np.sqrt(s), label=q)
              for s, q in zip(SIZES, ["Q1", "Q2", "Q3", "Q4"])]
    arm = "dynamic" if yk["prefix"] == "lma" else "fixed"
    # Below the axes, not on them: constrained_layout reserves no room for a
    # figure legend, so at -0.012 it sat on top of the H row's panels.
    # bbox_inches="tight" grows the saved canvas to include it.
    fig.legend(handles=pft_h + size_h, loc="upper center", ncol=6, frameon=False,
               fontsize=7, bbox_to_anchor=(0.5, 0.0),
               title=f"marker size: global quartiles of the {arm}-arm mean "
                     f"annual flux, per row",
               title_fontsize=6.5)
    fig.suptitle(yk["label"], fontsize=10.5, x=0.01, ha="left")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    note = f"   ({unsized} without a flux magnitude)" if unsized else ""
    print(f"  -> {out_png}   ({drawn} points){note}")
    return pd.DataFrame(records)


def common_fleet(d: pd.DataFrame) -> pd.DataFrame:
    """Keep only stations present in every dataset, so the columns are comparable."""
    n = d.groupby("station")["dataset"].nunique()
    keep = set(n[n == d["dataset"].nunique()].index)
    dropped = sorted(set(d["station"]) - keep)
    if dropped:
        print(f"  fleet: {len(keep)} stations in all {d['dataset'].nunique()} "
              f"datasets; dropped {len(dropped)} incomplete "
              f"({', '.join(dropped[:6])}{', ...' if len(dropped) > 6 else ''})")
    return d[d["station"].isin(keep)]


def load(root: Path, yk: dict, coords: pd.DataFrame,
         fleet: str = "common") -> pd.DataFrame:
    """One row per (dataset, station, pft, variable): value, magnitude, lat/lon."""
    s = y_values(root, yk).rename(columns={"slope": "value"})
    s = s[s["variable"].isin([c for _, c in VARS])].copy()
    if s.empty:
        raise Missing(f"{yk['prefix']}: no rows for the six flux variables")
    arm = "dyn" if yk["prefix"] == "lma" else "fixed"
    mag = flux_magnitude(root, arm)
    out = s.merge(mag, on=["dataset", "station", "variable"], how="left")
    out = out.merge(coords, on="station", how="inner")
    if out.empty:
        raise Missing("no station matched between the metric table and the site "
                      "lists -- check the station IDs")
    # pft comes from the metric table; site_pft from the lists. They are the
    # same field via two routes, so a disagreement means one of them is stale.
    out["pft"] = out["pft"].fillna(out["site_pft"])
    bad = out[(out["site_pft"] != "") & (out["pft"] != out["site_pft"])]
    if len(bad):
        print(f"  WARNING: {bad['station'].nunique()} station(s) disagree on "
              f"forest type between the metric table and the site lists: "
              f"{sorted(bad['station'].unique())[:5]}", file=sys.stderr)
    lost = s["station"].nunique() - out["station"].nunique()
    if lost:
        print(f"  note: {lost} station(s) have a metric but no coordinates")
    out = out.drop(columns=["site_pft"])
    if fleet == "common":
        out = common_fleet(out)
    else:
        print("  fleet: " + "  ".join(
            f"{ds}={g['station'].nunique()}" for ds, g in out.groupby("dataset")))
    if out.empty:
        raise Missing("no station survived the fleet filter")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--metrics", default=",".join(YAXIS),
                    help="comma-separated subset of " + ",".join(YAXIS))
    ap.add_argument("--basemap", type=Path, default=None,
                    help="us_eco_l3 shapefile for the CONUS outline; without it "
                         "the panels still plot but are much harder to read")
    ap.add_argument("--size", default="10x9", help="WxH inches")
    ap.add_argument("--vmax-pct", type=float, default=98.0,
                    help="percentile of |value| setting each row's colour limit")
    ap.add_argument("--fleet", default="common", choices=["common", "all"],
                    help="common (default): only stations present in all four "
                         "datasets, so the columns hold the same sample")
    ap.add_argument("--csv", default="map_stations.csv",
                    help="per-point table written beside the figures")
    a = ap.parse_args(argv)

    w, h = (float(x) for x in a.size.lower().split("x"))
    try:
        root = Path(a.results or resolve_out(".", create=False))
        out_dir = Path(a.out or resolve_figure("."))
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    out_dir.mkdir(parents=True, exist_ok=True)

    want = [m.strip() for m in a.metrics.split(",") if m.strip()]
    unknown = [m for m in want if m not in YAXIS]
    if unknown:
        print(f"ERROR: unknown metric(s) {unknown}; choose from "
              f"{list(YAXIS)}", file=sys.stderr)
        return 1

    # A WRONG PATH IS AN ERROR, NOT A NOTE. This used to print "no --basemap"
    # whether the flag was absent or pointed at a file that was not there, so a
    # typo and a deliberate omission produced the same line and the same
    # outline-less figures.
    if a.basemap is None:
        print("note: no --basemap given, so the panels have no CONUS outline.\n"
              "      The EPA shapefile is downloaded by "
              "slurm/submit_verify_pairing.sh to\n"
              "      $TC_INPUT_DATA/ecoregions/us_eco_l3.shp")
    elif not Path(a.basemap).is_file():
        print(f"ERROR: --basemap {a.basemap} does not exist.\n"
              f"       Run slurm/submit_verify_pairing.sh to download it, or "
              f"omit --basemap to draw without an outline.", file=sys.stderr)
        return 1

    try:
        coords = read_site_coords()
    except Missing as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    # A station outside the fixed frame would be clipped away with no trace on
    # the figure, which is the one failure a map cannot show you.
    cx, cy = albers(coords["lat"].to_numpy(float), coords["lon"].to_numpy(float))
    off = coords[(cx < XLIM[0]) | (cx > XLIM[1]) |
                 (cy < YLIM[0]) | (cy > YLIM[1])]
    print(f"station coordinates: {len(coords)}")
    if len(off):
        print(f"ERROR: {len(off)} station(s) fall outside the CONUS frame and "
              f"would be clipped invisibly: "
              f"{', '.join(off['station'].head(8))}", file=sys.stderr)
        return 1

    # Read the 41 MB shapefile ONCE, not once per metric: the rings are the same
    # in all five figures and reprojecting them five times is pure waste.
    try:
        outline = read_outline(a.basemap) if a.basemap else None
    except Missing as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    tables, failures = [], []
    for m in want:
        yk = YAXIS[m]
        print(f"\n{m}: {yk['label']}")
        try:
            d = load(root, yk, coords, a.fleet)
            print(f"  stations {d['station'].nunique()}   rows {len(d)}")
            tables.append(build(d, yk, out_dir / f"map_{yk['prefix']}.png",
                                (w, h), outline, a.vmax_pct))
        except Missing as e:
            failures.append(f"{m}: {e}")
            print(f"SKIP {m}: {e}")

    if tables:
        csv_path = out_dir / a.csv
        pd.concat(tables, ignore_index=True).to_csv(csv_path, index=False)
        print(f"\n-> {csv_path}")
    if failures:
        print(f"\n{len(failures)} figure(s) not produced:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
