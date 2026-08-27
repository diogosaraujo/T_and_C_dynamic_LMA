"""The three hourly figures: skill maps, Taylor diagram, RSR/SS error bars.

All read hourly_stats.csv and each writes ONE figure. There is no per-timestep
loop any more, because annual and seasonal skill were never defensible: the
model record is 1985-2020 and many towers start after 2015, so at US-HBK the
annual overlap is two years and a correlation from two points is +/-1 by
construction. Hourly gives tens of thousands of matched steps everywhere.

  maps    4 rows x 3 cols -- CONUS map, then the SS distribution for deciduous
          and for evergreen, so the map's spatial pattern and the fleet's
          distribution sit side by side.
  taylor  4 rows x 1 col.
  (the RSR/SS error bars are gone: the skill change is too small for a
   three-bar comparison to show anything the maps do not.)

THE SKILL SCORE IS THE SAME NUMBER EITHER WAY. 1 - RSR_dyn/RSR_fixed equals
1 - RMSE_dyn/RMSE_fixed identically, because the observed standard deviation
divides both arms and cancels. The colour bar says so, since the question keeps
coming up. RSR earns its place on the error-bar axis, where the raw RMSE would
carry units and make the four rows incomparable.

TAYLOR CIRCLES ARE LABELLED AS WHAT THEY ARE. The arcs centred on the tower are
contours of the centred RMS difference, normalised by the observed standard
deviation -- so they are RSR contours, and they are labelled with their values
rather than left to "closer is better". Note "centred": the mean is removed
first, so a Taylor diagram cannot show bias. That is what the bias column of
hourly_stats.csv is for.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_skill_maps import albers_xy, read_sites               # noqa: E402
from results_dir import NoResultsDir, resolve_figure, resolve_out  # noqa: E402

ROWS = [("GPP", "GPP"), ("ET", "ET"), ("H", "H"), ("LE", "LE")]
C_DYN, C_FIX, C_SS = "#b2182b", "#4d4d4d", "#2166ac"
RSR_VERY_GOOD = 0.5
# Natural Earth 110m state outlines: 48 kB, no authentication, and it does not
# depend on the EPA shapefile, whose old gaftp path now 404s. Fetched once and
# cached, so only the first run needs the network.
NE_URL = ("https://naciscdn.org/naturalearth/110m/cultural/"
          "ne_110m_admin_1_states_provinces.zip")
C_DEC, C_EVE = "#1b7837", "#2166ac"      # deciduous green, evergreen blue


def get_basemap(explicit, cache_dir: Path):
    """The CONUS outline, downloaded and cached if it is not already there.

    Earlier versions searched the filesystem for us_eco_l3.shp and silently drew
    no coastline when it was absent, which is how several rounds of maps shipped
    as bare scatter plots.
    """
    if explicit:
        return Path(explicit) if Path(explicit).exists() else None
    cache_dir.mkdir(parents=True, exist_ok=True)
    shp = cache_dir / "ne_110m_admin_1_states_provinces.shp"
    if not shp.exists():
        import io
        import urllib.request
        import zipfile
        try:
            req = urllib.request.Request(NE_URL, headers={"User-Agent": "research/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                blob = r.read()
            zipfile.ZipFile(io.BytesIO(blob)).extractall(cache_dir)
        except Exception as e:                                   # noqa: BLE001
            print(f"  basemap download failed: {type(e).__name__}: {e}")
            return None
    return shp if shp.exists() else None


def read_conus(shp: Path):
    """CONUS state rings in Albers metres, via pyshp and pyproj."""
    import shapefile
    from pyproj import Transformer
    tr = Transformer.from_crs("EPSG:4326", "EPSG:5070", always_xy=True)
    sf = shapefile.Reader(str(shp))
    flds = [f[0] for f in sf.fields[1:]]
    def col(rec, *names):
        for n in names:
            if n in flds:
                return str(rec[flds.index(n)])
        return ""
    rings = []
    for sr in sf.shapeRecords():
        who = col(sr.record, "admin", "sov_a3", "iso_a2")
        if not any(k in who for k in ("United States", "USA", "US")):
            continue
        # Skip Alaska and Hawaii: in frame they shrink the mainland to a strip.
        nm = col(sr.record, "name", "gn_name")
        if nm in ("Alaska", "Hawaii"):
            continue
        pts, parts = sr.shape.points, list(sr.shape.parts) + [len(sr.shape.points)]
        for a, b in zip(parts[:-1], parts[1:]):
            lon = [p[0] for p in pts[a:b]]
            lat = [p[1] for p in pts[a:b]]
            if len(lon) < 3:
                continue
            x, y = tr.transform(lon, lat)
            rings.append(list(zip(x, y)))
    return rings or None


def load(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path)
    need = {"station", "pft", "variable", "arm", "n", "rsr", "skill_score", "r",
            "sd_mod", "sd_obs"}
    miss = need - set(d.columns)
    if miss:
        raise SystemExit(f"ERROR: {path.name} lacks {', '.join(sorted(miss))}")
    d["pft"] = d["pft"].astype(str).str.strip().str.lower()
    return d


def wide(d: pd.DataFrame, var: str) -> pd.DataFrame:
    """One row per station for one variable, both arms side by side."""
    s = d[d["variable"] == var]
    f = s[s.arm == "fixed"].set_index("station")
    y = s[s.arm == "dyn"].set_index("station")
    out = pd.DataFrame({
        "pft": f["pft"], "n": f["n"],
        "rsr_fixed": f["rsr"], "rsr_dyn": y["rsr"].reindex(f.index),
        "r_fixed": f["r"], "r_dyn": y["r"].reindex(f.index),
        "sdr_fixed": f["sd_mod"] / f["sd_obs"],
        "sdr_dyn": (y["sd_mod"] / y["sd_obs"]).reindex(f.index),
        "ss": f["skill_score"]})
    return out.dropna(subset=["rsr_fixed"])


# --------------------------------------------------------------------- maps
def fig_maps(d, sites, out_png, basemap, figsize, linthresh, vmax):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm, LinearSegmentedColormap
    from matplotlib.lines import Line2D

    cmap = LinearSegmentedColormap.from_list(
        "ss_puor", ["#7f3b08", "#e08214", "#fee0b6", "#ffffff",
                    "#d8daeb", "#8073ac", "#2d004b"])
    norm = SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax, base=10)
    # PYSHP + PYPROJ, NOT GEOPANDAS. Job 39708 died on
    # "No module named 'geopandas'": it was deliberately left out of
    # requirements when the outline was optional, and then the outline became
    # core. pyshp and pyproj are already in the venv and do the whole job --
    # geopandas would drag in a stack for one polygon read.
    outline = read_conus(basemap) if basemap is not None else None

    fig, axes = plt.subplots(len(ROWS), 3, figsize=figsize,
                             constrained_layout=True,
                             gridspec_kw={"width_ratios": [2.1, 1, 1]})
    for i, (label, var) in enumerate(ROWS):
        w = wide(d, var)
        ax = axes[i, 0]
        if outline is not None:
            from matplotlib.patches import Polygon as MplPoly
            # One flat grey landmass, no state borders. The edge is drawn in the
            # SAME grey rather than left off: with edgecolor="none" the
            # antialiased seams between abutting states show as hairline gaps
            # across the fill.
            for ring in outline:
                ax.add_patch(MplPoly(ring, closed=True, facecolor="#ececec",
                                     edgecolor="#ececec", linewidth=0.6,
                                     zorder=0))
            # CONUS only: Alaska and Hawaii would shrink the mainland to a strip.
            ax.set_xlim(-2.4e6, 2.4e6); ax.set_ylim(2.5e5, 3.3e6)
        for mk, kind in (("o", "deciduous"), ("^", "evergreen")):
            sub = w[w.pft == kind]
            xs, ys, cs = [], [], []
            for sid, row in sub.iterrows():
                if sid not in sites or not np.isfinite(row["ss"]):
                    continue
                lat, lon, _ = sites[sid]
                x, y = albers_xy(lat, lon)
                xs.append(float(x)); ys.append(float(y)); cs.append(row["ss"])
            if xs:
                ax.scatter(xs, ys, c=cs, cmap=cmap, norm=norm, marker=mk, s=34,
                           edgecolor="0.25", linewidth=0.35, zorder=3)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.4); sp.set_color("0.6")
        ax.set_ylabel(label, fontsize=10)

        # Columns 2 and 3: the SS distribution, one forest type each.
        for j, kind in enumerate(("deciduous", "evergreen"), start=1):
            axh = axes[i, j]
            col = C_DEC if kind == "deciduous" else C_EVE
            v = w.loc[w.pft == kind, "ss"].to_numpy(float)
            v = v[np.isfinite(v)]
            if v.size:
                lim = max(0.05, np.nanpercentile(np.abs(v), 98) * 1.2)
                bins = np.linspace(-lim, lim, 19)
                # Transparent fill under a solid outline: the bars stay readable
                # where they cross the zero reference and the median line.
                axh.hist(v, bins=bins, color=col, alpha=.28, zorder=2)
                axh.hist(v, bins=bins, histtype="step", color=col,
                         linewidth=.9, zorder=3)
                axh.axvline(0, color="0.35", lw=.8, ls=":", zorder=1)
                axh.axvline(np.median(v), color=col, lw=1.6, zorder=4)
                axh.text(.03, .93, f"n={v.size}" + chr(10)
                         + f"med {np.median(v):+.3f}",
                         transform=axh.transAxes, fontsize=6, va="top",
                         color=col)
            axh.tick_params(labelsize=6)
            axh.set_yticks([])
            for sp in ("top", "right", "left"):
                axh.spines[sp].set_visible(False)
            if i == 0:
                axh.set_title(kind, fontsize=9, color=col)
            if i == len(ROWS) - 1:
                axh.set_xlabel("skill score", fontsize=7)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=axes[:, 0].tolist(), orientation="horizontal",
                      fraction=0.05, pad=0.02, extend="both")
    # The default SymLog ticks crowd -1e-2 and 1e-2 into each other around zero.
    ticks = [-1, -1e-1, -1e-2, 0, 1e-2, 1e-1, 1]
    cb.set_ticks(ticks)
    cb.set_ticklabels(["-1", "-0.1", "-0.01", "0", "0.01", "0.1", "1"])
    cb.ax.tick_params(labelsize=7)
    # Short enough not to be clipped by the colour bar's own width; the
    # RSR/RMSE equivalence goes underneath instead of on one long line.
    cb.set_label("skill score   1 - RMSE$_{dyn}$/RMSE$_{fixed}$", fontsize=8)
    handles = [Line2D([], [], marker="o", ls="", color="0.35", label="deciduous"),
               Line2D([], [], marker="^", ls="", color="0.35", label="evergreen")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=7, bbox_to_anchor=(0.28, -0.055))
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------- taylor
def _mini_legend(ax, C_FIX, C_DYN):
    """A small Taylor sketch used instead of a sentence of axis prose."""
    a = np.linspace(0, np.pi / 2, 120)
    for r in (0.5, 1.0, 1.5):
        ax.plot(r * np.cos(a), r * np.sin(a), color="0.85", lw=.4)
    for rr in (0.0, 0.6, 0.9, 0.99):
        th = np.arccos(rr)
        ax.plot([0, 1.6 * np.cos(th)], [0, 1.6 * np.sin(th)], color="0.85", lw=.4)
    th = np.linspace(0, np.pi, 120)
    for rms in (0.5, 1.0):
        ax.plot(1 + rms * np.cos(th), rms * np.sin(th), color="#9ecae1",
                lw=.5, ls=(0, (2, 2)))
    ax.plot([1], [0], marker="*", ms=8, color="0.15")
    ax.annotate("", xy=(1.05, 0.75), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color="0.3", lw=.8))
    ax.text(0.42, 0.48, "sd(model)/sd(tower)", fontsize=5.5, rotation=35,
            color="0.25", ha="center")
    ax.text(0.34, 1.30, "correlation", fontsize=5.5, color="0.25", ha="center")
    ax.text(1.42, 0.42, "RMSD", fontsize=5.5, color="#2166ac", ha="center")
    ax.text(1.02, -0.16, "tower", fontsize=5.5, color="0.15", ha="center")
    ax.set_xlim(-0.05, 1.85); ax.set_ylim(-0.28, 1.7)
    ax.set_aspect("equal"); ax.axis("off")


def fig_taylor(d, out_png, figsize):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import figure_taylor as T

    allr, allsd = [], []
    for _, v in ROWS:
        w = wide(d, v)
        allr += w[["r_fixed", "r_dyn"]].to_numpy(float).ravel().tolist()
        allsd += w[["sdr_fixed", "sdr_dyn"]].to_numpy(float).ravel().tolist()
    allr = np.array([x for x in allr if np.isfinite(x)])
    allsd = np.array([x for x in allsd if np.isfinite(x)])
    half = bool(allr.size and allr.min() < 0)
    # ONE STATION WAS STRETCHING THE AXIS. A single GPP site with a huge
    # sd(model)/sd(tower) pushed smax out and squashed everyone else into the
    # origin, in both arms. Scale to the 95th percentile and let the outliers
    # fall outside, counted and reported rather than silently accommodated.
    smax = float(np.clip(np.nanpercentile(allsd, 95) * 1.15, 1.6, 2.5))

    # 2x2 panels with the legend across the BOTTOM. Stacked in one column the
    # quarter-circles were tall and thin, wasting most of the page width; square
    # panels let each diagram use its full radius at 6.5in total width.
    ncol = 2
    nrow = -(-len(ROWS) // ncol)
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(nrow + 1, ncol,
                          height_ratios=[1] * nrow + [0.34])
    off = 0
    for i, (label, var) in enumerate(ROWS):
        r, c = divmod(i, ncol)
        ax = fig.add_subplot(gs[r, c])
        T.draw_axes(ax, smax, half, xlabel=(r == nrow - 1))
        th = np.linspace(0, 2 * np.pi, 400)
        for rms in (0.5, 1.0, 1.5):
            xc, yc = 1 + rms * np.cos(th), rms * np.sin(th)
            keep = (yc >= 0) & (np.hypot(xc, yc) <= smax)
            if not keep.any():
                continue
            k = int(np.argmax(keep & (yc > 0.15)))
            if keep[k]:
                ax.text(xc[k], yc[k], f"{rms:g}", fontsize=5.5, color="#2166ac",
                        ha="center", va="bottom")
        w = wide(d, var)
        for kind, mk, base in (("deciduous", "o", C_DEC), ("evergreen", "^", C_EVE)):
            sub = w[w.pft == kind]
            for arm, fc in (("fixed", base), ("dyn", "none")):
                r = sub[f"r_{arm}"].to_numpy(float)
                sd = sub[f"sdr_{arm}"].to_numpy(float)
                m = np.isfinite(r) & np.isfinite(sd)
                off += int((sd[m] > smax).sum())
                m &= sd <= smax
                if not m.any():
                    continue
                t = np.arccos(np.clip(r[m], -1, 1))
                ax.scatter(sd[m] * np.cos(t), sd[m] * np.sin(t), marker=mk, s=26,
                           facecolor=fc, edgecolor=base, linewidth=.8,
                           alpha=.85, zorder=4)
        ax.set_ylabel(label, fontsize=10)

    axl = fig.add_subplot(gs[nrow, 0]); _mini_legend(axl, C_DEC, C_EVE)
    axk = fig.add_subplot(gs[nrow, 1]); axk.axis("off")
    handles = [
        Line2D([], [], marker="o", ls="", color=C_DEC, label="deciduous, fixed LMA"),
        Line2D([], [], marker="o", ls="", mfc="none", mec=C_DEC, color="none",
               label="deciduous, dynamic LMA"),
        Line2D([], [], marker="^", ls="", color=C_EVE, label="evergreen, fixed LMA"),
        Line2D([], [], marker="^", ls="", mfc="none", mec=C_EVE, color="none",
               label="evergreen, dynamic LMA"),
        Line2D([], [], marker="*", ls="", color="0.15", label="AmeriFlux Tower")]
    axk.legend(handles=handles, loc="center", frameon=False, fontsize=7,
               ncol=2, handletextpad=.5, borderaxespad=0,
               columnspacing=1.0, labelspacing=.5)
    if off:
        axk.text(.5, .02, f"{off} point(s) beyond the axis", ha="center",
                 va="bottom", fontsize=6, color="0.45", transform=axk.transAxes)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stats", type=Path, default=None,
                    help="hourly_stats.csv (default $TC_RESULTS/hourly_stats.csv)")
    ap.add_argument("--basemap", type=Path, default=None)
    ap.add_argument("--linthresh", type=float, default=0.01)
    ap.add_argument("--vmax", type=float, default=1.0)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--map-size", default="7.5x10")
    ap.add_argument("--taylor-size", default="6.5x8")
    a = ap.parse_args(argv)

    try:
        stats = a.stats or resolve_out("hourly_stats.csv", create=False)
        out_dir = resolve_figure(a.out or ".")
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    if not Path(stats).is_file():
        print(f"ERROR: {stats} not found -- run submit_hourly_stats.sh first",
              file=sys.stderr)
        return 1

    d = load(Path(stats))
    sites = read_sites()
    print(f"stations: {d.station.nunique()}   variables: {sorted(d.variable.unique())}")
    print(f"matched hours per station: median {d.n.median():,.0f}  "
          f"min {d.n.min():,.0f}  max {d.n.max():,.0f}")

    bm = get_basemap(a.basemap, Path(stats).parent / "_basemap")
    print(f"basemap : {bm if bm else 'UNAVAILABLE -- maps will have no outline'}")

    dims = lambda s: tuple(float(x) for x in s.lower().split("x"))
    fig_maps(d, sites, out_dir / "hourly_skill_maps.png", bm, dims(a.map_size),
             a.linthresh, a.vmax)
    print(f"  -> {out_dir/'hourly_skill_maps.png'}")
    fig_taylor(d, out_dir / "hourly_taylor.png", dims(a.taylor_size))
    print(f"  -> {out_dir/'hourly_taylor.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
