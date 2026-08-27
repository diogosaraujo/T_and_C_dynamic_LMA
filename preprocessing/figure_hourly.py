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
  bars    4 rows x 2 cols, one per forest type.

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
BASEMAP_GUESSES = [
    "us_eco_l3/us_eco_l3.shp", "ecoregions/us_eco_l3.shp", "us_eco_l3.shp",
    "NA_CEC_Eco_Level3/NA_CEC_Eco_Level3.shp",
]


def find_basemap(explicit, roots) -> Path | None:
    """The CONUS outline, searched for rather than demanded on the command line.

    The maps shipped without a coastline for several rounds purely because
    --basemap was never passed, and nothing said so. Now it looks, and main()
    reports loudly when it comes up empty.
    """
    if explicit:
        p = Path(explicit)
        return p if p.exists() else None
    for r in roots:
        if not r:
            continue
        for g in BASEMAP_GUESSES:
            p = Path(r) / g
            if p.exists():
                return p
        hits = sorted(Path(r).glob("**/us_eco_l3*.shp"))
        if hits:
            return hits[0]
    return None


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
    outline = None
    if basemap is not None:
        import geopandas as gpd
        outline = gpd.read_file(basemap).to_crs("EPSG:5070").dissolve()

    fig, axes = plt.subplots(len(ROWS), 3, figsize=figsize,
                             constrained_layout=True,
                             gridspec_kw={"width_ratios": [2.1, 1, 1]})
    for i, (label, var) in enumerate(ROWS):
        w = wide(d, var)
        ax = axes[i, 0]
        if outline is not None:
            outline.boundary.plot(ax=ax, color="0.55", linewidth=0.45)
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
            v = w.loc[w.pft == kind, "ss"].to_numpy(float)
            v = v[np.isfinite(v)]
            if v.size:
                lim = max(0.05, np.nanpercentile(np.abs(v), 98) * 1.2)
                axh.hist(v, bins=np.linspace(-lim, lim, 19), color="0.55",
                         edgecolor="white", linewidth=.4)
                axh.axvline(0, color="0.3", lw=.8, ls=":")
                axh.axvline(np.median(v), color=C_SS, lw=1.4)
                axh.text(.03, .93, f"n={v.size}\nmed {np.median(v):+.3f}",
                         transform=axh.transAxes, fontsize=6, va="top")
            axh.tick_params(labelsize=6)
            axh.set_yticks([])
            for sp in ("top", "right", "left"):
                axh.spines[sp].set_visible(False)
            if i == 0:
                axh.set_title(kind, fontsize=9)
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
    cb.set_label("skill score   1 - RSR$_{dyn}$/RSR$_{fixed}$", fontsize=7.5)
    cb.ax.text(0.5, -2.6, "identical to 1 - RMSE$_{dyn}$/RMSE$_{fixed}$; "
               "> 0 means dynamic LMA is closer to the tower",
               transform=cb.ax.transAxes, ha="center", va="top",
               fontsize=6.5, color="0.3")
    handles = [Line2D([], [], marker="o", ls="", color="0.35", label="deciduous"),
               Line2D([], [], marker="^", ls="", color="0.35", label="evergreen")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=7, bbox_to_anchor=(0.28, -0.055))
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------- taylor
def fig_taylor(d, out_png, figsize):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import figure_taylor as T

    allr = pd.concat([wide(d, v)[["r_fixed", "r_dyn"]] for _, v in ROWS])
    half = bool(np.nanmin(allr.to_numpy(float)) < 0)
    sh = pd.concat([wide(d, v)[["sdr_fixed", "sdr_dyn"]] for _, v in ROWS])
    smax = float(min(3.0, max(1.6, np.nanpercentile(sh.to_numpy(float), 98) * 1.1)))

    fig, axes = plt.subplots(len(ROWS), 1, figsize=figsize,
                             constrained_layout=True, squeeze=False)
    for i, (label, var) in enumerate(ROWS):
        ax = axes[i, 0]
        T.draw_axes(ax, smax, half, xlabel=(i == len(ROWS) - 1))
        # Label the arcs centred on REF for what they are: contours of the
        # centred RMS difference normalised by sd(obs), i.e. RSR.
        th = np.linspace(0, 2 * np.pi, 400)
        for rms in (0.5, 1.0, 1.5):
            xc, yc = 1 + rms * np.cos(th), rms * np.sin(th)
            keep = (yc >= 0) & (np.hypot(xc, yc) <= smax)
            if not keep.any():
                continue
            k = np.argmax(keep & (yc > 0.12))
            if keep[k]:
                ax.text(xc[k], yc[k], f"RSR {rms:g}", fontsize=5.5,
                        color="#2166ac", ha="center", va="bottom")
        w = wide(d, var)
        for kind, mk in (("deciduous", "o"), ("evergreen", "^")):
            s = w[w.pft == kind]
            for arm, col, fc in (("fixed", C_FIX, C_FIX), ("dyn", C_DYN, "none")):
                r = s[f"r_{arm}"].to_numpy(float)
                sd = s[f"sdr_{arm}"].to_numpy(float)
                m = np.isfinite(r) & np.isfinite(sd)
                if not m.any():
                    continue
                t = np.arccos(np.clip(r[m], -1, 1))
                ax.scatter(sd[m] * np.cos(t), sd[m] * np.sin(t), marker=mk, s=24,
                           facecolor=fc, edgecolor=col, linewidth=.7,
                           alpha=.8, zorder=4)
        ax.set_ylabel(label, fontsize=10)

    handles = [
        Line2D([], [], marker="o", ls="", color=C_FIX, label="deciduous, fixed LMA"),
        Line2D([], [], marker="^", ls="", color=C_FIX, label="evergreen, fixed LMA"),
        Line2D([], [], marker="o", ls="", mfc="none", mec=C_DYN, color="none",
               label="deciduous, dynamic LMA"),
        Line2D([], [], marker="^", ls="", mfc="none", mec=C_DYN, color="none",
               label="evergreen, dynamic LMA"),
        Line2D([], [], marker="*", ls="", color="0.15", label="AmeriFlux Tower")]
    leg = fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
                     fontsize=7, bbox_to_anchor=(0.5, -0.030))
    fig.text(0.5, -0.062,
             "radius = sd(model)/sd(tower);  angle = correlation;  "
             "arcs centred on the tower are RSR = centred RMS difference / "
             "sd(tower)", ha="center", fontsize=6.5, color="0.3")
    fig.savefig(out_png, dpi=200, bbox_inches="tight",
                bbox_extra_artists=[leg])
    plt.close(fig)


# --------------------------------------------------------------------- bars
def fig_bars(d, out_png, figsize):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(len(ROWS), 2, figsize=figsize,
                             constrained_layout=True)
    rng = np.random.default_rng(0)
    for i, (label, var) in enumerate(ROWS):
        w = wide(d, var)
        vals = w[["rsr_fixed", "rsr_dyn", "ss"]].to_numpy(float).ravel()
        vals = vals[np.isfinite(vals)]
        hi = float(np.nanpercentile(vals, 98) * 1.15) if vals.size else 1.5
        lo = float(min(-0.25, np.nanpercentile(vals, 2) * 1.15)) if vals.size else -0.25
        for j, kind in enumerate(("evergreen", "deciduous")):
            ax = axes[i, j]
            s = w[w.pft == kind]
            ax.axhspan(lo, RSR_VERY_GOOD, color="#4daf4a", alpha=.05, zorder=0)
            ax.axhline(1.0, color="0.55", lw=.6, ls="--", zorder=1)
            ax.axhline(0.0, color="0.55", lw=.6, ls=":", zorder=1)
            for k, (col_, c, mk) in enumerate((("rsr_dyn", C_DYN, "o"),
                                               ("rsr_fixed", C_FIX, "o"),
                                               ("ss", C_SS, "D"))):
                v = s[col_].to_numpy(float); v = v[np.isfinite(v)]
                if not v.size:
                    continue
                ax.scatter(np.full(v.size, k) + rng.uniform(-.09, .09, v.size),
                           v, s=5, color=c, alpha=.28, edgecolor="none", zorder=2)
                ax.errorbar([k], [v.mean()],
                            yerr=[v.std(ddof=1) if v.size > 1 else 0], fmt=mk,
                            ms=5, color=c, capsize=3, lw=1.4, zorder=3)
            ax.set_ylim(lo, hi); ax.set_xlim(-.5, 2.5)
            ax.set_xticks([0, 1, 2])
            ax.set_xticklabels(["RSR$_{dyn}$", "RSR$_{fix}$", "SS"]
                               if i == len(ROWS) - 1 else [], fontsize=6)
            ax.tick_params(axis="y", labelsize=6, labelleft=(j == 0))
            ax.text(.03, .95, f"n={len(s)}", transform=ax.transAxes,
                    fontsize=5.5, va="top", color="0.35")
            if i == 0:
                ax.set_title(kind, fontsize=9)
            if j == 0:
                ax.set_ylabel(label, fontsize=10)
    handles = [Line2D([], [], marker="o", ls="", color=C_DYN, label="RSR dynamic LMA"),
               Line2D([], [], marker="o", ls="", color=C_FIX, label="RSR fixed LMA"),
               Line2D([], [], marker="D", ls="", color=C_SS, label="skill score"),
               Line2D([], [], ls="--", color="0.55", label="RSR = 1 (NSE = 0)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=7, bbox_to_anchor=(0.5, -0.035))
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
    ap.add_argument("--taylor-size", default="5x11")
    ap.add_argument("--bar-size", default="6.5x10")
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

    bm = find_basemap(a.basemap, [Path(stats).parent,
                                  Path(stats).parent.parent / "input_data",
                                  Path(__file__).resolve().parent.parent])
    print(f"basemap : {bm if bm else 'NOT FOUND -- maps will have no CONUS outline'}")

    dims = lambda s: tuple(float(x) for x in s.lower().split("x"))
    fig_maps(d, sites, out_dir / "hourly_skill_maps.png", bm, dims(a.map_size),
             a.linthresh, a.vmax)
    print(f"  -> {out_dir/'hourly_skill_maps.png'}")
    fig_taylor(d, out_dir / "hourly_taylor.png", dims(a.taylor_size))
    print(f"  -> {out_dir/'hourly_taylor.png'}")
    fig_bars(d, out_dir / "hourly_errorbars.png", dims(a.bar_size))
    print(f"  -> {out_dir/'hourly_errorbars.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
