#!/usr/bin/env python3
"""Flux change and interannual variability, ERA5 and the two SSPs.

Four panels stacked vertically -- ERA5-Land, GCM historical, ssp126, ssp585 --
each showing six variables as a violin with the median and interquartile range
drawn over it, deciduous above evergreen, sharing one x axis so the four are
comparable by position.

Absolute fluxes are a SEPARATE figure (--metrics absolute): six vertical
panels, one per variable, sixteen boxes each. Attached to every variable row
they came to 24 axes on one page and were unreadable.

MEDIAN ACROSS GCMS IS TAKEN AT THE METRIC LEVEL, NEVER AT THE FLUX LEVEL.
Each metric is computed inside one GCM first, and the median is taken across
the five resulting values. GCMs are not synchronised in time: their internal
variability is independent, so one model's drought year falls in a different
calendar year from another's. Taking a median of raw fluxes by calendar year
would average anomalies that are out of phase and would systematically damp
interannual variance -- which would corrupt the variability metric outright.
The metrics here are paired within-model differences: both arms see the same
model's weather, so that model's climate phase is common to both and cancels.
A median across five such estimates is a median of five estimates of the same
quantity, not an average across mismatched climates.

  flux         100*(dyn-fixed)/|fixed| per station-year, then median over GCMs
  variability  100*(sd_dyn/sd_fixed - 1) per station,    then median over GCMs

SYMLOG X, SCALED TO THE BOXES. Linear inside --linthresh, logarithmic
outside. The limits come from the interquartile ranges of the groups drawn,
not from the extremes, so violin tails extend past the view and are clipped
on purpose. Scaling to the extremes instead let a single station-year with a
near-zero denominator set the range and compress every median onto one line.
Every value is still computed and still shapes the box and the median; only
the view is bounded.

No station is excluded. Station-years whose denominator survives as
near-zero after the median are COUNTED AND REPORTED, so it is visible whether
the median removed them or not, rather than assumed either way.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from results_dir import NoResultsDir, resolve_figure, resolve_out  # noqa: E402
from station_metrics import DATASETS, load as load_effect, read_sites  # noqa: E402

VARS = [("GPP", "GPP"), ("LAI", "LAI_H"), ("TR", "T"),
        ("ET", "ET"), ("LE", "QE"), ("H", "H")]
UNITS = {"GPP": "gC m$^{-2}$", "LAI": "m$^2$ m$^{-2}$", "TR": "mm",
         "ET": "mm", "LE": "W m$^{-2}$", "H": "W m$^{-2}$"}
C_DEC, C_EVE = "#1b7837", "#2166ac"
PFTS = [("deciduous", C_DEC), ("evergreen", C_EVE)]
PANELS = [("era5", "ERA5-Land"), ("historical", "GCM historical"),
          ("ssp126", "GCM ssp126"), ("ssp585", "GCM ssp585")]
METRICS = {
    "flux": "dynamic - fixed LMA   (% of fixed-LMA annual flux)",
    "variability": "change in interannual SD   (%, dynamic vs fixed)",
    "absolute": "mean annual flux, native units",
}


def pct_ticks(ax, linthresh, lo, hi):
    """Percent tick labels, and limits taken from the IQR spread.

    The range comes from the boxes, not from the extremes: the whole point of
    the figure is the middle of the distribution, and a single station-year
    with a near-zero denominator would otherwise set the scale and compress
    every median into a line. Violin tails are clipped by this on purpose --
    they are still drawn, just outside the view.

    Labels carry a % sign. The default symlog formatter writes
    mantissa-and-exponent, which reads as a measurement rather than a percent
    change; a power of ten is written as a power rather than as 1000000%.
    """
    import numpy as _np
    from matplotlib.ticker import FuncFormatter, FixedLocator, MaxNLocator

    def fmt(v, _pos):
        if v == 0:
            return "0"
        a, sign = abs(v), ("-" if v < 0 else "")
        if a < 1000:
            return f"{sign}{a:g}%"
        return sign + f"$10^{{{int(round(_np.log10(a)))}}}$%"

    pad = 0.12 * (hi - lo) if hi > lo else max(abs(hi), 1.0) * 0.5
    lo_, hi_ = lo - pad, hi + pad
    if max(abs(lo_), abs(hi_)) <= linthresh:
        # Inside the linear region of the symlog axis, so ordinary ticks read
        # correctly and decade ticks would give one or two across the panel.
        ax.xaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5, 10]))
    else:
        e = int(_np.ceil(_np.log10(max(abs(lo_), abs(hi_)))))
        dec = [10.0 ** k for k in range(int(_np.log10(linthresh)), e + 1)]
        ticks = [-t for t in reversed(dec)] + [0.0] + dec
        ax.xaxis.set_major_locator(
            FixedLocator([t for t in ticks if lo_ <= t <= hi_]))
    ax.xaxis.set_major_formatter(FuncFormatter(fmt))
    ax.set_xlim(lo_, hi_)


class Missing(Exception):
    """A required table or column is absent. Never substituted."""


def flux_values(root: Path, near_zero: float) -> tuple[pd.DataFrame, list[str]]:
    """rel_pct per station-year, median across GCMs. Reports near-zero cases."""
    frames, notes = [], []
    for ds in DATASETS:
        d, why = load_effect(Path(root), ds, "annual")
        if d is None:
            notes.append(f"{ds}: {why}"); continue
        d = d[d["variable"].isin([c for _, c in VARS])].copy()
        d["dataset"] = ds
        # Flag the denominator BEFORE the median, so it can be said whether
        # the median removed the pathology or carried it through.
        f = pd.to_numeric(d["fixed"], errors="coerce")
        site = f.abs().groupby([d["gcm"], d["station"], d["variable"]]) \
                .transform("mean")
        d["_tiny"] = f.abs() < near_zero * site
        frames.append(d[["dataset", "gcm", "station", "year", "variable",
                         "rel_pct", "_tiny"]])
    if not frames:
        raise Missing("flux: no annual effect tables found -- " + "; ".join(notes))
    d = pd.concat(frames, ignore_index=True)
    key = ["dataset", "station", "year", "variable"]
    before = int(d["_tiny"].sum())
    g = d.groupby(key, as_index=False).agg(
        value=("rel_pct", "median"), tiny_frac=("_tiny", "mean"),
        n_gcm=("gcm", "nunique"))
    # A station-year whose denominator was tiny in MOST models still has a
    # pathological median; one tiny model out of five does not.
    after = int((g["tiny_frac"] > 0.5).sum())
    notes.append(f"near-zero denominators: {before} station-year-models "
                 f"before the median, {after} station-years still majority-tiny "
                 f"after it ({100*after/max(len(g),1):.3f}% of {len(g)})")
    if after:
        notes.append("  -> the median did NOT remove them; a denominator guard "
                     "is still required if these dominate the axis")
    else:
        notes.append("  -> the median removed them; no guard needed")
    g = g.merge(read_sites(), on="station", how="left")
    # An unmatched site list gives every row pft=NaN, which filters to nothing
    # and draws four empty panels with no error. Name it instead.
    miss = int(g["pft"].isna().sum())
    if miss == len(g):
        raise Missing("flux: no station matched the site lists, so every row "
                      "has no forest type; check that the effect tables and "
                      "deciduous/evergreen_ameriflux.csv use the same IDs")
    if miss:
        notes.append(f"  {miss} of {len(g)} rows have no forest type "
                     f"({g.loc[g['pft'].isna(), 'station'].nunique()} stations "
                     f"absent from the site lists)")
    return g, notes


def var_values(M: pd.DataFrame) -> pd.DataFrame:
    """100*(sd_ratio-1) per station, median across GCMs."""
    d = M[(M["freq"] == "annual") & (M["subset"] == "all")].copy()
    if d.empty:
        raise Missing("variability: no rows with freq=annual and subset=all")
    d["v"] = 100.0 * (pd.to_numeric(d["sd_ratio"], errors="coerce") - 1.0)
    return (d.groupby(["dataset", "station", "pft", "variable"],
                      as_index=False)["v"].median()
             .rename(columns={"v": "value"}))


def abs_values(M: pd.DataFrame) -> pd.DataFrame:
    """Mean annual flux per station and arm, median across GCMs."""
    d = M[(M["freq"] == "annual") & (M["subset"] == "all")].copy()
    out = []
    for arm, col in (("fixed", "mean_fixed"), ("dyn", "mean_dyn")):
        if col not in d.columns:
            raise Missing(f"absolute boxes: column {col!r} absent")
        g = (d.assign(v=pd.to_numeric(d[col], errors="coerce"))
               .groupby(["dataset", "station", "pft", "variable"],
                        as_index=False)["v"].median())
        g["arm"] = arm
        out.append(g)
    return pd.concat(out, ignore_index=True).rename(columns={"v": "value"})


def draw_rows(ax, d, linthresh):
    """Violin + median/IQR per variable row, deciduous above evergreen."""
    n = 0
    for i, (disp, col) in enumerate(VARS):
        for k, (pft, colr) in enumerate(PFTS):
            off = 0.20 if k == 0 else -0.20
            v = d.loc[(d["variable"] == col) & (d["pft"] == pft),
                      "value"].to_numpy(float)
            v = v[np.isfinite(v)]
            if v.size == 0:
                continue
            n += v.size
            y = -i + off
            if v.size > 2 and np.ptp(v) > 0:
                # Violin in SYMLOG space: a kernel density on raw values would
                # put essentially all its mass in one bin once a station-year
                # reaches 1e8, drawing a flat line. Estimate on the transformed
                # axis and map back, so the shape reflects what is displayed.
                t = np.arcsinh(v / linthresh)
                try:
                    pc = ax.violinplot([t], positions=[y], widths=0.34,
                                       orientation="horizontal",
                                       showextrema=False, showmedians=False)
                except TypeError:
                    pc = ax.violinplot([t], positions=[y], widths=0.34,
                                       vert=False, showextrema=False,
                                       showmedians=False)
                for b in pc["bodies"]:
                    p = b.get_paths()[0]
                    p.vertices[:, 0] = linthresh * np.sinh(p.vertices[:, 0])
                    b.set_facecolor(colr); b.set_alpha(0.28)
                    b.set_edgecolor(colr); b.set_linewidth(0.5)
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.plot([q1, q3], [y, y], color=colr, lw=2.2,
                    solid_capstyle="butt", zorder=4)
            ax.plot([med], [y], "o", mfc="white", mec=colr, mew=1.3, ms=4.6,
                    zorder=5)
    ax.axvline(0, color="0.15", lw=0.9, zorder=3)
    ax.set_xscale("symlog", linthresh=linthresh)
    ax.set_yticks([-i for i in range(len(VARS))])
    ax.set_yticklabels([d_ for d_, _ in VARS], fontsize=8)
    ax.set_ylim(-len(VARS) + 0.45, 0.55)
    ax.grid(axis="x", color="0.90", lw=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return n


def build_absolute(A, out_png, figsize):
    """Absolute fluxes: six vertical panels, sixteen boxes each.

    One panel per variable so each keeps its own units on y. Within a panel,
    four dataset groups of four boxes -- deciduous fixed, deciduous dynamic,
    evergreen fixed, evergreen dynamic -- so the treatment pair sits side by
    side and the scenario progression reads left to right.

    This is what makes the percent figures interpretable: a 5% change is a
    different quantity for GPP than for LAI, and only the absolute magnitude
    says which. Expect the fixed and dynamic boxes to overlap almost
    completely -- the station-to-station spread is far larger than the
    treatment effect, which is the honest picture and the reason the paired
    percent view exists at all.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    fig, axes = plt.subplots(len(VARS), 1, figsize=figsize,
                             constrained_layout=True)
    total = 0
    for ax, (disp, col) in zip(np.atleast_1d(axes), VARS):
        x, ticks, labels = 0.0, [], []
        for ds, lab in PANELS:
            first = x
            for pft, colr in PFTS:
                for arm in ("fixed", "dyn"):
                    v = A.loc[(A["dataset"] == ds) & (A["variable"] == col)
                              & (A["pft"] == pft) & (A["arm"] == arm),
                              "value"].to_numpy(float)
                    v = v[np.isfinite(v)]
                    if v.size:
                        total += v.size
                        bp = ax.boxplot([v], positions=[x], widths=0.68,
                                        showfliers=False, patch_artist=True,
                                        manage_ticks=False)
                        for b in bp["boxes"]:
                            b.set_facecolor(colr if arm == "fixed" else "white")
                            b.set_alpha(0.5 if arm == "fixed" else 1.0)
                            b.set_edgecolor(colr); b.set_linewidth(0.8)
                        for part in ("whiskers", "caps", "medians"):
                            for it in bp[part]:
                                it.set_color(colr); it.set_linewidth(0.8)
                    x += 1.0
            ticks.append((first + x - 1) / 2); labels.append(lab)
            x += 0.9                      # gap between dataset groups
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_ylabel(disp + "   " + UNITS.get(disp, ""), fontsize=7.5)
        ax.tick_params(axis="y", labelsize=7)
        ax.grid(axis="y", color="0.91", lw=0.5)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    if total == 0:
        raise Missing("absolute: no finite values in any panel")
    handles = [Patch(facecolor=c, alpha=.5, edgecolor=c, label=p)
               for p, c in PFTS]
    handles += [Patch(facecolor="0.55", alpha=.5, edgecolor="0.25",
                      label="fixed LMA (filled)"),
                Patch(facecolor="white", edgecolor="0.25",
                      label="dynamic LMA (open)")]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
               fontsize=7.5, bbox_to_anchor=(0.5, 0.0))
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {out_png}   ({total} station values)")


def build(d, metric, out_png, figsize, linthresh):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    have = [(k, lab) for k, lab in PANELS if (d["dataset"] == k).any()]
    if not have:
        raise Missing(f"{metric}: none of {[k for k,_ in PANELS]} present; "
                      f"table has {sorted(d['dataset'].unique())}")
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    outer = fig.add_gridspec(len(have), 1, hspace=0.06)
    total = 0
    # ONE shared x axis. Per-panel limits made the four impossible to compare:
    # the eye reads position, and the same position meant a different percent
    # in each panel.
    # Limits from the QUARTILES of every group drawn, not from the extremes.
    q = (d.groupby(["dataset", "variable", "pft"], observed=True)["value"]
           .quantile([0.25, 0.75]))
    q = q[np.isfinite(q)]
    lo = float(q.min()) if len(q) else -1.0
    hi = float(q.max()) if len(q) else 1.0
    ax0 = None
    for i, (ds, lab) in enumerate(have):
        ax = fig.add_subplot(outer[i], sharex=ax0)
        ax0 = ax0 or ax
        total += draw_rows(ax, d[d["dataset"] == ds], linthresh)
        pct_ticks(ax, linthresh, lo, hi)
        if i < len(have) - 1:
            ax.tick_params(labelbottom=False)
        ax.set_title(lab, fontsize=9, loc="left", pad=3)
    if total == 0:
        raise Missing(f"{metric}: every panel empty after filtering")
    handles = [Patch(facecolor=c, alpha=.45, edgecolor=c, label=p)
               for p, c in PFTS]
    fig.supxlabel(METRICS[metric], fontsize=8.5)
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               fontsize=7, bbox_to_anchor=(0.5, 0.0))
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {out_png}   ({total} values)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--metrics", default="flux,variability,absolute")
    ap.add_argument("--size", default="6.5x9")
    ap.add_argument("--linthresh", type=float, default=100.0,
                    help="symlog: linear inside +/- this many percent")
    ap.add_argument("--near-zero", type=float, default=0.01,
                    help="denominator below this fraction of the site mean "
                         "counts as near-zero (reported, not filtered)")
    a = ap.parse_args(argv)

    def dims(s):
        w, h = s.lower().split("x")
        return float(w), float(h)

    try:
        root = Path(a.results or resolve_out(".", create=False))
        out_dir = Path(a.out or resolve_figure("."))
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    out_dir.mkdir(parents=True, exist_ok=True)

    p = root / "station_metrics.csv"
    if not p.is_file():
        print(f"ERROR: {p} not found -- run station_metrics.py first",
              file=sys.stderr)
        return 1
    M = pd.read_csv(p, low_memory=False)

    failures = []
    for metric in [m.strip() for m in a.metrics.split(",") if m.strip()]:
        if metric not in METRICS:
            print(f"ERROR: unknown metric {metric!r}", file=sys.stderr)
            return 1
        try:
            if metric == "absolute":
                build_absolute(abs_values(M), out_dir / "effect_absolute.png",
                               dims(a.size))
                continue
            if metric == "flux":
                d, notes = flux_values(root, a.near_zero)
                for n in notes:
                    print(n)
            else:
                d = var_values(M)
            build(d, metric, out_dir / f"effect_{metric}.png",
                  dims(a.size), a.linthresh)
        except Missing as e:
            failures.append(str(e)); print(f"SKIP {metric}: {e}")
    if failures:
        print(f"\n{len(failures)} figure(s) not produced:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
