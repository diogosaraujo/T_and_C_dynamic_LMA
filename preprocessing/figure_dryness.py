#!/usr/bin/env python3
"""LMA sensitivity against dryness, one figure per dataset.

Six panels, one per flux, each plotting a station's dFlux/dLMA slope against
its dryness index phi = PET/P, split deciduous and evergreen, with an OLS line
per forest type and its slope, R2 and p in the corner.

    era5 | historical | ssp126 | ssp585        4 figures

WHY phi ON THE X AXIS. It places every station on the Budyko axis, so the
question becomes whether the LMA effect depends on where a site sits between
energy limitation (phi < 1, water to spare) and water limitation (phi > 1,
demand exceeds supply). The dashed line at phi = 1 is that boundary, and the
panels are shaded either side of it.

OLS, NOT LOESS. A straight line with a reported slope, R2 and p is a testable
claim; a smoother is a picture. The template's running lines are visibly
jagged, and those steps read as structure when they are the window crossing
individual points. If the relationship turns out to be genuinely non-monotonic
the residuals will say so.

MEDIAN ACROSS GCMS AT THE METRIC LEVEL, as everywhere else in this analysis.
Both the slope and phi are computed within one GCM and then reduced across
models, so a station is one point per figure rather than five. Averaging the
underlying fluxes or PET first would mix models whose internal variability is
out of phase.

Reads station_sensitivity.csv (predictor == LMA) and station_dryness.csv.
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
FLUX_UNITS = {"GPP": "gC m$^{-2}$ yr$^{-1}$", "LAI": "m$^2$ m$^{-2}$",
              "TR": "mm yr$^{-1}$", "ET": "mm yr$^{-1}$",
              "LE": "W m$^{-2}$", "H": "W m$^{-2}$"}
C_DEC, C_EVE = "#1b7837", "#2166ac"
PFTS = [("deciduous", C_DEC, "o"), ("evergreen", C_EVE, "^")]
PANELS = [("era5", "ERA5-Land"), ("historical", "GCM historical"),
          ("ssp126", "GCM ssp126"), ("ssp585", "GCM ssp585")]

# The x axis. phi is climate and varies by dataset; rooting depth is a static
# site parameter, so its four panels share one x and differ only in the slopes.
XAXIS = {
    "phi": dict(col="phi", label=r"dryness index  $\phi$ = PET / P",
                note="        (left of the dashed line: energy-limited)",
                vline=1.0, prefix="dryness"),
    "zr95": dict(col="zr95", label="rooting depth  ZR95$_H$  (mm)",
                 note="", vline=None, prefix="rootdepth"),
    "hc": dict(col="hc", label="canopy height  hc$_H$  (m)",
               note="", vline=None, prefix="height"),
}
# MOD_PARAM name and a sanity range for each static site parameter, so a
# value that cannot be right is caught rather than plotted.
SITE_PARAM = {"zr95": ("ZR95_H", 50.0, 6000.0),
              "hc": ("hc_H", 0.5, 80.0)}

# The y axis. Every metric already exists as a column; none is recomputed from
# model output. "how" says where the per-station value comes from:
#   sens  station_sensitivity.csv, one row per station already
#   metric station_metrics.csv, one row per station already
#   years  the annual effect tables, MEDIAN OF STATION-YEARS per station --
#          deliberately different from the violin figures, which pool the
#          station-years themselves, so the two sets do not report the same
#          number and their captions must say which.
YAXIS = {
    "lma":   dict(how="sens", pred="LMA", col="slope_dyn", prefix="lma",
                  label="dFlux/dLMA, dynamic arm", per_unit=" / (g m$^{-2}$)"),
    "flux":  dict(how="years", col="rel_pct", prefix="flux",
                  label="median change in annual flux", per_unit=" (%)"),
    "sd":    dict(how="metric", col="sd_ratio", prefix="sd",
                  label="change in interannual SD", per_unit=" (%)"),
    "spei":  dict(how="sens", pred="SPEI12", col="delta_slope", prefix="spei",
                  label=r"$\Delta$ dFlux/dSPEI12, dynamic - fixed",
                  per_unit=" / SPEI unit"),
    "ta":    dict(how="sens", pred="Ta", col="delta_slope", prefix="ta",
                  label=r"$\Delta$ dFlux/dTa, dynamic - fixed",
                  per_unit=" K$^{-1}$"),
}


# Where the fit text sits, per dataset and variable. Defaults to top-left;
# these move it clear of the points, which run down-right in most panels.
CORNER = {}
for _ds, _br, _tr in (
        ("era5",       ("GPP", "LAI", "ET", "TR", "LE"), ("H",)),
        ("historical", ("GPP", "LAI"),                   ()),
        ("ssp126",     ("GPP", "LAI", "TR", "ET", "LE"), ("H",)),
        ("ssp585",     ("GPP", "LAI", "TR", "ET", "LE"), ("H",))):
    for _v in _br:
        CORNER[(_ds, _v)] = "br"
    for _v in _tr:
        CORNER[(_ds, _v)] = "tr"

# x, y, ha, va and the direction successive lines stack.
ANCHOR = {"tl": (0.02, 0.97, "left", "top", -1),
          "tr": (0.98, 0.97, "right", "top", -1),
          "br": (0.98, 0.03, "right", "bottom", +1)}


class Missing(Exception):
    """A required table or column is absent. Never substituted."""


def ols(x: np.ndarray, y: np.ndarray) -> dict:
    """Slope, intercept, R2 and a two-sided p. NaNs where n is too small."""
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 3 or np.std(x[m]) == 0:
        return {"n": n, "b": np.nan, "a": np.nan, "r2": np.nan, "p": np.nan}
    xx, yy = x[m], y[m]
    b, a = np.polyfit(xx, yy, 1)
    res = yy - (a + b * xx)
    ss_res, ss_tot = float(np.sum(res ** 2)), float(np.sum((yy - yy.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    sxx = float(np.sum((xx - xx.mean()) ** 2))
    se = float(np.sqrt(ss_res / (n - 2) / sxx)) if n > 2 and sxx > 0 else np.nan
    p = np.nan
    if np.isfinite(se) and se > 0:
        from scipy import stats as st
        p = float(2 * st.t.sf(abs(b / se), n - 2))
    return {"n": n, "b": float(b), "a": float(a), "r2": r2, "p": p}


def site_param(root: Path, xkind: str) -> pd.DataFrame:
    """Station -> a static MOD_PARAM value, from the files the runs used.

    Read back out of MOD_PARAM rather than from the product it came from.
    ZR95_H is capped at the deepest soil layer -- T&C aborts in
    Root_Fraction_General otherwise -- so the fetched D95 and the depth the
    model actually rooted with differ wherever the cap bit, and only the
    second explains model behaviour. hc_H is taken the same way for the same
    reason: what the run used is what matters, not what was looked up.
    """
    import os
    import re
    mr = os.environ.get("MODEL_RUN") or os.environ.get("TC_MODEL_RUN")
    if not mr or not Path(mr).is_dir():
        raise Missing(f"MODEL_RUN is not set, so {SITE_PARAM[xkind][0]} "
                      f"cannot be read from the MOD_PARAM files the runs "
                      f"actually used")
    name, lo, hi = SITE_PARAM[xkind]
    rows = []
    for f in sorted(Path(mr).glob("*/**/MOD_PARAM_*.m")):
        st = f.parent
        while st.parent != Path(mr) and st.parent != st:
            st = st.parent
        # The name may carry a MATLAB index: ZR95_H is written "ZR95_H = [800];"
        # but canopy height is "hc_H(1,:) = [...];". Requiring "NAME =" matched
        # the first and silently missed the second, and the run failed with
        # "no hc_H found" against files that all contain it.
        m = re.search(rf"^\s*{name}\s*(?:\([^)]*\))?\s*=\s*\[?\s*([0-9.]+)",
                      f.read_text(errors="replace"), re.M)
        if m:
            rows.append({"station": st.name, xkind: float(m.group(1))})
    if not rows:
        raise Missing(f"no {name} found in any MOD_PARAM under {mr}")
    d = (pd.DataFrame(rows).groupby("station", as_index=False)[xkind]
           .median())
    bad = int(((d[xkind] < lo) | (d[xkind] > hi)).sum())
    print(f"  {name} from MOD_PARAM: {len(d)} stations, "
          f"{d[xkind].min():.1f}-{d[xkind].max():.1f}, "
          f"median {d[xkind].median():.1f}, {d[xkind].nunique()} distinct")
    if bad:
        print(f"  WARNING: {bad} station(s) outside the plausible range "
              f"{lo}-{hi} for {name}", file=sys.stderr)
    if d[xkind].nunique() < 4:
        print(f"  WARNING: only {d[xkind].nunique()} distinct {name} values -- "
              f"a regression across so few x positions is weak whatever its p",
              file=sys.stderr)
    return d


def y_values(root: Path, yk: dict) -> pd.DataFrame:
    """Per (dataset, station, pft, variable): the metric, median over GCMs."""
    how = yk["how"]
    if how == "years":
        frames = []
        for ds in DATASETS:
            d, why = load_effect(Path(root), ds, "annual")
            if d is None:
                print(f"  {ds}: {why}"); continue
            d = d[d["variable"].isin([c for _, c in VARS])].copy()
            d["dataset"] = ds
            frames.append(d[["dataset", "gcm", "station", "year", "variable",
                             "rel_pct"]])
        if not frames:
            raise Missing("flux: no annual effect tables found")
        d = pd.concat(frames, ignore_index=True)
        # Median of station-years, then over GCMs -- one value per station.
        g = (d.groupby(["dataset", "gcm", "station", "variable"],
                       as_index=False)["rel_pct"].median()
               .groupby(["dataset", "station", "variable"],
                        as_index=False)["rel_pct"].median()
               .rename(columns={"rel_pct": "slope"}))
        return g.merge(read_sites(), on="station", how="left")

    if how == "metric":
        p = Path(root) / "station_metrics.csv"
        if not p.is_file():
            raise Missing("station_metrics.csv not found -- run station_metrics.py")
        m = pd.read_csv(p, low_memory=False)
        m = m[(m["freq"] == "annual") & (m["subset"] == "all")].copy()
        m["slope"] = 100.0 * (pd.to_numeric(m["sd_ratio"], errors="coerce") - 1)
        return (m.groupby(["dataset", "station", "pft", "variable"],
                          as_index=False)["slope"].median())

    p = Path(root) / "station_sensitivity.csv"
    if not p.is_file():
        raise Missing("station_sensitivity.csv not found -- run station_metrics.py")
    S = pd.read_csv(p, low_memory=False)
    S = S[(S["freq"] == "annual") & (S["subset"] == "all")
          & (S["predictor"] == yk["pred"])].copy()
    if S.empty:
        raise Missing(f"no annual rows for predictor {yk['pred']!r}")
    S["slope"] = pd.to_numeric(S[yk["col"]], errors="coerce")
    return (S.groupby(["dataset", "station", "pft", "variable"],
                      as_index=False)["slope"].median())


def load(root: Path, xkind: str = "phi", yk: dict | None = None) -> pd.DataFrame:
    """One row per (dataset, station, pft, variable): slope and phi."""
    sp = Path(root) / "station_sensitivity.csv"
    dp = Path(root) / "station_dryness.csv"
    need = []
    if (yk or YAXIS["lma"])["how"] == "sens":
        need.append((sp, "run station_metrics.py"))
    if xkind == "phi":
        need.append((dp, "run build_dryness.py"))
    for p, why in need:
        if not p.is_file():
            raise Missing(f"{p.name} not found -- {why} first")
    s = y_values(root, yk or YAXIS["lma"])

    if xkind in SITE_PARAM:
        d = site_param(root, xkind)
        out = s.merge(d, on="station", how="inner")   # static: no dataset key
    else:
        D = pd.read_csv(dp, low_memory=False)
        if "phi" not in D.columns:
            raise Missing("station_dryness.csv has no phi column")
        D["phi"] = pd.to_numeric(D["phi"], errors="coerce")
        d = (D.groupby(["dataset", "station"], as_index=False)["phi"].median())
        out = s.merge(d, on=["dataset", "station"], how="inner")
    if out.empty:
        raise Missing(f"no station matched between the sensitivity table and "
                      f"the {xkind} values -- check the station IDs")
    lost = s["station"].nunique() - out["station"].nunique()
    if lost:
        print(f"  note: {lost} station(s) have a slope but no {xkind}")
    return out


def build(d, ds, label, out_png, figsize, xk, yk):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    sub = d[d["dataset"] == ds]
    if sub.empty:
        raise Missing(f"{ds}: no rows")
    fig, axes = plt.subplots(3, 2, figsize=figsize, constrained_layout=True,
                             sharex=True)
    xc = xk["col"]
    lo = float(np.nanmin(sub[xc])) * 0.95
    hi = float(np.nanmax(sub[xc])) * 1.05
    total = 0
    for ax, (disp, col) in zip(axes.ravel(), VARS):
        # Energy- and water-limited sides of the Budyko boundary.
        if xk["vline"] is not None:
            ax.axvspan(lo, xk["vline"], color="#eef2f7", zorder=0)
            ax.axvspan(xk["vline"], hi, color="#fbf0ea", zorder=0)
            ax.axvline(xk["vline"], color="0.35", lw=0.9, ls="--", zorder=2)
        ax.axhline(0.0, color="0.15", lw=0.8, zorder=2)
        txt = []
        for pft, colr, mk in PFTS:
            w = sub[(sub["variable"] == col) & (sub["pft"] == pft)]
            x = w[xc].to_numpy(float)
            y = w["slope"].to_numpy(float)
            m = np.isfinite(x) & np.isfinite(y)
            if not m.any():
                continue
            total += int(m.sum())
            ax.scatter(x[m], y[m], marker=mk, s=17, facecolor=colr,
                       edgecolor="white", linewidth=0.3, alpha=0.85, zorder=4)
            r = ols(x, y)
            if np.isfinite(r["b"]):
                xs = np.linspace(np.nanmin(x[m]), np.nanmax(x[m]), 50)
                ax.plot(xs, r["a"] + r["b"] * xs, color=colr, lw=1.6, zorder=5)
                star = ("***" if r["p"] < 0.001 else "**" if r["p"] < 0.01
                        else "*" if r["p"] < 0.05 else "")
                txt.append((colr, f"b={r['b']:.3g}  R$^2$={r['r2']:.2f}"
                                  f"{star}  n={r['n']}"))
        x0, y0, ha, va, step = ANCHOR[CORNER.get((ds, disp), "tl")]
        for i, (colr, t) in enumerate(txt):
            # Boxed: unboxed it sat on the zero line and the points, and the
            # numbers are the reportable part of the panel. Bottom-anchored
            # text stacks upward so the two PFT lines never overlap.
            j = (len(txt) - 1 - i) if step > 0 else i
            ax.text(x0, y0 + step * 0.085 * j, t, transform=ax.transAxes,
                    fontsize=5.6, color=colr, va=va, ha=ha, zorder=6,
                    bbox=dict(facecolor="white", alpha=0.78, edgecolor="none",
                              boxstyle="round,pad=0.18"))
        ax.set_title(disp, fontsize=9, loc="left", pad=2)
        ax.set_ylabel(("" if yk["per_unit"].startswith(" (")
                       else FLUX_UNITS.get(disp, "")) + yk["per_unit"],
                      fontsize=6)
        ax.tick_params(labelsize=7)
        ax.set_xlim(lo, hi)
        ax.grid(color="0.92", lw=0.5)
        ax.set_axisbelow(True)
        for sp_ in ("top", "right"):
            ax.spines[sp_].set_visible(False)
    if total == 0:
        raise Missing(f"{ds}: every panel empty")
    handles = [Line2D([], [], marker=mk, ls="", color=c, label=p)
               for p, c, mk in PFTS]
    fig.suptitle(f"{yk['label']}   —   {label}", fontsize=10, x=0.01,
                 ha="left")
    fig.supxlabel(xk["label"] + xk["note"], fontsize=8)
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, 0.0))
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {out_png}   ({total} points)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--datasets", default=",".join(k for k, _ in PANELS))
    ap.add_argument("--size", default="6.5x8")
    ap.add_argument("--y", default="lma", choices=list(YAXIS),
                    help="metric on the y axis")
    ap.add_argument("--x", default="phi", choices=list(XAXIS),
                    help="x axis: dryness index or rooting depth")
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

    try:
        d = load(root, a.x, YAXIS[a.y])
    except Missing as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    print(f"stations with slope and {a.x}: {d['station'].nunique()}  "
          f"rows: {len(d)}")

    xk, yk = XAXIS[a.x], YAXIS[a.y]
    failures = []
    want = [x.strip() for x in a.datasets.split(",")]
    for ds, label in PANELS:
        if ds not in want:
            continue
        try:
            build(d, ds, label, out_dir /
                  f"{yk['prefix']}_vs_{xk['prefix']}_{ds}.png",
                  dims(a.size), xk, yk)
        except Missing as e:
            failures.append(str(e)); print(f"SKIP {ds}: {e}")
    if failures:
        print(f"\n{len(failures)} figure(s) not produced:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
