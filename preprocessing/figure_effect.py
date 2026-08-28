#!/usr/bin/env python3
"""Set 2: horizontal violin + error bar figures of the LMA treatment effect.

MODEL vs MODEL. Nothing here touches tower data -- every x axis is a
dynamic-arm quantity expressed against the fixed arm, so the comparison is
between the two model runs and the tower never enters.

Five metrics, each its own figure:

  variability   100*(sd_dyn/sd_fixed - 1)          change in interannual SD
  flux          rel_pct, one point per STATION-YEAR % change in the annual flux
  sens_Ta       100*(b_dyn - b_fixed)/|b_fixed|    % change in dF/dTa
  sens_SPEI12   100*(b_dyn - b_fixed)/|b_fixed|    % change in dF/dSPEI12
  sens_LMA      slope_dyn_std                      standardised dF/dLMA

WHY sens_LMA IS NOT A PERCENT CHANGE. The fixed arm holds LMA constant, so
there is no LMA variation to regress against and station_metrics sets its
slope to NaN on purpose. A percent change against that is 0/0. The dynamic
arm's standardised slope is the honest quantity: it says how strongly the
flux tracks LMA when LMA is allowed to move, and the fixed arm's answer is
"not at all" by construction rather than by measurement.

POOLING DIFFERS BY METRIC, because the metrics are not the same shape.

  flux          one point per STATION-YEAR, straight from the annual effect
                tables. Every station is driven by the same forcing record,
                so they carry the same years and no site can outvote another
                by being longer. The script checks that and says so if the
                year counts actually differ.
  variability   one point per STATION. An interannual SD needs the whole
                series; there is no per-year value to pool.
  sensitivity   one point per STATION. One regression over all years yields
                one slope, for the same reason.

--flux-pool station falls back to the per-site mean if the year counts turn
out uneven after all. Note that the site mean is not robust: a year whose
fixed-arm flux is near zero gives a huge ratio that averaging does not damp,
and one station's summary went from 1.2% to 129.54% on exactly that earlier
in this project. Pooling station-years shows those years as individual
outliers instead of letting one of them move a site's whole point.

Datasets:
  era5        one panel
  gcm         3 panels (historical / ssp126 / ssp585), all GCMs pooled
  gcmmedian   3 panels, median across GCMs per station first
  gcmrows     GCM x scenario grid, one row per GCM

sens_LMA ignores those four and produces ONE figure:
ERA5-Land | GCM historical | GCM ssp126 | GCM ssp585.
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

# Display name -> name in the metrics tables. LE is stored as QE and LAI as
# LAI_H; using the display names against the CSV silently matches nothing.
VARS = [("GPP", "GPP"), ("LAI", "LAI_H"), ("T", "T"),
        ("ET", "ET"), ("LE", "QE"), ("H", "H")]

C_DEC, C_EVE = "#1b7837", "#2166ac"
PFTS = [("deciduous", C_DEC), ("evergreen", C_EVE)]
SCEN = ["historical", "ssp126", "ssp585"]

METRICS = {
    "variability": ("change in interannual SD  (%, dynamic vs fixed)",
                    "Change in interannual variability"),
    "flux": ("dynamic - fixed LMA  (% of fixed-LMA mean)",
             "Change in annual flux"),
    "sens_Ta": ("change in dFlux/dTa  (% of fixed-LMA slope)",
                "Change in temperature sensitivity"),
    "sens_SPEI12": ("change in dFlux/dSPEI12  (% of fixed-LMA slope)",
                    "Change in drought sensitivity"),
    "sens_LMA": ("standardised dFlux/dLMA, dynamic arm  (SD per SD)",
                 "Sensitivity to LMA"),
}


class Missing(Exception):
    """A required table or column is absent. Never substituted."""


def _pct_change(new: np.ndarray, old: np.ndarray) -> np.ndarray:
    """100*(new-old)/|old|, NaN where old is 0 or non-finite.

    A near-zero fixed-arm slope makes this explode, so the guard is on the
    denominator itself rather than on the result: a station with no fixed-arm
    sensitivity has no defined percent change in sensitivity.
    """
    out = np.full(new.shape, np.nan)
    ok = np.isfinite(new) & np.isfinite(old) & (np.abs(old) > 0)
    out[ok] = 100.0 * (new[ok] - old[ok]) / np.abs(old[ok])
    return out


def flux_station_years(root: Path) -> pd.DataFrame:
    """Per (dataset, gcm, scenario, station, pft, variable) x YEAR rel_pct."""
    frames, missing = [], []
    for ds in DATASETS:
        d, why = load_effect(Path(root), ds, "annual")
        if d is None:
            missing.append(why); continue
        d["dataset"] = ds
        frames.append(d[["dataset", "gcm", "scenario", "station", "year",
                         "variable", "rel_pct"]])
    if not frames:
        raise Missing("flux: no annual effect tables found -- " + "; ".join(missing))
    if missing:
        print("  note: " + "; ".join(missing))
    d = pd.concat(frames, ignore_index=True)
    d = d.merge(read_sites(), on="station", how="left")
    if "pft" not in d.columns:
        raise Missing("flux: the site table has no pft column")
    d = d.rename(columns={"rel_pct": "value"})

    # The claim that every station shares the forcing record is checkable, so
    # check it rather than assume it. Uneven year counts would mean a long
    # station carries more weight in the pooled violin than a short one.
    yr = d.groupby(["dataset", "station"])["year"].nunique()
    for ds, g in yr.groupby(level=0):
        if g.nunique() > 1:
            lo, hi = int(g.min()), int(g.max())
            print(f"  WARNING {ds}: station year counts are uneven "
                  f"({lo}-{hi}); pooled station-years weight the long "
                  f"stations more. --flux-pool station avoids that.")
    return d[["dataset", "gcm", "scenario", "station", "pft", "variable", "value"]]


def series(M: pd.DataFrame, S: pd.DataFrame, metric: str,
           flux_pool: str = "stationyear", root: Path | None = None) -> pd.DataFrame:
    """One row per (dataset, gcm, scenario, station, pft, variable, value)."""
    keep = ["dataset", "gcm", "scenario", "station", "pft", "variable"]
    if metric in ("variability", "flux"):
        if M is None:
            raise Missing("station_metrics.csv is required for " + metric)
        d = M[(M["freq"] == "annual") & (M["subset"] == "all")].copy()
        if d.empty:
            raise Missing(f"{metric}: no rows with freq=annual and subset=all")
        if metric == "variability":
            d["value"] = 100.0 * (pd.to_numeric(d["sd_ratio"],
                                                errors="coerce") - 1.0)
        else:
            if flux_pool == "stationyear":
                return flux_station_years(root)
            if "mean_rel_pct" not in d.columns:
                raise Missing("flux: mean_rel_pct absent from station_metrics.csv")
            d["value"] = pd.to_numeric(d["mean_rel_pct"], errors="coerce")
    else:
        if S is None:
            raise Missing("station_sensitivity.csv is required for " + metric)
        pred = metric.split("_", 1)[1]
        d = S[(S["freq"] == "annual") & (S["subset"] == "all")
              & (S["predictor"] == pred)].copy()
        if d.empty:
            have = sorted(S.get("predictor", pd.Series(dtype=str)).unique())
            raise Missing(f"{metric}: no rows for predictor={pred!r}; "
                          f"table has {have}")
        if pred == "LMA":
            d["value"] = pd.to_numeric(d["slope_dyn_std"], errors="coerce")
        else:
            d["value"] = _pct_change(
                pd.to_numeric(d["slope_dyn"], errors="coerce").to_numpy(float),
                pd.to_numeric(d["slope_fixed"], errors="coerce").to_numpy(float))
    for c in keep:
        if c not in d.columns:
            raise Missing(f"{metric}: column {c!r} absent from the metrics table")
    d["gcm"] = d["gcm"].fillna("").astype(str)
    d["scenario"] = d["scenario"].fillna("historical").astype(str)
    return d[keep + ["value"]]


def var_shares(d: pd.DataFrame) -> str | None:
    """Share of variance sitting between GCMs vs between stations.

    Type-I style: remove the GCM means, then the station means from what is
    left. With unequal cells these do not partition the variance exactly, so
    the numbers are read as 'roughly how much of the spread is which', which
    is what the panel note claims and no more.
    """
    v = d["value"].to_numpy(float)
    m = np.isfinite(v)
    if m.sum() < 8 or d["gcm"].nunique() < 2:
        return None
    x = d.loc[m, ["gcm", "station"]].copy()
    x["v"] = v[m]
    tot = float(np.var(x["v"], ddof=1))
    if not (tot > 0):
        return None
    g = x.groupby("gcm")["v"].transform("mean")
    res_g = x["v"] - g
    s = res_g.groupby(x["station"]).transform("mean")
    res_gs = res_g - s
    f_g = max(0.0, 1 - float(np.var(res_g, ddof=1)) / tot)
    f_s = max(0.0, (float(np.var(res_g, ddof=1))
                    - float(np.var(res_gs, ddof=1))) / tot)
    return (f"variance: GCM {100*f_g:.0f}%  site {100*f_s:.0f}%  "
            f"resid {100*max(0.0, 1-f_g-f_s):.0f}%")


def panel(ax, d: pd.DataFrame, note: str | None = None) -> int:
    """Horizontal violins, one variable per row, deciduous above evergreen."""
    import matplotlib.pyplot as plt      # noqa: F401  (backend already set)

    n_used = 0
    yticks, ylabels = [], []
    for i, (disp, col) in enumerate(VARS):
        y0 = -i
        yticks.append(y0); ylabels.append(disp)
        for k, (pft, colr) in enumerate(PFTS):
            off = 0.20 if k == 0 else -0.20
            v = d.loc[(d["variable"] == col) & (d["pft"] == pft),
                      "value"].to_numpy(float)
            v = v[np.isfinite(v)]
            if v.size == 0:
                continue
            n_used += v.size
            if v.size > 1 and np.ptp(v) > 0:
                kw = dict(positions=[y0 + off], widths=0.34,
                          showextrema=False, showmedians=False)
                try:        # vert= is deprecated in mpl 3.11, gone in 3.13
                    pc = ax.violinplot([v], orientation="horizontal", **kw)
                except TypeError:                       # mpl < 3.10
                    pc = ax.violinplot([v], vert=False, **kw)
                for b in pc["bodies"]:
                    b.set_facecolor(colr); b.set_alpha(0.30)
                    b.set_edgecolor(colr); b.set_linewidth(0.6)
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.plot([q1, q3], [y0 + off] * 2, color=colr, lw=2.4,
                    solid_capstyle="butt", zorder=3)
            ax.plot([med], [y0 + off], "o", mfc="white", mec=colr,
                    mew=1.4, ms=5, zorder=4)
    ax.axvline(0, color="0.15", lw=1.0, zorder=2)
    ax.set_yticks(yticks); ax.set_yticklabels(ylabels)
    ax.set_ylim(-len(VARS) + 0.4, 0.6)
    ax.grid(axis="x", color="0.88", lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if note:
        ax.text(0.99, 0.02, note, transform=ax.transAxes, ha="right",
                va="bottom", fontsize=6, color="0.40")
    return n_used


def _finish(fig, axes, xlabel, title, out_png, dpi=200):
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="s", ls="", color=c, alpha=.55,
                      mec=c, label=p) for p, c in PFTS]
    fig.suptitle(title, fontsize=10, x=0.01, ha="left")
    fig.supxlabel(xlabel, fontsize=9)
    # BELOW the x label, not on top of it. constrained_layout reserves room for
    # supxlabel but not for a figure legend, so "lower center" put the two in
    # the same place.
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, 0.0))
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    print(f"  -> {out_png}")


def fig_single(d, metric, out_png, figsize):
    import matplotlib.pyplot as plt
    xlabel, title = METRICS[metric]
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    n = panel(ax, d)
    if n == 0:
        raise Missing(f"{metric}: every variable empty after filtering")
    _finish(fig, ax, xlabel, title, out_png)
    plt.close(fig)


def fig_scen(d, metric, out_png, figsize, decompose=False):
    """Three panels, one per scenario, sharing an x axis."""
    import matplotlib.pyplot as plt
    xlabel, title = METRICS[metric]
    have = [s for s in SCEN if (d["scenario"] == s).any()]
    if not have:
        raise Missing(f"{metric}: none of {SCEN} present; "
                      f"table has {sorted(d['scenario'].unique())}")
    fig, axes = plt.subplots(1, len(have), figsize=figsize, sharex=True,
                             sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    for ax, sc in zip(axes, have):
        sub = d[d["scenario"] == sc]
        panel(ax, sub, var_shares(sub) if decompose else None)
        ax.set_title(sc, fontsize=9)
    _finish(fig, axes, xlabel, title, out_png)
    plt.close(fig)


def fig_grid(d, metric, out_png, figsize):
    """One row per GCM, one column per scenario."""
    import matplotlib.pyplot as plt
    xlabel, title = METRICS[metric]
    gcms = sorted(g for g in d["gcm"].unique() if g)
    have = [s for s in SCEN if (d["scenario"] == s).any()]
    if not gcms:
        raise Missing(f"{metric}: no named GCMs in the table")
    fig, axes = plt.subplots(len(gcms), len(have), figsize=figsize,
                             sharex=True, sharey=True,
                             constrained_layout=True, squeeze=False)
    for r, g in enumerate(gcms):
        for c, sc in enumerate(have):
            ax = axes[r][c]
            panel(ax, d[(d["gcm"] == g) & (d["scenario"] == sc)])
            if r == 0:
                ax.set_title(sc, fontsize=9)
            if c == 0:
                ax.annotate(g, xy=(-0.42, 0.5), xycoords="axes fraction",
                            rotation=90, ha="center", va="center", fontsize=8)
    _finish(fig, axes, xlabel, title, out_png)
    plt.close(fig)


def fig_lma(d, out_png, figsize):
    """sens_LMA on its own: ERA5-Land, then the three GCM scenarios.

    This metric gets one figure rather than the four framings the others get.
    There is no fixed-arm counterpart to difference against, so 'gcm_median'
    and 'gcm_by_model' would be re-cuts of a single quantity rather than
    comparisons, and the ERA5 panel belongs beside the scenarios instead of
    in a figure of its own.
    """
    import matplotlib.pyplot as plt
    xlabel, title = METRICS["sens_LMA"]
    era5 = d[d["dataset"] == "era5"]
    gcm = d[d["dataset"] != "era5"]
    cols = [("ERA5-Land", era5)]
    for sc in SCEN:
        sub = gcm[gcm["scenario"] == sc]
        if not sub.empty:
            cols.append((f"GCM {sc}", sub))
    if len(cols) == 1 and era5.empty:
        raise Missing("sens_LMA: neither ERA5 nor GCM rows present")
    fig, axes = plt.subplots(1, len(cols), figsize=figsize, sharex=True,
                             sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    for ax, (name, sub) in zip(axes, cols):
        panel(ax, sub)
        ax.set_title(name, fontsize=9)
    _finish(fig, axes, xlabel, title, out_png)
    plt.close(fig)


def gcm_median(d: pd.DataFrame) -> pd.DataFrame:
    """Median across GCMs, per station / scenario / variable."""
    g = (d.groupby(["scenario", "station", "pft", "variable"], observed=True)
           ["value"].median().reset_index())
    g["dataset"] = "gcm"; g["gcm"] = "median"
    return g


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--metrics", default=",".join(METRICS))
    ap.add_argument("--sets", default="era5,gcm,gcmmedian,gcmrows")
    ap.add_argument("--single-size", default="6.5x5.5")
    ap.add_argument("--scen-size", default="9.5x5.5")
    ap.add_argument("--grid-size", default="9.5x13")
    ap.add_argument("--lma-size", default="11x5.5")
    ap.add_argument("--flux-pool", default="stationyear",
                    choices=["stationyear", "station"],
                    help="pool the flux metric over station-years or sites")
    a = ap.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")

    def dims(s):
        w, h = s.lower().split("x")
        return float(w), float(h)

    try:
        root = a.results or resolve_out(".", create=False)
        # resolve_figure already IS $TC_FIGURES; passing "figures" nested it
        # one level deeper and wrote figures/figures/.
        out_dir = a.out or resolve_figure(".")
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    def read(name):
        p = Path(root) / name
        if not p.is_file():
            print(f"ERROR: {p} not found -- run station_metrics.py first",
                  file=sys.stderr)
            return None
        return pd.read_csv(p)

    M, S = read("station_metrics.csv"), read("station_sensitivity.csv")
    if M is None and S is None:
        return 1

    wanted = [m.strip() for m in a.metrics.split(",") if m.strip()]
    sets = [s.strip() for s in a.sets.split(",") if s.strip()]
    bad = [m for m in wanted if m not in METRICS]
    if bad:
        print(f"ERROR: unknown metric(s) {bad}; choose from {list(METRICS)}",
              file=sys.stderr)
        return 1

    failures = []
    for metric in wanted:
        try:
            d = series(M, S, metric, a.flux_pool, root)
        except Missing as e:
            failures.append(str(e)); print(f"SKIP {metric}: {e}"); continue
        era5 = d[d["dataset"] == "era5"]
        gcm = d[d["dataset"] != "era5"]
        print(f"{metric:<12} era5 {len(era5):>7} rows   gcm {len(gcm):>8} rows")
        if metric == "sens_LMA":
            png = out_dir / "effect_sens_LMA.png"
            try:
                fig_lma(d, png, dims(a.lma_size))
            except Missing as e:
                failures.append(f"sens_LMA: {e}"); print(f"  SKIP: {e}")
            continue
        jobs = []
        if "era5" in sets:
            jobs.append(("era5", lambda p, dd=era5: fig_single(
                dd, metric, p, dims(a.single_size))))
        if "gcm" in sets:
            jobs.append(("gcm", lambda p, dd=gcm: fig_scen(
                dd, metric, p, dims(a.scen_size), decompose=True)))
        if "gcmmedian" in sets:
            jobs.append(("gcm_median", lambda p, dd=gcm: fig_scen(
                gcm_median(dd), metric, p, dims(a.scen_size))))
        if "gcmrows" in sets:
            jobs.append(("gcm_by_model", lambda p, dd=gcm: fig_grid(
                dd, metric, p, dims(a.grid_size))))
        for tag, fn in jobs:
            png = out_dir / f"effect_{metric}_{tag}.png"
            try:
                fn(png)
            except Missing as e:
                failures.append(f"{metric}/{tag}: {e}")
                print(f"  SKIP {tag}: {e}")

    if failures:
        print(f"\n{len(failures)} figure(s) not produced:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
