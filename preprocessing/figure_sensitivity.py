#!/usr/bin/env python3
"""Regression slopes per flux, in their own units. No ratios anywhere.

Six panels, one per variable, each with its OWN y axis in that variable's
slope units -- so nothing needs standardising and nothing gets divided.
The x axis is categorical: deciduous and evergreen side by side within a
dataset, a wider gap, then the next dataset.

    era5   historical   ssp126   ssp585
    D E      D E         D E      D E

WHY ABSOLUTE AND NOT PERCENT. A slope is a difference quantity: it passes
through zero legitimately and changes sign across stations, so
100*(b_dyn-b_fixed)/|b_fixed| runs to +/-inf on both sides of zero and two
stations with nearly identical physics land at +8000% and -8000%. A symlog
axis would render that legibly, which is worse than not rendering it --
the artefacts would look like findings. The difference b_dyn - b_fixed has
no denominator and cannot do this.

It also answers what a ratio cannot: "-20% of GPP sensitivity to
temperature" is a different physical statement depending on the sign of
the fixed slope. The supplementary figures show each arm's slope on its
own, so reduced-versus-grew-less is readable directly.

    --metric ta        delta_slope, dFlux/dTa        MAIN
    --metric spei      delta_slope, dFlux/dSPEI12    MAIN
    --metric lma       slope_dyn, dFlux/dLMA         MAIN
    --metric ta_fixed / ta_dyn / spei_fixed / spei_dyn   SUPPLEMENTARY

LMA has no fixed-arm slope to difference against -- that arm holds LMA
constant, so station_metrics sets it NaN on purpose -- which is why the
LMA figure shows the dynamic slope itself. In its own units it needs no
standardising, so the station-variability contamination of a standardised
slope does not arise.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from results_dir import NoResultsDir, resolve_figure, resolve_out  # noqa: E402

VARS = [("GPP", "GPP"), ("LAI", "LAI_H"), ("TR", "T"),
        ("ET", "ET"), ("LE", "QE"), ("H", "H")]
C_DEC, C_EVE = "#1b7837", "#2166ac"
PFTS = [("deciduous", C_DEC), ("evergreen", C_EVE)]
DATASETS = ["era5", "historical", "ssp126", "ssp585"]
LABELS = {"era5": "ERA5-Land", "historical": "GCM hist",
          "ssp126": "GCM ssp126", "ssp585": "GCM ssp585"}

# metric -> (predictor, column, y-axis label)
METRICS = {
    "ta":         ("Ta", "delta_slope", "dFlux/dTa   dynamic - fixed"),
    "spei":       ("SPEI12", "delta_slope", "dFlux/dSPEI12   dynamic - fixed"),
    "lma":        ("LMA", "slope_dyn", "dFlux/dLMA   dynamic arm"),
    "ta_fixed":   ("Ta", "slope_fixed", "dFlux/dTa   fixed LMA"),
    "ta_dyn":     ("Ta", "slope_dyn", "dFlux/dTa   dynamic LMA"),
    "spei_fixed": ("SPEI12", "slope_fixed", "dFlux/dSPEI12   fixed LMA"),
    "spei_dyn":   ("SPEI12", "slope_dyn", "dFlux/dSPEI12   dynamic LMA"),
}
# Numerator: the flux's own units. These panels plot SLOPES, so the axis unit
# is flux per unit of predictor -- the denominator below. Labelling them with
# the bare flux unit said the panel showed a flux, which it does not.
FLUX_UNITS = {"GPP": "gC m$^{-2}$ yr$^{-1}$", "LAI": "m$^2$ m$^{-2}$",
              "TR": "mm yr$^{-1}$", "ET": "mm yr$^{-1}$",
              "LE": "W m$^{-2}$", "H": "W m$^{-2}$"}
# Denominator: the predictor's units. SPEI is a standardised index and is
# therefore dimensionless, so "per SPEI unit" is the honest phrasing rather
# than inventing one. LMA is dry mass per area.
PRED_UNITS = {"Ta": " K$^{-1}$", "SPEI12": " / SPEI unit",
              "LMA": " / (g m$^{-2}$)"}


def yunit(disp: str, pred: str) -> str:
    return FLUX_UNITS.get(disp, "") + PRED_UNITS.get(pred, "")


class Missing(Exception):
    """A required column or subset is absent. Never substituted."""


def positions():
    """x positions: two per dataset, a wider gap between datasets."""
    pos, ticks = {}, []
    x = 0.0
    for ds in DATASETS:
        for k, (pft, _) in enumerate(PFTS):
            pos[(ds, pft)] = x + k * 0.62
        ticks.append(x + 0.31)
        x += 2.05
    return pos, ticks


def panel(ax, d, var_col, pos, clip_pct):
    """Boxes with a translucent violin over each. Returns points drawn."""
    n = 0
    for ds in DATASETS:
        for pft, colr in PFTS:
            v = d.loc[(d["variable"] == var_col) & (d["dataset"] == ds)
                      & (d["pft"] == pft), "value"].to_numpy(float)
            v = v[np.isfinite(v)]
            if v.size == 0:
                continue
            n += v.size
            x = pos[(ds, pft)]
            if v.size > 2 and np.ptp(v) > 0:
                try:
                    pc = ax.violinplot([v], positions=[x], widths=0.55,
                                       showextrema=False, showmedians=False)
                except TypeError:                       # older matplotlib
                    pc = ax.violinplot([v], positions=[x], widths=0.55)
                for b in pc["bodies"]:
                    b.set_facecolor(colr); b.set_alpha(0.22)
                    b.set_edgecolor("none"); b.set_zorder(3)
            bp = ax.boxplot([v], positions=[x], widths=0.30, showfliers=False,
                            patch_artist=True, manage_ticks=False)
            for b in bp["boxes"]:
                b.set_facecolor("white"); b.set_edgecolor(colr)
                b.set_linewidth(1.0); b.set_zorder(4)
            for part in ("whiskers", "caps"):
                for it in bp[part]:
                    it.set_color(colr); it.set_linewidth(0.9); it.set_zorder(4)
            for m in bp["medians"]:
                m.set_color(colr); m.set_linewidth(1.8); m.set_zorder(5)
    ax.axhline(0, color="0.15", lw=0.9, zorder=2)
    if clip_pct:
        # A handful of stations can stretch a native-unit axis the way one
        # station stretched the Taylor radius. Scale to a percentile and SAY
        # how many fall outside rather than dropping them from the boxes --
        # the boxes are computed on everything either way.
        v = d.loc[d["variable"] == var_col, "value"].to_numpy(float)
        v = v[np.isfinite(v)]
        if v.size > 4:
            lo, hi = np.nanpercentile(v, [clip_pct, 100 - clip_pct])
            pad = 0.12 * (hi - lo) if hi > lo else 1.0
            off = int(((v < lo - pad) | (v > hi + pad)).sum())
            ax.set_ylim(lo - pad, hi + pad)
            if off:
                ax.text(0.99, 0.02, f"{off} off scale", transform=ax.transAxes,
                        ha="right", va="bottom", fontsize=5.5, color="0.45")
    return n


def build(S, metric, out_png, figsize, clip_pct):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    pred, col, ylab = METRICS[metric]
    d = S[(S["freq"] == "annual") & (S["subset"] == "all")
          & (S["predictor"] == pred)].copy()
    if d.empty:
        have = sorted(S["predictor"].dropna().unique())
        raise Missing(f"{metric}: no rows for predictor={pred!r}; "
                      f"table has {have}")
    if col not in d.columns:
        raise Missing(f"{metric}: column {col!r} absent from the table")
    for c in ("dataset", "pft", "variable"):
        if c not in d.columns:
            raise Missing(f"{metric}: column {c!r} absent from the table")
    d["value"] = pd.to_numeric(d[col], errors="coerce")
    if not np.isfinite(d["value"]).any():
        raise Missing(f"{metric}: every {col} value is NaN for predictor "
                      f"{pred!r} -- nothing to plot")

    pos, ticks = positions()
    fig, axes = plt.subplots(3, 2, figsize=figsize, constrained_layout=True)
    total = 0
    for ax, (disp, vcol) in zip(axes.ravel(), VARS):
        total += panel(ax, d, vcol, pos, clip_pct)
        ax.set_title(disp, fontsize=9, loc="left")
        ax.set_ylabel(yunit(disp, pred), fontsize=6.5)
        ax.set_xticks(ticks)
        ax.set_xticklabels([LABELS[x] for x in DATASETS], fontsize=7)
        ax.tick_params(axis="y", labelsize=7)
        ax.grid(axis="y", color="0.90", lw=0.6)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    if total == 0:
        raise Missing(f"{metric}: every panel empty after filtering")
    handles = [Patch(facecolor=c, alpha=.35, edgecolor=c, label=p)
               for p, c in PFTS]
    fig.suptitle(ylab, fontsize=10, x=0.01, ha="left")
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, 0.0))
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {out_png}   ({total} station values)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--metrics", default=",".join(METRICS))
    ap.add_argument("--size", default="6.5x9")
    ap.add_argument("--clip-pct", type=float, default=1.0,
                    help="y axis spans this to 100-this percentile; 0 = no clip")
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

    p = root / "station_sensitivity.csv"
    if not p.is_file():
        print(f"ERROR: {p} not found -- run station_metrics.py first",
              file=sys.stderr)
        return 1
    S = pd.read_csv(p, low_memory=False)

    wanted = [m.strip() for m in a.metrics.split(",") if m.strip()]
    bad = [m for m in wanted if m not in METRICS]
    if bad:
        print(f"ERROR: unknown metric(s) {bad}; choose from {list(METRICS)}",
              file=sys.stderr)
        return 1

    failures = []
    for metric in wanted:
        try:
            build(S, metric, out_dir / f"sens_{metric}.png",
                  dims(a.size), a.clip_pct)
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
