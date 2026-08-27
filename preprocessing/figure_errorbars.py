"""RMSE and skill score by forest type: 4x4 error-bar panels, ERA5-Land only.

Rows are GPP, ET, H, LE. Columns pair each forest type with its drought subset:

    evergreen / all | evergreen / drought | deciduous / all | deciduous / drought

Every panel carries three error bars across the stations of that type -- RMSE of
the DYNAMIC arm, RMSE of the FIXED arm, and the skill score. Mean with +/-1 SD,
with the individual stations drawn behind so the spread is visible rather than
merely summarised: a mean that hides a bimodal fleet is worse than no mean.

TWO Y-AXES, ON PURPOSE. RMSE carries the flux's units (gC/m2/d, mm/d, W/m2) and
the skill score is dimensionless. Sharing one axis would either squash the skill
score against zero or force the RMSEs off-scale, and it would invite the two to
be read as comparable numbers. RMSE is on the left in its own units; SS is on the
right, centred on zero with its own scale, and marked in a different colour.

WHY SS IS NOT JUST A FUNCTION OF THE TWO RMSEs HERE. It is per station -- each
station's own 1 - RMSE_dyn/RMSE_fixed -- and the mean of a ratio is not the
ratio of the means. The SS bar answers "how does a typical station fare"; the
two RMSE bars answer "how large are the errors". They are related but not
redundant, and the difference between them is itself informative.

ERA5-LAND ONLY, by design. A GCM does not reproduce the actual weather a tower
measured, so a station-by-station RMSE against observations is only meaningful
for the reanalysis-driven runs.

NO MINIMUM-n FILTER. A station with three drought periods is drawn, and n is
printed on the panel, so the drought columns can be judged rather than silently
thinned. Correlation is not used here, so RMSE needs only one point to exist.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drought_labels as DL                                      # noqa: E402
from figure_skill_maps import (SEASON_FIGS, Missing, read_model,  # noqa: E402
                               read_sites, read_tower, rmse)
from results_dir import NoResultsDir, resolve_figure             # noqa: E402

# Row order as requested: GPP, ET, H, LE. Units are for the LEFT axis, which
# really is an RMSE in those units -- unlike the skill maps, where printing them
# would have been misleading.
ROWS = [("GPP", "GPP", "gC m$^{-2}$ d$^{-1}$"),
        ("ET",  "ET",  "mm d$^{-1}$"),
        ("H",   "H",   "W m$^{-2}$"),
        ("LE",  "QE",  "W m$^{-2}$")]

COLS = [("evergreen", "all"), ("evergreen", "drought"),
        ("deciduous", "all"), ("deciduous", "drought")]

C_RMSE_DYN, C_RMSE_FIX, C_SS = "#b2182b", "#4d4d4d", "#2166ac"


def per_station(model, tower, spei, var, threshold):
    """{station: {(subset): (rmse_fixed, rmse_dyn, ss, n)}}."""
    from collections import defaultdict
    rec = defaultdict(list)
    for (sid, year, period, v), (f, d) in model.items():
        if v != var:
            continue
        obs = tower.get((sid, year, period, var))
        if obs is None or not np.isfinite(obs):
            continue
        rec[sid].append((f, d, obs, spei.get((sid, year, period), np.nan)))

    out = {}
    for sid, r in rec.items():
        f = np.array([x[0] for x in r]); d = np.array([x[1] for x in r])
        o = np.array([x[2] for x in r]); s = np.array([x[3] for x in r])
        dry = np.isfinite(s) & (s <= threshold)
        got = {}
        for sub, m in (("all", np.ones(f.size, bool)), ("drought", dry)):
            if not m.any():
                continue
            rf, rd = rmse(f[m], o[m]), rmse(d[m], o[m])
            ss = (1 - rd / rf) if rf > 0 else np.nan
            got[sub] = (rf, rd, ss, int(m.sum()))
        if got:
            out[sid] = got
    return out


def make_figure(tables, sites, title, out_png, figsize):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(len(ROWS), len(COLS), figsize=figsize,
                             constrained_layout=True)
    for i, (label, mvar, unit) in enumerate(ROWS):
        per = tables.get(mvar, {})
        # ONE SCALE PER ROW, both axes. With per-panel autoscaling the four
        # columns of a row came out on different limits -- GPP topped at 1.6 for
        # evergreen and 1.7 for deciduous -- so comparing forest types by eye,
        # which is the entire point of this layout, silently compared different
        # rulers. Limits are computed across the whole row first.
        rr, sv = [], []
        for sid, g in per.items():
            if sid not in sites:
                continue
            for sub, vals in g.items():
                rr += [vals[0], vals[1]]
                sv.append(vals[2])
        rr = np.array([x for x in rr if np.isfinite(x)], float)
        sv = np.array([x for x in sv if np.isfinite(x)], float)
        rlim = (0, float(np.nanpercentile(rr, 99) * 1.15)) if rr.size else (0, 1)
        slim = float(max(0.05, np.nanpercentile(np.abs(sv), 95) * 1.25))             if sv.size else 0.5
        for j, (pft, sub) in enumerate(COLS):
            ax = axes[i, j]
            rf, rd, ss = [], [], []
            for sid, g in per.items():
                if sid not in sites or sites[sid][2] != pft or sub not in g:
                    continue
                a, b, c, _ = g[sub]
                rf.append(a); rd.append(b); ss.append(c)
            rf, rd = np.array(rf, float), np.array(rd, float)
            ss = np.array(ss, float)
            n = len(rf)

            axr = ax.twinx()
            for k, (vals, col) in enumerate(((rd, C_RMSE_DYN), (rf, C_RMSE_FIX))):
                v = vals[np.isfinite(vals)]
                if not v.size:
                    continue
                ax.scatter(np.full(v.size, k) + np.random.default_rng(0)
                           .uniform(-.09, .09, v.size), v, s=5, color=col,
                           alpha=.28, edgecolor="none", zorder=2)
                ax.errorbar([k], [v.mean()], yerr=[v.std(ddof=1) if v.size > 1 else 0],
                            fmt="o", ms=5, color=col, capsize=3, lw=1.4, zorder=3)
            v = ss[np.isfinite(ss)]
            if v.size:
                axr.scatter(np.full(v.size, 2) + np.random.default_rng(1)
                            .uniform(-.09, .09, v.size), v, s=5, color=C_SS,
                            alpha=.28, edgecolor="none", zorder=2)
                axr.errorbar([2], [v.mean()],
                             yerr=[v.std(ddof=1) if v.size > 1 else 0],
                             fmt="D", ms=5, color=C_SS, capsize=3, lw=1.4, zorder=3)
            axr.axhline(0, color=C_SS, lw=.5, ls=":", zorder=1)
            ax.set_ylim(*rlim)
            axr.set_ylim(-slim, slim)

            ax.set_xlim(-0.5, 2.5)
            ax.set_xticks([0, 1, 2])
            ax.set_xticklabels(["RMSE$_{dyn}$", "RMSE$_{fix}$", "SS"]
                               if i == len(ROWS) - 1 else [], fontsize=6)
            ax.tick_params(axis="y", labelsize=6,
                           labelleft=(j == 0))
            axr.tick_params(axis="y", labelsize=6, colors=C_SS)
            ax.text(.03, .95, f"n={n}", transform=ax.transAxes, fontsize=5.5,
                    va="top", color="0.35")
            if i == 0:
                ax.set_title(f"{pft}\n{sub} steps", fontsize=8)
            if j == 0:
                ax.set_ylabel(f"{label}\nRMSE ({unit})", fontsize=7)
            if j == len(COLS) - 1:
                axr.set_ylabel("skill score", fontsize=7, color=C_SS)
            else:
                axr.set_yticklabels([])

    handles = [
        Line2D([], [], marker="o", ls="", color=C_RMSE_DYN, label="RMSE dynamic"),
        Line2D([], [], marker="o", ls="", color=C_RMSE_FIX, label="RMSE fixed"),
        Line2D([], [], marker="D", ls="", color=C_SS,
               label="skill score (right axis)")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=7, bbox_to_anchor=(0.5, -0.03))
    fig.suptitle(f"{title}\nmean $\\pm$ 1 SD across stations; "
                 f"points are individual stations", fontsize=9)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model-dir", type=Path, required=True)
    ap.add_argument("--tower-dir", type=Path, required=True)
    ap.add_argument("--step", default="all",
                    choices=["all", "annual", "monthly"] + SEASON_FIGS)
    ap.add_argument("--gpp", default="NT", choices=["NT", "DT"])
    ap.add_argument("--threshold", type=float, default=-1.0)
    ap.add_argument("--figsize", default="6.5x10")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)

    w, h = (float(v) for v in a.figsize.lower().split("x"))
    try:
        out_dir = resolve_figure(a.out or ".")
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    sites = read_sites()
    steps = ([("annual", "annual", ["ANN"]), ("monthly", "monthly", None)] +
             [(s, "seasonal", [s]) for s in SEASON_FIGS])
    if a.step != "all":
        steps = [s for s in steps if s[0] == a.step]

    for name, freq, keep in steps:
        try:
            model = read_model(a.model_dir / f"era5_{freq}.csv", freq)
            if keep is not None:
                model = {k: v for k, v in model.items() if k[2] in keep}
            spei = DL.station_spei({s: (v[0], v[1]) for s, v in sites.items()}, freq)
            tower = read_tower(a.tower_dir, freq, sites, a.gpp)
        except (Missing, DL.NoLabel) as e:
            print(f"ERROR [{name}]: {e}", file=sys.stderr)
            return 1
        tables = {mv: per_station(model, tower, spei, mv, a.threshold)
                  for _, mv, _ in ROWS}
        for _, mv, _ in ROWS:
            t = tables[mv]
            nd = sum(1 for g in t.values() if "drought" in g)
            print(f"  {name:<9}{mv:<5} {len(t):>3} stations, {nd:>3} with a "
                  f"drought subset", flush=True)
        png = out_dir / f"errorbars_{name}.png"
        make_figure(tables, sites, f"RMSE and skill score, {name} steps",
                    png, (w, h))
        print(f"  -> {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
