"""Taylor diagrams: how each arm stands against the towers, station by station.

Same 4x2 portrait shape as figure_skill_maps.py and the same six time steps
(annual, monthly, DJF, MAM, JJA, SON). Rows are GPP, ET, LE, H; the left column
uses every step and the right only the drought ones.

    marker shape   circle = deciduous, triangle = evergreen
    marker fill    FILLED = fixed LMA, OUTLINED = dynamic LMA

so each station contributes two markers per panel and the question the figure
answers is whether the outlined marker sits closer to REF than the filled one.

WHAT A TAYLOR DIAGRAM SHOWS. Radius is the model's standard deviation, angle is
its correlation with the observations, and -- because those two plus the
observed standard deviation fix it -- the distance from the REF point is the
centred RMS difference. A point at REF matches the tower in both variability and
phase. Radius > 1 means the model is too variable, radius < 1 too flat, and
swinging away from the x-axis means the timing is off.

NORMALISED BY THE OBSERVED SIGMA, so REF sits at 1.0 in every panel. Without
that a single panel would mix stations whose fluxes differ by an order of
magnitude and the spread would be about site productivity rather than about
model skill. It also lets the four rows share one axis convention despite
carrying gC/m2/d, mm/d and W/m2.

NEGATIVE CORRELATIONS ARE NOT DISCARDED. A drought-only subset can be short
enough for the correlation to go negative; those points are real and the axes
extend to a half circle when any appear, rather than dropping the station and
quietly shrinking the sample. The log says when this happens.

STATIONS NEED ENOUGH STEPS. A correlation over four months is not a measurement,
so --min-n applies here as it does to the skill maps, and a station whose
observed sigma is zero is dropped by name.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drought_labels as DL                                      # noqa: E402
from figure_skill_maps import (ROWS, SEASON_FIGS, Missing,       # noqa: E402
                               read_model, read_sites, read_tower)
from results_dir import NoResultsDir, resolve_figure             # noqa: E402

# Correlation rays to draw. Crowded near 1 on purpose: that is where a good
# model sits, and evenly spaced ticks would leave the interesting region blank.
CORR_TICKS = [0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
CORR_TICKS_NEG = [-0.9, -0.6, -0.3] + CORR_TICKS


def stats(model: dict, tower: dict, spei: dict, var: str,
          threshold: float, min_n: int):
    """{station: {(cond, arm): (sigma_hat, corr, n)}} for one variable."""
    per_st = defaultdict(list)
    for (sid, year, period, v), (f, d) in model.items():
        if v != var:
            continue
        obs = tower.get((sid, year, period, var))
        if obs is None or not np.isfinite(obs):
            continue
        per_st[sid].append((f, d, obs, spei.get((sid, year, period), np.nan)))

    out, dropped = {}, []
    for sid, recs in per_st.items():
        f = np.array([r[0] for r in recs]); d = np.array([r[1] for r in recs])
        o = np.array([r[2] for r in recs]); s = np.array([r[3] for r in recs])
        dry = np.isfinite(s) & (s <= threshold)
        got = {}
        for cond, m in (("all", np.ones(f.size, bool)), ("drought", dry)):
            if m.sum() < min_n:
                continue
            so = float(np.std(o[m]))
            if so <= 0:
                dropped.append(f"{sid}/{cond}: observed sigma is zero")
                continue
            for arm, series in (("fixed", f), ("dyn", d)):
                sm = float(np.std(series[m]))
                if sm <= 0:
                    dropped.append(f"{sid}/{cond}/{arm}: model sigma is zero")
                    continue
                r = float(np.corrcoef(series[m], o[m])[0, 1])
                if not np.isfinite(r):
                    dropped.append(f"{sid}/{cond}/{arm}: correlation undefined")
                    continue
                got[(cond, arm)] = (sm / so, r, int(m.sum()))
        if got:
            out[sid] = got
    return out, dropped


# --------------------------------------------------------------------- drawing
def draw_axes(ax, smax: float, half: bool, xlabel: bool = False):
    """Correlation rays, sigma arcs and centred-RMSD contours about REF."""
    ticks = CORR_TICKS_NEG if half else CORR_TICKS
    th_max = np.pi if half else np.pi / 2
    a = np.linspace(0, th_max, 300)
    sig_ticks = [s for s in np.arange(0.5, smax + 0.01, 0.5)]

    # sigma arcs, centred on the origin
    for s in sig_ticks:
        ax.plot(s * np.cos(a), s * np.sin(a), color="0.86", lw=0.4, zorder=0)
    # correlation rays
    for r in ticks:
        th = np.arccos(r)
        ax.plot([0, smax * np.cos(th)], [0, smax * np.sin(th)],
                color="0.86", lw=0.4, zorder=0)
        ax.text(smax * 1.05 * np.cos(th), smax * 1.05 * np.sin(th), f"{r:g}",
                fontsize=5, ha="center", va="center", color="0.35")
    # Centred-RMSD contours: circles about REF at (1, 0), clipped to the wedge.
    # These are the whole point of the diagram -- distance from REF IS the
    # centred RMS difference -- so they are drawn dark enough to see, which the
    # first version at colour "0.9" was not.
    th = np.linspace(0, 2 * np.pi, 400)
    for rms in np.arange(0.25, smax + 0.5, 0.25):
        xc, yc = 1 + rms * np.cos(th), rms * np.sin(th)
        keep = (yc >= -1e-9) & (np.hypot(xc, yc) <= smax)
        if not half:
            keep &= (xc >= -1e-9)
        if not keep.any():
            continue
        ax.plot(np.where(keep, xc, np.nan), np.where(keep, yc, np.nan),
                color="#9ecae1", lw=0.45, ls=(0, (2, 2)), zorder=0)

    # RADIAL TICK LABELS. Without them the radius is unreadable -- the reader
    # can see a point is far out but not that it is 1.5x the tower's variability.
    ax.plot([0, smax], [0, 0], color="0.5", lw=0.5, zorder=1)
    for s in sig_ticks:
        ax.plot([s, s], [0, -smax * 0.018], color="0.5", lw=0.5, zorder=1)
        ax.text(s, -smax * 0.055, f"{s:g}", fontsize=5.5, ha="center",
                va="top", color="0.3")
    if xlabel:
        ax.text(smax * 0.5, -smax * 0.17, "sd(model) / sd(tower)",
                fontsize=6, ha="center", va="top", color="0.3")

    ax.plot([1], [0], marker="*", ms=8, color="0.15", zorder=5)
    ax.set_xlim(-smax * 1.10 if half else -smax * 0.06, smax * 1.14)
    ax.set_ylim(-smax * 0.20, smax * 1.14)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def make_figure(tables: dict, sites: dict, title: str, out_png: Path,
                figsize, min_n: int) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    # Does anything need the left half-plane?
    allr = [v[1] for t in tables.values() for st in t.values() for v in st.values()]
    half = bool(allr) and min(allr) < 0
    smax = 1.6
    if allr:
        sh = [v[0] for t in tables.values() for st in t.values() for v in st.values()]
        smax = float(min(3.0, max(1.6, np.percentile(sh, 98) * 1.1)))

    fig, axes = plt.subplots(len(ROWS), 2, figsize=figsize,
                             constrained_layout=True)
    for i, (label, mvar, unit) in enumerate(ROWS):
        for j, cond in enumerate(("all", "drought")):
            ax = axes[i, j]
            draw_axes(ax, smax, half, xlabel=(i == len(ROWS) - 1))
            got = tables.get(mvar, {})
            for sid, entries in got.items():
                if sid not in sites:
                    continue
                _, _, ft = sites[sid]
                mk = "o" if ft == "deciduous" else "^"
                for arm in ("fixed", "dyn"):
                    v = entries.get((cond, arm))
                    if v is None:
                        continue
                    sh, r, _ = v
                    th = np.arccos(np.clip(r, -1, 1))
                    x, y = sh * np.cos(th), sh * np.sin(th)
                    if arm == "fixed":                       # FILLED
                        ax.scatter([x], [y], marker=mk, s=22, c="#4d4d4d",
                                   edgecolor="none", alpha=0.75, zorder=4)
                    else:                                    # OUTLINED
                        ax.scatter([x], [y], marker=mk, s=26, facecolor="none",
                                   edgecolor="#b2182b", linewidth=0.7, zorder=4)
            if i == 0:
                ax.set_title("all steps" if cond == "all" else "drought only",
                             fontsize=9)
            if j == 0:
                ax.set_ylabel(f"{label}\n({unit})", fontsize=8)

    handles = [
        Line2D([], [], marker="o", ls="", color="#4d4d4d", label="deciduous, fixed"),
        Line2D([], [], marker="^", ls="", color="#4d4d4d", label="evergreen, fixed"),
        Line2D([], [], marker="o", ls="", mfc="none", mec="#b2182b",
               color="none", label="deciduous, dynamic"),
        Line2D([], [], marker="^", ls="", mfc="none", mec="#b2182b",
               color="none", label="evergreen, dynamic"),
        Line2D([], [], marker="*", ls="", color="0.15", label="tower (REF)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=7, bbox_to_anchor=(0.5, -0.045))
    fig.suptitle(f"{title}\nradius = sd(model)/sd(tower), angle = correlation; "
                 f"closer to REF is better", fontsize=9)
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
    ap.add_argument("--min-n", type=int, default=6)
    ap.add_argument("--figsize", default="6x10")
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

        tables, dropped = {}, []
        for _, mvar, _ in ROWS:
            t, dr = stats(model, tower, spei, mvar, a.threshold, a.min_n)
            tables[mvar] = t
            dropped += [f"{mvar} {x}" for x in dr]
        n = {mv: len(t) for mv, t in tables.items()}
        neg = sum(1 for t in tables.values() for st in t.values()
                  for v in st.values() if v[1] < 0)
        print(f"  {name}: stations per row {n}"
              + (f"; {neg} negative correlation(s) -- axes extended" if neg else ""))
        if dropped:
            print(f"    dropped {len(dropped)}: {dropped[0]}"
                  + (f" ... and {len(dropped)-1} more" if len(dropped) > 1 else ""))
        png = out_dir / f"taylor_{name}.png"
        make_figure(tables, sites, f"Taylor diagram, {name} steps", png, (w, h),
                    a.min_n)
        print(f"  -> {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
