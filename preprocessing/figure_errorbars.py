"""RMSE and skill score by forest type: 4x4 error-bar panels, ERA5-Land only.

Rows are GPP, ET, H, LE. Columns pair each forest type with its drought subset:

    evergreen / all | evergreen / drought | deciduous / all | deciduous / drought

Every panel carries three error bars across the stations of that type -- RMSE of
the DYNAMIC arm, RMSE of the FIXED arm, and the skill score. Mean with +/-1 SD,
with the individual stations drawn behind so the spread is visible rather than
merely summarised: a mean that hides a bimodal fleet is worse than no mean.

ONE DIMENSIONLESS AXIS. The error bars are RSR = RMSE / SD(observed), not raw
RMSE. RMSE carries the flux's units, which forced a twin axis and left the four
rows on incomparable scales -- an RMSE of 0.8 gC/m2/d and one of 0.8 W/m2 sat at
the same height meaning nothing alike. RSR is unitless (both terms rescale
together, verified exactly across three unit systems), so GPP, ET, H and LE now
share one ruler and "which flux does dynamic LMA help most" is readable off the
figure rather than inferred.

RSR = 1 is drawn as a reference: there the model's error equals the observed
standard deviation, i.e. no better than predicting the observed mean, which is
NSE = 0. Moriasi et al. (2007) rate RSR <= 0.5 very good and > 0.7
unsatisfactory; the 0.5 line is shaded.

THE SKILL SCORE IS UNCHANGED BY THIS. The observed SD cancels between the arms,
so 1 - RSR_dyn/RSR_fixed equals 1 - RMSE_dyn/RMSE_fixed identically. Only the
bars are renormalised.

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

# Row order as requested: GPP, ET, H, LE. The third field is kept only for the
# CSV, where the raw RMSE is still dimensional; the axis itself is unitless.
ROWS = [("GPP", "GPP", "gC m$^{-2}$ d$^{-1}$"),
        ("ET",  "ET",  "mm d$^{-1}$"),
        ("H",   "H",   "W m$^{-2}$"),
        ("LE",  "QE",  "W m$^{-2}$")]

# TWO COLUMNS, one per forest type, all steps. The drought columns are gone:
# a drought subset rests on far fewer periods than the all-steps one, so the
# pair could be read as a drought effect when part of it is sample size.
# Drought is still computed and written to the CSV beside these.
COLS = [("evergreen", "all"), ("deciduous", "all")]

C_RMSE_DYN, C_RMSE_FIX, C_SS = "#b2182b", "#4d4d4d", "#2166ac"
# Moriasi et al. (2007) performance bands for RSR.
RSR_VERY_GOOD, RSR_UNSATISFACTORY = 0.5, 0.7


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
            sdo = float(np.std(o[m], ddof=1)) if m.sum() > 1 else np.nan
            # RSR = RMSE / SD(obs). Raw RMSE is kept alongside so the
            # dimensional values remain available in the CSV.
            rsr_f = rf / sdo if (np.isfinite(sdo) and sdo > 0) else np.nan
            rsr_d = rd / sdo if (np.isfinite(sdo) and sdo > 0) else np.nan
            got[sub] = (rf, rd, ss, int(m.sum()), sdo, rsr_f, rsr_d)
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
        # One scale per row. RSR and SS are both dimensionless and both O(1),
        # so they share a single axis -- no twin spine, and the four columns of
        # a row are directly comparable.
        vals_all = []
        for sid, g in per.items():
            if sid in sites:
                for v in g.values():
                    vals_all += [v[5], v[6], v[2]]
        vals_all = np.array([x for x in vals_all if np.isfinite(x)], float)
        hi = float(np.nanpercentile(vals_all, 98) * 1.15) if vals_all.size else 1.5
        lo = float(min(-0.25, np.nanpercentile(vals_all, 2) * 1.15))             if vals_all.size else -0.25
        for j, (pft, sub) in enumerate(COLS):
            ax = axes[i, j]
            rsr_f, rsr_d, ss = [], [], []
            for sid, g in per.items():
                if sid not in sites or sites[sid][2] != pft or sub not in g:
                    continue
                v = g[sub]
                rsr_f.append(v[5]); rsr_d.append(v[6]); ss.append(v[2])
            series = [(np.array(rsr_d, float), C_RMSE_DYN, "o"),
                      (np.array(rsr_f, float), C_RMSE_FIX, "o"),
                      (np.array(ss, float), C_SS, "D")]
            n = int(np.isfinite(series[0][0]).sum())

            # RSR = 1: error equal to the observed SD, i.e. no better than the
            # observed mean (NSE = 0). Below 0.5 is Moriasi's "very good".
            ax.axhspan(lo, RSR_VERY_GOOD, color="#4daf4a", alpha=.05, zorder=0)
            ax.axhline(1.0, color="0.55", lw=.6, ls="--", zorder=1)
            ax.axhline(0.0, color="0.55", lw=.6, ls=":", zorder=1)
            for k, (v, col, mk) in enumerate(series):
                v = v[np.isfinite(v)]
                if not v.size:
                    continue
                ax.scatter(np.full(v.size, k) + np.random.default_rng(k)
                           .uniform(-.09, .09, v.size), v, s=5, color=col,
                           alpha=.28, edgecolor="none", zorder=2)
                ax.errorbar([k], [v.mean()],
                            yerr=[v.std(ddof=1) if v.size > 1 else 0],
                            fmt=mk, ms=5, color=col, capsize=3, lw=1.4, zorder=3)
            ax.set_ylim(lo, hi)
            ax.set_xlim(-0.5, 2.5)
            ax.set_xticks([0, 1, 2])
            ax.set_xticklabels(["RSR$_{dyn}$", "RSR$_{fix}$", "SS"]
                               if i == len(ROWS) - 1 else [], fontsize=6)
            ax.tick_params(axis="y", labelsize=6, labelleft=(j == 0))
            ax.text(.03, .95, f"n={n}", transform=ax.transAxes, fontsize=5.5,
                    va="top", color="0.35")
            if i == 0:
                ax.set_title(pft + "\n" + sub + " steps", fontsize=8)
            if j == 0:
                ax.set_ylabel(label, fontsize=9)

    handles = [
        Line2D([], [], marker="o", ls="", color=C_RMSE_DYN, label="RSR dynamic"),
        Line2D([], [], marker="o", ls="", color=C_RMSE_FIX, label="RSR fixed"),
        Line2D([], [], marker="D", ls="", color=C_SS, label="skill score"),
        Line2D([], [], ls="--", color="0.55", label="RSR = 1 (NSE = 0)")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=7, bbox_to_anchor=(0.5, -0.03))
    sub2 = ("RSR = RMSE/SD(obs), dimensionless; mean $\\pm$ 1 SD across "
            "stations, points are individual stations")
    fig.suptitle(title + chr(10) + sub2, fontsize=9)
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
        import csv as _csv
        csv_path = out_dir / f"errorbars_{name}.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as fh:
            w = _csv.writer(fh)
            w.writerow(["station", "pft", "variable", "subset", "n",
                        "sd_obs", "rmse_fixed", "rmse_dyn",
                        "rsr_fixed", "rsr_dyn", "skill_score"])
            for _, mv, _u in ROWS:
                for sid, g in sorted(tables[mv].items()):
                    pft = sites[sid][2] if sid in sites else ""
                    for sub, v in sorted(g.items()):
                        w.writerow([sid, pft, mv, sub, v[3],
                                    f"{v[4]:.6g}", f"{v[0]:.6g}", f"{v[1]:.6g}",
                                    f"{v[5]:.6g}", f"{v[6]:.6g}", f"{v[2]:.6g}"])
        print(f"  -> {csv_path}")
        png = out_dir / f"errorbars_{name}.png"
        make_figure(tables, sites, f"RMSE and skill score, {name} steps",
                    png, (w, h))
        print(f"  -> {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
