"""Skill-score maps: does dynamic LMA beat fixed LMA against the towers?

One figure per time step (annual, monthly, DJF, MAM, JJA, SON). Each is a 4x2
portrait panel: rows are GPP, ET, LE, H; the left column uses every time step and
the right only the drought ones.

    SS = 1 - RMSE_dyn / RMSE_fixed

positive where the dynamic arm is closer to the tower, negative where it is
worse, zero where the two are indistinguishable -- which is why the colour scale
is white in the middle. This is the Murphy skill score, with the FIXED arm as the
reference forecast, and it is the standard way to ask "did the added mechanism
help". Note the request read "1 - RMSE dyn - RMSE fixed"; a plain difference of
two RMSEs is not a skill score (it carries the variable's units, is unbounded,
and has no natural centre), so it is read here as the ratio. --skill diff selects
a symmetric alternative, (RMSE_fixed - RMSE_dyn) / (RMSE_fixed + RMSE_dyn), which
is bounded [-1, 1] and also white at zero.

UNITS ARE NORMALISED TO PER-DAY INTENSITIES before any RMSE is taken -- gC/m2/d,
mm/d, W/m2. The model table stores monthly GPP as a SUM and monthly LE as a MEAN,
while the towers report daily rates throughout, so without this an annual RMSE
and a monthly RMSE would not be the same quantity and the six figures could not
be read against each other. n_days in the model table is what makes the
conversion exact rather than assuming 30.

DROUGHT USES THE MATCHING ACCUMULATION -- SPEI-3 for a month, SPEI-3 at the
season's last month, SPEI-12 at the water-year end for a year. See
drought_labels.py.

STATIONS. Only those with BOTH an ERA5-Land pair and tower data are drawn;
deciduous are circles, evergreen triangles. A station with too few overlapping
steps is dropped and named rather than plotted from a two-point RMSE.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drought_labels as DL                                      # noqa: E402
from results_dir import NoResultsDir, resolve_figure            # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
SITE_LISTS = [REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
              REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv"]

# Row order of the figure, and the model variable behind each.
ROWS = [("GPP", "GPP", "gC m$^{-2}$ d$^{-1}$"),
        ("ET",  "ET",  "mm d$^{-1}$"),
        ("LE",  "QE",  "W m$^{-2}$"),
        ("H",   "H",   "W m$^{-2}$")]

# Latent heat of vaporisation, J/kg. LE [W/m2] -> ET [mm/d] is LE*86400/(LAMBDA).
LAMBDA = 2.45e6
W_TO_MM_D = 86400.0 / LAMBDA          # 0.03527 mm/d per W/m2

# Model variables stored as a SUM over the period; the rest are means. Dividing
# the sums by n_days puts everything on a per-day footing.
MODEL_SUMS = {"GPP", "ET"}

SEASON_FIGS = ["DJF", "MAM", "JJA", "SON"]


class Missing(Exception):
    """Something needed is absent, named so it is never silently skipped."""


# ---------------------------------------------------------------- station table
def read_sites() -> dict:
    """{station: (lat, lon, 'deciduous'|'evergreen')}."""
    out = {}
    for path in SITE_LISTS:
        if not path.exists():
            raise Missing(f"site list not found: {path}")
        with path.open(newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                sid = (row.get("StationID") or "").strip()
                ft = (row.get("ForestType") or "").strip().lower()
                try:
                    lat, lon = float(row["Lat"]), float(row["Lon"])
                except (KeyError, TypeError, ValueError):
                    continue
                if sid:
                    out.setdefault(sid, (lat, lon, ft))
    return out


# ------------------------------------------------------------------ model side
def read_model(path: Path, freq: str) -> dict:
    """{(station, year, period, var): (fixed, dyn)} as per-day intensities."""
    if not path.is_file():
        raise Missing(f"model table not found: {path}")
    want = {m for _, m, _ in ROWS}
    out = {}
    with path.open(newline="", encoding="utf-8-sig") as fh:
        rd = csv.DictReader(fh)
        need = {"station", "year", "period", "variable", "fixed", "dyn", "n_days"}
        miss = need - set(rd.fieldnames or [])
        if miss:
            raise Missing(f"{path.name} lacks {', '.join(sorted(miss))}")
        for r in rd:
            var = r["variable"]
            if var not in want:
                continue
            try:
                f, d, n = float(r["fixed"]), float(r["dyn"]), int(r["n_days"])
                year = int(r["year"])
            except (TypeError, ValueError):
                continue
            if n <= 0 or not (np.isfinite(f) and np.isfinite(d)):
                continue
            period = int(r["period"]) if freq == "monthly" else r["period"]
            if var in MODEL_SUMS:            # sum over the period -> per day
                f, d = f / n, d / n
            out[(r["station"], year, period, var)] = (f, d)
    return out


# ------------------------------------------------------------------ tower side
# FLUXNET (ONEFlux) column names. Each entry is the ordered list of candidates
# tried, so a site missing the preferred name falls back to a documented
# alternative rather than to silence. Verified names are confirmed once a real
# archive lands; until then a file that matches none of these is reported with
# the columns it DOES have, which is the only useful thing to say.
TOWER_COLS = {
    "GPP_NT": ["GPP_NT_VUT_REF", "GPP_NT_VUT_MEAN", "GPP_NT_CUT_REF"],
    "GPP_DT": ["GPP_DT_VUT_REF", "GPP_DT_VUT_MEAN", "GPP_DT_CUT_REF"],
    "LE":     ["LE_F_MDS", "LE_CORR", "LE"],
    "H":      ["H_F_MDS", "H_CORR", "H"],
}
# Quality flags, where present: ONEFlux uses 0 = measured, 1 = good gap-fill.
QC_OF = {"LE": "LE_F_MDS_QC", "H": "H_F_MDS_QC",
         "GPP_NT": "NEE_VUT_REF_QC", "GPP_DT": "NEE_VUT_REF_QC"}
FLUXNET_NA = {-9999.0, -9999}


def _pick(header, cands):
    for c in cands:
        if c in header:
            return c
    return None


def _fluxnet_file(root: Path, sid: str, res: str):
    """The <res> (MM/YY) FLUXNET csv for one site, wherever the unpack put it."""
    pats = [f"**/*{sid}*FLUXNET*{res}*.csv", f"**/*{sid}*_{res}_*.csv"]
    for pat in pats:
        hits = sorted(Path(root).glob(pat))
        if hits:
            return hits[0]
    return None


def read_tower(root: Path, freq: str, sites, gpp: str = "NT",
               max_qc: float = 1.0) -> dict:
    """{(station, year, period, var): observed} as per-day intensities.

    GPP IS NOT A MEASUREMENT. Eddy covariance measures net exchange; GPP is
    partitioned from it by nighttime extrapolation (NT) or a daytime
    light-response fit (DT), and the two disagree. --gpp picks one; running both
    and reading their spread as the observational uncertainty is the honest use.

    ET IS LE. The tower has one latent-heat measurement, converted here with
    LE*86400/LAMBDA. The ET and LE rows therefore share an observation and differ
    only in which MODEL pathway they test -- ET from T&C's water balance, QE from
    its energy balance -- which makes the pair an internal consistency check.

    MONTHLY AND ANNUAL COME FROM THE ARCHIVE'S OWN MM AND YY FILES rather than
    from aggregating half-hourly data: ONEFlux applies its gap-filling and
    u*-filtering before aggregating, and re-deriving them here would silently
    use a different method than the GPP partitioning already assumes. SEASONS
    are built from MM, since no seasonal file exists -- a season needs all three
    of its months present or it is dropped.
    """
    res = {"monthly": "MM", "seasonal": "MM", "annual": "YY"}[freq]
    key = f"GPP_{gpp}"
    out, missing = {}, []
    monthly: dict = defaultdict(dict)

    for sid in sorted(sites):
        path = _fluxnet_file(root, sid, res)
        if path is None:
            missing.append(f"{sid}: no {res} file")
            continue
        with path.open(newline="", encoding="utf-8-sig") as fh:
            rd = csv.DictReader(fh)
            header = set(rd.fieldnames or [])
            tstamp = _pick(header, ["TIMESTAMP", "TIMESTAMP_START"])
            cols = {k: _pick(header, v) for k, v in TOWER_COLS.items()}
            if tstamp is None or cols.get(key) is None or cols["LE"] is None:
                missing.append(
                    f"{sid}: {path.name} has no usable "
                    f"{'timestamp' if tstamp is None else key if not cols.get(key) else 'LE'} "
                    f"column; saw {sorted(header)[:12]}")
                continue
            for r in rd:
                ts = (r.get(tstamp) or "").strip()
                if not ts:
                    continue
                try:
                    year = int(ts[:4])
                    month = int(ts[4:6]) if len(ts) >= 6 else None
                except ValueError:
                    continue

                def val(name):
                    c = cols.get(name)
                    if not c:
                        return np.nan
                    try:
                        v = float(r[c])
                    except (TypeError, ValueError, KeyError):
                        return np.nan
                    if v in FLUXNET_NA or v <= -9990:
                        return np.nan
                    qc = QC_OF.get(name)
                    if qc and qc in header:
                        try:
                            if float(r[qc]) > max_qc:
                                return np.nan
                        except (TypeError, ValueError):
                            pass
                    return v

                g, le, h = val(key), val("LE"), val("H")
                rec = {"GPP": g, "QE": le, "H": h,
                       "ET": le * W_TO_MM_D if np.isfinite(le) else np.nan}
                if freq == "annual":
                    for var, v in rec.items():
                        if np.isfinite(v):
                            out[(sid, year, "ANN", var)] = v
                elif month:
                    monthly[(sid, year, month)] = rec

    if freq == "monthly":
        for (sid, year, month), rec in monthly.items():
            for var, v in rec.items():
                if np.isfinite(v):
                    out[(sid, year, month, var)] = v
    elif freq == "seasonal":
        # DJF is filed under the year of its January, matching the flux table.
        want = {"DJF": [(-1, 12), (0, 1), (0, 2)], "MAM": [(0, 3), (0, 4), (0, 5)],
                "JJA": [(0, 6), (0, 7), (0, 8)], "SON": [(0, 9), (0, 10), (0, 11)]}
        years = {y for _, y, _ in monthly}
        for sid in sorted(sites):
            for year in sorted(years):
                for season, months in want.items():
                    recs = [monthly.get((sid, year + dy, m)) for dy, m in months]
                    if any(r is None for r in recs):
                        continue          # a partial season is not a season
                    for var in ("GPP", "QE", "H", "ET"):
                        vs = [r[var] for r in recs if np.isfinite(r[var])]
                        if len(vs) == 3:
                            out[(sid, year, season, var)] = float(np.mean(vs))

    if not out:
        head = ["no tower data could be read."] + missing[:10]
        if len(missing) > 10:
            head.append(f"... and {len(missing) - 10} more")
        raise Missing("\n  ".join(head))
    if missing:
        print(f"  tower: {len(missing)} station(s) unusable, e.g. {missing[0]}")
    return out


# ------------------------------------------------------------------ skill score
def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))


def skill(model: dict, tower: dict, spei: dict, var: str, freq: str,
          how: str, threshold: float, min_n: int):
    """{station: (SS_all, SS_drought, n_all, n_drought)}."""
    per_st = defaultdict(list)
    for (sid, year, period, v), (f, d) in model.items():
        if v != var:
            continue
        obs = tower.get((sid, year, period, var))
        if obs is None or not np.isfinite(obs):
            continue
        s = spei.get((sid, year, period), np.nan)
        per_st[sid].append((f, d, obs, s))

    out = {}
    for sid, recs in per_st.items():
        f = np.array([r[0] for r in recs]); d = np.array([r[1] for r in recs])
        o = np.array([r[2] for r in recs]); s = np.array([r[3] for r in recs])
        dry = np.isfinite(s) & (s <= threshold)
        res = []
        for m in (np.ones(f.size, bool), dry):
            if m.sum() < min_n:
                res.append((np.nan, int(m.sum()))); continue
            rf, rd = rmse(f[m], o[m]), rmse(d[m], o[m])
            if how == "diff":
                ss = ((rf - rd) / (rf + rd)) if (rf + rd) > 0 else np.nan
            else:
                ss = (1.0 - rd / rf) if rf > 0 else np.nan
            res.append((ss, int(m.sum())))
        out[sid] = (res[0][0], res[1][0], res[0][1], res[1][1])
    return out


# ---------------------------------------------------------------------- drawing
def albers_xy(lat, lon):
    """CONUS Albers Equal Area (EPSG:5070). pyproj, so no cartopy dependency."""
    from pyproj import Transformer
    tr = Transformer.from_crs("EPSG:4326", "EPSG:5070", always_xy=True)
    x, y = tr.transform(np.asarray(lon, float), np.asarray(lat, float))
    return np.asarray(x), np.asarray(y)


def make_figure(table: dict, sites: dict, title: str, out_png: Path,
                basemap: Path | None, linthresh: float, vmax: float,
                figsize) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm, LinearSegmentedColormap
    from matplotlib.lines import Line2D

    # Colour-blind-safe diverging pair, white at the centre. Orange and purple
    # stay distinguishable under deuteranopia and protanopia, where a red-green
    # pair does not; this is the ColorBrewer PuOr anchor set.
    cmap = LinearSegmentedColormap.from_list(
        "ss_puor", ["#7f3b08", "#e08214", "#fee0b6", "#ffffff",
                    "#d8daeb", "#8073ac", "#2d004b"])
    norm = SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax, base=10)

    outline = None
    if basemap is not None and Path(basemap).exists():
        import geopandas as gpd
        g = gpd.read_file(basemap)
        outline = g.to_crs("EPSG:5070").dissolve()

    fig, axes = plt.subplots(len(ROWS), 2, figsize=figsize,
                             constrained_layout=True)
    for i, (label, mvar, unit) in enumerate(ROWS):
        for j, (col, which) in enumerate((("all steps", 0), ("drought only", 1))):
            ax = axes[i, j]
            if outline is not None:
                outline.boundary.plot(ax=ax, color="0.6", linewidth=0.4)
            got = table.get(mvar, {})
            for marker, kind in (("o", "deciduous"), ("^", "evergreen")):
                xs, ys, cs = [], [], []
                for sid, vals in got.items():
                    if sid not in sites:
                        continue
                    lat, lon, ft = sites[sid]
                    if ft != kind:
                        continue
                    ss = vals[which]
                    if not np.isfinite(ss):
                        continue
                    x, y = albers_xy(lat, lon)
                    xs.append(float(x)); ys.append(float(y)); cs.append(ss)
                if xs:
                    ax.scatter(xs, ys, c=cs, cmap=cmap, norm=norm, marker=marker,
                               s=34, edgecolor="0.25", linewidth=0.35, zorder=3)
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_linewidth(0.4); sp.set_color("0.6")
            if i == 0:
                ax.set_title(col, fontsize=9)
            if j == 0:
                ax.set_ylabel(f"{label}\n({unit})", fontsize=8)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=axes, orientation="horizontal",
                      fraction=0.04, pad=0.01, extend="both")
    cb.set_label("skill score   1 - RMSE$_{dyn}$/RMSE$_{fixed}$   "
                 "(> 0: dynamic LMA closer to tower)", fontsize=8)
    # Legend under the colour bar, not "upper right": at 6 inches wide it
    # collided with the right-hand column title.
    handles = [Line2D([], [], marker="o", ls="", color="0.35", label="deciduous"),
               Line2D([], [], marker="^", ls="", color="0.35", label="evergreen")]
    cb.ax.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
                 fontsize=8, bbox_to_anchor=(0.5, -1.6))
    fig.suptitle(title, fontsize=10)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model-dir", type=Path, required=True,
                    help="directory holding era5_{annual,monthly,seasonal}.csv")
    ap.add_argument("--tower-dir", type=Path, required=True,
                    help="FLUXNET archives from submit_fluxnet_download.sh")
    ap.add_argument("--step", default="all",
                    choices=["all", "annual", "monthly"] + SEASON_FIGS)
    ap.add_argument("--skill", default="murphy", choices=["murphy", "diff"])
    ap.add_argument("--gpp", default="NT", choices=["NT", "DT"])
    ap.add_argument("--threshold", type=float, default=-1.0)
    ap.add_argument("--min-n", type=int, default=6,
                    help="fewest overlapping steps an RMSE may rest on")
    ap.add_argument("--basemap", type=Path, default=None,
                    help="us_eco_l3 shapefile for the CONUS outline")
    ap.add_argument("--linthresh", type=float, default=0.05)
    ap.add_argument("--vmax", type=float, default=1.0)
    ap.add_argument("--figsize", default="6x10",
                    help="WxH inches; portrait by default")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)

    w, h = (float(v) for v in a.figsize.lower().split("x"))
    try:
        # $TC_FIGURES, a sibling of model_run -- figures are their own product,
        # regenerated freely and never wanted in the repo.
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
        model_csv = a.model_dir / f"era5_{freq}.csv"
        try:
            model = read_model(model_csv, freq)
            if keep is not None:
                model = {k: v for k, v in model.items() if k[2] in keep}
            spei = DL.station_spei({s: (v[0], v[1]) for s, v in sites.items()}, freq)
            tower = read_tower(a.tower_dir, freq, sites, a.gpp)
        except (Missing, DL.NoLabel) as e:
            print(f"ERROR [{name}]: {e}", file=sys.stderr)
            return 1

        table = {mvar: skill(model, tower, spei, mvar, freq, a.skill,
                             a.threshold, a.min_n)
                 for _, mvar, _ in ROWS}
        png = out_dir / f"skill_{name}.png"
        make_figure(table, sites, f"Skill score, {name} steps", png,
                    a.basemap, a.linthresh, a.vmax, (w, h))
        print(f"  -> {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
