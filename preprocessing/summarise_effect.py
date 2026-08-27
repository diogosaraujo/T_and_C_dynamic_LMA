"""Read every effect table and say what the experiment found.

Walks the daily, monthly, seasonal and annual tables in $TC_RESULTS -- ERA5 and
all three GCM scenarios -- and writes one markdown report plus a tidy CSV of the
per-station numbers behind it.

WHAT IT REPORTS, and why each is the defensible form of the question:

  SIZE        median across stations of |rel_ann_pct|. Median, not mean: with
              ~92 stations one pathological site turned a fleet doing ~1.2% into
              a printed 129.54% earlier in this project. rel_ann, not rel_pct:
              rel_pct divides by the same period's own value, which vanishes out
              of season and produced a "50465967%" leakage figure that was a
              February denominator rather than an effect.

  DIRECTION   the signed median, plus the share of stations whose own median has
              that sign. A 3% effect that is +3% at half the stations and -3% at
              the other half is not the same finding as one that is +3%
              everywhere, and |rel| alone cannot tell them apart.

  SEASONALITY which period carries the largest effect. The mechanism is
              LMA -> SLA -> leaf area, so the effect should concentrate around
              leaf-out and peak season; a peak in midwinter is the signature of
              a vanishing denominator, not of ecology.

  DROUGHT     effect in drought periods against normal ones, joined from
              drought_periods_*.csv on (station, freq, year, period) so each
              step carries the accumulation that matches it: SPEI-12 at the
              water-year end for a year, SPEI-3 for a month, SPEI-3 at the
              season's last month for a season. The same definition the figures
              use. Steps the index cannot label -- the first N-1 months of an
              N-month accumulation -- are "unlabelled", never folded into
              "normal".

  FOREST TYPE deciduous against evergreen, because the LMA retrieval is noisier
              for evergreen (25 of 53 evergreen stations meet both noise
              criteria against 8 of 39 deciduous), so an evergreen-driven result
              deserves more caution than a deciduous one.

STALE TABLES ARE REFUSED, NOT SUMMARISED. A file without rel_ann_pct predates
the fix that made the relative measure meaningful, and summarising it would
produce numbers that look fine and are not. Such files are named and skipped.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from results_dir import NoResultsDir, resolve_out                 # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
SITE_LISTS = [REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
              REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv"]

# The fluxes worth a headline. Ratios are deliberately excluded from the size
# table: they are dimensionless and their rel_ann is a different kind of number
# from a flux's, so putting them in one ranked column would invite a comparison
# that does not mean anything.
FLUXES = ["LAI_H", "GPP", "NPP", "ET", "T", "EG", "EIn", "Lk", "QE", "H", "Rn"]
RATIOS = ["Tfrac", "Bowen", "WUE"]

DATASETS = [("era5", "ERA5-Land"), ("historical", "GCM historical"),
            ("ssp126", "SSP1-2.6"), ("ssp585", "SSP5-8.5")]
# Daily is deliberately NOT summarised. Those tables are day-of-year
# CLIMATOLOGIES: every value is already a multi-year mean, so a "drought" class
# built from ~4.8 years cancels less noise than an "all" class built from ~36
# and the contrast is not like for like. Monthly, seasonal and annual keep one
# row per period per year, which is what these summaries need.
FREQS = ["monthly", "seasonal", "annual"]

SEASON_ORDER = {"DJF": 1, "MAM": 2, "JJA": 3, "SON": 4}

# A period whose fixed-arm value is below this fraction of the station's record
# mean is excluded from the relative measure: the ratio there is not a change,
# it is a division by nothing.
TINY_DEN_FRAC = 0.01


def table_path(root: Path, ds: str, freq: str) -> Path:
    return root / (f"era5_{freq}.csv" if ds == "era5"
                   else f"gcm_{freq}_{ds}.csv")


def read_sites() -> pd.DataFrame:
    rows = []
    for p in SITE_LISTS:
        if p.exists():
            d = pd.read_csv(p, encoding="utf-8-sig")
            if {"StationID", "ForestType"} <= set(d.columns):
                rows.append(d[["StationID", "ForestType", "US_L3NAME"]]
                            .rename(columns={"StationID": "station",
                                             "ForestType": "pft",
                                             "US_L3NAME": "ecoregion"}))
    if not rows:
        raise SystemExit("no site list found")
    s = pd.concat(rows).drop_duplicates("station")
    s["pft"] = s["pft"].str.strip().str.lower()
    return s


def read_drought(root: Path, ds: str) -> pd.DataFrame | None:
    """Per-period drought labels for one dataset, from drought_periods_*.csv.

    REPLACES drought_years*.csv, which held the ANNUAL MEAN of monthly SPEI-12
    -- an average over twelve overlapping 12-month windows, which is not a
    quantity anyone wants and is not the definition the figures use. The new
    table labels each step with the accumulation that matches it: SPEI-12 at the
    water-year end for a year, SPEI-3 for a month, SPEI-3 at the season's last
    month for a season. Verified at ACCESS-CM2/US-Bar/2000, where September
    SPEI-12 is -0.7268 against an annual mean of -0.6143.

    Returned keyed on (station, freq, year, period) -- and on gcm too for the
    GCM runs, because a model's dry years are its own.
    """
    p = root / ("drought_periods_era5.csv" if ds == "era5"
                else "drought_periods_gcm.csv")
    if not p.is_file():
        return None
    d = pd.read_csv(p)
    need = {"station", "freq", "year", "period", "class"}
    if not need <= set(d.columns):
        return None
    if ds == "era5":
        keep = ["station", "freq", "year", "period", "class"]
    else:
        d = d[d.get("scenario", "").astype(str) == ds]
        keep = ["station", "gcm", "freq", "year", "period", "class"]
    d = d[keep].copy()
    d["period"] = d["period"].astype(str)
    d["year"] = pd.to_numeric(d["year"], errors="coerce").astype("Int64")
    return d.drop_duplicates()


def load(path: Path, freq: str) -> tuple[pd.DataFrame | None, str]:
    """The table as a tidy frame, or (None, reason)."""
    if not path.is_file():
        return None, "not generated"
    head = pd.read_csv(path, nrows=0)
    cols = set(head.columns)
    if "rel_ann_pct" not in cols:
        return None, ("STALE: no rel_ann_pct column -- predates the fix that "
                      "made the relative measure meaningful; regenerate it")
    period_col = "doy" if freq == "daily" else "period"
    need = {"station", "key", "variable", "diff", "rel_pct", "rel_ann_pct",
            "fixed", "dyn", period_col}
    miss = need - cols
    if miss:
        return None, f"missing column(s): {', '.join(sorted(miss))}"

    use = sorted(need | ({"class"} if "class" in cols else set())
                 | ({"year"} if "year" in cols else set()))
    df = pd.read_csv(path, usecols=use,
                     dtype={"station": "category", "variable": "category",
                            "key": "category"})
    df = df.rename(columns={period_col: "period"})
    df["freq"] = freq                      # the join matches window to step
    # THE DENOMINATOR IS THE FIXED ARM'S OWN VALUE FOR THAT PERIOD (rel_pct),
    # which is what "how much did dynamic LMA change this flux" actually means.
    # rel_ann divides by the record mean instead, and so understates the change
    # in a low period and overstates it in a high one -- it answers a different
    # question and is kept only as a secondary column.
    #
    # The cost is a denominator that can approach zero: 14.8% of station-months
    # have monthly GPP below 1% of that station's record mean, almost all
    # deciduous winters, and the untamed tail reaches 279037% for GPP and 17.7
    # million % for EG. Those rows are excluded and counted. The MEDIANS this
    # script reports would survive them anyway -- a median ignores a tail -- but
    # a number that cannot be defended should not be in the pool at all.
    scale = (df.groupby(["station", "variable"], observed=True)["fixed"]
               .transform(lambda x: np.mean(np.abs(x))))
    df["tiny_den"] = df["fixed"].abs() < TINY_DEN_FRAC * scale
    df["rel"] = df["rel_pct"].where(~df["tiny_den"])
    df["rel"] = df["rel"].replace([np.inf, -np.inf], np.nan)
    if "class" not in df.columns:
        df["class"] = "all"
    # The GCM key is "<scenario>/<gcm>:<arm>"; the drought labels are per GCM.
    if df["key"].astype(str).str.contains("/").any():
        df["gcm"] = (df["key"].astype(str).str.split("/").str[1]
                     .str.split(":").str[0])
    return df, "ok"


def per_station(df: pd.DataFrame) -> pd.DataFrame:
    """One row per (station, variable, class): the station's own medians.

    Two stages on purpose. Collapsing every row at once would let a station with
    more periods or more years weigh more than one with fewer, and the fleet
    statement is about stations, not about rows.
    """
    # TWO WAYS AN "all" CLASS ARISES, and conflating them double-counts.
    # The daily tables store a SEPARATE climatology per class, so their "all"
    # rows already cover the whole record and the drought rows are additional
    # rows over the same days -- relabelling everything "all" would count those
    # days twice. The monthly/seasonal/annual tables have one row per period,
    # labelled drought or normal by the join, and no "all" at all; there, "all"
    # has to be synthesised by taking every row.
    cls = set(df["class"].astype(str))
    if "all" in cls:
        stacked = df
    else:
        base = df.copy()
        base["class"] = "all"
        extra = df[df["class"].astype(str).isin(("drought", "normal"))]
        stacked = pd.concat([base, extra], ignore_index=True) if len(extra) else base
    g = stacked.groupby(["station", "variable", "class"], observed=True)
    out = g.agg(abs_rel=("rel", lambda s: np.nanmedian(np.abs(s))),
                signed_rel=("rel", "median"),
                abs_rel_ann=("rel_ann_pct", lambda s: np.nanmedian(np.abs(s))),
                signed_diff=("diff", "median"),
                mean_fixed=("fixed", "mean"),
                mean_dyn=("dyn", "mean"),
                dropped=("tiny_den", "mean"),
                n=("rel", "size")).reset_index()
    # The MEAN SHIFT: how far the dynamic arm's long-run mean sits from the
    # fixed arm's. Distinct from the median of the per-period changes, and the
    # distinction matters -- historical and ERA5 share a mean by construction
    # while the SSP arms do not, because the projected LMA trends away from the
    # historical baseline it was referenced to.
    with np.errstate(divide="ignore", invalid="ignore"):
        out["mean_shift_pct"] = np.where(
            out["mean_fixed"].abs() > 0,
            100 * (out["mean_dyn"] - out["mean_fixed"]) / out["mean_fixed"].abs(),
            np.nan)
    return out.dropna(subset=["abs_rel"])


def fleet(ps: pd.DataFrame, cls: str = "all") -> pd.DataFrame:
    """Across stations: median size, direction, and how many stations agree."""
    d = ps[ps["class"] == cls]
    rows = []
    for var, grp in d.groupby("variable", observed=True):
        sgn = grp["signed_rel"].dropna()
        if sgn.empty:
            continue
        pos = float((sgn > 0).mean() * 100)
        rows.append({"variable": str(var),
                     "median_abs_pct": float(np.nanmedian(grp["abs_rel"])),
                     "median_signed_pct": float(np.nanmedian(sgn)),
                     "stations": int(len(grp)),
                     "pct_stations_positive": pos,
                     "agreement_pct": max(pos, 100 - pos)})
    return pd.DataFrame(rows)


def seasonality(df: pd.DataFrame, freq: str) -> pd.DataFrame:
    """Which period carries the biggest effect, per variable."""
    # Same trap as per_station: after the drought join every row is
    # drought/normal and there is no "all" class, which left the monthly
    # seasonality section silently EMPTY in the first run.
    d = df[df["class"].astype(str) == "all"] if         "all" in set(df["class"].astype(str)) else df
    g = (d.groupby(["variable", "period"], observed=True)["rel"]
          .apply(lambda s: np.nanmedian(np.abs(s))).reset_index(name="abs_rel"))
    if g.empty:
        return g
    idx = g.groupby("variable", observed=True)["abs_rel"].idxmax()
    return g.loc[idx].reset_index(drop=True)


def fmt_table(df: pd.DataFrame, cols: list, heads: list, prec=2) -> list[str]:
    L = ["| " + " | ".join(heads) + " |",
         "|" + "|".join(["---"] * len(heads)) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, bool) or not isinstance(v, (int, float, np.number)):
                cells.append(str(v))
            elif isinstance(v, (int, np.integer)):
                cells.append(str(int(v)))       # station COUNTS are not 5.00
            else:
                cells.append(f"{v:.{prec}f}")
        L.append("| " + " | ".join(cells) + " |")
    return L


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None,
                    help="directory holding the effect tables (default $TC_RESULTS)")
    ap.add_argument("--out", type=Path, default=None,
                    help="markdown report; a bare name lands in $TC_RESULTS")
    ap.add_argument("--freqs", default=",".join(FREQS))
    a = ap.parse_args(argv)

    try:
        root = a.results or resolve_out(".", create=False)
        out_md = resolve_out(a.out or "effect_summary.md")
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    sites = read_sites()
    freqs = [f.strip() for f in a.freqs.split(",") if f.strip()]
    L = ["# Fixed vs dynamic LMA — summary of effect tables", "",
         "Percentages are **`100 x (dyn - fixed) / |fixed|` for that same "
         "period**, aggregated as a median over periods within a station and "
         "then a median across stations. Medians throughout, because one "
         "pathological station otherwise sets the whole number.", "",
         f"Periods whose fixed-arm value falls below {TINY_DEN_FRAC:.0%} of the "
         "station's record mean are excluded: there the ratio is a division by "
         "nothing rather than a change. The share dropped is reported per "
         "table. Day-of-year climatologies are not summarised here at all -- "
         "see the note in the source.", ""]
    tidy, skipped, found = [], [], []

    for freq in freqs:
        blocks = []
        for ds, label in DATASETS:
            p = table_path(root, ds, freq)
            df, why = load(p, freq)
            if df is None:
                skipped.append((p.name, why))
                continue
            found.append(p.name)
            dro = read_drought(root, ds)
            if dro is not None and "year" in df.columns:
                # Join on the PERIOD as well as the year: a July step takes
                # July's SPEI-3, not the whole year's label. Joining on year
                # alone was what let the coarser annual definition leak into
                # the monthly and seasonal sections.
                keys = ["station", "freq", "year", "period"]
                if "gcm" in dro.columns and "gcm" in df.columns:
                    keys.append("gcm")
                elif "gcm" in dro.columns:
                    dro = dro.drop(columns=["gcm"]).drop_duplicates()
                df = df.drop(columns=["class"])
                df["period"] = df["period"].astype(str)
                df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
                df = df.merge(dro, on=keys, how="left")
                # Unlabelled is NOT normal: the first N-1 months of an N-month
                # index have no value, so those steps must not silently join
                # the non-drought pool.
                df["class"] = df["class"].fillna("unlabelled")
            ps = per_station(df)
            ps.insert(0, "freq", freq); ps.insert(1, "dataset", ds)
            ps = ps.merge(sites, on="station", how="left")
            tidy.append(ps)
            blocks.append((label, df, ps))

        if not blocks:
            continue
        L += [f"## {freq.capitalize()}", ""]

        # -- size and direction, fluxes only
        L += ["### Effect size and direction (all steps)", ""]
        head = ["variable"] + [lab for lab, _, _ in blocks]
        L += ["| " + " | ".join(head) + " |",
              "|" + "|".join(["---"] * len(head)) + "|"]
        per_ds = {lab: fleet(ps).set_index("variable") for lab, _, ps in blocks}
        for var in FLUXES + RATIOS:
            cells = []
            for lab, _, _ in blocks:
                t = per_ds[lab]
                if var in t.index:
                    r = t.loc[var]
                    cells.append(f"{r['median_abs_pct']:.2f}% "
                                 f"({r['median_signed_pct']:+.2f}, "
                                 f"{r['agreement_pct']:.0f}% agree)")
                else:
                    cells.append("—")
            if any(c != "—" for c in cells):
                L.append(f"| {var} | " + " | ".join(cells) + " |")
        L += ["", "Each cell: median |effect|, then the signed median and the "
                  "share of stations agreeing on that sign.", ""]

        # -- mean shift, the thing a median of per-period changes cannot show
        L += ["### Long-run mean shift (dynamic minus fixed)", "",
              "Median across stations of `100 x (mean(dyn) - mean(fixed)) / "
              "|mean(fixed)|`. Near zero means the two arms share a long-run "
              "mean and the treatment only redistributes; a non-zero value "
              "means dynamic LMA moves the mean itself.", ""]
        head2 = ["variable"] + [lab for lab, _, _ in blocks]
        L += ["| " + " | ".join(head2) + " |",
              "|" + "|".join(["---"] * len(head2)) + "|"]
        for var in FLUXES:
            cells = []
            for lab, _, ps in blocks:
                g = ps[(ps["class"] == "all") & (ps["variable"].astype(str) == var)]
                v = np.nanmedian(g["mean_shift_pct"]) if len(g) else np.nan
                cells.append("—" if not np.isfinite(v) else f"{v:+.2f}%")
            if any(c != "—" for c in cells):
                L.append(f"| {var} | " + " | ".join(cells) + " |")
        L += [""]

        # -- how much was excluded by the denominator guard
        drops = []
        for lab, _, ps in blocks:
            g = ps[ps["class"] == "all"]
            if len(g):
                drops.append(f"{lab} {100*np.nanmean(g['dropped']):.1f}%")
        if drops:
            L += ["Rows excluded by the small-denominator guard: "
                  + "; ".join(drops) + ".", ""]

        # -- seasonality
        L += ["### Where in the year the effect peaks", ""]
        for lab, df, _ in blocks:
            s = seasonality(df, freq)
            if s.empty:
                continue
            s = s[s["variable"].astype(str).isin(FLUXES)]
            s = s.sort_values("abs_rel", ascending=False).head(6)
            peaks = ", ".join(f"{r['variable']} at {r['period']} "
                              f"({r['abs_rel']:.2f}%)" for _, r in s.iterrows())
            L += [f"- **{lab}**: {peaks}"]
        L += [""]

        # -- drought contrast
        # DROUGHT AGAINST NORMAL, NOT AGAINST "ALL" -- and never at daily
        # resolution. The daily tables are CLIMATOLOGIES: the "all" class
        # averages ~36 years and the "drought" class ~4.8, and averaging fewer
        # years cancels less year-to-year noise, so |effect| comes out larger in
        # the drought composite whether or not drought has anything to do with
        # it. That inflation is severe -- the first version of this report put
        # GPP at a 11.9x drought ratio, while the same contrast computed
        # monthly, where no cross-year averaging happens, gives 1.18x.
        #
        # Monthly, seasonal and annual keep one row per period per year, so
        # drought and normal are built the same way and the comparison is fair.
        rows = []
        if freq == "daily":
            L += ["### Drought years", "",
                  "**Not reported at daily resolution.** These tables are "
                  "climatologies: the *all* class averages ~36 years and "
                  "*drought* ~4.8, and a shorter average cancels less "
                  "year-to-year noise, which inflates the drought effect for "
                  "reasons unrelated to drought. Read the monthly, seasonal and "
                  "annual sections instead, where every row is one period of "
                  "one year and the two classes are built alike.", ""]
        else:
            for lab, _, ps in blocks:
                cls = set(ps["class"].astype(str).unique())
                if not {"drought", "normal"} <= cls:
                    continue
                n_ = fleet(ps, "normal").set_index("variable")
                d_ = fleet(ps, "drought").set_index("variable")
                for var in FLUXES:
                    if var in n_.index and var in d_.index:
                        base = n_.loc[var, "median_abs_pct"]
                        rows.append({"dataset": lab, "variable": var,
                                     "normal": base,
                                     "drought": d_.loc[var, "median_abs_pct"],
                                     "ratio": (d_.loc[var, "median_abs_pct"] / base
                                               if base else np.nan),
                                     "n_dry": int(d_.loc[var, "stations"])})
        if rows:
            r = pd.DataFrame(rows).sort_values("ratio", ascending=False).head(12)
            L += ["### Drought years against normal years", "",
                  "Like for like: both classes are built from single "
                  "period-years, so the two are directly comparable. Ratio > 1 "
                  "means the treatment bites harder in dry years.", ""]
            L += fmt_table(r, ["dataset", "variable", "normal", "drought",
                               "ratio", "n_dry"],
                           ["dataset", "variable", "normal %", "drought %",
                            "ratio", "stations"])
            L += [""]

        # -- forest type
        L += ["### Deciduous against evergreen (all steps)", ""]
        rows = []
        for lab, _, ps in blocks:
            d = ps[(ps["class"] == "all") & ps["variable"].astype(str).isin(FLUXES)]
            for pft, g in d.groupby("pft", observed=True):
                if not isinstance(pft, str):
                    continue
                rows.append({"dataset": lab, "pft": pft,
                             "stations": int(g["station"].nunique()),
                             "median_abs_pct": float(np.nanmedian(g["abs_rel"]))})
        if rows:
            L += fmt_table(pd.DataFrame(rows),
                           ["dataset", "pft", "stations", "median_abs_pct"],
                           ["dataset", "forest type", "stations", "median abs effect %"])
        L += [""]

    if skipped:
        L += ["## Inputs not summarised", ""]
        for name, why in skipped:
            L.append(f"- `{name}` — {why}")
        L += [""]

    if not tidy:
        print("ERROR: no usable table found. Skipped:", file=sys.stderr)
        for name, why in skipped:
            print(f"  {name}: {why}", file=sys.stderr)
        return 1

    out_md.write_text("\n".join(L), encoding="utf-8")
    tidy_csv = out_md.with_name(out_md.stem + "_per_station.csv")
    pd.concat(tidy).to_csv(tidy_csv, index=False)
    print(f"read {len(found)} table(s); skipped {len(skipped)}")
    for name, why in skipped:
        print(f"  ! {name}: {why}")
    print(f"-> {out_md}")
    print(f"-> {tidy_csv}")
    return 1 if skipped else 0


if __name__ == "__main__":
    sys.exit(main())
