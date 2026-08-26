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

  DROUGHT     effect in drought years against normal ones. Daily tables carry a
              class column already; the others join the drought CSVs on the
              year, so the label is the annual SPEI-12 classification applied to
              every period inside that year -- coarser than the per-period
              SPEI-3 the figures use, and labelled as such in the report.

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
FREQS = ["daily", "monthly", "seasonal", "annual"]

SEASON_ORDER = {"DJF": 1, "MAM": 2, "JJA": 3, "SON": 4}


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
    """{station, year} -> class, from whichever drought CSV covers this dataset."""
    p = root / ("drought_years.csv" if ds == "era5" else "drought_years_gcm.csv")
    if not p.is_file():
        return None
    d = pd.read_csv(p)
    if "class" not in d.columns or "year" not in d.columns:
        return None
    if ds == "era5":
        d = d[d.get("scenario", "era5_land").astype(str) == "era5_land"]
        return d[["station", "year", "class"]].drop_duplicates()
    d = d[d.get("scenario", "").astype(str) == ds]
    # Keyed by GCM as well: a model's dry years are its own.
    return d[["station", "gcm", "year", "class"]].drop_duplicates()


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
    need = {"station", "key", "variable", "diff", "rel_ann_pct", period_col}
    miss = need - cols
    if miss:
        return None, f"missing column(s): {', '.join(sorted(miss))}"

    use = sorted(need | ({"class"} if "class" in cols else set())
                 | ({"year"} if "year" in cols else set()))
    df = pd.read_csv(path, usecols=use,
                     dtype={"station": "category", "variable": "category",
                            "key": "category"})
    df = df.rename(columns={period_col: "period"})
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
    out = g.agg(abs_rel=("rel_ann_pct", lambda s: np.nanmedian(np.abs(s))),
                signed_rel=("rel_ann_pct", "median"),
                signed_diff=("diff", "median"),
                n=("rel_ann_pct", "size")).reset_index()
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
    g = (d.groupby(["variable", "period"], observed=True)["rel_ann_pct"]
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
         "Sizes are the **median across stations** of `rel_ann_pct`, the "
         "difference as a percentage of that variable's own mean magnitude over "
         "the record. Median because one pathological station otherwise sets the "
         "number; `rel_ann` because `rel_pct` divides by a denominator that "
         "vanishes out of season.", ""]
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
                keys = ["station", "year"] + (["gcm"] if "gcm" in dro.columns
                                              and "gcm" in df.columns else [])
                df = df.drop(columns=["class"]).merge(dro, on=keys, how="left")
                df["class"] = df["class"].fillna("all")
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
