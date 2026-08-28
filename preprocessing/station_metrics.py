"""Per-station fixed-vs-dynamic metrics and climate sensitivities.

Two tidy tables, both keyed per station AND per GCM -- the five members are
never averaged together here, so any pooling downstream is visible in the code
that does it. That matters: aggregating silently once made the SSP5-8.5 LAI mean
shift read -4.75% instead of -2.98%, because it counted one station's five
members as five stations.

  station_metrics.csv
    dataset,gcm,scenario,station,pft,variable,freq,period,subset,
    n,sd_fixed,sd_dyn,sd_ratio,mean_rel_pct,median_rel_pct,
    mean_fixed,mean_dyn,mean_shift_pct

    freq=annual   period=ANN            interannual variability            (1)
    freq=monthly  period=ALL            SD over every monthly value        (2)
    freq=monthly  period=1..12          that month's series across years   (9)
    freq=seasonal period=DJF|MAM|JJA|SON that season's series across years (3,10)
    sd_ratio = sd_dyn/sd_fixed                                             (7)
    mean_rel_pct = mean of (dyn-fixed)/|fixed| over the periods in the cell (8)
    subset in {all, drought, normal}                                      (11)

  station_sensitivity.csv
    dataset,gcm,scenario,station,pft,variable,freq,period,predictor,subset,n,
    slope_fixed,slope_dyn,delta_slope,se_fixed,se_dyn,
    r2_fixed,r2_dyn,p_fixed,p_dyn,
    slope_fixed_std,slope_dyn_std,delta_slope_std,corr_ta_spei

    predictor Ta       flux on that period's mean air temperature         (13)
    predictor SPEI12   annual: SPEI-12 at the water-year end              (14)
    predictor SPEI3    monthly: that month's SPEI-3
                       seasonal: SPEI-3 at the season's last month
    predictor LMA      on that year's LMA -- DYNAMIC ARM ONLY             (15)

    AT EVERY FREQUENCY. Each fit runs ACROSS YEARS within one period, so a
    monthly table holds twelve separate slopes per station and variable -- July
    GPP against July SPEI-3 across years -- and a seasonal table holds four.
    Pooling periods before fitting would average a July sensitivity with a
    January one and hide the seasonality the mechanism predicts.

DELTA_SLOPE IS THE RESULT; the per-arm slopes are context. Both arms see
identical forcing, so any moisture confounding inside the temperature slope --
the concern Jung et al. (2017) and Humphrey et al. (2018) raise about apparent
temperature sensitivity -- sits in both and cancels in the difference. The
absolute slopes carry that ambiguity; delta_slope does not. This mirrors the
C4MIP gamma framework, where the carbon-climate feedback is also read off paired
factorial runs rather than from one.

corr_ta_spei is carried as a diagnostic, not a correction: under SSP5-8.5 the
drought share doubles across the century (16.8% -> 34.5%) while temperature
rises monotonically, so the two axes are far less separable there than under
SSP1-2.6, where the SPEI trend is flat.

NO MINIMUM-n GUARD, BY REQUEST. A cell with one drought year yields NaN for the
SD and the slope; n is always written so it can be filtered afterwards. The
drought columns are therefore patchier than the rest and any fleet median over
them must apply n.

ITEMS 4-6 ARE NOT DUPLICATED HERE. The annual, monthly and seasonal timeseries
of (dyn-fixed)/fixed are already the rel_pct column of the tables this reads.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as _st

sys.path.insert(0, str(Path(__file__).resolve().parent))
from results_dir import NoResultsDir, resolve_out                 # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
SITE_LISTS = [REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
              REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv"]

DATASETS = ["era5", "historical", "ssp126", "ssp585"]
FREQS = ["annual", "monthly", "seasonal"]
# Forcing, carried as predictors rather than results: identical in both arms.
PREDICTOR_VARS = {"Ta", "Pr"}


def read_sites() -> pd.DataFrame:
    rows = []
    for p in SITE_LISTS:
        if p.exists():
            d = pd.read_csv(p, encoding="utf-8-sig")
            rows.append(d[["StationID", "ForestType"]].rename(
                columns={"StationID": "station", "ForestType": "pft"}))
    s = pd.concat(rows).drop_duplicates("station")
    s["pft"] = s["pft"].str.strip().str.lower()
    return s


def table_path(root: Path, ds: str, freq: str) -> Path:
    return root / (f"era5_{freq}.csv" if ds == "era5" else f"gcm_{freq}_{ds}.csv")


def load(root: Path, ds: str, freq: str):
    p = table_path(root, ds, freq)
    if not p.is_file():
        return None, f"{p.name}: not generated"
    d = pd.read_csv(p, usecols=lambda c: c in {
        "station", "key", "year", "period", "variable", "fixed", "dyn", "rel_pct"})
    for c in ("station", "key", "year", "period", "variable", "fixed", "dyn"):
        if c not in d.columns:
            return None, f"{p.name}: missing {c}"
    # key is "<scenario>/<gcm>:<arm>" for GCM runs, "era5_land:<arm>" otherwise.
    k = d["key"].astype(str)
    # "<scenario>/<gcm>:<arm>" for GCM runs, "era5_land:<arm>" otherwise. A
    # regex handles both: chaining .str.split("/").str[1].str.split(":") blew up
    # on the ERA5 tables, where no key contains "/" so the intermediate column
    # was all-NaN and had float dtype by the time the second .str ran.
    d["gcm"] = k.str.extract(r"^[^/]+/([^:]+):", expand=False).fillna("")
    d["scenario"] = ds
    d["freq"] = freq
    d = d.drop(columns=["key"])
    # rel_pct is (dyn-fixed)/|fixed| for that period; infinities come from a
    # fixed arm of exactly zero and are not numbers.
    if "rel_pct" not in d.columns:
        with np.errstate(divide="ignore", invalid="ignore"):
            d["rel_pct"] = 100 * (d["dyn"] - d["fixed"]) / d["fixed"].abs()
    d["rel_pct"] = pd.to_numeric(d["rel_pct"], errors="coerce")
    d["rel_pct"] = d["rel_pct"].replace([np.inf, -np.inf], np.nan)
    return d, None


def attach_drought(d: pd.DataFrame, lab: pd.DataFrame | None) -> pd.DataFrame:
    """Add class and spei, joining on the window that matches the step."""
    if lab is None:
        d["class"], d["spei"] = "unlabelled", np.nan
        return d
    L = lab[["gcm", "scenario", "station", "freq", "year", "period", "spei",
             "class"]].copy()
    for f in (d, L):
        f["period"] = f["period"].astype(str)
        f["year"] = pd.to_numeric(f["year"], errors="coerce").astype("Int64")
    # The ERA5 label table writes scenario as "era5_land"; load() sets it to
    # the dataset name, "era5". They never matched, so every ERA5 row joined
    # to nothing: spei all NaN (no SPEI sensitivity) and class all
    # "unlabelled" (empty drought and normal subsets). Normalise before the
    # join rather than after.
    L["scenario"] = L["scenario"].astype(str).replace({"era5_land": "era5"})
    keys = ["station", "freq", "year", "period", "scenario"]
    if d["gcm"].astype(str).str.len().gt(0).any():
        keys.append("gcm")
    else:
        L = L.drop(columns=["gcm"])
    out = d.merge(L, on=keys, how="left")
    # A join that matches NOTHING is a bug, not a dataset without droughts.
    # Silently filling "unlabelled" is what hid this for a whole run.
    hit = int(out["class"].notna().sum())
    if hit == 0 and len(L):
        print(f"  ERROR: drought labels matched 0 of {len(out)} rows on "
              f"{keys}; label scenarios={sorted(L['scenario'].unique())[:4]}, "
              f"data scenarios={sorted(d['scenario'].astype(str).unique())[:4]}",
              file=sys.stderr)
    elif hit < len(out):
        print(f"  labels matched {hit}/{len(out)} rows", flush=True)
    out["class"] = out["class"].fillna("unlabelled")
    return out


def _cell(g: pd.DataFrame) -> dict:
    f, v, r = g["fixed"].to_numpy(float), g["dyn"].to_numpy(float), \
              g["rel_pct"].to_numpy(float)
    n = int(np.isfinite(f).sum())
    sd_f = float(np.nanstd(f, ddof=1)) if n > 1 else np.nan
    sd_d = float(np.nanstd(v, ddof=1)) if n > 1 else np.nan
    mf, md = float(np.nanmean(f)), float(np.nanmean(v))
    return {"n": n, "sd_fixed": sd_f, "sd_dyn": sd_d,
            "sd_ratio": (sd_d / sd_f) if (sd_f and np.isfinite(sd_f) and sd_f > 0)
                        else np.nan,
            "mean_rel_pct": float(np.nanmean(r)) if np.isfinite(r).any() else np.nan,
            "median_rel_pct": float(np.nanmedian(r)) if np.isfinite(r).any() else np.nan,
            "mean_fixed": mf, "mean_dyn": md,
            "mean_shift_pct": (100 * (md - mf) / abs(mf)) if mf else np.nan}


def metrics(d: pd.DataFrame) -> pd.DataFrame:
    """One row per (gcm, station, variable, freq, period, subset)."""
    out = []
    base = ["gcm", "scenario", "station", "variable", "freq"]
    for subset in ("all", "drought", "normal"):
        s = d if subset == "all" else d[d["class"] == subset]
        if s.empty:
            continue
        # per period
        for keys, g in s.groupby(base + ["period"], observed=True, dropna=False):
            rec = dict(zip(base + ["period"], keys))
            out.append({**rec, "subset": subset, **_cell(g)})
        # monthly ALL: the SD over every monthly value in the series (item 2),
        # which is a different quantity from the SD of one month across years.
        if (s["freq"] == "monthly").any():
            m = s[s["freq"] == "monthly"]
            for keys, g in m.groupby(base, observed=True, dropna=False):
                rec = dict(zip(base, keys))
                out.append({**rec, "period": "ALL", "subset": subset, **_cell(g)})
    return pd.DataFrame(out)


def _ols(x: np.ndarray, y: np.ndarray) -> dict:
    """Slope with standard error, R2 and a two-sided p, or NaNs."""
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 3 or np.nanstd(x[m]) == 0:
        return {"n": n, "slope": np.nan, "se": np.nan, "r2": np.nan,
                "p": np.nan, "slope_std": np.nan}
    xx, yy = x[m], y[m]
    b, a = np.polyfit(xx, yy, 1)
    pred = a + b * xx
    ss_res = float(np.sum((yy - pred) ** 2))
    ss_tot = float(np.sum((yy - yy.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    sxx = float(np.sum((xx - xx.mean()) ** 2))
    se = float(np.sqrt(ss_res / (n - 2) / sxx)) if n > 2 and sxx > 0 else np.nan
    p = np.nan
    if np.isfinite(se) and se > 0:
        p = float(2 * _st.t.sf(abs(b / se), n - 2))
    sx, sy = np.std(xx, ddof=1), np.std(yy, ddof=1)
    return {"n": n, "slope": float(b), "se": se, "r2": r2, "p": p,
            "slope_std": float(b * sx / sy) if sy > 0 else np.nan}


def sensitivity(d: pd.DataFrame, lma: pd.DataFrame | None,
                freq: str = "annual") -> pd.DataFrame:
    """Flux regressed on Ta, SPEI and LMA, per arm, per subset, PER PERIOD.

    AT EVERY FREQUENCY, NOT ONLY ANNUAL. For monthly the regression is that
    month's flux against that month's SPEI-3, across years -- twelve separate
    fits per station and variable. For seasonal it is the season's flux against
    SPEI-3 at the season's last month, four fits. Annual is the flux against
    SPEI-12 at the water-year end, one fit.

    Each fit is therefore a genuine interannual sensitivity for that period, not
    a pooled slope across periods: a July slope and a January slope are
    different quantities, and averaging them before fitting would hide the
    seasonality the mechanism predicts.

    LMA is annual whatever the frequency -- one value per year -- so a monthly
    fit asks "does July GPP respond to that year's LMA", which is a meaningful
    question and the one the dynamic arm is built to answer.
    """
    out = []
    idx = ["gcm", "scenario", "station"]
    # Predictors live as their own 'variable' rows; pivot them out per station.
    # Predictors are per (station, year, PERIOD): Ta for that period and the
    # SPEI whose window matches it.
    key_p = idx + ["year", "period"]
    pred = (d[d["variable"].isin(PREDICTOR_VARS)]
            .pivot_table(index=key_p, columns="variable",
                         values="fixed", aggfunc="mean").reset_index())
    spei = d[key_p + ["spei", "class"]].drop_duplicates(subset=key_p)
    pred = pred.merge(spei, on=key_p, how="left")
    if lma is not None:                     # LMA is annual at every frequency
        pred = pred.merge(lma, on=idx + ["year"], how="left")

    # MERGE ONCE, THEN GROUP -- not a merge per group. The first version put
    # this join inside the loop, so historical monthly ran 71,400 merges against
    # a 153,000-row predictor table and job 39688 was still going after two
    # hours with the two SSP scenarios, each ~3x larger, still ahead of it.
    flux = d[~d["variable"].isin(PREDICTOR_VARS)]
    flux = flux.merge(pred, on=key_p, how="left", suffixes=("", "_p"))
    cls_all = (flux["class_p"] if "class_p" in flux.columns
               else flux.get("class"))
    flux = flux.assign(_cls=cls_all)
    for keys, j in flux.groupby(idx + ["variable", "period"], observed=True,
                                dropna=False):
        rec = dict(zip(idx + ["variable", "period"], keys))
        rec["freq"] = freq
        cls = j["_cls"]
        for subset in ("all", "drought", "normal"):
            s = j if subset == "all" else j[cls == subset]
            if s.empty:
                continue
            spei_name = "SPEI12" if freq == "annual" else "SPEI3"
            for name, col in (("Ta", "Ta"), (spei_name, "spei"), ("LMA", "LMA")):
                if col not in s.columns:
                    continue
                x = pd.to_numeric(s[col], errors="coerce").to_numpy(float)
                rf = _ols(x, s["fixed"].to_numpy(float))
                rd = _ols(x, s["dyn"].to_numpy(float))
                if name == "LMA":       # fixed arm has no LMA variation
                    rf = {k: np.nan for k in rf}
                    rf["n"] = rd["n"]
                ta = pd.to_numeric(s.get("Ta"), errors="coerce")
                sp = pd.to_numeric(s.get("spei"), errors="coerce")
                c = np.nan
                if ta is not None and sp is not None:
                    mm = ta.notna() & sp.notna()
                    if mm.sum() > 2:
                        c = float(np.corrcoef(ta[mm], sp[mm])[0, 1])
                out.append({**rec, "predictor": name, "subset": subset,
                            "n": rd["n"],
                            "slope_fixed": rf["slope"], "slope_dyn": rd["slope"],
                            "delta_slope": rd["slope"] - rf["slope"]
                                           if np.isfinite(rf["slope"]) else np.nan,
                            "se_fixed": rf["se"], "se_dyn": rd["se"],
                            "r2_fixed": rf["r2"], "r2_dyn": rd["r2"],
                            "p_fixed": rf["p"], "p_dyn": rd["p"],
                            "slope_fixed_std": rf["slope_std"],
                            "slope_dyn_std": rd["slope_std"],
                            "delta_slope_std": rd["slope_std"] - rf["slope_std"]
                                               if np.isfinite(rf["slope_std"]) else np.nan,
                            "corr_ta_spei": c})
    return pd.DataFrame(out)


def read_lma(root: Path, ds: str) -> pd.DataFrame | None:
    """Yearly LMA per station from the lma_effect annual tables."""
    tag = "era5_land_ic" if ds == "era5" else ds
    p = root / f"lma_effect_{tag}_annual.csv"
    if not p.is_file():
        return None
    d = pd.read_csv(p)
    cols = {c.lower(): c for c in d.columns}
    if "lma" not in cols or "year" not in cols:
        return None
    sid = cols.get("stationid") or cols.get("station")
    out = d[[sid, cols["year"], cols["lma"]]].rename(
        columns={sid: "station", cols["year"]: "year", cols["lma"]: "LMA"})
    out["scenario"] = ds
    # USE THE MODEL COLUMN WHEN IT IS THERE. This used to select only
    # station/year/LMA, deduplicate on (station, year) -- silently keeping
    # whichever of the five GCMs happened to come first and discarding the
    # other four -- and then hardcode gcm="". That empty string can never
    # match a GCM run's real model name, so the join in sensitivity() found
    # nothing and every GCM LMA regression was skipped without a word.
    if "gcm" in cols:
        out["gcm"] = d[cols["gcm"]].astype(str).fillna("")
        keys = ["gcm", "station", "year"]
    elif ds != "era5":
        # HARD STOP, not a note. A GCM dataset whose LMA table has no model
        # column cannot produce per-model LMA rows, and the run would take 35
        # minutes to finish looking healthy while silently reproducing the
        # previous output. That happened three times. Fail in one second
        # instead, and say exactly which command regenerates the file --
        # analyze_lma_effect.py has an extract stage and a report stage, and
        # only the REPORT stage writes this table.
        raise SystemExit(
            f"ERROR: {p.name} has no gcm column, so LMA cannot be attributed "
            f"per model for {ds}.
"
            f"       Regenerate it first:
"
            f"         sbatch -p SOE_legacy -A efthymios "
            f"slurm/submit_lma_effect.sh --report
"
            f"       then confirm:  head -1 {p}
"
            f"       Re-running station_metrics before that reproduces the "
            f"previous output exactly.")
    else:
        out["gcm"] = ""
        keys = ["station", "year"]
    n0 = len(out)
    out = out.drop_duplicates(keys)
    print(f"  lma {ds:<11}{n0:>8} rows -> {len(out):>8} unique on {keys}, "
          f"{out['gcm'].nunique()} model(s)", flush=True)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--out-prefix", default="station")
    ap.add_argument("--datasets", default=",".join(DATASETS))
    ap.add_argument("--update", action="store_true",
                    help="keep rows of OTHER datasets already in the output "
                         "and replace only the ones recomputed here")
    a = ap.parse_args(argv)
    try:
        root = a.results or resolve_out(".", create=False)
        out_m = resolve_out(f"{a.out_prefix}_metrics.csv")
        out_s = resolve_out(f"{a.out_prefix}_sensitivity.csv")
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    sites = read_sites()
    labs = {}
    for src, f in (("era5", "drought_periods_era5.csv"),
                   ("gcm", "drought_periods_gcm.csv")):
        p = root / f
        labs[src] = pd.read_csv(p) if p.is_file() else None
        print(f"labels {src:<5}: {'ok' if labs[src] is not None else 'MISSING ' + f}")

    all_m, all_s, skipped = [], [], []
    for ds in [x.strip() for x in a.datasets.split(",")]:
        lab = labs["era5"] if ds == "era5" else labs["gcm"]
        per_freq = {}
        for freq in FREQS:
            d, why = load(root, ds, freq)
            if d is None:
                skipped.append(why); continue
            d = attach_drought(d, lab)
            m = metrics(d)
            m["dataset"] = ds
            all_m.append(m)
            per_freq[freq] = d
            print(f"  {ds:<11}{freq:<9}{len(d):>9} rows -> {len(m):>7} metric rows",
                  flush=True)
        for freq, dd in per_freq.items():
            s = sensitivity(dd, read_lma(root, ds), freq)
            if s.empty:
                continue
            s["dataset"] = ds
            all_s.append(s)
            print(f"  {ds:<11}{'sens ' + freq:<14}{len(s):>7} regression rows",
                  flush=True)

    if not all_m:
        print("ERROR: nothing computed", file=sys.stderr)
        for w in skipped:
            print(f"  {w}", file=sys.stderr)
        return 1

    done = [x.strip() for x in a.datasets.split(",")]

    def splice(new_rows: pd.DataFrame, path: Path) -> pd.DataFrame:
        """Replace only the recomputed datasets, keep the rest as they were.

        Recomputing one dataset costs a fraction of the full run, but the
        writer overwrote the whole file, so a targeted rerun would silently
        delete every other dataset's rows. --update keeps them.
        """
        if not a.update or not path.is_file():
            return new_rows
        old_rows = pd.read_csv(path, low_memory=False)
        if "dataset" not in old_rows.columns:
            print(f"  WARNING: {path.name} has no dataset column; cannot "
                  f"splice, writing only {done}", file=sys.stderr)
            return new_rows
        keep = old_rows[~old_rows["dataset"].isin(done)]
        print(f"  {path.name}: kept {len(keep)} rows from "
              f"{sorted(set(keep['dataset']))}, replaced "
              f"{len(old_rows) - len(keep)} with {len(new_rows)} for {done}")
        return pd.concat([keep, new_rows], ignore_index=True)

    M = pd.concat(all_m).merge(sites, on="station", how="left")
    M = splice(M, out_m)
    M.to_csv(out_m, index=False)
    print(f"\n-> {out_m}  ({len(M)} rows)")
    if all_s:
        S = pd.concat(all_s).merge(sites, on="station", how="left")
        S = splice(S, out_s)
        S.to_csv(out_s, index=False)
        print(f"-> {out_s}  ({len(S)} rows)")
    if skipped:
        print(f"\nSKIPPED {len(skipped)}:")
        for w in skipped:
            print(f"  ! {w}")
    return 1 if skipped else 0


if __name__ == "__main__":
    sys.exit(main())
