#!/usr/bin/env python3
"""The numbers behind every figure, as CSV.

Each figure computes its statistics inside the drawing code and only renders
them, so nothing is quotable without reading a value off a panel -- and the
dryness panels print b, R2 and p at 5.6 pt, which is fine for seeing a pattern
and useless for writing a sentence. This writes the same quantities to disk.

THE AGGREGATION RULES ARE THE FIGURES' OWN, not a second implementation.
Median across GCMs at the metric level, quartiles from the same percentiles the
error bars span, station-years for the flux metric and stations for everything
else. A table that disagreed with the figure it describes would be worse than
no table.

    table_flux.csv         dataset x variable x pft
    table_variability.csv  dataset x variable x pft
    table_absolute.csv     dataset x variable x pft x arm
    table_sensitivity.csv  dataset x variable x pft x predictor  (3 slopes)
    table_fit.csv          dataset x variable x pft x predictor x arm
    table_dryness.csv      dataset x variable x pft   (OLS b, R2, p, n)

Read-only over station_metrics.csv, station_sensitivity.csv and
station_dryness.csv. Nothing is recomputed from the model output.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from results_dir import NoResultsDir, resolve_out                 # noqa: E402
from station_metrics import DATASETS, load as load_effect, read_sites  # noqa: E402
from figure_dryness import ols                                    # noqa: E402

VARS = {"GPP": "GPP", "LAI_H": "LAI", "T": "TR",
        "ET": "ET", "QE": "LE", "H": "H"}
DISP = list(VARS.values())


def spread(g: pd.Series) -> dict:
    """n, median and the quartiles the error bars actually span."""
    v = pd.to_numeric(g, errors="coerce").to_numpy(float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {"n": 0}
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    return {"n": int(v.size), "median": med, "q1": q1, "q3": q3,
            "mean": float(v.mean()), "min": float(v.min()),
            "max": float(v.max()),
            "pct_negative": float(100.0 * (v < 0).mean())}


def by(d: pd.DataFrame, keys: list[str], col: str) -> pd.DataFrame:
    out = []
    for k, g in d.groupby(keys, observed=True, dropna=False):
        rec = dict(zip(keys, k if isinstance(k, tuple) else (k,)))
        out.append({**rec, **spread(g[col])})
    t = pd.DataFrame(out)
    if "variable" in t.columns:
        t["variable"] = t["variable"].map(VARS).fillna(t["variable"])
    return t


def flux_table(root: Path) -> pd.DataFrame:
    """Station-years, median across GCMs first -- figure_flux's rule."""
    frames = []
    for ds in DATASETS:
        d, why = load_effect(Path(root), ds, "annual")
        if d is None:
            print(f"  {ds}: {why}"); continue
        d = d[d["variable"].isin(VARS)].copy()
        d["dataset"] = ds
        frames.append(d[["dataset", "gcm", "station", "year", "variable",
                         "rel_pct"]])
    if not frames:
        return pd.DataFrame()
    d = pd.concat(frames, ignore_index=True)
    g = (d.groupby(["dataset", "station", "year", "variable"], as_index=False)
           ["rel_pct"].median()
           .merge(read_sites(), on="station", how="left"))
    return by(g, ["dataset", "variable", "pft"], "rel_pct")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--prefix", default="table_")
    a = ap.parse_args(argv)
    try:
        root = Path(a.results or resolve_out(".", create=False))
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    def put(t: pd.DataFrame, name: str) -> None:
        if t is None or t.empty:
            print(f"  SKIP {name}: nothing to write"); return
        p = resolve_out(f"{a.prefix}{name}.csv")
        t.to_csv(p, index=False)
        print(f"  -> {p}  ({len(t)} rows)")

    print("flux (station-years):")
    put(flux_table(root), "flux")

    M = root / "station_metrics.csv"
    if M.is_file():
        m = pd.read_csv(M, low_memory=False)
        m = m[(m["freq"] == "annual") & (m["subset"] == "all")
              & (m["variable"].isin(VARS))].copy()
        # Variability: 100*(sd_ratio-1), median over GCMs, one point per station.
        m["v"] = 100.0 * (pd.to_numeric(m["sd_ratio"], errors="coerce") - 1.0)
        v = (m.groupby(["dataset", "station", "pft", "variable"],
                       as_index=False)["v"].median())
        print("variability (stations):")
        put(by(v, ["dataset", "variable", "pft"], "v"), "variability")

        abs_rows = []
        for arm, col in (("fixed", "mean_fixed"), ("dyn", "mean_dyn")):
            if col not in m.columns:
                continue
            g = (m.assign(x=pd.to_numeric(m[col], errors="coerce"))
                   .groupby(["dataset", "station", "pft", "variable"],
                            as_index=False)["x"].median())
            t = by(g, ["dataset", "variable", "pft"], "x")
            t["arm"] = arm
            abs_rows.append(t)
        print("absolute fluxes (stations, native units):")
        put(pd.concat(abs_rows, ignore_index=True) if abs_rows else None,
            "absolute")
    else:
        print(f"  SKIP: {M.name} not found")

    S = root / "station_sensitivity.csv"
    if S.is_file():
        s = pd.read_csv(S, low_memory=False)
        s = s[(s["freq"] == "annual") & (s["subset"] == "all")
              & (s["variable"].isin(VARS))].copy()
        parts = []
        for col, tag in (("delta_slope", "delta"), ("slope_fixed", "fixed"),
                         ("slope_dyn", "dyn")):
            if col not in s.columns:
                continue
            g = (s.assign(x=pd.to_numeric(s[col], errors="coerce"))
                   .groupby(["dataset", "station", "pft", "variable",
                             "predictor"], as_index=False)["x"].median())
            t = by(g, ["dataset", "variable", "pft", "predictor"], "x")
            t["quantity"] = tag
            parts.append(t)
        print("sensitivity slopes (stations, flux per predictor unit):")
        put(pd.concat(parts, ignore_index=True) if parts else None,
            "sensitivity")

        # Goodness of fit: r keeps the sign, p is summarised by how often it
        # clears 0.05 -- a median p is hard to read and easy to misquote.
        fits = []
        for arm in ("fixed", "dyn"):
            need = [f"r2_{arm}", f"p_{arm}", f"slope_{arm}"]
            if any(c not in s.columns for c in need):
                continue
            w = s.copy()
            r2 = pd.to_numeric(w[f"r2_{arm}"], errors="coerce").clip(lower=0)
            sg = np.sign(pd.to_numeric(w[f"slope_{arm}"], errors="coerce"))
            w["r"] = sg * np.sqrt(r2)
            w["p"] = pd.to_numeric(w[f"p_{arm}"], errors="coerce")
            g = (w.groupby(["dataset", "station", "pft", "variable",
                            "predictor"], as_index=False)
                   .agg(r=("r", "median"), p=("p", "median")))
            t = by(g, ["dataset", "variable", "pft", "predictor"], "r")
            sig = (g.assign(ok=g["p"] < 0.05)
                     .groupby(["dataset", "variable", "pft", "predictor"],
                              as_index=False)["ok"].mean())
            sig["variable"] = sig["variable"].map(VARS).fillna(sig["variable"])
            sig["pct_p_lt_05"] = 100.0 * sig.pop("ok")
            t = t.merge(sig, on=["dataset", "variable", "pft", "predictor"],
                        how="left")
            t["arm"] = arm
            fits.append(t)
        print("goodness of fit (median Pearson r, share significant):")
        put(pd.concat(fits, ignore_index=True) if fits else None, "fit")
    else:
        print(f"  SKIP: {S.name} not found")

    D = root / "station_dryness.csv"
    if S.is_file() and D.is_file():
        try:
            from figure_dryness import load as load_dry
            d = load_dry(root)
        except Exception as e:                                   # noqa: BLE001
            print(f"  SKIP dryness: {type(e).__name__}: {e}")
            d = None
        if d is not None:
            rows = []
            for k, g in d.groupby(["dataset", "variable", "pft"],
                                  observed=True, dropna=False):
                r = ols(g["phi"].to_numpy(float), g["slope"].to_numpy(float))
                rows.append(dict(zip(["dataset", "variable", "pft"], k)) | r)
            t = pd.DataFrame(rows)
            t["variable"] = t["variable"].map(VARS).fillna(t["variable"])
            t = t.rename(columns={"b": "ols_slope", "a": "ols_intercept"})
            print("dryness (the OLS printed on each panel):")
            put(t, "dryness")
    else:
        print("  SKIP dryness: needs station_dryness.csv (run build_dryness.py)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
