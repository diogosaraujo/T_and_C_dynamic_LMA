#!/usr/bin/env python3
"""Is the collapse driven by LMA the PLSR never saw?

Both arms are driven by IDENTICAL forcing -- same meteorology, same soil,
same parameters. The only difference is Sl. So a forcing problem on its own
cannot produce a dynamic-only collapse; it would hit both arms and show up
as "both". The difference must enter through the SLA channel, and LMA is
derived from climate, so a forcing problem can only reach the canopy by
making LMA extreme.

That is the testable question, and it is the whole point of this script:

  EXTRAPOLATION ARTEFACT  collapse-year LMA lies outside anything the
                          station saw historically. Future SSP climate is
                          outside the observational envelope the PLSR was
                          trained on, so the prediction may simply be
                          running off the end of its calibration.
  REAL FEEDBACK           collapse-year LMA sits inside the historical
                          range, and the canopy fails anyway -- LMA up,
                          Sl down, LAI = Sl*B(1) down, less absorbed
                          light, less carbon to B(1), LAI down again.

Also reports the lead-lag. If LMA jumps in the years BEFORE onset, the
input is driving the collapse. If LAI falls first and LMA follows, the
feedback loop is running. These need different fixes and look the same in
a count.

No model reruns. Joins dieoff_events.csv against the annual LMA series
already on disk.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from results_dir import NoResultsDir, resolve_out                 # noqa: E402

HIST = ["era5", "historical"]
FUT = ["ssp126", "ssp585"]


def find_lma(root: Path) -> tuple[pd.DataFrame | None, list[str]]:
    """Yearly LMA per (dataset, gcm, station). Says what it looked for.

    read_lma() in station_metrics builds the name as
    lma_effect_<ds>_annual.csv, which is a different convention from the
    effect tables themselves (gcm_annual_<scenario>.csv). That mismatch is
    the standing suspect for the empty GCM LMA sensitivity, so this lists
    what actually exists instead of failing silently.
    """
    notes, frames = [], []
    cands = sorted(Path(root).glob("*lma*annual*.csv"))
    notes.append("files matching *lma*annual*.csv: "
                 + (", ".join(p.name for p in cands) if cands else "NONE"))
    for p in cands:
        try:
            d = pd.read_csv(p, low_memory=False)
        except Exception as e:                                   # noqa: BLE001
            notes.append(f"  {p.name}: unreadable ({type(e).__name__})")
            continue
        low = {c.lower(): c for c in d.columns}
        sid = low.get("station") or low.get("stationid")
        yr, lma = low.get("year"), low.get("lma")
        if not (sid and yr and lma):
            notes.append(f"  {p.name}: no station/year/LMA columns "
                         f"({list(d.columns)[:6]})")
            continue
        g = d[[sid, yr, lma]].rename(columns={sid: "station", yr: "year",
                                              lma: "LMA"})
        g["gcm"] = d[low["gcm"]].astype(str) if "gcm" in low else ""
        # Dataset from the filename, since that is the only place it lives.
        stem = p.stem.replace("lma_effect_", "").replace("_annual", "")
        g["dataset"] = ("era5" if stem.startswith("era5") else stem)
        frames.append(g)
        notes.append(f"  {p.name}: {len(g)} rows -> dataset={g['dataset'].iloc[0]}")
    if not frames:
        return None, notes
    L = pd.concat(frames, ignore_index=True)
    L["year"] = pd.to_numeric(L["year"], errors="coerce")
    L["LMA"] = pd.to_numeric(L["LMA"], errors="coerce")
    return L.dropna(subset=["year", "LMA"]), notes


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--events", default="dieoff_events.csv")
    ap.add_argument("--out", default="dieoff_lma_check.csv")
    ap.add_argument("--lead", type=int, default=5,
                    help="years before onset to inspect")
    a = ap.parse_args(argv)
    try:
        root = Path(a.results or resolve_out(".", create=False))
        out_p = resolve_out(a.out)
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    p = root / a.events
    if not p.is_file():
        print(f"ERROR: {p} not found -- run dieoff_summary.py first",
              file=sys.stderr)
        return 1
    E = pd.read_csv(p, low_memory=False)
    E = E[E["variable"].isin(["LAI_H", "GPP", "T"])].copy()
    E["year"] = pd.to_numeric(E["year"], errors="coerce")
    E["gcm"] = E["gcm"].fillna("").astype(str)

    L, notes = find_lma(root)
    for n in notes:
        print(n)
    if L is None:
        print("\nERROR: no usable LMA series found. Without it this question "
              "cannot be answered -- not guessing.", file=sys.stderr)
        return 1

    # Historical envelope per station: what the PLSR actually produced in a
    # climate it was calibrated for.
    h = L[L["dataset"].isin(HIST)]
    if h.empty:
        print("\nERROR: no historical LMA rows; cannot define an envelope",
              file=sys.stderr)
        return 1
    env = (h.groupby("station")["LMA"]
             .agg(hist_min="min", hist_max="max", hist_med="median",
                  hist_n="size").reset_index())
    print(f"\nhistorical envelope from {len(h)} rows, "
          f"{env['station'].nunique()} stations")

    f = L[L["dataset"].isin(FUT)].merge(env, on="station", how="left")
    if f.empty:
        print("\nERROR: no future LMA rows; the SSP LMA series is missing, "
              "so the extrapolation question is untestable", file=sys.stderr)
        return 1
    f["above_hist"] = f["LMA"] > f["hist_max"]
    f["below_hist"] = f["LMA"] < f["hist_min"]
    f["outside"] = f["above_hist"] | f["below_hist"]
    # How far outside, in units of the station's own historical spread.
    span = (f["hist_max"] - f["hist_min"]).replace(0, np.nan)
    f["excess"] = np.where(f["above_hist"], (f["LMA"] - f["hist_max"]) / span,
                  np.where(f["below_hist"], (f["hist_min"] - f["LMA"]) / span,
                           0.0))

    key = ["dataset", "gcm", "station", "year"]
    ev = E[E["arm"].isin(["dynamic_only", "both"])][key].drop_duplicates()
    ev["collapse"] = True
    j = f.merge(ev, on=key, how="left")
    j["collapse"] = j["collapse"].fillna(False)

    print("\nIS COLLAPSE-YEAR LMA OUTSIDE THE HISTORICAL RANGE?")
    t = (j.groupby(["dataset", "collapse"], observed=True)
           .agg(station_years=("LMA", "size"),
                pct_outside=("outside", lambda s: round(100 * s.mean(), 2)),
                pct_above=("above_hist", lambda s: round(100 * s.mean(), 2)),
                med_LMA=("LMA", "median"),
                med_excess_spans=("excess", "median"),
                max_excess_spans=("excess", "max")))
    print(t.to_string())
    print("\nRead: if pct_outside is far higher for collapse=True, the "
          "collapses sit where the PLSR is extrapolating.\nIf the two rows "
          "are similar, collapse-year LMA is ordinary and the feedback is "
          "real.")

    # Lead-lag: LMA in the years before the first collapse at each station.
    onset = (E[E["arm"].isin(["dynamic_only", "both"])]
             .groupby(["dataset", "gcm", "station"], observed=True)["year"]
             .min().reset_index().rename(columns={"year": "onset"}))
    q = f.merge(onset, on=["dataset", "gcm", "station"], how="inner")
    q["lag"] = q["year"] - q["onset"]
    w = q[q["lag"].between(-a.lead, 2)]
    if not w.empty:
        print(f"\nLMA RELATIVE TO ONSET (lag 0 = first collapse year), "
              f"median across stations:")
        print(w.groupby(["dataset", "lag"], observed=True)
               .agg(n=("LMA", "size"), med_LMA=("LMA", "median"),
                    pct_outside=("outside",
                                 lambda s: round(100 * s.mean(), 1))).to_string())
        print("\nRead: LMA rising in the years BEFORE lag 0 means the input "
              "drives the collapse.\nLMA flat until lag 0 means the feedback "
              "is running.")

    j.to_csv(out_p, index=False)
    print(f"\n-> {out_p}  ({len(j)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
