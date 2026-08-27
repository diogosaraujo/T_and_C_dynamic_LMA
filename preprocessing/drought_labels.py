"""Which SPEI value labels a given model time step, at each frequency.

The accumulation window has to match the step being labelled, otherwise a
"drought month" is being judged by a year of antecedent conditions and the
composite means something else entirely. So:

    monthly   -> SPEI-3 for THAT month
    seasonal  -> SPEI-3 for the LAST month of the season (DJF->Feb, MAM->May,
                 JJA->Aug, SON->Nov). A 3-month accumulation ending in the
                 season's final month covers exactly that season.
    annual    -> SPEI-12 for SEPTEMBER, the end of the water year falling in
                 that calendar year. A 12-month accumulation ending in September
                 spans Oct-Sep, so it carries the snowpack and the growing
                 season that actually drove that year's fluxes -- a December
                 SPEI-12 would instead be dominated by the autumn after them.

SEASONAL YEARS FOLLOW THE FLUX TABLE. analyze_period_effect files DJF under the
year of its January, so DJF 2002 is Dec 2001 + Jan/Feb 2002 and its label is
SPEI-3 at February 2002 -- the same convention on both sides, which is what
makes the join correct.

ACCUMULATION OFFSETS ARE REAL. An N-month index has no value for the first N-1
months, and the stacks are trimmed rather than padded: SPEI-3 starts 1980-03 and
SPEI-12 starts 1980-12. A missing key is therefore normal near the start of the
record and is reported as absent, never as zero.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from era5_predictors import DEFAULT_ERA5_ROOT, Era5Monthly      # noqa: E402

INDEX_FOR = {"monthly": "SPEI3_ts", "seasonal": "SPEI3_ts", "annual": "SPEI12_ts"}
SEASON_END_MONTH = {"DJF": 2, "MAM": 5, "JJA": 8, "SON": 11}
WATER_YEAR_END = 9          # September closes the Oct-Sep water year


class NoLabel(Exception):
    """No SPEI value exists for this step, and why."""


def key_for(freq: str, year: int, period) -> int:
    """The year*100+month key whose SPEI labels this (year, period) step."""
    if freq == "monthly":
        m = int(period)
        if not 1 <= m <= 12:
            raise NoLabel(f"month {period} out of range")
        return year * 100 + m
    if freq == "seasonal":
        p = str(period)
        if p not in SEASON_END_MONTH:
            raise NoLabel(f"unknown season {period!r}")
        return year * 100 + SEASON_END_MONTH[p]
    if freq == "annual":
        return year * 100 + WATER_YEAR_END
    raise NoLabel(f"unknown freq {freq!r}")


def station_spei(sites: dict, freq: str, root: Path = DEFAULT_ERA5_ROOT) -> dict:
    """{(station, year, period): spei} for every step the index can label.

    sites is {station_id: (lat, lon)}. Stations whose pixel is all-NaN are
    omitted; the caller reports them rather than treating them as non-drought.
    """
    index = INDEX_FOR.get(freq)
    if index is None:
        raise NoLabel(f"unknown freq {freq!r}")
    store = Era5Monthly(root)
    out: dict = {}
    try:
        for sid, (lat, lon) in sorted(sites.items()):
            ser = store.pixel_series(lat, lon)
            vals = np.asarray(ser["si"][index], dtype=float)
            keys = np.asarray(ser["si_time"][index], dtype=int)
            if not np.isfinite(vals).any():
                continue
            by_key = {int(k): float(v) for k, v in zip(keys, vals)
                      if np.isfinite(v)}
            for k, v in by_key.items():
                y, m = divmod(k, 100)
                if freq == "monthly":
                    out[(sid, y, m)] = v
                elif freq == "seasonal":
                    for p, em in SEASON_END_MONTH.items():
                        if m == em:
                            out[(sid, y, p)] = v
                elif m == WATER_YEAR_END:
                    out[(sid, y, "ANN")] = v
    finally:
        store.close()
    return out


def is_drought(spei: float, threshold: float = -1.0) -> bool:
    return np.isfinite(spei) and spei <= threshold


# ---------------------------------------------------------------- label table
def _rows(by_key: dict, freq: str, sid: str, tag: tuple, threshold: float):
    """One row per period the index can label, for one station."""
    out = []
    for k, v in by_key.items():
        y, m = divmod(int(k), 100)
        if freq == "monthly":
            periods = [m]
        elif freq == "seasonal":
            periods = [p for p, em in SEASON_END_MONTH.items() if m == em]
        else:
            periods = ["ANN"] if m == WATER_YEAR_END else []
        for p in periods:
            out.append(tag + (sid, freq, y, p, round(float(v), 4),
                              "drought" if v <= threshold else "normal"))
    return out


def build_table(sites: dict, source: str, threshold: float = -1.0,
                gcms=None, scenarios=None, era5_root: Path = DEFAULT_ERA5_ROOT,
                freqs=("annual", "monthly", "seasonal")):
    """Rows of (source, gcm, scenario, station, freq, year, period, spei, class).

    ONE TABLE, ONE DEFINITION. Everything downstream joins on this rather than
    re-deriving drought, so the figures and the metrics cannot drift apart --
    which they had: classify_drought was writing the ANNUAL MEAN of monthly
    SPEI-12, an average over twelve overlapping 12-month windows, while the
    figures used SPEI-12 at the water-year end. Only the latter is meant.
    """
    import classify_drought as CD
    rows, skipped = [], []
    need = {INDEX_FOR[f] for f in freqs}

    if source == "era5":
        store = Era5Monthly(era5_root)
        try:
            for sid, (lat, lon) in sorted(sites.items()):
                ser = store.pixel_series(lat, lon)
                for freq in freqs:
                    idx = INDEX_FOR[freq]
                    vals = np.asarray(ser["si"][idx], float)
                    keys = np.asarray(ser["si_time"][idx], int)
                    if not np.isfinite(vals).any():
                        skipped.append(f"{sid} {idx}: all-NaN pixel"); continue
                    by = {int(k): float(v) for k, v in zip(keys, vals)
                          if np.isfinite(v)}
                    rows += _rows(by, freq, sid, ("era5", "", "era5_land"),
                                  threshold)
        finally:
            store.close()
        return rows, skipped

    from gcm_variables import GCMS, SCENARIOS
    gcms = gcms or list(GCMS)
    scenarios = scenarios or list(SCENARIOS)
    for idx in sorted(need):
        acc = int("".join(c for c in idx if c.isdigit()) or 12)
        freqs_here = [f for f in freqs if INDEX_FOR[f] == idx]
        for g in gcms:
            for sc in scenarios:
                f = CD.gcm_spei_file(g, sc, acc)
                if f is None:
                    skipped.append(f"{g}/{sc}: no spei_{acc}_{g}_{sc}_*.nc")
                    continue
                print(f"  {g:<15}{sc:<12}{f.name}", flush=True)
                for sid, (vals, keys) in CD.gcm_series(f, sites).items():
                    vals = np.asarray(vals, float)
                    if not np.isfinite(vals).any():
                        skipped.append(f"{sid} {g}/{sc} {idx}: all-NaN")
                        continue
                    by = {int(k): float(v) for k, v in zip(keys, vals)
                          if np.isfinite(v)}
                    for freq in freqs_here:
                        rows += _rows(by, freq, sid, ("gcm", g, sc), threshold)
    return rows, skipped


def main(argv=None) -> int:
    import argparse
    import csv
    from results_dir import NoResultsDir, resolve_out
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", choices=["era5", "gcm"], default="era5")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--threshold", type=float, default=-1.0)
    ap.add_argument("--gcms", default=None)
    ap.add_argument("--scenarios", default=None)
    ap.add_argument("--stations", default=None)
    ap.add_argument("--era5-root", type=Path, default=DEFAULT_ERA5_ROOT)
    a = ap.parse_args(argv)

    import classify_drought as CD
    want = {s.strip() for s in a.stations.split(",")} if a.stations else None
    sites = CD.read_sites(want)
    try:
        out = resolve_out(a.out or f"drought_periods_{a.source}.csv")
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    print(f"source    : {a.source}")
    print(f"stations  : {len(sites)}")
    print(f"threshold : {a.threshold}")
    print("windows   : annual = SPEI-12 at September (water-year end); "
          "monthly = SPEI-3 that month; seasonal = SPEI-3 at the season's "
          "last month")
    print("", flush=True)
    rows, skipped = build_table(
        sites, a.source, a.threshold,
        [g.strip() for g in a.gcms.split(",")] if a.gcms else None,
        [s.strip() for s in a.scenarios.split(",")] if a.scenarios else None,
        a.era5_root)
    if not rows:
        print("ERROR: nothing labelled", file=sys.stderr)
        for s_ in skipped[:10]:
            print(f"  {s_}", file=sys.stderr)
        return 1
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["source", "gcm", "scenario", "station", "freq", "year",
                    "period", "spei", "class"])
        w.writerows(rows)
    n_dry = sum(1 for r in rows if r[-1] == "drought")
    print("")
    print(f"{len(rows)} row(s), {n_dry} drought ({100*n_dry/len(rows):.1f}%)")
    for freq in ("annual", "monthly", "seasonal"):
        sub = [r for r in rows if r[4] == freq]
        if sub:
            d = sum(1 for r in sub if r[-1] == "drought")
            print(f"  {freq:<9} {len(sub):>8} rows, {100*d/len(sub):5.1f}% drought")
    if skipped:
        print("")
        print(f"SKIPPED {len(skipped)}:")
        for s_ in skipped[:8]:
            print(f"  ! {s_}")
    print("")
    print(f"-> {out}")
    return 1 if skipped else 0


if __name__ == "__main__":
    sys.exit(main())
