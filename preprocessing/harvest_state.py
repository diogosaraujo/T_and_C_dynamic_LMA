#!/usr/bin/env python3
"""Take the final vegetation state out of finished runs, to restart from.

    python harvest_state.py --from '*/era5_land/fixed_lma'          # round 1
    python harvest_state.py --from '*/historical/*/spinup'          # round 2
    python harvest_state.py --from '*/era5_land/fixed_lma' --drift  # equilibrium check

WHY. The evergreen sites start with a leaf pool far below what the climate
supports, so LAI ramps up over the first years and that ramp is in the analysis
window. Restarting from an already-equilibrated state removes it without throwing
data away (Dr. Paschalis, 2026-08-19).

WHAT IS CARRIED, AND WHAT IS NOT. MOD_PARAM sets four initial-condition fields and
MAIN_FRAME's `run(PARAM_IC)` (line 222) executes after the arrays are preallocated
(lines 112-113), so those four assignments survive:

    LAI_H(1,:)   B_H(1,:,:)   PHE_S_H(1,:)   AgeL_H(1,:)

Everything else VEGGIE_UNIT carries between days -- dflo_H, AgeDL_H, and the
running integrators NPPI_H, TdpI_H, Bfac_weekH, NupI_H, PARI_H, NBLI_H, NBLeaf_H,
and the N/P/K reserves -- starts at its default. Those are short-memory quantities
(days to weeks) and re-equilibrate quickly, but this is a four-field restart, not a
full one, and it is not described as more than that. Extending it means adding new
assignments to MOD_PARAM, which is a change to the substitution table, not here.

ALL EIGHT POOLS ARE CARRIED, B(6) INCLUDED. Heartwood never leaves -- Wm = 0 for
all 8 PFTs -- so that pool grows without bound and US-Ha2 ends a 36-year run at
176% of the observed biomass mean (CLAUDE.md 9). Modelled biomass is therefore not
a validation target. It is still carried, for three reasons:

  * Zeroing it would make the pool depend on where the simulation happens to be
    cut into runs. Restart at 2014 instead of 2020 and the same simulated years
    give a different heartwood. A monotonic drift is at least reproducible.
  * Heartwood is produced from sapwood (Ss = ds*B(2)). Carrying the sapwood that
    generated it while discarding the product leaves the state self-inconsistent.
  * "Restart from the final state" is one rule; "except pool 6" is a rule plus an
    exception every later reader has to be told about.

It changes nothing mechanically either way here: B(6) reaches the model only
through Vegetation_Structural_Attributes, which MAIN_FRAME_SLA calls inside
`if OPT_VCA == 1` (line 416), and OPT_VCA = 0. If OPT_VCA is ever switched on, the
accumulated value needs reviewing THEN, deliberately -- not pre-emptily discarded
now by a rule nobody will remember.

NO FALLBACKS. A missing RES, a missing field, a NaN, an ambiguous array shape or a
non-finite pool is an error that names the station and stops. A restart built from
a state that was quietly patched is worse than no restart, because nothing
downstream can tell the difference.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

POOLS = 8
# The four fields MOD_PARAM can set. Order matters only for the CSV header.
FIELDS = ["LAI_H", "PHE_S_H", "AgeL_H"]


def res_file(d: Path) -> Path:
    hits = sorted(d.glob("RES_*.mat"))
    if len(hits) != 1:
        raise SystemExit(f"ERROR: {d} has {len(hits)} RES_*.mat, expected exactly 1")
    return hits[0]


def daily_series(f, key: str, station: str) -> np.ndarray:
    """A (ndays,) daily state series, with the canopy axis squeezed out.

    MATLAB v7.3 hands h5py the transpose, and cc == 1 for these runs, so the
    stored shape is (1, ndays) rather than (ndays, 1). Squeezing is safe; ravel
    is not, which is why this does not reuse analyze_lma_effect._read -- that
    helper ravels anything with a size-1 axis and would flatten B_H's pools.
    """
    if key not in f:
        raise SystemExit(f"ERROR: {station}: RES has no '{key}'")
    a = np.squeeze(np.asarray(f[key][()], dtype=float))
    if a.ndim != 1:
        raise SystemExit(f"ERROR: {station}: '{key}' is {a.shape} after squeeze, "
                         f"expected one series (is cc > 1?)")
    return a


def pools_series(f, station: str) -> np.ndarray:
    """B_H as (ndays, 8). Identifies the pool axis by length, never by position."""
    if "B_H" not in f:
        raise SystemExit(f"ERROR: {station}: RES has no 'B_H'")
    a = np.squeeze(np.asarray(f["B_H"][()], dtype=float))
    if a.ndim != 2:
        raise SystemExit(f"ERROR: {station}: B_H is {a.shape} after squeeze, "
                         f"expected 2 axes (days x pools)")
    if a.shape[0] == POOLS and a.shape[1] != POOLS:
        a = a.T                                   # (8, ndays) -> (ndays, 8)
    elif a.shape[1] != POOLS:
        raise SystemExit(f"ERROR: {station}: B_H is {a.shape}, neither axis is "
                         f"{POOLS} pools -- cannot tell days from pools")
    return a


def harvest(d: Path, station: str, key: str, drift: bool) -> dict:
    import h5py
    with h5py.File(res_file(d), "r") as f:
        B = pools_series(f, station)
        rec = {"station": station, "key": key, "ndays": B.shape[0]}

        last = B[-1, :]
        if not np.isfinite(last).all():
            raise SystemExit(f"ERROR: {station}: final B_H has non-finite pools: {last}")
        if (last < 0).any():
            raise SystemExit(f"ERROR: {station}: final B_H has negative pools: {last}")
        for i in range(POOLS):
            rec[f"B{i + 1}"] = float(last[i])

        for name in FIELDS:
            s = daily_series(f, name, station)
            if s.size != B.shape[0]:
                raise SystemExit(f"ERROR: {station}: '{name}' has {s.size} days, "
                                 f"B_H has {B.shape[0]}")
            v = float(s[-1])
            if not np.isfinite(v):
                raise SystemExit(f"ERROR: {station}: final {name} is not finite")
            rec[name] = v

        if drift:
            # Is it actually spun up? Evergreen wood pools need decades
            # (CLAUDE.md 5). Compare the live pools against a decade earlier;
            # a flat answer is evidence, a large one says run another cycle.
            back = min(B.shape[0] - 1, 10 * 365)
            prev = B[-1 - back, :]
            live = [0, 1, 2, 3]                   # leaf, sapwood, fine root, reserve
            num = float(np.abs(last[live] - prev[live]).sum())
            den = float(np.abs(last[live]).sum()) or 1.0
            rec["drift_10yr_pct"] = round(100.0 * num / den, 2)
            rec["drift_years"] = round(back / 365.0, 1)
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True,
                    help="model_run root")
    ap.add_argument("--from", dest="pattern", required=True,
                    help="glob under the root selecting arm directories, e.g. "
                         "'*/era5_land/fixed_lma' or '*/historical/*/spinup'")
    ap.add_argument("--out", type=Path, default=None,
                    help="CSV to write (default <root>/initial_state.csv); "
                         "existing rows with the same (station, key) are replaced")
    ap.add_argument("--drift", action="store_true",
                    help="add a 10-year drift column: how far the live pools moved "
                         "over the last decade of the run")
    a = ap.parse_args(argv)

    if not a.root.is_dir():
        print(f"ERROR: model_run root not found: {a.root}", file=sys.stderr)
        return 1
    out = a.out or (a.root / "initial_state.csv")

    dirs = sorted(p for p in a.root.glob(a.pattern) if p.is_dir())
    if not dirs:
        print(f"ERROR: '{a.pattern}' matched no directory under {a.root}. Nothing "
              f"was harvested, which is not the same as an empty result.",
              file=sys.stderr)
        return 1
    print(f"model_run : {a.root}\npattern   : {a.pattern}\nmatched   : {len(dirs)}\n")

    rows = []
    for d in dirs:
        rel = d.relative_to(a.root)
        station, key = rel.parts[0], Path(*rel.parts[1:]).as_posix()
        rows.append(harvest(d, station, key, a.drift))

    # Merge with whatever is already there, replacing by (station, key) so a later
    # round does not lose an earlier one.
    old = []
    if out.is_file():
        with open(out, newline="", encoding="utf-8-sig") as fh:
            old = [r for r in csv.DictReader(fh)]
    fresh = {(r["station"], r["key"]) for r in rows}
    merged = [r for r in old if (r["station"], r["key"]) not in fresh] + rows

    cols = list(rows[0])
    for r in merged:                              # older rounds may lack --drift
        for c in r:
            if c not in cols:
                cols.append(c)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in merged:
            w.writerow(r)

    print(f"{'station':<10}{'B1 leaf':>10}{'B2 sap':>10}{'B4 res':>10}"
          f"{'LAI_H':>8}{'PHE':>5}" + (f"{'drift10y':>10}" if a.drift else ""))
    for r in rows[:12]:
        line = (f"{r['station']:<10}{r['B1']:>10.1f}{r['B2']:>10.1f}{r['B4']:>10.1f}"
                f"{r['LAI_H']:>8.2f}{r['PHE_S_H']:>5.0f}")
        print(line + (f"{r['drift_10yr_pct']:>9.1f}%" if a.drift else ""))
    if len(rows) > 12:
        print(f"  ... and {len(rows) - 12} more")

    if a.drift:
        worst = sorted(rows, key=lambda r: -r["drift_10yr_pct"])[:5]
        print(f"\nlargest 10-year drift in the live pools (B1-B4):")
        for r in worst:
            print(f"  {r['station']:<10}{r['key']:<28}{r['drift_10yr_pct']:>7.1f}%")
        print("  A few percent means equilibrated. Tens of percent means the wood\n"
              "  pools are still moving and one more spin-up cycle is warranted.")

    print(f"\n{len(rows)} state(s) harvested, {len(merged)} row(s) in {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
