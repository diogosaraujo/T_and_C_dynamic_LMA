#!/usr/bin/env python3
"""Is the fixed-vs-dynamic LMA treatment actually doing anything?

The first full array produced dyn_lma runs that were BIT-IDENTICAL to fixed_lma at
every station: MAIN_FRAME_SLA updated Sl_H each year, but VEGETATION_DYNAMIC reads
Sl from VegH_Param_Dyn, which Restating_parameters filled once before the time loop.
Every flux, pool and state matched exactly and the experiment measured nothing.

That failure is invisible in any single-run figure -- both arms look perfectly
healthy on their own. It only shows up when the two are differenced, so this is the
check to run before spending any time on analysis.

    python check_treatment_effect.py US-Ha2 US-HBK
    python check_treatment_effect.py --all

Exit 1 if any station's arms are identical, so a SLURM wrapper can gate on it.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

MODEL_RUN = Path(os.environ.get(
    "MODEL_RUN",
    Path(os.environ.get("TC_INPUT_DATA",
                        "/vol_efthymios/NFS07/dd1136/T_and_C/input_data")).parent
    / "model_run"))

# LAI first: LMA enters T&C only as Sl in LAI = Sl*B(1), so if LAI does not move,
# nothing downstream can. The rest are the fluxes the study reports.
VARS = ["LAI_H", "An_H", "B_H", "T_H", "EIn_H", "EG", "QE", "H", "Rn", "SWE", "Lk"]


def arm_file(root: Path, station: str, arm: str) -> Path | None:
    d = root / station / "era5_land" / arm
    hits = sorted(d.glob("RES_*.mat"))
    return hits[0] if len(hits) == 1 else None


def compare(root: Path, station: str) -> int:
    import h5py
    fx, dy = arm_file(root, station, "fixed_lma"), arm_file(root, station, "dyn_lma")
    if fx is None or dy is None:
        print(f"{station:<9} SKIP  missing RES for "
              f"{'fixed_lma' if fx is None else ''}{'dyn_lma' if dy is None else ''}")
        return 0

    with h5py.File(fx, "r") as a, h5py.File(dy, "r") as b:
        sl_a = float(np.array(a["Sl_H"]).ravel()[0])
        sl_b = float(np.array(b["Sl_H"]).ravel()[0])
        print(f"\n{station}")
        print(f"  Sl_H  fixed={sl_a:.6f}  dyn(final)={sl_b:.6f}")

        rows, identical = [], []
        for k in VARS:
            if k not in a or k not in b:
                continue
            x, y = np.array(a[k]), np.array(b[k])
            if x.shape != y.shape:
                rows.append((k, "shape mismatch", "", ""))
                continue
            same = np.array_equal(x, y)
            if same:
                identical.append(k)
            d = np.abs(x - y)
            denom = np.nanmean(np.abs(x)) or 1.0
            rows.append((k, "IDENTICAL" if same else "differs",
                         f"{np.nanmax(d):.6g}", f"{100*np.nanmean(d)/denom:.3f}%"))

        print(f"  {'variable':<8} {'verdict':<10} {'max|diff|':>12} {'mean|diff| as % of |fixed|':>28}")
        for k, verdict, mx, pc in rows:
            flag = "  <--" if verdict == "IDENTICAL" else ""
            print(f"  {k:<8} {verdict:<10} {mx:>12} {pc:>28}{flag}")

        if "LAI_H" in identical:
            print("\n  FAIL: LAI_H is identical between arms. LMA enters T&C only via\n"
                  "        LAI = Sl*B(1), so the treatment is not propagating. Check that\n"
                  "        MAIN_FRAME_SLA refreshes VegH_Param_Dyn.Sl inside the loop and\n"
                  "        that model_run/Code is not a stale copy.")
            return 1
        print(f"\n  OK: the arms differ. LAI_H mean |diff| = "
              f"{100*np.nanmean(np.abs(np.array(a['LAI_H'])-np.array(b['LAI_H'])))/(np.nanmean(np.abs(np.array(a['LAI_H']))) or 1):.2f}% "
              f"of the fixed-arm mean.")
        if identical:
            print(f"  note: still identical for {identical} -- expected only if that\n"
                  f"        variable is genuinely insensitive to leaf area.")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stations", nargs="*")
    ap.add_argument("--root", type=Path, default=MODEL_RUN)
    ap.add_argument("--all", action="store_true",
                    help="every station under the model_run root")
    a = ap.parse_args(argv)

    if not a.root.is_dir():
        print(f"ERROR: model_run root not found: {a.root}", file=sys.stderr)
        return 1
    stations = a.stations
    if a.all or not stations:
        stations = sorted(p.name for p in a.root.iterdir()
                          if p.is_dir() and (p / "era5_land").is_dir())
    print(f"model_run : {a.root}\nstations  : {len(stations)}")

    bad = sum(compare(a.root, s) for s in stations)
    print(f"\n{'=' * 60}")
    if bad:
        print(f"{bad} station(s) with NO treatment effect -- do not analyse these runs")
    else:
        print("every station's arms differ: the fixed-vs-dynamic treatment is live")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
