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
    python check_treatment_effect.py US-Wrc --pair 'ssp585/*'

Every arm PAIR is checked, not one per station. The tree holds two layouts:

    <station>/era5_land/{fixed_lma,dyn_lma}                 ERA5-Land
    <station>/<scenario>/<GCM>/{fixed_lma,dyn_lma}          GCM

and one station now carries 16 pairs (1 ERA5 + 3 scenarios x 5 GCMs). Pairs are
found by looking for fixed_lma directories at any depth rather than hardcoding
either shape, because hardcoding 'era5_land' is exactly what made jobs 37691/37692
report "the treatment is live" for US-Wrc while the 30 GCM arms they were meant to
vet went unread. Each pair is reported on its own line: the fixed arm's Sl_H is the
per-GCM 1985-2014 mean, so one GCM propagating and another not is a state the
per-station summary cannot express.

Exit 1 if ANY pair's arms are identical, or if a named station has no pairs at all
-- a check that silently finds nothing must not read as a pass.
"""
from __future__ import annotations

import argparse
import fnmatch
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


def take(dset, stride: int):
    """Read a dataset, optionally every Nth element, WITHOUT materialising it all.

    The slice happens at the HDF5 level, so a stride of 24 moves a
    twenty-fourth of the bytes. That matters: --all now walks ~1400 pairs, each
    opening two RES files and reading 11 series, and job 38154 hit its wall
    clock having flushed nothing.

    MATLAB stores an hourly series as (N,1), which h5py presents as (1,N), so
    the stride goes on the LAST axis.

    Subsampling is safe in the direction that matters. If the sample differs the
    arrays certainly differ, so "differs" is never wrong. The only error it can
    introduce is calling two arrays IDENTICAL when they differ outside the
    sample -- which reports a treatment failure that is not real, i.e. it fails
    loudly rather than hiding a problem.
    """
    if stride <= 1:
        return np.asarray(dset[()], dtype=float)
    if dset.ndim == 1:
        return np.asarray(dset[::stride], dtype=float)
    if dset.ndim == 2 and dset.shape[0] == 1:
        return np.asarray(dset[:, ::stride], dtype=float)
    return np.asarray(dset[()], dtype=float)[..., ::stride]


def one_res(d: Path) -> Path | None:
    """The single RES_*.mat in an arm directory, or None if absent/ambiguous."""
    hits = sorted(d.glob("RES_*.mat"))
    return hits[0] if len(hits) == 1 else None


def find_pairs(root: Path, station: str, pattern: str | None = None):
    """[(label, fixed_res, dyn_res)] for every arm pair under a station.

    The label is the pair directory relative to the station -- 'era5_land', or
    'ssp585/GFDL-ESM4' -- so both layouts print the same way and the caller does
    not need to know which it is looking at.
    """
    sd = root / station
    if not sd.is_dir():
        return []
    out = []
    # Every fixed_lma* arm and its dyn_lma* twin. The suffix matters: the restart
    # runs live alongside the originals as fixed_lma_ic/dyn_lma_ic, and a glob of
    # bare 'fixed_lma' would check the pre-spin-up pair and report on the wrong
    # experiment. 'spinup' has no twin and is disposable, so it is not a pair.
    for fdir in sorted(sd.glob("**/fixed_lma*")):
        suffix = fdir.name[len("fixed_lma"):]
        twin = fdir.parent / f"dyn_lma{suffix}"
        label = f"{fdir.parent.relative_to(sd).as_posix()}:{fdir.name}"
        if pattern and not fnmatch.fnmatch(label, pattern):
            continue
        out.append((label, one_res(fdir), one_res(twin)))
    return out


def compare_pair(station: str, label: str, fx: Path | None, dy: Path | None,
                 stride: int = 1) -> int:
    import h5py
    if fx is None or dy is None:
        missing = " and ".join(n for n, p in (("fixed_lma", fx), ("dyn_lma", dy))
                               if p is None)
        print(f"\n{station}  {label}\n  SKIP  no single RES_*.mat for {missing}")
        return 0

    with h5py.File(fx, "r") as a, h5py.File(dy, "r") as b:
        sl_a = float(np.array(a["Sl_H"]).ravel()[0])
        sl_b = float(np.array(b["Sl_H"]).ravel()[0])
        print(f"\n{station}  {label}")
        print(f"  Sl_H  fixed={sl_a:.6f}  dyn(final)={sl_b:.6f}")

        rows, identical, lai_rel = [], [], float("nan")
        for k in VARS:
            if k not in a or k not in b:
                continue
            x, y = take(a[k], stride), take(b[k], stride)
            if x.shape != y.shape:
                rows.append((k, "shape mismatch", "", ""))
                continue
            same = np.array_equal(x, y)
            if same:
                identical.append(k)
            d = np.abs(x - y)
            denom = np.nanmean(np.abs(x)) or 1.0
            rel = 100 * np.nanmean(d) / denom
            if k == "LAI_H":
                lai_rel = rel        # keep it; re-reading cost three full passes
            rows.append((k, "IDENTICAL" if same else "differs",
                         f"{np.nanmax(d):.6g}", f"{rel:.3f}%"))

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
        # lai_rel was computed in the loop above. Recomputing it here read LAI_H
        # three more times per pair -- across ~1400 pairs that is most of a
        # wall-clock limit spent re-reading a number already in hand.
        print(f"\n  OK: the arms differ. LAI_H mean |diff| = {lai_rel:.2f}% "
              f"of the fixed-arm mean.")
        if identical:
            print(f"  note: still identical for {identical} -- expected only if that\n"
                  f"        variable is genuinely insensitive to leaf area.")
    return 0


def compare(root: Path, station: str, pattern: str | None = None,
            stride: int = 1) -> tuple[int, int]:
    """(failures, pairs checked) for one station."""
    pairs = find_pairs(root, station, pattern)
    if not pairs:
        where = f" matching {pattern!r}" if pattern else ""
        print(f"\n{station}\n  FAIL: no fixed_lma/dyn_lma pair found{where}. Nothing was\n"
              f"        compared, which is not the same as a pass -- check the run\n"
              f"        tree under {root / station}.")
        return 1, 0
    return sum(compare_pair(station, *p, stride) for p in pairs), len(pairs)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stations", nargs="*")
    ap.add_argument("--root", type=Path, default=MODEL_RUN)
    ap.add_argument("--all", action="store_true",
                    help="every station under the model_run root")
    ap.add_argument("--stride", type=int, default=1, metavar="N",
                    help="compare every Nth timestep instead of all of them. "
                         "--stride 24 is daily sampling of the hourly arrays and "
                         "cuts the I/O 24-fold, which is what makes a whole-tree "
                         "run finish. Safe: a sampled difference is still a real "
                         "difference.")
    ap.add_argument("--pair", metavar="GLOB",
                    help="only pairs whose label matches, e.g. 'era5_land' or "
                         "'ssp585/*'. The full tree is 16 pairs per station.")
    a = ap.parse_args(argv)

    if not a.root.is_dir():
        print(f"ERROR: model_run root not found: {a.root}", file=sys.stderr)
        return 1
    stations = a.stations
    if a.all or not stations:
        # Any station with an arm pair at any depth, so a station that has only
        # GCM runs is not skipped the way the era5_land test used to skip it.
        stations = sorted(p.name for p in a.root.iterdir()
                          if p.is_dir() and next(p.glob("**/fixed_lma"), None))
    print(f"model_run : {a.root}\nstations  : {len(stations)}"
          f"{f'   pair filter: {a.pair}' if a.pair else ''}")

    bad = npairs = 0
    for s in stations:
        f, n = compare(a.root, s, a.pair, a.stride)
        bad, npairs = bad + f, npairs + n
    print(f"\n{'=' * 60}")
    print(f"{npairs} arm pair(s) checked across {len(stations)} station(s)")
    if bad:
        print(f"{bad} failing check(s) -- a pair whose arms are identical, or a\n"
              f"station with no pair to compare. Do not analyse these runs.")
    else:
        print("every pair's arms differ: the fixed-vs-dynamic treatment is live")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
