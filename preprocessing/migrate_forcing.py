#!/usr/bin/env python3
"""Move already-built forcing out of input_data and into the run tree.

    python migrate_forcing.py --dry-run     # what would move, moves nothing
    python migrate_forcing.py               # do it
    python migrate_forcing.py --link        # leave a symlink behind (rollback aid)

ONE-TIME. New runs land in the right place already: build_gcm_meteo.py and
build_meteo_input.py stamp a dest_dir into each raw file and finish_meteo.m writes
there. This is only for the 1,515 GCM files and ~101 ERA5 files built before that.

    input_data/gcm_meteo/<scen>/<GCM_us>/Meteo_<ST>_<GCM_us>_<scen>_<years>.mat
        -> model_run/<ST>/<scen>/<GCM>/Meteo_<ST>_<GCM_us>_<scen>_<years>.mat
    input_data/meteo/Meteo_<ST>_<years>.mat
        -> model_run/<ST>/era5_land/Meteo_<ST>_<years>.mat

WHY NOT PARSE THE FILENAME. 'Meteo_US_Wrc_GFDL_ESM4_historical_1985_2014.mat'
splits on underscores that belong to the station, the GCM, the scenario and the
year tag alike, and no rule separates them. The station list and the GCM list are
known, so the mapping is built from those and the filename is only matched, never
dissected. A file that matches no known (station, GCM, scenario) is reported and
left alone rather than guessed at.

The move is os.rename within one filesystem -- instant, no copy, no second 54 GB.
It falls back to copy+unlink across filesystems, which is slow but correct.

The directory name differs on purpose between the two trees: the forcing tree uses
mat_name's underscored GCM (GFDL_ESM4) and the run tree uses the plain one
(GFDL-ESM4). The run tree wins; the underscored form survives only inside the
FILE name, where it belongs, since that is what GO's load() refers to.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_gcm_model_run import (GCM_METEO, MODEL_RUN, SCENARIOS,  # noqa: E402
                                 YEAR_TAG, read_stations)
from build_model_run import mat_name                               # noqa: E402
from gcm_variables import GCMS                                     # noqa: E402

ERA5_METEO = GCM_METEO.parent / "meteo"


def gcm_moves(meteo_root: Path, run_root: Path, stations):
    """[(src, dst)] for the GCM forcing, built from the known station/GCM lists."""
    out, unmatched = [], []
    for scen in SCENARIOS:
        for gcm in GCMS:
            d = meteo_root / scen / mat_name(gcm)
            if not d.is_dir():
                continue
            wanted = {}
            for st in stations:
                sid = st["station_id"]
                fn = (f"Meteo_{mat_name(sid)}_{mat_name(gcm)}_{scen}"
                      f"_{YEAR_TAG[scen]}.mat")
                wanted[fn] = run_root / sid / scen / gcm / fn
            for p in sorted(d.glob("Meteo_*.mat")):
                if p.name in wanted:
                    out.append((p, wanted[p.name]))
                else:
                    unmatched.append(p)
    return out, unmatched


def era5_moves(meteo_root: Path, run_root: Path, stations):
    """[(src, dst)] for the ERA5-Land forcing. The year tag is read off the file."""
    out, unmatched = [], []
    if not meteo_root.is_dir():
        return out, unmatched
    by_prefix = {f"Meteo_{mat_name(st['station_id'])}_": st["station_id"]
                 for st in stations}
    for p in sorted(meteo_root.glob("Meteo_*.mat")):
        if p.name.endswith("_raw.mat"):
            continue                      # intermediate, not forcing
        # Longest prefix wins: US_Ha2 must not claim a US_Ha2x file.
        hit = max((k for k in by_prefix if p.name.startswith(k)),
                  key=len, default=None)
        if hit is None:
            unmatched.append(p)
            continue
        out.append((p, run_root / by_prefix[hit] / "era5_land" / p.name))
    return out, unmatched


def move(src: Path, dst: Path, link: bool) -> str:
    if dst.exists():
        return "dst exists"
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.rename(src, dst)                       # same filesystem: instant
    except OSError:
        shutil.copy2(src, dst)                    # across filesystems
        src.unlink()
    if link:
        try:
            os.symlink(dst, src)
        except OSError:
            pass
    return "moved"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meteo", type=Path, default=GCM_METEO)
    ap.add_argument("--era5-meteo", type=Path, default=ERA5_METEO)
    ap.add_argument("--root", type=Path, default=MODEL_RUN,
                    help="model_run root the forcing moves into")
    ap.add_argument("--link", action="store_true",
                    help="leave a symlink at the old path; lets a half-migrated "
                         "tree keep working while you check the new one")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    if not a.root.is_dir():
        print(f"ERROR: model_run root not found: {a.root}", file=sys.stderr)
        return 1
    stations = read_stations()
    print(f"model_run : {a.root}\ngcm meteo : {a.meteo}\nera5 meteo: {a.era5_meteo}")
    print(f"stations  : {len(stations)}\n")

    g, g_un = gcm_moves(a.meteo, a.root, stations)
    e, e_un = era5_moves(a.era5_meteo, a.root, stations)
    print(f"GCM  forcing : {len(g)} file(s) to move")
    print(f"ERA5 forcing : {len(e)} file(s) to move")
    for p in (g_un + e_un)[:10]:
        print(f"  ? unmatched, left alone: {p}")
    if len(g_un) + len(e_un) > 10:
        print(f"  ... and {len(g_un) + len(e_un) - 10} more unmatched")

    if a.dry_run:
        for src, dst in (g + e)[:3]:
            print(f"\n  would move {src}\n          -> {dst}")
        print(f"\nDRY RUN -- nothing moved. {len(g) + len(e)} file(s) would move.")
        return 0

    tally: dict[str, int] = {}
    for src, dst in g + e:
        outcome = move(src, dst, a.link)          # called ONCE -- it is not idempotent
        tally[outcome] = tally.get(outcome, 0) + 1
    print()
    for k in sorted(tally):
        print(f"  {k:<12} {tally[k]}")
    moved = tally.get("moved", 0)
    print(f"\n{moved}/{len(g) + len(e)} file(s) moved into {a.root}")
    if tally.get("dst exists"):
        print("  'dst exists' means the run tree already had that file -- the source "
              "was left in place rather than overwriting it. Compare and delete by hand.")
    # Non-zero only if something was left behind, so a wrapper can gate on it.
    return 0 if moved == len(g) + len(e) and not (g_un or e_un) else 1


if __name__ == "__main__":
    sys.exit(main())
