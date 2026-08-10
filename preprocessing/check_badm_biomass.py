#!/usr/bin/env python3
"""How many stations report biomass in AmeriFlux BADM, and in what units?

Motivation: the T&C initial carbon pools B_H(1:8) are currently transplanted from
US_xRM, a single subalpine conifer site. The heartwood pool B(6) is the one that
matters most, because Vegetation_Structural_Attributes.m computes

    TBio = 0.02 * (B(1)+B(2)+B(3)+B(4)+B(6))    [ton DM / ha]

and TBio drives Allocation_Coefficients. With B(6) = 0 a mature forest is presented
to the model as a ~21 t DM/ha sapling. Setting B(6) from a site's real biomass,

    B(6) = max(0, 50*TBio_target - (B1+B2+B3+B4))

replaces that guess -- IF the sites report biomass. This script answers that, and
deliberately answers it BEFORE any conversion logic is written: it reports the raw
variables, values and units it finds rather than assuming a unit convention, since
BADM mixes gC/m2, kgDM/m2 and tDM/ha across sites.

    python check_badm_biomass.py                    # coverage over the site list
    python check_badm_biomass.py --csv out.csv      # also dump every value found

Reads the BADM/BIF files already downloaded by download_ameriflux.py; downloads
nothing itself.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

from inspect_ameriflux_badm import find_badm_files, read_badm_rows, MISSING_VALUES

REPO_ROOT = Path(__file__).resolve().parents[1]
PREPROC = Path(__file__).resolve().parent
INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA",
                                 "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
DEFAULT_DIR = INPUT_ROOT / "ameriflux"
DEFAULT_SITE_LISTS = [
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv",
]

# What counts as biomass, grouped by which T&C pool it could seed. Matching is on
# the variable name with a word boundary, not a bare substring: a plain "in"
# test on "LAI" also hits "LAI_STATISTIC" and, worse, "REPLAI"-style names in
# other groups.
CATEGORIES = {
    "aboveground": (r"^AG_BIOMASS", r"^BIOMASS_AG", r"^AGB\b"),
    "root":        (r"^ROOT_BIOMASS", r"^BG_BIOMASS", r"^BIOMASS_BG"),
    "lai":         (r"^LAI\b", r"^LAI_TOT"),
    "other_bio":   (r"BIOMASS",),      # anything else carrying the word
}
# Fields that qualify a value rather than being one.
QUALIFIER = re.compile(r"_(UNIT|DATE|ORGAN|APPROACH|METHOD|COMMENT|STATISTIC|"
                       r"SPP|SPECIES|LOCATION|PUBLICATION|DEPTH)$", re.I)


def classify(variable: str) -> str | None:
    v = variable.upper()
    for cat, pats in CATEGORIES.items():
        if any(re.search(p, v) for p in pats):
            return cat
    return None


def is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


def read_site_list(paths) -> dict[str, str]:
    """StationID -> forest type, for the stations we actually model."""
    out = {}
    for p in paths:
        if not Path(p).is_file():
            print(f"  ! site list not found: {p}", file=sys.stderr)
            continue
        for r in csv.DictReader(open(p, newline="", encoding="utf-8-sig")):
            sid = (r.get("StationID") or r.get("SITE_ID") or "").strip()
            if sid:
                out[sid] = (r.get("ForestType") or "").strip().lower()
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, default=DEFAULT_DIR,
                    help="root holding one directory per station")
    ap.add_argument("--site-list", type=Path, nargs="*", default=None)
    ap.add_argument("--csv", type=Path, default=None,
                    help="write every biomass value found to this CSV")
    a = ap.parse_args(argv)

    if not a.dir.is_dir():
        print(f"ERROR: BADM directory not found: {a.dir}\n"
              f"       Set TC_INPUT_DATA or pass --dir.", file=sys.stderr)
        return 1

    wanted = read_site_list(a.site_list or DEFAULT_SITE_LISTS)
    print(f"BADM root  : {a.dir}")
    print(f"site list  : {len(wanted)} stations we model\n")

    found: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    qualifiers: dict[str, dict[str, str]] = defaultdict(dict)
    no_badm, scanned = [], 0

    for station_dir in sorted(p for p in a.dir.iterdir() if p.is_dir()):
        sid = station_dir.name
        if wanted and sid not in wanted:
            continue
        scanned += 1
        files = find_badm_files(station_dir)
        if not files:
            no_badm.append(sid)
            continue
        for path in files:
            try:
                rows = read_badm_rows(path)
            except Exception as exc:                       # noqa: BLE001
                print(f"  ! {sid}: cannot read {path.name}: {exc}", file=sys.stderr)
                continue
            for row in rows:
                var = (row.get("VARIABLE") or "").strip()
                val = str(row.get("DATAVALUE") or "").strip()
                if not var or val in MISSING_VALUES:
                    continue
                cat = classify(var)
                if cat is None:
                    continue
                if QUALIFIER.search(var):
                    qualifiers[sid][var.upper()] = val       # units, dates, organ
                elif is_number(val):
                    found[sid][cat].append(
                        {"variable": var, "value": val,
                         "group": (row.get("VARIABLE_GROUP") or "").strip()})

    # ---------------------------------------------------------------- report
    print(f"{'=' * 70}\nCOVERAGE over {scanned} station directories\n{'=' * 70}")
    for cat in ("aboveground", "root", "lai"):
        have = [s for s in found if found[s].get(cat)]
        by_type = defaultdict(int)
        for s in have:
            by_type[wanted.get(s, "?")] += 1
        detail = ", ".join(f"{k} {v}" for k, v in sorted(by_type.items()))
        print(f"  {cat:<12} {len(have):>4} / {scanned} stations   ({detail})")
    if no_badm:
        print(f"\n  {len(no_badm)} station(s) with no BADM file at all: "
              f"{', '.join(no_badm[:8])}{' ...' if len(no_badm) > 8 else ''}")

    agb = sorted(s for s in found if found[s].get("aboveground"))
    if agb:
        print(f"\n{'=' * 70}\nABOVEGROUND BIOMASS -- the one that sets B(6)\n{'=' * 70}")
        print(f"  {'station':<9} {'type':<10} {'variable':<26} {'value':>12}  unit")
        for sid in agb:
            for e in found[sid]["aboveground"][:4]:
                unit = next((v for k, v in qualifiers[sid].items()
                             if k.startswith(e["variable"].upper()[:12])
                             and k.endswith("UNIT")), "")
                if not unit:
                    unit = next((v for k, v in qualifiers[sid].items()
                                 if k.endswith("UNIT")), "?")
                print(f"  {sid:<9} {wanted.get(sid,'?'):<10} {e['variable'][:26]:<26} "
                      f"{e['value']:>12}  {unit}")
        units = {u for s in agb for k, u in qualifiers[s].items() if k.endswith("UNIT")}
        print(f"\n  distinct units reported: {sorted(units) if units else 'NONE STATED'}")
        print("  -> convert per station, not globally: BADM mixes gC/m2, kgDM/m2 "
              "and tDM/ha.")
    else:
        print("\n  No aboveground biomass anywhere. B(6) cannot be set from BADM; "
              "fall back to a gridded product (NBCD/GEDI) or a PFT constant.")

    if a.csv:
        a.csv.parent.mkdir(parents=True, exist_ok=True)
        with open(a.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["StationID", "ForestType", "category", "variable_group",
                        "variable", "value", "qualifiers"])
            for sid in sorted(found):
                for cat, entries in found[sid].items():
                    for e in entries:
                        q = "; ".join(f"{k}={v}" for k, v in
                                      sorted(qualifiers[sid].items())
                                      if k.startswith(e["variable"].upper()[:12]))
                        w.writerow([sid, wanted.get(sid, ""), cat, e["group"],
                                    e["variable"], e["value"], q])
        print(f"\nfull dump -> {a.csv}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
