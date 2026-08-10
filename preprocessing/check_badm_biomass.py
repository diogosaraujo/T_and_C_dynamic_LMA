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
    # Only TREE seeds B(6): it is the standing stem biomass the heartwood pool
    # represents. SHRUB/CROP/OTHER are understory or non-woody components and
    # were lumped in here at first, which made a 70 gC/m2 shrub layer look like a
    # forest's aboveground biomass.
    "ag_tree":     (r"^AG_BIOMASS_TREE",),
    "ag_other":    (r"^AG_BIOMASS_(SHRUB|CROP|OTHER)",),
    "root":        (r"^ROOT_BIOMASS", r"^BG_BIOMASS", r"^BIOMASS_BG"),
    "lai":         (r"^LAI\b", r"^LAI_TOT"),
    # Discovery bucket: anything carrying BIOMASS/AGB/BA/STOCK that the specific
    # patterns above did not claim. It exists to answer "is there biomass we are
    # not looking at?" -- its distinct variable names are printed, so an
    # unanticipated field shows up as a name to add rather than staying invisible.
    "unclaimed":   (r"BIOMASS", r"\bAGB\b", r"^BASAL_AREA", r"^STAND_", r"_STOCK"),
}
# Fields that qualify a value rather than being one. The suffixes are NOT anchored
# with $ any more: BADM appends further tokens, so AG_BIOMASS_DATE_UNC,
# AG_BIOMASS_OTHER_SPATIAL_VARIABILITY and AG_BIOMASS_TREE_SPATIAL_REPLICATES all
# escaped an end-anchored pattern and were counted as measurements.
QUALIFIER = re.compile(r"_(UNIT|DATE|ORGAN|APPROACH|METHOD|COMMENT|STATISTIC|"
                       r"SPP|SPECIES|LOCATION|PUBLICATION|DEPTH|UNC|SPATIAL|"
                       r"VARIABILITY|REPLICATE)", re.I)

# Standing tree AGB below this is not a mature forest -- it is an annual
# increment, a subplot, or an understory component filed under the wrong
# variable. Flagged rather than dropped, so the report shows what was rejected.
AGB_PLAUSIBLE_MIN_TDMHA = 20.0


def to_t_dm_ha(value: float, unit: str) -> float | None:
    """Convert a BADM biomass value to [ton DM / ha], the unit TBio uses.

    Returns None when the unit is unrecognised rather than guessing: the two units
    seen in this dataset differ by a factor of 20, so a wrong assumption is not a
    rounding error.
    """
    u = (unit or "").strip().lower().replace(" ", "")
    if u in ("gcm-2", "gc/m2", "gcm^-2"):
        return value * 2.0 / 100.0          # gC/m2 -> gDM/m2 (f_C=0.5) -> tDM/ha
    if u in ("kgdmm-2", "kgdm/m2", "kgdmm^-2"):
        return value * 10.0                 # kgDM/m2 -> tDM/ha
    if u in ("tdmha-1", "tdm/ha", "mgdmha-1"):
        return value
    return None


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
    for cat in ("ag_tree", "ag_other", "root", "lai", "unclaimed"):
        have = [s for s in found if found[s].get(cat)]
        by_type = defaultdict(int)
        for s in have:
            by_type[wanted.get(s, "?")] += 1
        detail = ", ".join(f"{k} {v}" for k, v in sorted(by_type.items()))
        print(f"  {cat:<10} {len(have):>4} / {scanned} stations   ({detail})")
    if no_badm:
        print(f"\n  {len(no_badm)} station(s) with no BADM file at all: "
              f"{', '.join(no_badm[:8])}{' ...' if len(no_badm) > 8 else ''}")

    # Name every variable the specific patterns did not claim, so a biomass field
    # we have not thought of is visible instead of silently uncounted.
    unclaimed = sorted({e["variable"] for s in found
                        for e in found[s].get("unclaimed", [])})
    if unclaimed:
        print(f"\n  UNCLAIMED biomass-like variables ({len(unclaimed)} distinct) -- "
              f"add to CATEGORIES if any of these are usable:")
        for v in unclaimed[:25]:
            n = sum(1 for s in found for e in found[s].get("unclaimed", [])
                    if e["variable"] == v)
            print(f"    {v:<40} {n:>4} value(s)")
        if len(unclaimed) > 25:
            print(f"    ... and {len(unclaimed) - 25} more")
    else:
        print("\n  No unclaimed biomass-like variables: the patterns cover "
              "everything BADM reports here.")

    def unit_for(sid, var):
        pref = var.upper()
        for k, v in qualifiers[sid].items():
            if k.startswith(pref) and "UNIT" in k:
                return v
        for k, v in qualifiers[sid].items():          # fall back to the group unit
            if k.startswith("AG_BIOMASS") and "UNIT" in k:
                return v
        return ""

    print(f"\n{'=' * 70}\nTREE AGB -> TBio_target -> B(6)\n{'=' * 70}")
    print(f"  {'station':<9} {'type':<10} {'raw':>12} {'unit':<11} {'t DM/ha':>9}  use?")
    usable, rejected = {}, []
    for sid in sorted(s for s in found if found[s].get("ag_tree")):
        best = None
        for e in found[sid]["ag_tree"]:
            unit = unit_for(sid, e["variable"])
            t = to_t_dm_ha(float(e["value"]), unit)
            ok = t is not None and t >= AGB_PLAUSIBLE_MIN_TDMHA
            # Several records per site (different years/plots): keep the largest
            # plausible one, which is the standing stock rather than an increment.
            if ok and (best is None or t > best[0]):
                best = (t, e, unit)
            print(f"  {sid:<9} {wanted.get(sid,'?'):<10} {e['value']:>12} {unit:<11} "
                  f"{('%.1f' % t) if t is not None else '   ?':>9}  "
                  f"{'yes' if ok else ('unit?' if t is None else 'too small')}")
        if best:
            usable[sid] = best[0]
        else:
            rejected.append(sid)

    print(f"\n  {len(usable)} station(s) with a usable standing tree AGB "
          f"(>= {AGB_PLAUSIBLE_MIN_TDMHA:g} t DM/ha)")
    if rejected:
        print(f"  {len(rejected)} rejected (increment, subplot or unknown unit): "
              f"{', '.join(rejected)}")
    if usable:
        print(f"\n  {'station':<9} {'type':<10} {'t DM/ha':>9} {'implied B(6) gC/m2':>20}")
        for sid, t in sorted(usable.items()):
            # B(6) = 50*TBio - (B1+B2+B3+B4); the active pools are ~1054 gC/m2 in
            # the current initialisation. Shown to size the correction, not to
            # prescribe it -- the builder recomputes with that station's own B1.
            print(f"  {sid:<9} {wanted.get(sid,'?'):<10} {t:>9.1f} "
                  f"{max(0.0, 50*t - 1054):>20.0f}")

        # ---- PFT means: one IC vector per forest type, not per station -------
        # Coverage is far too sparse to initialise 91 stations individually, and
        # per-site tuning is not the goal. A single deciduous vector, built the
        # same way US_xRM's evergreen vector was, is the deliverable.
        print(f"\n{'=' * 70}\nPFT MEANS -- one IC vector per forest type\n{'=' * 70}")
        by_pft = defaultdict(list)
        for sid, t in usable.items():
            by_pft[wanted.get(sid, "?")].append(t)
        for ft in sorted(by_pft):
            v = sorted(by_pft[ft])
            mean = sum(v) / len(v)
            med = v[len(v) // 2] if len(v) % 2 else (v[len(v)//2 - 1] + v[len(v)//2]) / 2
            print(f"  {ft:<10} n={len(v):<3} mean {mean:7.1f}  median {med:7.1f}  "
                  f"range {min(v):.0f}-{max(v):.0f} t DM/ha")
            print(f"  {'':<10} -> TBio_target {mean:.0f} implies "
                  f"B(6) = 50*{mean:.0f} - (B1+B2+B3+B4)")
            if len(v) < 5:
                print(f"  {'':<10} !! n = {len(v)} is a thin sample; report the "
                      f"range, and cross-check against the literature value for "
                      f"the biome before adopting the mean.")
    else:
        print("\n  Nothing usable. B(6) cannot be set per site from BADM; the "
              "fallback is a gridded product (GEDI L4A/L4B or NBCD).")

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
