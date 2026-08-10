#!/usr/bin/env python3
"""Check every AmeriFlux station against the EPA Level III ecoregion it is paired to.

The pairing drives which ecoregion's PLSR fit supplies a station's LMA series, so a
station filed under the wrong ecoregion is silently modelled with another region's
leaf traits. Self-consistency checks cannot catch that -- the pairing file can be
perfectly consistent and still put a Maryland site in the Southeastern Plains -- so
this does the only test that settles it: point-in-polygon against the EPA shapefile.

The EPA distributes us_eco_l3 in USA Contiguous Albers Equal Area (GRS80), so the
station coordinates have to be projected before the polygon test. The projection
parameters are READ FROM THE .prj rather than hard-coded, and the result is checked
against control stations whose ecoregion is not in doubt; if any control fails the
script aborts instead of reporting, because a wrong projection would otherwise
produce a confident, entirely fictional set of mismatches.

    python verify_ecoregion_pairing.py --shapefile $TC_INPUT_DATA/ecoregions/us_eco_l3.shp
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PAIRING = REPO_ROOT / "T&C" / "dynamic_lma_test" / "ecoregion_ameriflux_pairing.csv"

# Stations whose Level III ecoregion is unambiguous. These exist to prove the
# projection and polygon test are right before any mismatch is believed. They are
# deliberately spread across the country and across ecoregion shapes.
CONTROLS = {
    "US-xRM": "Southern Rockies",
    "US-UMB": "Northern Lakes and Forests",
    "US-Ha1": "Northeastern Highlands",
    "US-Me2": "Eastern Cascades Slopes and Foothills",
    "US-MMS": "Interior Plateau",
    "US-Blk": "Middle Rockies",
}


# --------------------------------------------------------------- projection

class AlbersFromPrj:
    """Ellipsoidal Albers equal-area conic, parameterised from a .prj file.

    Snyder, Map Projections -- A Working Manual (1987), pp. 101-102. Only the
    forward direction is needed: stations go to the shapefile's coordinate system,
    never the other way.
    """

    def __init__(self, prj_text: str):
        def num(key, default=None):
            m = re.search(r'PARAMETER\["' + key + r'",\s*(-?[0-9.eE+]+)\]',
                          prj_text, re.I)
            if m:
                return float(m.group(1))
            if default is None:
                raise ValueError(f"{key} not found in .prj")
            return default

        m = re.search(r'SPHEROID\["[^"]*",\s*(-?[0-9.eE+]+),\s*(-?[0-9.eE+]+)',
                      prj_text, re.I)
        if not m:
            raise ValueError("SPHEROID not found in .prj")
        self.a = float(m.group(1))
        inv_f = float(m.group(2))
        f = 1.0 / inv_f
        self.e2 = 2 * f - f * f
        self.e = math.sqrt(self.e2)

        if not re.search(r'Albers', prj_text, re.I):
            raise ValueError("the .prj is not an Albers projection; "
                             "this script cannot project to it")

        phi1 = math.radians(num("Standard_Parallel_1"))
        phi2 = math.radians(num("Standard_Parallel_2"))
        self.lon0 = math.radians(num("Central_Meridian"))
        phi0 = math.radians(num("Latitude_Of_Origin"))
        self.fe = num("False_Easting", 0.0)
        self.fn = num("False_Northing", 0.0)

        m1, m2 = self._m(phi1), self._m(phi2)
        q0, q1, q2 = self._q(phi0), self._q(phi1), self._q(phi2)
        self.n = ((m1 * m1 - m2 * m2) / (q2 - q1) if abs(phi1 - phi2) > 1e-12
                  else math.sin(phi1))
        self.C = m1 * m1 + self.n * q1
        self.rho0 = self.a * math.sqrt(self.C - self.n * q0) / self.n

    def _m(self, phi):
        s = math.sin(phi)
        return math.cos(phi) / math.sqrt(1 - self.e2 * s * s)

    def _q(self, phi):
        s = math.sin(phi)
        if self.e < 1e-12:
            return 2 * s
        return (1 - self.e2) * (s / (1 - self.e2 * s * s)
                                - (1 / (2 * self.e))
                                * math.log((1 - self.e * s) / (1 + self.e * s)))

    def __call__(self, lat, lon):
        phi, lam = math.radians(lat), math.radians(lon)
        rho = self.a * math.sqrt(self.C - self.n * self._q(phi)) / self.n
        theta = self.n * (lam - self.lon0)
        return (self.fe + rho * math.sin(theta),
                self.fn + self.rho0 - rho * math.cos(theta))


# ------------------------------------------------------------ polygon test

def point_in_shape(x, y, points, parts):
    """Even-odd ray casting over every ring of a shapefile polygon.

    Shapefile holes are inner rings with opposite winding; the even-odd rule
    handles them without needing to know which ring is which.
    """
    inside = False
    bounds = list(parts) + [len(points)]
    for pi in range(len(parts)):
        ring = points[bounds[pi]:bounds[pi + 1]]
        n = len(ring)
        j = n - 1
        for i in range(n):
            xi, yi = ring[i]
            xj, yj = ring[j]
            if (yi > y) != (yj > y):
                if x < (xj - xi) * (y - yi) / (yj - yi) + xi:
                    inside = not inside
            j = i
    return inside


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairing", type=Path, default=DEFAULT_PAIRING)
    ap.add_argument("--shapefile", type=Path, required=True,
                    help="us_eco_l3.shp from the EPA Ecoregions download")
    ap.add_argument("--out", type=Path, default=None,
                    help="write the full per-station result to this CSV")
    a = ap.parse_args()

    try:
        import shapefile  # pyshp
    except ImportError:
        print("ERROR: pyshp is not installed. Add 'pyshp' to preprocessing/"
              "requirements.txt and re-run slurm/setup_env.sh", file=sys.stderr)
        return 1
    if not a.shapefile.is_file():
        print(f"ERROR: shapefile not found: {a.shapefile}", file=sys.stderr)
        return 1
    prj = a.shapefile.with_suffix(".prj")
    if not prj.is_file():
        print(f"ERROR: no .prj beside the shapefile ({prj}); the projection cannot "
              "be determined and will not be guessed", file=sys.stderr)
        return 1

    proj = AlbersFromPrj(prj.read_text(encoding="utf-8", errors="replace"))
    print(f"projection : Albers a={proj.a:.1f} n={proj.n:.6f} "
          f"lon0={math.degrees(proj.lon0):.1f} (read from {prj.name})")

    sf = shapefile.Reader(str(a.shapefile))
    names = [f[0] for f in sf.fields[1:]]
    try:
        i_code = names.index("US_L3CODE")
        i_name = names.index("US_L3NAME")
    except ValueError:
        print(f"ERROR: US_L3CODE/US_L3NAME not among fields: {names}", file=sys.stderr)
        return 1

    polys = []
    for sr in sf.iterShapeRecords():
        s = sr.shape
        if not s.points:
            continue
        polys.append((s.bbox, s.points, list(s.parts),
                      str(sr.record[i_code]).strip(), str(sr.record[i_name]).strip()))
    print(f"shapefile  : {len(polys)} Level III polygons\n")

    def lookup(lat, lon):
        x, y = proj(lat, lon)
        for bbox, pts, parts, code, name in polys:
            if bbox[0] <= x <= bbox[2] and bbox[1] <= y <= bbox[3]:
                if point_in_shape(x, y, pts, parts):
                    return code, name
        return None, None

    rows = [r for r in csv.DictReader(open(a.pairing, newline="", encoding="utf-8-sig"))
            if (r.get("StationID") or "").strip()]

    # Controls first: a wrong projection must stop the run, not colour the report.
    print("control stations (projection + polygon test must reproduce these):")
    bad_ctrl = 0
    for r in rows:
        sid = r["StationID"].strip()
        if sid not in CONTROLS:
            continue
        _, got = lookup(float(r["Lat"]), float(r["Lon"]))
        ok = (got == CONTROLS[sid])
        bad_ctrl += (not ok)
        print(f"  {'ok  ' if ok else 'FAIL'} {sid:<8} expected {CONTROLS[sid]:<38} got {got}")
    if bad_ctrl:
        print(f"\nERROR: {bad_ctrl} control station(s) failed -- the projection or the "
              "polygon test is wrong, so no mismatch below would be trustworthy. "
              "Not reporting.", file=sys.stderr)
        return 1
    print()

    out, mismatch, outside = [], [], []
    for r in rows:
        sid = r["StationID"].strip()
        lat, lon = float(r["Lat"]), float(r["Lon"])
        code, name = lookup(lat, lon)
        assigned_code = (r.get("US_L3CODE") or "").strip()
        assigned = (r.get("US_L3NAME") or "").strip()
        status = ("outside" if name is None else
                  "ok" if name == assigned else "MISMATCH")
        out.append({"StationID": sid, "StationName": r.get("StationName", ""),
                    "Lat": lat, "Lon": lon, "ForestType": r.get("ForestType", ""),
                    "assigned_L3CODE": assigned_code, "assigned_L3NAME": assigned,
                    "actual_L3CODE": code or "", "actual_L3NAME": name or "",
                    "status": status})
        if status == "MISMATCH":
            mismatch.append(out[-1])
        elif status == "outside":
            outside.append(out[-1])

    print(f"{'=' * 74}\n{len(rows) - len(mismatch) - len(outside)}/{len(rows)} "
          f"stations fall inside the ecoregion they are paired to\n{'=' * 74}")
    if mismatch:
        print(f"\n{len(mismatch)} MISMATCH(ES) -- paired to an ecoregion that does not "
              f"contain the station:")
        for m in sorted(mismatch, key=lambda d: d["StationID"]):
            print(f"  {m['StationID']:<8} {m['StationName'][:38]:<38}")
            print(f"           paired to [{m['assigned_L3CODE']:>2}] {m['assigned_L3NAME']}")
            print(f"           actually  [{m['actual_L3CODE']:>2}] {m['actual_L3NAME']}")
    if outside:
        print(f"\n{len(outside)} station(s) fell outside every polygon (coastal "
              f"sliver, or outside CONUS -- check individually):")
        for m in outside:
            print(f"  {m['StationID']:<8} {m['Lat']:.3f}, {m['Lon']:.3f}  "
                  f"paired to {m['assigned_L3NAME']}")

    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        with open(a.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
            w.writeheader()
            w.writerows(out)
        print(f"\nfull result -> {a.out}")

    return 1 if mismatch else 0


if __name__ == "__main__":
    sys.exit(main())
