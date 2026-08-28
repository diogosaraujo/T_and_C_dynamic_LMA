#!/usr/bin/env python3
"""Dryness index phi = PET / P per station, per GCM, per dataset.

    phi = sum(PET) / sum(P)   over the whole period of each dataset

PET comes from the SAME products the PLSR predictors used -- ERA5-Land for the
ERA5 arm, NEX-GDDP for the GCM arms -- so the x axis of the dryness figures
sits on the same footing as the LMA on the y axis. A Priestley-Taylor PET from
T&C's own Rn would have been free (Rn is already in the annual tables) but it
is a different quantity, and the LMA slopes are themselves downstream of the
PLSR's PET.

P NEEDS NO EXTRACTION. Precipitation is already a variable row in the annual
effect tables, in mm/yr, identical in both arms because it is forcing. Reading
it from there rather than re-extracting keeps P and the LMA slopes on one
vintage and removes a whole class of unit error.

    era5_annual.csv / gcm_annual_<scenario>.csv   variable == "Pr"

TWO UNIT CONVERSIONS, AND THEY DIFFER. Getting either wrong shifts every
station relative to phi = 1, which is the line the figure is read against.

    ERA5-Land   data_all is potential_evaporation in m/day and NEGATIVE
                (downward flux). mm/day = -1000 * value, then * days in month.
    NEX-GDDP    pet is labelled "mm" on a monthly file but the values are
                mm/DAY: mean 2.9, max 13.7. Monthly totals would be 50-150.
                Summing the twelve values directly gives ~35 mm/yr, which
                against P ~ 1000 puts every station at phi ~ 0.035 -- uniformly
                hyper-humid, the whole sample crushed against the left edge,
                and phi = 1 never reached. Multiply by days in month.

READ CONTIGUOUS LATITUDE ROWS FROM THE .mat. It is MATLAB v7.3, so column
major: h5py reports (time, lat, lon) but lon is the fastest-varying axis.
Slicing axis 0 gathers millions of scattered elements out of a 27 GB file and
stalls. One latitude row is 3600 contiguous values and returns instantly, so
stations are grouped by latitude row and every station on a row is read at
once.

THE LONGITUDE CONVENTION IS DETERMINED, NOT ASSUMED. ERA5-Land grids appear
both as 0..359.9 and as -180..179.9, and at 40N a probe lands on land under
either. Both are tried, all stations sampled under each, and the one yielding
more finite cells wins -- a wrong convention puts every US station in Asia and
returns NaN or nonsense at most of them, so the test is decisive. The choice
and both counts are printed.

Output: $TC_RESULTS/station_dryness.csv
    dataset,gcm,station,pft,lat,lon,years,pet_mm_yr,p_mm_yr,phi
"""
from __future__ import annotations

import argparse
import calendar
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from results_dir import NoResultsDir, resolve_out                 # noqa: E402
from station_metrics import table_path, read_sites                # noqa: E402

GCMS = ["GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
        "MRI-ESM2-0", "UKESM1-0-LL"]
SCEN = {"historical": "historical", "ssp126": "ssp126", "ssp585": "ssp585"}
NEX = Path("/vol_efthymios/NFS07/Data/CMIP6/NEXGDDP")
ERA5_PET = Path("/vol_efthymios/NFS07/Data/ERA5_Land/monthly/"
                "e_1980_2022_monthly.mat")
ERA5_YEAR0 = 1980          # data_all starts at January 1980, 516 months
# ERA5-Land 0.1 deg: 1801 latitudes from +90 descending, 3600 longitudes.
ERA5_DLAT = 0.1
ERA5_NLON = 3600


def days_per_month(year: int) -> np.ndarray:
    return np.array([calendar.monthrange(year, m)[1] for m in range(1, 13)],
                    dtype=float)


def station_coords() -> pd.DataFrame:
    """Station lat/lon/pft. Coordinates come from the pairing used elsewhere."""
    from figure_skill_maps import read_sites as read_xy
    xy = read_xy()                       # {station: (lat, lon, pft)}
    rows = [{"station": s, "lat": v[0], "lon": v[1]} for s, v in xy.items()]
    d = pd.DataFrame(rows)
    if d.empty:
        raise SystemExit("ERROR: no station coordinates available")
    return d.merge(read_sites(), on="station", how="left")


# ----------------------------------------------------------------- ERA5
def era5_pet(st: pd.DataFrame, y0: int, y1: int) -> pd.DataFrame:
    """Annual PET (mm) per station from the ERA5-Land monthly .mat."""
    import h5py

    lat_i = np.clip(np.round((90.0 - st["lat"].to_numpy(float)) / ERA5_DLAT),
                    0, 1800).astype(int)
    lon = st["lon"].to_numpy(float)
    # The two candidate conventions, resolved below by which one lands on land.
    cand = {"0-360": np.round(np.mod(lon, 360.0) / ERA5_DLAT).astype(int)
                     % ERA5_NLON,
            "-180-180": np.round((lon + 180.0) / ERA5_DLAT).astype(int)
                        % ERA5_NLON}

    months = [(y, m) for y in range(y0, y1 + 1) for m in range(12)]
    idx = [(y - ERA5_YEAR0) * 12 + m for y, m in months]
    with h5py.File(ERA5_PET, "r") as f:
        d = f["data_all"]
        n_t = d.shape[0]
        idx = [i for i in idx if 0 <= i < n_t]
        if not idx:
            raise SystemExit(f"ERROR: {y0}-{y1} outside the ERA5 PET record")
        probe = int(idx[len(idx) // 2])
        rows = sorted(set(lat_i.tolist()))
        # ONE READ PER LATITUDE ROW, not per station: contiguous and instant,
        # where slicing the time axis gathers scattered elements and stalls.
        strip = {r: np.asarray(d[probe, r, :], dtype=float) for r in rows}
        # DECIDE ON THE LAND MASK, NOT ON THE STATIONS. Counting finite
        # stations cannot separate the two conventions: a site at -105 maps to
        # 255E (Colorado) under one and 75E (central Asia) under the other, and
        # BOTH are land in ERA5-Land. The first version scored 116 against 118
        # and picked the wrong one on two stations of noise, which would have
        # sampled Asia and produced entirely plausible-looking numbers.
        # The Atlantic at 40N is ocean, hence NaN, and that is unambiguous.
        row40 = np.asarray(d[probe, 500, :], dtype=float)
        fin = np.isfinite(row40)
        frac = {"0-360": float(fin[2900:3500].mean()),      # lon -70..-10
                "-180-180": float(fin[1100:1700].mean())}
        best = min(frac, key=frac.get)
        print(f"  longitude convention: {best}   (finite fraction where the "
              f"Atlantic must be: " +
              ", ".join(f"{k}={v:.3f}" for k, v in frac.items()) + ")")
        if frac[best] > 0.05 or (max(frac.values()) - min(frac.values())) < 0.1:
            raise SystemExit(
                "ERROR: the longitude convention is not resolved -- neither "
                "candidate puts ocean where the Atlantic is at 40N. Refusing "
                "to sample on a guess; check the grid before rerunning.")
        lon_i = cand[best]
        n_ok = int(np.isfinite([strip[r][j]
                                for r, j in zip(lat_i, lon_i)]).sum())
        print(f"  stations on finite cells: {n_ok} of {len(lat_i)}")

        out = []
        for r in rows:
            sel = np.where(lat_i == r)[0]
            cols = lon_i[sel]
            series = np.full((len(idx), len(sel)), np.nan)
            for k, t in enumerate(idx):
                series[k] = np.asarray(d[t, r, :], dtype=float)[cols]
            # m/day, negative -> mm/day, positive.
            series = -1000.0 * series
            dd = np.array([days_per_month(y)[m] for y, m in months
                           if 0 <= (y - ERA5_YEAR0) * 12 + m < n_t])
            annual = (series * dd[:, None]).reshape(-1, 12, len(sel)).sum(axis=1)
            for j, s in enumerate(sel):
                out.append({"station": st["station"].iloc[s],
                            "pet_mm_yr": float(np.nanmean(annual[:, j]))})
    return pd.DataFrame(out)


# ----------------------------------------------------------------- GCM
def gcm_pet(st: pd.DataFrame, gcm: str, scen: str,
            y0: int, y1: int) -> pd.DataFrame | None:
    """Annual PET (mm) per station from the NEX-GDDP monthly files."""
    import netCDF4 as nc

    files = {}
    for y in range(y0, y1 + 1):
        hits = sorted((NEX / gcm / scen / "pet").glob(
            f"pet_month_{gcm}_{scen}_*_{y}.nc"))
        if hits:
            files[y] = hits[0]
    if not files:
        return None

    lat_i = lon_i = None
    tot, cnt = np.zeros(len(st)), 0
    for y, p in sorted(files.items()):
        with nc.Dataset(p) as d:
            if lat_i is None:
                la = np.asarray(d.variables["latitude"][:], dtype=float)
                lo = np.asarray(d.variables["longitude"][:], dtype=float)
                # Stored 2D but a plain meshgrid; take the vectors.
                lav = la[:, 0] if la.ndim == 2 else la
                lov = lo[0, :] if lo.ndim == 2 else lo
                slon = st["lon"].to_numpy(float)
                if np.nanmax(lov) > 180.0:
                    slon = np.mod(slon, 360.0)
                lat_i = np.abs(lav[None, :]
                               - st["lat"].to_numpy(float)[:, None]).argmin(1)
                lon_i = np.abs(lov[None, :] - slon[:, None]).argmin(1)
            a = np.ma.filled(d.variables["pet"][:], np.nan).astype(float)
        # mm/DAY despite the "mm" label -- see the module docstring.
        yr = (a[:, lat_i, lon_i] * days_per_month(y)[:, None]).sum(axis=0)
        tot += yr
        cnt += 1
    return pd.DataFrame({"station": st["station"],
                         "pet_mm_yr": tot / max(cnt, 1)})


# ----------------------------------------------------------------- P
def precip(root: Path, ds: str) -> pd.DataFrame:
    """Mean annual precipitation per station and GCM, from the annual table."""
    p = table_path(Path(root), ds, "annual")
    if not Path(p).is_file():
        raise SystemExit(f"ERROR: {Path(p).name} not found; P comes from the "
                         f"annual effect tables, not from a re-extraction")
    d = pd.read_csv(p, usecols=lambda c: c in
                    {"station", "key", "year", "variable", "fixed"},
                    low_memory=False)
    d = d[d["variable"] == "Pr"]
    if d.empty:
        raise SystemExit(f"ERROR: no Pr rows in {Path(p).name}")
    d["gcm"] = d["key"].astype(str).str.extract(r"^[^/]+/([^:]+):",
                                                expand=False).fillna("")
    g = (d.groupby(["gcm", "station"], as_index=False)["fixed"].mean()
           .rename(columns={"fixed": "p_mm_yr"}))
    yrs = (int(d["year"].min()), int(d["year"].max()))
    g["years"] = f"{yrs[0]}-{yrs[1]}"
    return g


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--out", default="station_dryness.csv")
    ap.add_argument("--datasets", default="era5,historical,ssp126,ssp585")
    ap.add_argument("--gcms", default=",".join(GCMS))
    a = ap.parse_args(argv)
    try:
        root = Path(a.results or resolve_out(".", create=False))
        out_p = resolve_out(a.out)
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    st = station_coords()
    print(f"stations with coordinates: {len(st)}")
    gcms = [g.strip() for g in a.gcms.split(",") if g.strip()]
    rows, missing = [], []

    for ds in [x.strip() for x in a.datasets.split(",")]:
        P = precip(root, ds)
        y0, y1 = (int(v) for v in P["years"].iloc[0].split("-"))
        print(f"\n{ds}: P from the annual table, {y0}-{y1}, "
              f"{P['station'].nunique()} stations")
        if ds == "era5":
            pet = era5_pet(st, y0, y1)
            pet["gcm"] = ""
            got = [("", pet)]
        else:
            got = []
            for g in gcms:
                pg = gcm_pet(st, g, SCEN[ds], y0, y1)
                if pg is None:
                    missing.append(f"{ds}/{g}: no pet files"); continue
                pg["gcm"] = g
                got.append((g, pg))
                print(f"  {g:<16} PET mean {np.nanmean(pg['pet_mm_yr']):8.1f} mm/yr")
        for g, pg in got:
            m = (pg.merge(P[P["gcm"] == g], on=["station", "gcm"], how="inner")
                   .merge(st, on="station", how="left"))
            m["dataset"] = ds
            m["phi"] = m["pet_mm_yr"] / m["p_mm_yr"].replace(0, np.nan)
            rows.append(m)

    if not rows:
        print("ERROR: nothing computed", file=sys.stderr)
        for w in missing:
            print(f"  {w}", file=sys.stderr)
        return 1
    D = pd.concat(rows, ignore_index=True)
    cols = ["dataset", "gcm", "station", "pft", "lat", "lon", "years",
            "pet_mm_yr", "p_mm_yr", "phi"]
    D = D[[c for c in cols if c in D.columns]]
    D.to_csv(out_p, index=False)
    print(f"\n-> {out_p}  ({len(D)} rows)")
    print("\nphi by dataset (median across GCMs first, then across stations):")
    med = (D.groupby(["dataset", "station"], as_index=False)["phi"].median())
    print(med.groupby("dataset")["phi"]
             .describe()[["count", "min", "25%", "50%", "75%", "max"]]
             .to_string(float_format=lambda x: f"{x:.3f}"))
    # A CONUS forest sample should straddle phi = 1, roughly 0.4-1.5. Every
    # station on one side means a unit or convention error, not a wet decade.
    frac = float((med["phi"] > 1).mean())
    print(f"\nfraction water-limited (phi > 1): {100*frac:.1f}%")
    if frac in (0.0, 1.0):
        print("  WARNING: every station falls on one side of phi = 1, which is "
              "what a unit or longitude error looks like -- check the PET "
              "conversions before using this", file=sys.stderr)
    if missing:
        print("\nNOT READ:", file=sys.stderr)
        for w in missing:
            print(f"  {w}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
