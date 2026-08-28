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

READ ONE CHUNK PER STATION. The .mat is chunked (516, 15, 1) with gzip, so
every chunk holds the COMPLETE time series for 15 latitudes at ONE longitude.
d[:, lat, lon] is therefore a single chunk, while d[t, lat, :] touches 3600
chunks and decompresses 516 months of each to keep one value. Reading by
latitude row -- the obvious layout given h5py reports (time, lat, lon) -- did
the second ~51,000 times and had to be cancelled after the fifteen GCM jobs
had all finished in minutes.

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
# pev, NOT e. In ERA5 "e" is total (actual) evaporation and "pev" is
# potential. Reading e_1980_2022_monthly.mat put every station below phi = 1
# with a maximum of 0.934 -- the signature of AET/P, which the water balance
# bounds at about 1, where PET/P has no such bound. Arizona ponderosa sites
# cannot be energy-limited.
ERA5_PET = Path("/vol_efthymios/NFS07/Data/ERA5_Land/monthly/"
                "pev_1980_2022_monthly.mat")
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
    """Annual PET (mm) per station from the ERA5-Land monthly .mat.

    ONE CHUNK PER STATION. The dataset is chunked (516, 15, 1) with gzip:
    every chunk holds the COMPLETE time series for 15 latitudes at a single
    longitude. So d[:, lat, lon] -- a whole time series at one point -- costs
    exactly one chunk read, while d[t, lat, :] -- one month of one latitude --
    touches 3600 separate compressed chunks and decodes all 516 months of each
    to keep one value. The first version did the latter ~51,000 times and had
    to be cancelled; this does the former 118 times.
    """
    import h5py

    lat_i = np.clip(np.round((90.0 - st["lat"].to_numpy(float)) / ERA5_DLAT),
                    0, 1800).astype(int)
    lon = st["lon"].to_numpy(float)
    cand = {"0-360": np.round(np.mod(lon, 360.0) / ERA5_DLAT).astype(int)
                     % ERA5_NLON,
            "-180-180": np.round((lon + 180.0) / ERA5_DLAT).astype(int)
                        % ERA5_NLON}

    months = [(y, m) for y in range(y0, y1 + 1) for m in range(12)]
    with h5py.File(ERA5_PET, "r") as f:
        d = f["data_all"]
        n_t = d.shape[0]
        idx = np.array([(y - ERA5_YEAR0) * 12 + m for y, m in months])
        keep = (idx >= 0) & (idx < n_t)
        if not keep.any():
            raise SystemExit(f"ERROR: {y0}-{y1} outside the ERA5 PET record")
        idx, months = idx[keep], [mm for mm, k in zip(months, keep) if k]
        probe = int(idx[len(idx) // 2])

        # Convention decided on the LAND MASK, using point reads rather than a
        # whole row: the Atlantic at 40N is ocean and therefore NaN, while a
        # US station maps to land under BOTH conventions, so counting finite
        # stations cannot separate them (job 40398 scored 116 vs 118 and chose
        # wrong). Ten probes per candidate, one chunk each.
        probes = np.linspace(0, 599, 10).astype(int)
        frac = {}
        for k, (lo_j, hi_j) in (("0-360", (2900, 3500)),
                                ("-180-180", (1100, 1700))):
            cols = lo_j + probes
            v = np.array([d[probe, 500, int(c)] for c in cols], dtype=float)
            frac[k] = float(np.isfinite(v).mean())
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

        dd = np.array([days_per_month(y)[m] for y, m in months], dtype=float)
        out, n_ok = [], 0
        for s in range(len(st)):
            ts = np.asarray(d[:, int(lat_i[s]), int(lon_i[s])], dtype=float)
            v = -1000.0 * ts[idx]                 # m/day, negative -> mm/day
            if np.isfinite(v).any():
                n_ok += 1
            annual = (v * dd).reshape(-1, 12).sum(axis=1)
            out.append({"station": st["station"].iloc[s],
                        "pet_mm_yr": float(np.nanmean(annual))})
        print(f"  stations on finite cells: {n_ok} of {len(st)}")
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
    r0 = r1 = c0 = c1 = None
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
                gi = np.abs(lav[None, :]
                            - st["lat"].to_numpy(float)[:, None]).argmin(1)
                gj = np.abs(lov[None, :] - slon[:, None]).argmin(1)
                # BOUNDING BOX. Reading the whole global field costs 83 MB per
                # file and 1035 files is ~86 GB to extract 92 points each. The
                # stations span a small part of the grid, so read that window
                # and shift the indices into it: same numbers, ~3% of the I/O.
                pad = 2
                r0, r1 = max(int(gi.min()) - pad, 0), int(gi.max()) + pad + 1
                c0, c1 = max(int(gj.min()) - pad, 0), int(gj.max()) + pad + 1
                lat_i, lon_i = gi - r0, gj - c0
                print(f"    window rows {r0}:{r1} cols {c0}:{c1} "
                      f"({100.0*(r1-r0)*(c1-c0)/(len(lav)*len(lov)):.1f}% "
                      f"of the grid)", flush=True)
            a = np.ma.filled(d.variables["pet"][:, r0:r1, c0:c1],
                             np.nan).astype(float)
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
