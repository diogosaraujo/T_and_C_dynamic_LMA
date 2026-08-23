#!/usr/bin/env python3
"""Seasonal signature of the LMA treatment: the dyn-fixed difference by day of year.

    python analyze_daily_effect.py --root $MODEL_RUN --pair 'era5_land:*_ic'
    python analyze_daily_effect.py --root $MODEL_RUN --pair 'historical/*:*' \
        --drought drought_years.csv

WHY DAILY, WHEN THE ANNUAL NUMBERS ALREADY EXIST

analyze_lma_effect.py reduces each run to annual series, which is right for the
fleet-wide headline and for future trends. But an annual mean hides compensating
changes: +8% in spring against -4% in summer averages to +1% and looks like
nothing happened.

The mechanism here is PHENOLOGICAL. LMA enters T&C only as Sl in LAI = Sl*B(1)
and never touches Vmax, so the treatment acts on leaf area -- which means the
effect should concentrate around leaf-out, peak LAI and senescence, and should
look different for deciduous and evergreen canopies. A day-of-year climatology
shows that directly instead of leaving it to be inferred.

WHAT IT WRITES

Long format, one row per (station, key, doy, variable):

    station,key,doy,variable,fixed,dyn,diff,rel,n_years[,drought]

'fixed' and 'dyn' are the multi-year mean for that day of year, 'diff' is
dyn - fixed, and 'rel' is diff as a percentage of |fixed|. Long format because
the downstream questions -- group by PFT, by ecoregion, by drought class -- are
all group-bys, and a wide table would need reshaping for every one of them.

DAY OF YEAR, ACROSS THREE CALENDARS. The GCM forcing is mapped onto real dates by
build_gcm_meteo.real_dates before it ever reaches T&C, so every run carries real
(year, month, day) in Datam regardless of whether its GCM used 365-day, 360-day
or standard. Day of year is therefore computed from month/day against a fixed
NON-LEAP table, so 1 March is always day 60 and the climatology bins line up
across calendars and across leap and common years. 29 February is DROPPED -- one
day in 1460, and keeping it would either misalign every subsequent day in leap
years or need a 366th bin fed by a quarter of the sample.

NO FALLBACKS. A pair missing a RES, a RES missing a required field, or a
fixed/dyn pair whose arrays disagree in length is reported by name and skipped;
the exit code is non-zero if anything was skipped. A seasonal cycle quietly
built from three variables instead of twelve is worse than no seasonal cycle.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_treatment_effect import find_pairs                    # noqa: E402

# Cumulative days before each month, non-leap. DOY = CUM[month-1] + day.
CUM = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]

# Daily arrays T&C already stores, and how a day's value is formed.
DAILY_MEAN = ["LAI_H"]
DAILY_SUM = ["NPP_H", "ANPP_H", "RA_H"]
# Hourly arrays, aggregated to a day. Water fluxes are mm/h -> sum; energy and
# state are means, because a daily mean W/m2 is the quantity anyone plots.
HOURLY_SUM = ["T_H", "T_L", "EG", "EIn_H", "EIn_L", "EIn_urb", "EIn_rock",
              "ESN", "ESN_In", "ELitter", "Lk", "Pr"]
HOURLY_MEAN = ["QE", "H", "Rn", "G", "Ta", "Ds"]
# What actually gets reported. Ordered by distance from the mechanism: LAI is the
# mediator, everything else is downstream of it.
REPORT = ["LAI_H", "GPP", "NPP", "ET", "T", "EG", "EIn", "Lk",
          "QE", "H", "Rn", "Tfrac", "Bowen", "WUE"]


class Unusable(Exception):
    """This pair cannot produce a climatology, and why."""


def read_run(path: Path) -> dict:
    """Daily series for one run, plus the (month, day, year) of each day."""
    import h5py
    with h5py.File(path, "r") as f:
        dm = np.asarray(f["Datam"][()], dtype=float)
        if dm.shape[0] != 4:
            dm = dm.T                       # v7.3 hands h5py the transpose
        need = DAILY_MEAN + DAILY_SUM + ["QE"]
        miss = [k for k in need if k not in f]
        if miss:
            raise Unusable(f"RES has no {', '.join(miss)}")

        def flat(k):
            a = np.asarray(f[k][()], dtype=float)
            return a.ravel() if (a.ndim == 1 or 1 in a.shape) else a

        nh = flat("QE").size
        nd = flat("LAI_H").size
        if nh < 24 or nd < 1:
            raise Unusable(f"{nh} hourly and {nd} daily steps")
        # The daily arrays run one step per day; Datam is hourly. Take every
        # 24th hourly timestamp as that day's date.
        yr, mo, da = (dm[i][:nh][::24][:nd].astype(int) for i in range(3))
        if yr.size != nd:
            raise Unusable(f"{yr.size} dates for {nd} daily values")

        out = {"year": yr, "mo": mo, "da": da}
        for k in DAILY_MEAN + DAILY_SUM:
            v = flat(k)
            if v.size != nd:
                raise Unusable(f"'{k}' has {v.size} days, LAI_H has {nd}")
            out[k] = v
        # Hourly -> daily. Trim to whole days so the reshape is exact.
        nfull = (nh // 24) if (nh // 24) <= nd else nd
        for k in HOURLY_SUM + HOURLY_MEAN:
            if k not in f:
                continue
            v = flat(k)[:nfull * 24].reshape(nfull, 24)
            agg = v.sum(axis=1) if k in HOURLY_SUM else v.mean(axis=1)
            out[k] = np.pad(agg, (0, nd - nfull), constant_values=np.nan)
    return out


def derive(d: dict) -> dict:
    """The reported quantities. Partitions matter more here than totals."""
    # A variable absent from RES contributes zero -- but as a zero ARRAY, not a
    # scalar. A scalar propagates through the sums and only fails later, at the
    # point where the result is indexed by the day mask, with an error that
    # names neither the variable nor the run.
    nd = np.asarray(d["LAI_H"]).size
    have = lambda k: d[k] if k in d else np.zeros(nd)
    ET = sum(have(k) for k in HOURLY_SUM if k not in ("Lk", "Pr"))
    T = have("T_H") + have("T_L")
    EIn = sum(have(k) for k in ("EIn_H", "EIn_L", "EIn_urb", "EIn_rock", "ELitter"))
    GPP = d["NPP_H"] + d["RA_H"]
    with np.errstate(divide="ignore", invalid="ignore"):
        out = {
            "LAI_H": d["LAI_H"], "GPP": GPP, "NPP": d["NPP_H"],
            "ET": ET, "T": T, "EG": have("EG"), "EIn": EIn, "Lk": have("Lk"),
            "QE": have("QE"), "H": have("H"), "Rn": have("Rn"),
            "Tfrac": np.where(ET > 0, T / ET, np.nan),
            "Bowen": np.where(np.abs(have("QE")) > 1e-6, have("H") / have("QE"), np.nan),
            "WUE": np.where(ET > 0, GPP / ET, np.nan),
        }
    return out


def doy_of(mo, da):
    """Day of year on a fixed non-leap calendar; 29 February -> 0 (dropped)."""
    d = np.array([0 if (m == 2 and x == 29) else CUM[m - 1] + x
                  for m, x in zip(mo, da)], dtype=int)
    return d


def climatology(fx: dict, dy: dict, years=None):
    """Per-DOY multi-year mean of each arm, for the days both arms share."""
    n = min(fx["LAI_H"].size, dy["LAI_H"].size)
    if n == 0:
        raise Unusable("no overlapping days")
    doy = doy_of(fx["mo"][:n], fx["da"][:n])
    keep = doy > 0
    if years is not None:
        keep &= np.isin(fx["year"][:n], list(years))
    if not keep.any():
        raise Unusable("no days left after filtering")

    F, D = derive({k: v[:n] for k, v in fx.items() if k != "year"}), \
           derive({k: v[:n] for k, v in dy.items() if k != "year"})
    doy, yr = doy[keep], fx["year"][:n][keep]
    rows = []
    for var in REPORT:
        if var not in F or var not in D:
            continue
        f_all, d_all = np.asarray(F[var])[keep], np.asarray(D[var])[keep]
        for j in range(1, 366):
            m = doy == j
            if not m.any():
                continue
            f = np.nanmean(f_all[m])
            v = np.nanmean(d_all[m])
            if not (np.isfinite(f) and np.isfinite(v)):
                continue
            rel = 100.0 * (v - f) / abs(f) if abs(f) > 1e-12 else np.nan
            rows.append((j, var, f, v, v - f, rel, int(np.unique(yr[m]).size)))
    return rows


def read_drought(path: Path):
    """{(station, year): class} from a CSV with station,year,class columns."""
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            out[(r["station"].strip(), int(r["year"]))] = r["class"].strip()
    if not out:
        raise SystemExit(f"ERROR: {path} has no rows")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--pair", default=None,
                    help="glob over pair labels, e.g. 'era5_land:*_ic' or "
                         "'historical/*:*'. Labels are <scenario>[/<GCM>]:<arm>.")
    ap.add_argument("--stations", default=None, help="comma-separated subset")
    ap.add_argument("--out", type=Path, default=None,
                    help="output CSV (default <root>/daily_effect.csv)")
    ap.add_argument("--drought", type=Path, default=None,
                    help="CSV of station,year,class -- writes one climatology per "
                         "class so drought and non-drought can be compared")
    a = ap.parse_args(argv)

    if not a.root.is_dir():
        print(f"ERROR: model_run root not found: {a.root}", file=sys.stderr)
        return 1
    out = a.out or (a.root / "daily_effect.csv")
    dro = read_drought(a.drought) if a.drought else None

    stations = ([s.strip() for s in a.stations.split(",")] if a.stations else
                sorted(p.name for p in a.root.iterdir()
                       if p.is_dir() and next(p.glob("**/fixed_lma*"), None)))
    print(f"model_run : {a.root}\npair      : {a.pair or '(all)'}\n"
          f"stations  : {len(stations)}"
          f"{'   drought classes from ' + a.drought.name if dro else ''}\n")

    rows, skipped, npairs = [], [], 0
    for st in stations:
        for label, fxp, dyp in find_pairs(a.root, st, a.pair):
            if fxp is None or dyp is None:
                skipped.append((st, label, "one arm has no RES"))
                continue
            npairs += 1
            try:
                fx, dy = read_run(fxp), read_run(dyp)
                groups = {"all": None}
                if dro is not None:
                    yrs = set(fx["year"].tolist())
                    for cls in sorted({dro.get((st, y)) for y in yrs} - {None}):
                        groups[cls] = {y for y in yrs if dro.get((st, y)) == cls}
                for cls, yrsel in groups.items():
                    for r in climatology(fx, dy, yrsel):
                        rows.append((st, label, cls) + r)
            except Unusable as e:
                skipped.append((st, label, str(e)))
            except Exception as e:                       # noqa: BLE001
                skipped.append((st, label, f"{type(e).__name__}: {e}"))

    if skipped:
        print(f"SKIPPED -- {len(skipped)} pair(s):")
        for st, label, why in skipped[:20]:
            print(f"  ! {st:<10}{label:<28}{why}")
        if len(skipped) > 20:
            print(f"  ... and {len(skipped) - 20} more")
        print()
    if not rows:
        print("ERROR: nothing was analysed.", file=sys.stderr)
        return 1

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["station", "key", "class", "doy", "variable",
                    "fixed", "dyn", "diff", "rel_pct", "n_years"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], r[3], r[4],
                        f"{r[5]:.6g}", f"{r[6]:.6g}", f"{r[7]:.6g}",
                        "" if not np.isfinite(r[8]) else f"{r[8]:.3f}", r[9]])

    # A one-screen summary: which variables move, and when in the year.
    print(f"{'variable':<9}{'mean |rel|':>11}{'peak |rel|':>11}{'peak doy':>10}")
    for var in REPORT:
        sel = [r for r in rows if r[4] == var and r[2] == "all"
               and np.isfinite(r[8])]
        if not sel:
            continue
        by_doy = {}
        for r in sel:
            by_doy.setdefault(r[3], []).append(abs(r[8]))
        means = {d: float(np.mean(v)) for d, v in by_doy.items()}
        peak = max(means, key=means.get)
        print(f"{var:<9}{np.mean(list(means.values())):>10.2f}%"
              f"{means[peak]:>10.2f}%{peak:>10}")

    print(f"\n{npairs} pair(s), {len(rows)} row(s) -> {out}")
    return 1 if skipped else 0


if __name__ == "__main__":
    sys.exit(main())
