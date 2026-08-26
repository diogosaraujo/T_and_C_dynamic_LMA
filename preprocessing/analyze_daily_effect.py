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

    station,key,class,doy,variable,fixed,dyn,diff,rel_pct,rel_ann_pct,n_years

'fixed' and 'dyn' are the multi-year mean for that day of year and 'diff' is
dyn - fixed. Long format because the downstream questions -- group by PFT, by
ecoregion, by drought class -- are all group-bys, and a wide table would need
reshaping for every one of them.

TWO RELATIVE COLUMNS, AND WHICH ONE IS VALID DEPENDS ON THE VARIABLE.

  FLUXES (LAI_H, GPP, NPP, ET, T, EG, EIn, Lk, QE, H, Rn) -> read rel_ann_pct.
  These are extensive and go to zero out of season, so rel_pct's same-day
  denominator explodes: job 38604 reported ground evaporation at a "32348190%
  peak" on doy 43, which is a February EG of ~1e-9 mm/d underneath, not an
  effect. rel_ann_pct divides by the fixed arm's mean magnitude over the whole
  record -- a percentage of a typical day, finite year-round.

  RATIOS (Tfrac, Bowen, WUE) -> read rel_pct, the opposite way round. They are
  intensive and dimensionless, so an annual mean magnitude is meaningless: a
  winter Bowen near 160 sits beside a summer 0.7 and swamps the scale. Their
  denominator is healthy because climatology() rebuilds them from aggregated
  fluxes rather than averaging daily ratios.

Measured on a clean +2% change in H with QE untouched: rel_pct returns 2.000%
on every day of the year while rel_ann_pct returns 0.017% in summer and 3.967%
in winter. The printed summary picks the right column per variable and names
which one it used.

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
from results_dir import NoResultsDir, resolve_out                # noqa: E402

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
# Ratios, and the two fluxes each is built from. Kept out of the per-day
# averaging path in climatology() -- see the note there.
RATIOS = {"Tfrac": ("T", "ET"), "Bowen": ("H", "QE"), "WUE": ("GPP", "ET")}

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
        # The daily arrays run one step per day; Datam is hourly, so every 24th
        # hourly timestamp is that day's date. T&C's daily arrays are ONE STEP
        # LONGER than the hourly record supports -- 13150 values against 13149
        # dated days on a 1985-2020 run -- which analyze_lma_effect.py also
        # handles ("daily array runs one step past the hours"). Demanding exact
        # equality here made all 92 pairs fail in job 38485.
        #
        # Tolerate that one step and nothing more: a larger gap is a truncated
        # or mismatched run and must still fail loudly.
        yr, mo, da = (dm[i][:nh][::24].astype(int) for i in range(3))
        ndates = yr.size
        if not 0 <= nd - ndates <= 1:
            raise Unusable(f"{ndates} dated days against {nd} daily values")
        # Trim to the dates rather than padding them. analyze_lma_effect pads,
        # repeating the last year so its annual sums keep every value; here a
        # padded day would land as a duplicate sample on one day of year, and an
        # invented date is worse than a dropped one at 1 in 13149. Both scripts
        # take daily index k to be hourly day k.
        nd = min(nd, ndates)
        yr, mo, da = yr[:nd], mo[:nd], da[:nd]

        out = {"year": yr, "mo": mo, "da": da}
        for k in DAILY_MEAN + DAILY_SUM:
            v = flat(k)
            if v.size < nd:
                raise Unusable(f"'{k}' has {v.size} days, short of the {nd} dated")
            out[k] = v[:nd]
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

    # RATIOS ARE REBUILT FROM THE AGGREGATED FLUXES, NOT AVERAGED.
    # derive() forms Bowen, Tfrac and WUE per DAY, and averaging a daily ratio
    # across years is not the ratio of the averages. Out of season the
    # denominator flux is near zero, single days reach Bowen ~ 160 against a
    # summer 0.7, and those days dominate any mean they enter -- including the
    # annual scale rel_ann divides by. Tested: on a clean +2% change in H, the
    # mean-of-ratios route reports the true summer effect as 0.017% while
    # ratio-of-means recovers 2.00% exactly.
    #
    # So per day of year, sum the components first and divide once. That is the
    # quantity anyone means by "the Bowen ratio on day j" anyway.
    for var, (num, den) in RATIOS.items():
        if not all(k in F and k in D for k in (num, den)):
            continue
        fn, fd = np.asarray(F[num])[keep], np.asarray(F[den])[keep]
        dn, dd = np.asarray(D[num])[keep], np.asarray(D[den])[keep]
        per_doy = {}
        for j in range(1, 366):
            m = doy == j
            if not m.any():
                continue
            fdm, ddm = np.nansum(fd[m]), np.nansum(dd[m])
            if not (abs(fdm) > 1e-9 and abs(ddm) > 1e-9):
                continue                      # no denominator flux on this day
            per_doy[j] = (np.nansum(fn[m]) / fdm, np.nansum(dn[m]) / ddm,
                          int(np.unique(yr[m]).size))
        if not per_doy:
            continue
        scale = float(np.mean([abs(f) for f, _, _ in per_doy.values()]))
        for j, (f, v, ny) in sorted(per_doy.items()):
            rel = 100.0 * (v - f) / abs(f) if abs(f) > 1e-12 else np.nan
            rel_ann = (100.0 * (v - f) / scale if scale > 1e-12 else np.nan)
            rows.append((j, var, f, v, v - f, rel, rel_ann, ny))

    for var in REPORT:
        if var in RATIOS:
            continue                          # handled above
        if var not in F or var not in D:
            continue
        f_all, d_all = np.asarray(F[var])[keep], np.asarray(D[var])[keep]
        # Dividing by the SAME DAY's fixed-arm value explodes wherever that
        # value is near zero, which for a seasonal flux is most of winter. Job
        # 38604 reported ground evaporation at "89817% mean, 32348190% peak",
        # peaking on doy 43 -- that is a February EG of ~1e-9 mm/d in the
        # denominator, not an effect. Every variable whose peak landed in
        # midwinter (T 362, EIn 351, Lk 342, Rn 333, Tfrac 8, WUE 24) was
        # reporting the same artefact, while the ones that stayed physical
        # peaked at leaf-out (LAI 108, GPP 109, ET 99, QE 113).
        #
        # So also normalise by a denominator that does not vanish: the fixed
        # arm's mean magnitude across the whole record. rel_ann is then "this
        # day's difference as a percentage of a typical day", comparable across
        # variables and stations, and finite whenever the variable is not
        # identically zero. rel_pct is kept because it is the right measure
        # where the denominator is healthy -- near peak season -- but it is the
        # one to distrust, and the summary no longer uses it.
        if not np.isfinite(f_all).any():
            continue
        scale = float(np.nanmean(np.abs(f_all)))
        for j in range(1, 366):
            m = doy == j
            if not m.any():
                continue
            # Check for any finite value BEFORE averaging. nanmean over an
            # all-NaN day is what raised the "Mean of empty slice" warnings in
            # job 38604's stderr; the result was discarded a line later anyway,
            # so this skips the day instead of computing a NaN noisily.
            fj, dj = f_all[m], d_all[m]
            if not (np.isfinite(fj).any() and np.isfinite(dj).any()):
                continue
            f = np.nanmean(fj)
            v = np.nanmean(dj)
            if not (np.isfinite(f) and np.isfinite(v)):
                continue
            rel = 100.0 * (v - f) / abs(f) if abs(f) > 1e-12 else np.nan
            rel_ann = (100.0 * (v - f) / scale
                       if np.isfinite(scale) and scale > 1e-12 else np.nan)
            rows.append((j, var, f, v, v - f, rel, rel_ann,
                         int(np.unique(yr[m]).size)))
    return rows


def read_drought(path: Path):
    """{(station, gcm, scenario, year): class} from classify_drought.py's CSV.

    Keyed on the run family, not just the station, because a GCM's dry years are
    its own. ACCESS-CM2's driest years at US-NR1 are 2013/1984/2009 against an
    observed 2012/2013/2002/2006 -- applying ERA5 labels to a GCM run would
    mislabel almost every year. An era5 table carries gcm="" and
    scenario="era5_land"; a gcm table carries both.
    """
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rd = csv.DictReader(fh)
        if "gcm" not in (rd.fieldnames or []):
            raise SystemExit(f"ERROR: {path} predates the gcm/scenario columns. "
                             f"Re-run classify_drought.py -- an old table cannot "
                             f"say which run family its years belong to.")
        for r in rd:
            out[(r["station"].strip(), r["gcm"].strip(),
                 r["scenario"].strip(), int(r["year"]))] = r["class"].strip()
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
    # A bare "--out era5_daily.csv" used to land in the working directory, i.e.
    # preprocessing/, putting a 108 MB table inside the repo. Relative names now
    # resolve under $TC_RESULTS; absolute ones are honoured as given.
    try:
        out = resolve_out(a.out or "daily_effect.csv")
        # --drought is an INPUT, but it comes from classify_drought.py, which
        # writes to the same place -- so a bare name resolves there too. Do not
        # create anything for it; a missing label file must be reported, not
        # papered over with an empty directory.
        dro_path = resolve_out(a.drought, create=False) if a.drought else None
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    if dro_path is not None and not dro_path.is_file():
        print(f"ERROR: drought labels not found: {dro_path}\n"
              f"       Generate them with slurm/submit_classify_drought.sh.",
              file=sys.stderr)
        return 1
    dro = read_drought(dro_path) if dro_path else None

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
                    # The pair label carries the run family: 'era5_land:arm' or
                    # 'historical/GFDL-ESM4:arm'.
                    scen, _, _ = label.partition(":")
                    gcm = ""
                    if "/" in scen:
                        scen, _, gcm = scen.partition("/")
                    yrs = set(fx["year"].tolist())
                    cls_of = {y: dro.get((st, gcm, scen, y)) for y in yrs}
                    if not any(cls_of.values()):
                        skipped.append((st, label, f"no drought labels for "
                                                   f"({gcm or '-'}, {scen})"))
                    for cls in sorted({c for c in cls_of.values() if c}):
                        groups[cls] = {y for y in yrs if cls_of[y] == cls}
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
                    "fixed", "dyn", "diff", "rel_pct", "rel_ann_pct",
                    "n_years"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], r[3], r[4],
                        f"{r[5]:.6g}", f"{r[6]:.6g}", f"{r[7]:.6g}",
                        "" if not np.isfinite(r[8]) else f"{r[8]:.3f}",
                        "" if not np.isfinite(r[9]) else f"{r[9]:.3f}", r[10]])

    # A one-screen summary: which variables move, and when in the year. Read off
    # rel_ann (r[9]), not rel_pct (r[8]) -- see the note in climatology(). The
    # numbers are percentages of a typical day for that variable, so "peak doy"
    # is where the treatment bites hardest rather than where the denominator
    # happened to be smallest.
    print(f"{'variable':<9}{'mean |rel|':>11}{'peak |rel|':>11}{'peak doy':>10}"
          f"   basis")
    for var in REPORT:
        # WHICH COLUMN IS TRUSTWORTHY DEPENDS ON THE VARIABLE.
        # A flux is extensive and vanishes out of season, so its same-day
        # denominator explodes -- read rel_ann, the percentage of a typical day.
        # A ratio is intensive and, now that it is built from aggregated fluxes,
        # has a healthy denominator year-round -- but its annual mean magnitude
        # is meaningless, because winter Bowen values of ~160 sit alongside a
        # summer 0.7. Tested on a clean +2% change in H: rel_pct gives 2.000%
        # every day, rel_ann gives 0.017% in summer and 3.967% in winter.
        col, basis = ((8, "% of same day") if var in RATIOS else
                      (9, "% of annual mean"))
        sel = [r for r in rows if r[4] == var and r[2] == "all"
               and np.isfinite(r[col])]
        if not sel:
            continue
        by_doy = {}
        for r in sel:
            by_doy.setdefault(r[3], []).append(abs(r[col]))
        means = {d: float(np.mean(v)) for d, v in by_doy.items()}
        peak = max(means, key=means.get)
        print(f"{var:<9}{np.mean(list(means.values())):>10.2f}%"
              f"{means[peak]:>10.2f}%{peak:>10}   {basis}")

    print(f"\n{npairs} pair(s), {len(rows)} row(s) -> {out}")
    return 1 if skipped else 0


if __name__ == "__main__":
    sys.exit(main())
