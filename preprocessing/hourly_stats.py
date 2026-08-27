"""Model-versus-tower statistics at HOURLY resolution, per station and arm.

WHY HOURLY IS THE ONLY DEFENSIBLE RESOLUTION FOR SKILL HERE. The annual
comparison rests on the overlap between the model record (1985-2020) and the
tower record, and at many sites that is almost nothing: US-HBK runs 2016-2024,
so the annual overlap is TWO YEARS. A correlation from two points is +/-1 by
construction. At n = 3 no correlation is significant at all, at n = 10 the
threshold is |r| = 0.63, and only past n ~ 30 does |r| = 0.36 become
distinguishable from zero. Hourly gives thousands of matched steps at every
station, so the sampling problem disappears.

    hourly_stats.csv
      station,pft,variable,arm,n,mean_obs,sd_obs,mean_mod,sd_mod,
      rmse,rsr,bias,r,skill_score

TWO THINGS THAT MUST BE RIGHT, AND ARE EASY TO GET WRONG:

TIME ZONE. AmeriFlux timestamps are LOCAL STANDARD TIME; T&C runs in UTC with
DeltaGMT = 0. Across CONUS that is a 5-8 hour offset, so an unshifted comparison
would misalign the diurnal cycle by most of a working day -- and it would look
plausible, because daily and coarser aggregates absorb it. The per-site offset
comes from the FLUXNET BIF (UTC_OFFSET; -5 at US-HBK), not from longitude, which
is wrong near time-zone boundaries.

CARBON UNITS. At HH the tower reports GPP in umol CO2 m-2 s-1, not the
gC m-2 d-1 of its own DD/MM/YY files -- confirmed on the data: a p99 of 34 and
27% negative values, which is a night-time respiration signal in an
instantaneous rate. T&C's hourly An_H and Rdark_H are in the same unit, since
VEG_DYN_RES forms GPP = 1.0368*(An + Rdark) in gC/m2/d and
1.0368 = 12.0 g/mol x 86400 s/d x 1e-6. So hourly model GPP is An_H + Rdark_H
with NO conversion. The DAILY GPP used elsewhere in this project is a different
quantity and must not be mixed in here.

NOTHING IS STORED PER TIMESTEP. 92 stations x 2 arms x 36 years of hourly values
is far too large to write out; each station is reduced to its statistics as it is
read, and only those are kept.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_treatment_effect import find_pairs                     # noqa: E402
from figure_skill_maps import Missing, read_sites                 # noqa: E402
from results_dir import NoResultsDir, resolve_out                 # noqa: E402

LAMBDA = 2.45e6                       # J/kg
W_TO_MM_H = 3600.0 / LAMBDA           # W/m2 -> mm/h
FLUXNET_NA = -9990.0                  # anything at or below is missing

# Model hourly arrays. ET is assembled from the water-flux components, exactly
# as the daily pipeline does, minus Lk and Pr which are not evaporation.
ET_PARTS = ["T_H", "T_L", "EG", "EIn_H", "EIn_L", "EIn_urb", "EIn_rock",
            "ESN", "ESN_In", "ELitter"]
NEED = ["QE", "H", "An_H", "Rdark_H"] + ET_PARTS

TOWER_COLS = {"GPP_NT": ["GPP_NT_VUT_REF"], "GPP_DT": ["GPP_DT_VUT_REF"],
              "LE": ["LE_F_MDS"], "H": ["H_F_MDS"]}
QC_OF = {"LE": "LE_F_MDS_QC", "H": "H_F_MDS_QC",
         "GPP_NT": "NEE_VUT_REF_QC", "GPP_DT": "NEE_VUT_REF_QC"}


def utc_offset(root: Path, sid: str, lon: float | None = None):
    """Hours to ADD to local standard time to reach UTC.

    THREE SOURCES, IN ORDER, AND THE ONE USED IS REPORTED. Job 39700 skipped
    every station because this demanded a UTC_OFFSET row at a fixed column
    index in the BIF. US-HBK has one; most sites apparently do not, and losing
    the entire fleet over a missing metadata field is far worse than a possible
    one-hour error.

    The error is also benign for the headline result: an offset that is wrong by
    an hour degrades BOTH arms identically, so the fixed-versus-dynamic skill
    score survives it. Only the absolute RMSE and correlation suffer, which is
    why the source is written into the output rather than hidden.
    """
    # "*BIF*" also matches the five BIFVARINFO files, and BIFVARINFO sorts
    # BEFORE BIF_ ("V" < "_"), so hits[0] was always the wrong file and every
    # site fell through to the longitude fallback -- job 39703 reported
    # "utc offset from longitude" everywhere while US-Bar's BIF held
    # "US-Bar,4579,UTC_OFFSET,UTC_OFFSET,-5" all along. Same collision as the
    # FLUXMET/ERA5/BIFVARINFO one in the data files.
    hits = sorted(Path(root).glob(f"**/*{sid}*_BIF_*.csv"))
    if hits:
        # The BIF is not UTF-8 at every site: job 39698 died on byte 0xb5, a
        # micro sign in a units string. Only one row is wanted, so decode
        # leniently. Columns are named, so find them rather than assume index 3.
        with hits[0].open(newline="", encoding="utf-8-sig", errors="replace") as fh:
            rd = csv.reader(fh)
            head = next(rd, [])
            up = [c.strip().upper() for c in head]
            iv = up.index("VARIABLE") if "VARIABLE" in up else 3
            idv = up.index("DATAVALUE") if "DATAVALUE" in up else 4
            for row in rd:
                if not row:
                    continue
                # Match on ANY field: the variable name has appeared in both the
                # VARIABLE and VARIABLE_GROUP columns.
                if not any(c.strip().upper() == "UTC_OFFSET" for c in row):
                    continue
                for j in (idv, len(row) - 1):
                    if j < len(row):
                        try:
                            return -float(row[j]), "BIF"
                        except ValueError:
                            continue
    if lon is not None and np.isfinite(lon):
        # Standard time zones follow 15-degree bands closely enough across
        # CONUS; boundaries are irregular, so this can be an hour out.
        return -float(round(lon / 15.0)), "longitude"
    return None, None


def read_model_hourly(path: Path) -> dict:
    """Hourly model series keyed by UTC hour (year*1000000+month*10000+day*100+h)."""
    import h5py
    with h5py.File(path, "r") as f:
        dm = np.asarray(f["Datam"][()], dtype=float)
        if dm.shape[0] != 4:
            dm = dm.T
        if dm.shape[0] < 4:
            raise Missing(f"{path.name}: Datam has no hour row")

        def flat(k):
            a = np.asarray(f[k][()], dtype=float)
            if a.ndim > 1 and 1 not in a.shape:
                a = a.sum(axis=1) if a.shape[0] < a.shape[1] else a.sum(axis=1)
            return np.asarray(a).ravel()

        if "QE" not in f:
            raise Missing(f"{path.name}: no QE")
        nh = flat("QE").size
        yr, mo, da, hh = (dm[i][:nh].astype(int) for i in range(4))
        key = yr * 1000000 + mo * 10000 + da * 100 + hh

        have = lambda k: flat(k)[:nh] if k in f else np.zeros(nh)
        out = {"key": key,
               "LE": have("QE"), "H": have("H"),
               "GPP": have("An_H") + have("Rdark_H")}
        et = np.zeros(nh)
        for k in ET_PARTS:
            et = et + have(k)
        out["ET"] = et
    return out


def read_tower_hourly(root: Path, sid: str, gpp: str, max_qc: float,
                      lon: float | None = None):
    """Tower hourly series keyed by UTC hour, or (None, reason)."""
    # HALF-HOURLY *OR* HOURLY. The resolution varies by site and is carried in
    # the filename -- US-Ha1, US-MMS, US-Ho1 and US-Cwt publish HR, and a
    # HH-only glob dropped four of the longest records in the network. The
    # groupby below aggregates either one to whole hours unchanged: HH averages
    # two records per hour, HR passes one through.
    hits0 = sorted(Path(root).glob(f"**/*{sid}*FLUXMET_HH_*.csv"))
    hits0 += sorted(Path(root).glob(f"**/*{sid}*FLUXMET_HR_*.csv"))
    if not hits0:
        # No FLUXNET archive at all -- one of the 32 sites ONEFlux has not
        # processed. An expected absence, not a parsing failure.
        return None, f"{sid}: no FLUXNET archive (site has no ONEFlux product)"
    off, src = utc_offset(root, sid, lon)
    if off is None:
        return None, f"{sid}: no UTC offset from BIF or longitude"
    hits = hits0

    import pandas as pd
    key = f"GPP_{gpp}"
    want = {"TIMESTAMP_START"} | {c for v in (TOWER_COLS[key], TOWER_COLS["LE"],
                                              TOWER_COLS["H"]) for c in v}
    want |= set(QC_OF.values())
    d = pd.read_csv(hits[0], usecols=lambda c: c in want, na_values=[-9999])
    if "TIMESTAMP_START" not in d.columns:
        return None, f"{sid}: {hits[0].name} has no TIMESTAMP_START"

    ts = d["TIMESTAMP_START"].astype("Int64").astype(str).str.zfill(12)
    # Half-hourly -> hourly, and LOCAL STANDARD TIME -> UTC in one step.
    lst_hour = (ts.str[:4].astype(int) * 1000000 + ts.str[4:6].astype(int) * 10000
                + ts.str[6:8].astype(int) * 100 + ts.str[8:10].astype(int))
    d["_h"] = lst_hour
    for name, cands in (("GPP", TOWER_COLS[key]), ("LE", TOWER_COLS["LE"]),
                        ("H", TOWER_COLS["H"])):
        col = next((c for c in cands if c in d.columns), None)
        d[name] = np.nan if col is None else d[col]
        qc = QC_OF.get(name if name != "GPP" else key)
        if qc and qc in d.columns:
            d.loc[d[qc] > max_qc, name] = np.nan
    g = d.groupby("_h")[["GPP", "LE", "H"]].mean()

    # Shift to UTC by whole hours. datetime arithmetic on a packed integer is
    # wrong across day and month boundaries, so decode, shift, re-encode.
    import datetime as dt
    keys = g.index.to_numpy()
    y, rem = np.divmod(keys, 1000000)
    mth, rem = np.divmod(rem, 10000)
    day, hr = np.divmod(rem, 100)
    utc = np.array([
        int((dt.datetime(int(a), int(b), int(c), int(e))
             + dt.timedelta(hours=off)).strftime("%Y%m%d%H"))
        for a, b, c, e in zip(y, mth, day, hr)], dtype=np.int64)
    out = {"key": utc, "GPP": g["GPP"].to_numpy(float),
           "LE": g["LE"].to_numpy(float), "H": g["H"].to_numpy(float)}
    out["ET"] = out["LE"] * W_TO_MM_H
    out["offset_source"] = src
    return out, None


def stats(mod: np.ndarray, obs: np.ndarray) -> dict:
    m = np.isfinite(mod) & np.isfinite(obs)
    n = int(m.sum())
    if n < 3:
        return {"n": n}
    a, b = mod[m], obs[m]
    sdo = float(np.std(b, ddof=1))
    rmse = float(np.sqrt(np.mean((a - b) ** 2)))
    r = float(np.corrcoef(a, b)[0, 1]) if np.std(a) > 0 and sdo > 0 else np.nan
    return {"n": n, "mean_obs": float(b.mean()), "sd_obs": sdo,
            "mean_mod": float(a.mean()), "sd_mod": float(np.std(a, ddof=1)),
            "rmse": rmse, "rsr": rmse / sdo if sdo > 0 else np.nan,
            "bias": float(a.mean() - b.mean()), "r": r}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True, help="$MODEL_RUN")
    ap.add_argument("--tower-dir", type=Path, required=True)
    ap.add_argument("--pair", default="era5_land:*_ic")
    ap.add_argument("--gpp", default="NT", choices=["NT", "DT"])
    ap.add_argument("--max-qc", type=float, default=1.0)
    ap.add_argument("--stations", default=None)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)

    try:
        out = resolve_out(a.out or "hourly_stats.csv")
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    sites = read_sites()
    want = ([s.strip() for s in a.stations.split(",")] if a.stations
            else sorted(p.name for p in a.root.iterdir()
                        if p.is_dir() and next(p.glob("**/fixed_lma*"), None)))

    rows, skipped, n_src = [], [], {}
    for sid in want:
        lon = sites[sid][1] if sid in sites else None
        tw, why = read_tower_hourly(a.tower_dir, sid, a.gpp, a.max_qc, lon)
        if tw is None:
            skipped.append(why); continue
        pairs = list(find_pairs(a.root, sid, a.pair))
        if not pairs:
            skipped.append(f"{sid}: no {a.pair} pair"); continue
        label, fxp, dyp = pairs[0]
        if fxp is None or dyp is None:
            skipped.append(f"{sid}: one arm has no RES"); continue
        try:
            arms = {"fixed": read_model_hourly(fxp), "dyn": read_model_hourly(dyp)}
        except Exception as e:                                   # noqa: BLE001
            skipped.append(f"{sid}: {type(e).__name__}: {e}"); continue

        pft = sites[sid][2] if sid in sites else ""
        got = {}
        for arm, md in arms.items():
            idx = {k: i for i, k in enumerate(md["key"])}
            sel = np.array([idx.get(k, -1) for k in tw["key"]])
            ok = sel >= 0
            for var in ("GPP", "ET", "LE", "H"):
                mv = np.full(tw["key"].size, np.nan)
                mv[ok] = md[var][sel[ok]]
                st = stats(mv, tw[var])
                got[(arm, var)] = st
        for var in ("GPP", "ET", "LE", "H"):
            f_, d_ = got[("fixed", var)], got[("dyn", var)]
            ss = (1 - d_["rmse"] / f_["rmse"]) if (f_.get("rmse", 0) or 0) > 0 else np.nan
            for arm, st in (("fixed", f_), ("dyn", d_)):
                rows.append([sid, pft, var, arm, st.get("n", 0),
                             *[f"{st.get(k, float('nan')):.6g}" for k in
                               ("mean_obs", "sd_obs", "mean_mod", "sd_mod",
                                "rmse", "rsr", "bias", "r")],
                             f"{ss:.6g}"])
        n0 = got[("fixed", "LE")].get("n", 0)
        src = tw.get("offset_source", "?")
        n_src[src] = n_src.get(src, 0) + 1
        print(f"  {sid:<9}{pft:<10} n={n0:>7} matched hours   "
              f"utc offset from {src}", flush=True)

    if not rows:
        print("ERROR: nothing computed", file=sys.stderr)
        for w in skipped[:15]:
            print(f"  {w}", file=sys.stderr)
        return 1
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["station", "pft", "variable", "arm", "n", "mean_obs",
                    "sd_obs", "mean_mod", "sd_mod", "rmse", "rsr", "bias",
                    "r", "skill_score"])
        w.writerows(rows)
    print(f"\n{len(rows)} row(s) over {len(rows)//8} station(s) -> {out}")
    if skipped:
        print(f"\nSKIPPED {len(skipped)}:")
        for w in skipped[:15]:
            print(f"  ! {w}")
        if len(skipped) > 15:
            print(f"  ... and {len(skipped)-15} more")
    return 1 if skipped else 0


if __name__ == "__main__":
    sys.exit(main())
