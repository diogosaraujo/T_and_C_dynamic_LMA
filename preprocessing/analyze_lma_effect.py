#!/usr/bin/env python3
"""Fleet-wide assessment of the fixed-vs-dynamic LMA treatment.

Reads BOTH arms' RES_*.mat for every station, reduces each to annual series, and
writes a report plus the tables behind it.

Two stages, because the RES files are ~300 MB each and there are ~200 of them:

  EXTRACT   one station -> a small JSON of annual series in $TC_INPUT_DATA/lma_effect/
            Resumable and array-parallel; re-reading 60 GB happens once.
  REPORT    every cached JSON -> lma_effect_report.md + three CSVs.

    python analyze_lma_effect.py --station US-Ha2      # one station
    python analyze_lma_effect.py --index 7             # array form, 1-based
    python analyze_lma_effect.py --all                 # every station, serially
    python analyze_lma_effect.py --report              # synthesise the document

WHAT IT MEASURES, and why these and not the raw means.

The fixed arm uses the MEAN of each station's PLSR LMA series, so both arms share
a long-term mean by construction. Any metric based on the mean is guaranteed to
find nothing -- that is the design, not a result. The treatment lives entirely in
the interannual signal, so the report scores:

  diff_pct     mean difference. Expected ~0. Non-zero only through Jensen's
               inequality, since the flux response to SLA is nonlinear and
               mean(f(SLA)) != f(mean(SLA)). Reported to show it IS ~0.
  diff_sd_pct  sd of the year-by-year difference -- the actual size of the effect.
  r_sla        correlation of that difference with the SLA anomaly. Near 1 means
               the model is tracking the treatment rather than solver noise.
  slope        d(ln flux)/d(ln SLA). The transferable number: model sensitivity,
               independent of how variable a given station's LMA input happens to
               be. SLA is a prescribed input, not a measured predictor, so there
               is no regression dilution and the slope is unbiased even where the
               LMA series is noisy.
  var_ratio    sd(dynamic)/sd(fixed) for the flux itself. >1 dynamic LMA ADDS
               interannual variance, <1 it DAMPS it.
  r_clim       correlation of the fixed-arm (climate-driven) anomaly with the SLA
               anomaly. Explains var_ratio: independent -> variances add;
               anti-correlated -> LMA opposes the climate signal and damps it.

It also scores the LMA INPUT itself, because a flux response is only meaningful
if the input carries signal. At US-Ha2 the series is statistically white (no
trend, lag-1 autocorrelation -0.02, year-to-year jumps matching pure noise), and
a canopy whose needles persist 3-7 years cannot physically do that. Where
jump_ratio ~ 1 and trend_r2 ~ 0, the flux response is the model faithfully
tracking retrieval noise, and the variance-based numbers should not be read as
ecology. The sensitivity slopes stay valid.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA",
                                 "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
MODEL_RUN = Path(os.environ.get("MODEL_RUN", INPUT_ROOT.parent / "model_run"))
CACHE = INPUT_ROOT / "lma_effect"
SITE_LISTS = [REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
              REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv"]
ARMS = ["fixed_lma", "dyn_lma"]
F_C = 0.5                       # LMA is dry mass; 0.5 converts to gC
GROW = [5, 6, 7, 8, 9]          # growing-season months for the energy diagnostics

# Fluxes carried through to the report. Order is the order they appear.
FLUXES = ["GPP", "NPP", "ANPP", "LAI_mean", "LAI_max", "ET", "T", "EG", "EIn",
          "Lk", "QE", "H", "Rn", "Bowen", "Tfrac", "WUE"]
UNITS = {"GPP": "gC/m2/yr", "NPP": "gC/m2/yr", "ANPP": "gC/m2/yr",
         "LAI_mean": "-", "LAI_max": "-", "ET": "mm/yr", "T": "mm/yr",
         "EG": "mm/yr", "EIn": "mm/yr", "Lk": "mm/yr", "QE": "W/m2",
         "H": "W/m2", "Rn": "W/m2", "Bowen": "-", "Tfrac": "-", "WUE": "gC/mm"}


# --------------------------------------------------------------- site metadata
def site_table() -> dict[str, dict]:
    out = {}
    for p in SITE_LISTS:
        if not p.is_file():
            print(f"  ! site list missing: {p}", file=sys.stderr)
            continue
        for r in csv.DictReader(open(p, newline="", encoding="utf-8-sig")):
            sid = (r.get("StationID") or "").strip()
            if not sid:
                continue
            out[sid] = {"forest_type": (r.get("ForestType") or "").strip().lower(),
                        "eco_code": (r.get("US_L3CODE") or "").strip(),
                        "eco_name": (r.get("US_L3NAME") or "").strip(),
                        "name": (r.get("StationName") or "").strip(),
                        "lat": r.get("Lat"), "lon": r.get("Lon")}
    return out


# ------------------------------------------------------------------- extraction
def _read(f, key):
    a = np.asarray(f[key][()], dtype=float)
    return a.ravel() if (a.ndim == 1 or 1 in a.shape) else a


def annual_from_res(path: Path) -> dict:
    """Reduce one RES file to annual series. Mirrors GRAPH_MOD's aggregation."""
    import h5py
    HOURLY = ["QE", "H", "Rn", "G", "T_H", "T_L", "EG", "EIn_H", "EIn_L", "EIn_urb",
              "EIn_rock", "ESN", "ESN_In", "EICE", "ELitter", "Lk", "Rh", "Rd",
              "Pr", "Ds", "Ta"]
    DAILY = ["LAI_H", "NPP_H", "RA_H", "ANPP_H"]
    with h5py.File(path, "r") as f:
        dm = np.asarray(f["Datam"][()], dtype=float)
        if dm.shape[0] != 4:
            dm = dm.T
        h = {k: _read(f, k) for k in HOURLY if k in f}
        d = {k: _read(f, k) for k in DAILY if k in f}
        # 'f' is infiltration; guard the name clash with the file handle
        infil = _read(f, "f") if "f" in f else None
    missing = [k for k in HOURLY + DAILY if k not in h and k not in d]
    n = len(h["QE"])
    yr_h, mo_h = dm[0][:n].astype(int), dm[1][:n].astype(int)
    nd = len(d["LAI_H"])
    yr_d = yr_h[::24][:nd]
    if len(yr_d) < nd:                       # daily array runs one step past the hours
        yr_d = np.concatenate([yr_d, np.repeat(yr_d[-1], nd - len(yr_d))])

    ET = sum(h[k] for k in ["T_H", "T_L", "EG", "EIn_H", "EIn_L", "EIn_urb",
                            "EIn_rock", "ESN", "ESN_In"])
    EIN = sum(h[k] for k in ["ELitter", "EIn_H", "EIn_L", "EIn_urb", "EIn_rock"])
    GPP_d = d["NPP_H"] + d["RA_H"]
    grow = np.isin(mo_h, GROW)
    day = grow & (h["Rn"] > 20.0)

    years = [int(y) for y in range(yr_h.min(), yr_h.max() + 1)
             if (yr_h == y).sum() > 350 * 24 and (yr_d == y).sum() > 350]
    o = {"years": years, "missing_vars": missing}
    kh = {y: yr_h == y for y in years}
    kd = {y: yr_d == y for y in years}

    def S(v):                                             # annual sum
        return [float(np.nansum(v[kh[y]])) for y in years]

    def M(v, mask=None):                                  # annual mean
        return [float(np.nanmean(v[kh[y] & mask])) if mask is not None
                else float(np.nanmean(v[kh[y]])) for y in years]

    o["ET"] = S(ET); o["T"] = S(h["T_H"] + h["T_L"]); o["EG"] = S(h["EG"])
    o["EIn"] = S(EIN); o["Lk"] = S(h["Lk"]); o["Rh"] = S(h["Rh"]); o["Rd"] = S(h["Rd"])
    o["Pr"] = S(h["Pr"])
    o["Inf"] = S(infil) if infil is not None else [float("nan")] * len(years)
    o["QE"] = M(h["QE"]); o["H"] = M(h["H"]); o["Rn"] = M(h["Rn"]); o["G"] = M(h["G"])
    o["Ta"] = M(h["Ta"]); o["VPD"] = M(h["Ds"])
    o["QE_gs"] = M(h["QE"], day); o["H_gs"] = M(h["H"], day); o["Rn_gs"] = M(h["Rn"], day)
    o["GPP"] = [float(np.nansum(GPP_d[kd[y]])) for y in years]
    o["NPP"] = [float(np.nansum(d["NPP_H"][kd[y]])) for y in years]
    o["ANPP"] = [float(np.nansum(d["ANPP_H"][kd[y]])) for y in years]
    o["LAI_mean"] = [float(np.nanmean(d["LAI_H"][kd[y]])) for y in years]
    o["LAI_max"] = [float(np.nanmax(d["LAI_H"][kd[y]])) for y in years]
    # Bowen from growing-season DAYTIME sums, not from the annual means. A 24-h
    # or all-season ratio is dominated by winter hours where LE -> 0 and the
    # ratio explodes; that version is not interpretable.
    bh = np.array([float(np.nansum(h["H"][kh[y] & day])) for y in years])
    bq = np.array([float(np.nansum(h["QE"][kh[y] & day])) for y in years])
    o["Bowen"] = list(np.where(bq != 0, bh / np.where(bq == 0, np.nan, bq), np.nan))
    o["Tfrac"] = list(np.array(o["T"]) / np.array(o["ET"]))
    o["WUE"] = list(np.array(o["GPP"]) / np.array(o["ET"]))
    return o


def lma_series(d: Path, station: str):
    """(years, LMA) for an arm directory. LMA_<ST>.mat is v7, not v7.3."""
    p = d / f"LMA_{station.replace('-', '_')}.mat"
    if not p.is_file():
        return None, None
    try:
        from scipy.io import loadmat
        m = loadmat(p)
        return (np.asarray(m["years"]).ravel().astype(int),
                np.asarray(m["LMA"], dtype=float).ravel())
    except (NotImplementedError, ValueError):
        import h5py
        with h5py.File(p, "r") as f:
            return (_read(f, "years").astype(int), _read(f, "LMA"))


def extract(station: str, meta: dict, force: bool = False) -> str:
    out = CACHE / f"{station}.json"
    if out.is_file() and not force:
        return "cached"
    rec = {"station": station, **meta, "arms": {}}
    for arm in ARMS:
        d = MODEL_RUN / station / "era5_land" / arm
        res = d / f"RES_{station.replace('-', '_')}.mat"
        if not res.is_file():
            rec["arms"][arm] = None
            continue
        rec["arms"][arm] = annual_from_res(res)
        rec["arms"][arm]["res_mtime"] = res.stat().st_mtime
    ly, lv = lma_series(MODEL_RUN / station / "era5_land" / "dyn_lma", station)
    if ly is not None:
        rec["lma_years"] = [int(x) for x in ly]
        rec["lma"] = [float(x) for x in lv]
    if rec["arms"]["fixed_lma"] is None or rec["arms"]["dyn_lma"] is None:
        rec["status"] = "incomplete"
    else:
        rec["status"] = "ok"
    CACHE.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec), encoding="utf-8")
    return rec["status"]


# ---------------------------------------------------------------- diagnostics
def input_quality(years, lma):
    """Does the LMA series carry signal, or is it year-to-year noise?

    jump_ratio compares the observed mean |year-to-year change| against what an
    independent series of the same variance would give (sqrt(2/pi)*sqrt(2)*sd).
    ~1.0 is indistinguishable from white noise; <1 means the series is smoother
    than noise, i.e. it has temporal structure.
    """
    y, v = np.asarray(years, float), np.asarray(lma, float)
    g = np.isfinite(v)
    y, v = y[g], v[g]
    if len(v) < 6:
        return {}
    sl, ic = np.polyfit(y, v, 1)
    res = v - (sl * y + ic)
    r2 = float(np.corrcoef(y, v)[0, 1] ** 2)
    ac = float(np.corrcoef(res[:-1], res[1:])[0, 1]) if len(res) > 3 else float("nan")
    sd = float(v.std(ddof=1))
    exp_jump = np.sqrt(2 / np.pi) * sd * np.sqrt(2)
    return {"n": int(len(v)), "lma_mean": float(v.mean()), "lma_sd": sd,
            "lma_cv": float(100 * sd / v.mean()),
            "sla_mean": float(np.mean(1 / (v * F_C))),
            "sla_cv": float(100 * np.std(1 / (v * F_C), ddof=1) / np.mean(1 / (v * F_C))),
            "trend_pct_record": float(100 * sl * len(v) / v.mean()), "trend_r2": r2,
            "detrended_cv": float(100 * res.std(ddof=1) / v.mean()), "lag1_ac": ac,
            "mean_jump": float(np.abs(np.diff(v)).mean()),
            "jump_ratio": float(np.abs(np.diff(v)).mean() / exp_jump) if exp_jump else float("nan")}


def _r(a, b):
    g = np.isfinite(a) & np.isfinite(b)
    if g.sum() < 4 or np.std(a[g]) == 0 or np.std(b[g]) == 0:
        return float("nan")
    return float(np.corrcoef(a[g], b[g])[0, 1])


def flux_metrics(rec) -> dict:
    Fx, Dy = rec["arms"]["fixed_lma"], rec["arms"]["dyn_lma"]
    yrs = [y for y in Fx["years"] if y in set(Dy["years"])]
    if len(yrs) < 6:
        return {}
    fi = [Fx["years"].index(y) for y in yrs]
    di = [Dy["years"].index(y) for y in yrs]
    lm = dict(zip(rec.get("lma_years", []), rec.get("lma", [])))
    sla = np.array([1 / (lm[y] * F_C) if y in lm else np.nan for y in yrs])
    if not np.isfinite(sla).sum() >= 6:
        return {}
    x = np.log(sla / np.nanmean(sla))                # SLA anomaly, log space
    sa = 100 * (sla / np.nanmean(sla) - 1)
    out = {"n_years": len(yrs), "year_first": yrs[0], "year_last": yrs[-1]}
    for k in FLUXES:
        if k not in Fx or k not in Dy:
            continue
        a = np.array(Fx[k], float)[fi]
        b = np.array(Dy[k], float)[di]
        ma = float(np.nanmean(a))
        if not np.isfinite(ma) or ma == 0:
            continue
        dp = 100 * (b - a) / ma                      # annual difference, % of mean
        m = {"fixed_mean": ma, "dyn_mean": float(np.nanmean(b)),
             "diff_pct": float(100 * (np.nanmean(b) / ma - 1)),
             "diff_sd_pct": float(np.nanstd(dp, ddof=1)),
             "diff_max_abs_pct": float(np.nanmax(np.abs(dp))),
             "r_sla": _r(dp, sa),
             "sd_fixed": float(np.nanstd(a, ddof=1)), "sd_dyn": float(np.nanstd(b, ddof=1)),
             "var_ratio": float(np.nanstd(b, ddof=1) / np.nanstd(a, ddof=1))
                          if np.nanstd(a, ddof=1) > 0 else float("nan"),
             "r_clim": _r(100 * (a / ma - 1), sa)}
        # sensitivity: only where both arms are strictly positive every year, so
        # the log is defined. Rh/Rd are ~0 at forest sites and drop out here.
        if np.all(np.isfinite(a)) and np.all(np.isfinite(b)) and a.min() > 0 and b.min() > 0:
            yv = np.log(b / a)
            g = np.isfinite(x) & np.isfinite(yv)
            if g.sum() >= 4 and np.std(x[g]) > 0:
                m["slope"] = float(np.polyfit(x[g], yv[g], 1)[0])
                m["slope_r"] = float(np.corrcoef(x[g], yv[g])[0, 1])
        out[k] = m
    return out


# -------------------------------------------------------------------- report
def fmt(v, nd=2, w=0):
    s = "n/a" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{nd}f}"
    return s.rjust(w) if w else s


def group_stat(rows, key, metric):
    v = np.array([r["flux"][key][metric] for r in rows
                  if key in r["flux"] and metric in r["flux"][key]
                  and np.isfinite(r["flux"][key][metric])], float)
    return v


def report(out_dir: Path):
    recs = []
    for p in sorted(CACHE.glob("*.json")):
        r = json.loads(p.read_text(encoding="utf-8"))
        if r.get("status") != "ok":
            recs.append(r); continue
        r["flux"] = flux_metrics(r)
        r["input"] = input_quality(r.get("lma_years", []), r.get("lma", []))
        if r["flux"]:
            recs.append(r)
        else:
            r["status"] = "too few overlapping years"; recs.append(r)
    ok = [r for r in recs if r.get("status") == "ok" and r.get("flux")]
    bad = [r for r in recs if r not in ok]
    types = sorted({r["forest_type"] for r in ok})

    L = []
    A = L.append
    A("# Fixed vs dynamic LMA: fleet assessment\n")
    A(f"Stations with both arms: **{len(ok)}**"
      + (f" — {', '.join(f'{t} {sum(1 for r in ok if r['forest_type']==t)}' for t in types)}"
         if types else "") + f". Excluded/incomplete: {len(bad)}.\n")
    if ok:
        yy = [r["flux"]["n_years"] for r in ok]
        A(f"Record length {min(yy)}–{max(yy)} years (median {int(np.median(yy))}).\n")
    A("\nThe fixed arm is driven by the **mean** of each station's PLSR LMA series, so both "
      "arms share a long-term mean by construction. A near-zero `diff_pct` is therefore the "
      "expected result and not evidence of a null effect — the treatment lives in the "
      "interannual signal, which is what everything below scores.\n")

    # ---- 1. is the input real?
    A("\n## 1. Does the LMA input carry signal?\n")
    A("A flux response is only meaningful if the driver is. `trend_r2` is the fraction of "
      "LMA variance in a linear trend; `lag1_ac` the autocorrelation of the detrended "
      "residual; `jump_ratio` the observed mean year-to-year change divided by what pure "
      "white noise of the same variance would give. **jump_ratio ≈ 1 with trend_r2 ≈ 0 and "
      "lag1_ac ≈ 0 is statistically indistinguishable from noise.** That matters most for "
      "evergreen, whose needles persist 3–7 years: canopy-mean LMA is a multi-year running "
      "average by construction and physically cannot jump year to year.\n")
    A("\n| forest type | n | LMA mean | LMA CV % | detrended CV % | trend r² | lag-1 AC | jump ratio |")
    A("|---|---|---|---|---|---|---|---|")
    for t in types:
        g = [r["input"] for r in ok if r["forest_type"] == t and r.get("input")]
        if not g:
            continue
        col = lambda k: np.array([x[k] for x in g if np.isfinite(x.get(k, np.nan))], float)
        A(f"| {t} | {len(g)} | {fmt(col('lma_mean').mean(),1)} | {fmt(col('lma_cv').mean(),1)} "
          f"| {fmt(col('detrended_cv').mean(),1)} | {fmt(col('trend_r2').mean(),2)} "
          f"| {fmt(col('lag1_ac').mean(),2)} | {fmt(col('jump_ratio').mean(),2)} |")
    nz = [r for r in ok if r.get("input", {}).get("jump_ratio", 0) > 0.9
          and r["input"].get("trend_r2", 1) < 0.1]
    A(f"\nStations whose LMA series is indistinguishable from white noise "
      f"(jump_ratio > 0.9 and trend_r² < 0.1): **{len(nz)} of {len(ok)}**"
      + (f" — {sum(1 for r in nz if r['forest_type']=='evergreen')} evergreen, "
         f"{sum(1 for r in nz if r['forest_type']=='deciduous')} deciduous." if nz else "."))
    A("\nAt those stations the flux response is the model faithfully tracking retrieval "
      "noise. The **variance** numbers in §3 should not be read as ecology there; the "
      "**sensitivity slopes** in §4 remain valid, because SLA is a prescribed model input "
      "rather than a measured predictor, so there is no regression dilution.\n")

    # ---- 2. mean effect
    A("\n## 2. Mean effect (expected ~0)\n")
    A("\n| flux | unit | " + " | ".join(f"{t} fixed | {t} dyn | {t} Δ%" for t in types) + " |")
    A("|---|---|" + "---|" * (3 * len(types)))
    for k in FLUXES:
        cells = []
        for t in types:
            g = [r for r in ok if r["forest_type"] == t]
            f_, d_, p_ = (group_stat(g, k, "fixed_mean"), group_stat(g, k, "dyn_mean"),
                          group_stat(g, k, "diff_pct"))
            nd = 3 if k in ("Bowen", "Tfrac", "WUE") else (2 if "LAI" in k else 1)
            cells += [fmt(f_.mean(), nd) if f_.size else "n/a",
                      fmt(d_.mean(), nd) if d_.size else "n/a",
                      (f"{p_.mean():+.3f}" if p_.size else "n/a")]
        A(f"| {k} | {UNITS.get(k,'')} | " + " | ".join(cells) + " |")

    # ---- 3. size of the interannual effect
    A("\n## 3. Size of the interannual effect\n")
    A("`diff_sd_pct` is the standard deviation of the year-by-year difference, as a "
      "percentage of the flux mean — the real magnitude of the treatment. `max` is the "
      "largest single-year excursion. `r_sla` says how much of it is the treatment rather "
      "than solver noise. `var_ratio` = sd(dynamic)/sd(fixed): above 1 dynamic LMA **adds** "
      "interannual variance, below 1 it **damps** it.\n")
    for t in types:
        g = [r for r in ok if r["forest_type"] == t]
        A(f"\n**{t}** (n={len(g)})\n")
        A("\n| flux | diff sd % | max \\|diff\\| % | r(SLA) | var ratio | r_clim |")
        A("|---|---|---|---|---|---|")
        for k in FLUXES:
            s, mx = group_stat(g, k, "diff_sd_pct"), group_stat(g, k, "diff_max_abs_pct")
            rs, vr = group_stat(g, k, "r_sla"), group_stat(g, k, "var_ratio")
            rc = group_stat(g, k, "r_clim")
            if not s.size:
                continue
            A(f"| {k} | {fmt(s.mean(),2)} | {fmt(mx.mean(),2)} | {fmt(rs.mean(),2)} "
              f"| {fmt(vr.mean(),3)} | {fmt(rc.mean(),2)} |")

    # ---- 4. sensitivity
    A("\n## 4. Sensitivity — d(ln flux)/d(ln SLA)\n")
    A("The transferable result: how strongly each flux responds per unit of LMA change, "
      "independent of how variable a given station's input happens to be. A LAI slope above "
      "1 is superlinear — leaf area feeds back into leaf carbon; below 1 means something is "
      "capping it, and the candidate is the leaf-carbon limit `B(1) < LtR*B(3)`.\n")
    A("\n| flux | " + " | ".join(f"{t} slope | {t} sd | {t} r" for t in types) + " |")
    A("|---|" + "---|" * (3 * len(types)))
    for k in FLUXES:
        cells = []
        for t in types:
            g = [r for r in ok if r["forest_type"] == t]
            s, rr = group_stat(g, k, "slope"), group_stat(g, k, "slope_r")
            cells += [fmt(s.mean(), 3) if s.size else "n/a",
                      fmt(s.std(), 3) if s.size else "n/a",
                      fmt(rr.mean(), 2) if rr.size else "n/a"]
        if any(c != "n/a" for c in cells):
            A(f"| {k} | " + " | ".join(cells) + " |")
    A("\nWater-vs-carbon asymmetry (transpiration slope / GPP slope):\n")
    for t in types:
        g = [r for r in ok if r["forest_type"] == t]
        tt, gg = group_stat(g, "T", "slope"), group_stat(g, "GPP", "slope")
        if tt.size and gg.size and gg.mean() != 0:
            A(f"  * **{t}**: {tt.mean()/gg.mean():.2f}×")

    # ---- 5. by ecoregion
    A("\n## 5. By ecoregion\n")
    A("\n| ecoregion | type | n | GPP diff sd % | ET diff sd % | GPP slope | T slope | LAI slope |")
    A("|---|---|---|---|---|---|---|---|")
    seen = {}
    for r in ok:
        seen.setdefault((r["eco_name"], r["forest_type"]), []).append(r)
    for (eco, t), g in sorted(seen.items()):
        A(f"| {eco or '?'} | {t} | {len(g)} "
          f"| {fmt(group_stat(g,'GPP','diff_sd_pct').mean(),2)} "
          f"| {fmt(group_stat(g,'ET','diff_sd_pct').mean(),2)} "
          f"| {fmt(group_stat(g,'GPP','slope').mean(),3)} "
          f"| {fmt(group_stat(g,'T','slope').mean(),3)} "
          f"| {fmt(group_stat(g,'LAI_mean','slope').mean(),3)} |")

    # ---- 6. extremes
    A("\n## 6. Stations at the extremes\n")
    for k in ["GPP", "ET"]:
        s = sorted((r for r in ok if k in r["flux"]),
                   key=lambda r: r["flux"][k]["diff_sd_pct"], reverse=True)
        A(f"\n**Largest {k} response**\n")
        A("\n| station | type | ecoregion | diff sd % | max % | slope | SLA CV % | jump ratio |")
        A("|---|---|---|---|---|---|---|---|")
        for r in s[:8]:
            m = r["flux"][k]
            A(f"| {r['station']} | {r['forest_type'][:4]} | {r['eco_name'][:28]} "
              f"| {fmt(m['diff_sd_pct'],2)} | {fmt(m['diff_max_abs_pct'],2)} "
              f"| {fmt(m.get('slope'),3)} | {fmt(r.get('input',{}).get('sla_cv'),1)} "
              f"| {fmt(r.get('input',{}).get('jump_ratio'),2)} |")

    # ---- 7. what did not run
    if bad:
        A("\n## 7. Stations not assessed\n")
        A("\n| station | type | reason |")
        A("|---|---|---|")
        for r in sorted(bad, key=lambda x: x["station"]):
            arms = r.get("arms", {})
            why = r.get("status", "?")
            if arms:
                miss = [a for a in ARMS if not arms.get(a)]
                if miss:
                    why = f"no RES for {', '.join(miss)}"
            A(f"| {r['station']} | {r.get('forest_type','')} | {why} |")

    A("\n## Caveats\n")
    A("* The mean effect is ~0 **by design**. Do not report it as a null result.\n"
      "* Where `jump_ratio` ≈ 1, the LMA series is noise and the variance numbers "
      "reflect that, not ecology. Sensitivity slopes are unaffected.\n"
      "* `Rh`/`Rd` runoff is ~0 mm/yr at free-draining forest sites and its percentage "
      "differences are meaningless; they are excluded from the sensitivity table.\n"
      "* Tower validation cannot separate the arms: the treatment is 2–7% of the "
      "model–observation error, and bias is identical by construction.\n")

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "lma_effect_report.md").write_text("\n".join(L), encoding="utf-8")

    # ---- machine-readable companions
    with open(out_dir / "lma_effect_metrics.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        cols = ["fixed_mean", "dyn_mean", "diff_pct", "diff_sd_pct", "diff_max_abs_pct",
                "r_sla", "sd_fixed", "sd_dyn", "var_ratio", "r_clim", "slope", "slope_r"]
        w.writerow(["StationID", "ForestType", "EcoCode", "EcoName", "n_years",
                    "flux", "unit"] + cols)
        for r in ok:
            for k in FLUXES:
                if k not in r["flux"]:
                    continue
                m = r["flux"][k]
                w.writerow([r["station"], r["forest_type"], r["eco_code"], r["eco_name"],
                            r["flux"]["n_years"], k, UNITS.get(k, "")]
                           + [f"{m[c]:.6g}" if isinstance(m.get(c), float)
                              and np.isfinite(m[c]) else m.get(c, "") for c in cols])
    with open(out_dir / "lma_effect_input_quality.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        keys = ["n", "lma_mean", "lma_sd", "lma_cv", "sla_mean", "sla_cv",
                "trend_pct_record", "trend_r2", "detrended_cv", "lag1_ac",
                "mean_jump", "jump_ratio"]
        w.writerow(["StationID", "ForestType", "EcoName"] + keys)
        for r in ok:
            i = r.get("input", {})
            w.writerow([r["station"], r["forest_type"], r["eco_name"]]
                       + [i.get(k, "") for k in keys])
    with open(out_dir / "lma_effect_annual.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["StationID", "ForestType", "year", "flux", "fixed", "dynamic",
                    "diff", "LMA", "SLA"])
        for r in ok:
            lm = dict(zip(r.get("lma_years", []), r.get("lma", [])))
            Fx, Dy = r["arms"]["fixed_lma"], r["arms"]["dyn_lma"]
            for k in FLUXES:
                if k not in Fx or k not in Dy:
                    continue
                for y in Fx["years"]:
                    if y not in Dy["years"]:
                        continue
                    a = Fx[k][Fx["years"].index(y)]
                    b = Dy[k][Dy["years"].index(y)]
                    L_ = lm.get(y)
                    w.writerow([r["station"], r["forest_type"], y, k,
                                f"{a:.6g}", f"{b:.6g}", f"{b-a:.6g}",
                                "" if L_ is None else f"{L_:.4f}",
                                "" if L_ is None else f"{1/(L_*F_C):.6f}"])

    print(f"\nreport  -> {out_dir/'lma_effect_report.md'}")
    print(f"metrics -> {out_dir/'lma_effect_metrics.csv'}")
    print(f"input   -> {out_dir/'lma_effect_input_quality.csv'}")
    print(f"annual  -> {out_dir/'lma_effect_annual.csv'}")
    print(f"\n{len(ok)} stations assessed, {len(bad)} not")
    return 0 if ok else 1


# ------------------------------------------------------------------------ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--station", action="append", help="extract these stations")
    ap.add_argument("--index", type=int, help="1-based index into the station list (array form)")
    ap.add_argument("--all", action="store_true", help="extract every station serially")
    ap.add_argument("--report", action="store_true", help="synthesise from the cache")
    ap.add_argument("--model-run", type=Path, default=None)
    ap.add_argument("--cache", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None, help="where the report is written")
    ap.add_argument("--force", action="store_true", help="re-extract even if cached")
    a = ap.parse_args(argv)

    global MODEL_RUN, CACHE
    if a.model_run:
        MODEL_RUN = a.model_run
    if a.cache:
        CACHE = a.cache
    out_dir = a.out or CACHE

    sites = site_table()
    have = sorted(s for s in sites if (MODEL_RUN / s).is_dir())
    print(f"model_run : {MODEL_RUN}")
    print(f"cache     : {CACHE}")
    print(f"site list : {len(sites)} stations, {len(have)} present in model_run\n")

    if a.report:
        return report(out_dir)

    if a.index is not None:
        if a.index < 1 or a.index > len(have):
            print(f"index {a.index} outside 1..{len(have)} -- nothing to do")
            return 0
        targets = [have[a.index - 1]]
    elif a.station:
        targets = a.station
    elif a.all:
        targets = have
    else:
        ap.error("give --station / --index / --all, or --report")

    rc = 0
    for i, s in enumerate(targets, 1):
        if s not in sites:
            print(f"  {s}: not in the site lists -- skipped"); continue
        try:
            st = extract(s, sites[s], force=a.force)
        except Exception as e:                                   # noqa: BLE001
            print(f"  [{i}/{len(targets)}] {s}: FAILED -- {type(e).__name__}: {e}")
            rc = 1
            continue
        print(f"  [{i}/{len(targets)}] {s}: {st}", flush=True)
    print(f"\ncache now holds {len(list(CACHE.glob('*.json')))} station file(s)")
    print("next: python analyze_lma_effect.py --report")
    return rc


if __name__ == "__main__":
    sys.exit(main())
