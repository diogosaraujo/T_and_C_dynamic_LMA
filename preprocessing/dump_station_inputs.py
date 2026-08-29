#!/usr/bin/env python3
"""Extract everything the station tables need, so they can be built off-cluster.

The authoritative inputs live on the cluster and nowhere else: MOD_PARAM in each
run directory, the forcing .mat beside it, the fetched root-depth table, and the
FLUXNET archives. Formatting a Word table is iterative and does not belong in a
SLURM job. So this runs once on the HPC and writes three small CSVs that carry
every number the tables need; build_station_table.py then works from those,
locally, as many times as the layout takes.

    station_inputs.csv   one row per station: every station-specific value,
                         both the FETCHED and the APPLIED rooting depth, the
                         forcing .mat's own Lat/Lon/Zbas, and the tower-overlap
                         counts that decide the flux-comparison subset
    mod_param_values.csv long: station x parameter x value, EVERY assignment in
                         MOD_PARAM. This is what lets the per-PFT table be built
                         by finding the parameters that are actually constant
                         within a vegetation type, rather than by assuming a list
    tower_overlap.csv    station x variable x n, the usable model/tower pairs

WHY THE FORCING .mat IS READ. Zbas is the elevation the radiation partition
actually used, and it is not necessarily the AmeriFlux registry value that
build_station_table currently reports, nor the number written in CLAUDE.md
(US-xRM: 2753 there against 2743 from the API). The .mat settles it, so the
table cites what the model ran with. Lat and Lon are read the same way for the
same reason. These files are '-v7.3', i.e. HDF5, so h5py reads them; scipy is
tried first in case an older run produced a v7 file.

WHY BOTH ROOTING DEPTHS. ZR95_H in MOD_PARAM is the Schenk & Jackson D95 CAPPED
at the soil column depth -- T&C aborts in Root_Fraction_General otherwise. The
cap bit at several stations, and a table that shows only the applied value hides
that. Both are written, and the table marks the capped ones.

THE FLUX-COMPARISON SUBSET IS COMPUTED, NOT ASSUMED. A station enters the
GPP/ET/LE/H comparison only if it has a FLUXNET archive AND enough annual steps
that overlap the model run. That test is figure_skill_maps' own read_model and
read_tower, imported rather than reimplemented, so the table's marked subset and
the skill figures' station set cannot disagree.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from results_dir import NoResultsDir, resolve_out                 # noqa: E402

PREPROC = Path(__file__).resolve().parent
REPO_ROOT = PREPROC.parent

# A MATLAB assignment whose right-hand side is a literal: a number, a bracketed
# list, NaN/Inf, or a simple arithmetic expression of those. Function calls and
# control flow do not match, which is the point -- they are not parameters.
ASSIGN = re.compile(
    r"^\s*([A-Za-z_]\w*)\s*(?:\(([^)]*)\))?\s*=\s*([^;]+);", re.M)
LITERAL = re.compile(r"^[\s\[\]0-9.,eE+\-*/()]*$|^\s*\[?\s*(?:NaN|Inf|-Inf)\s*\]?\s*$",
                     re.I)


def strip_comments(text: str) -> str:
    """Drop MATLAB comments. Lines are handled one at a time so a '%' inside a
    bracketed list on a later line cannot swallow an earlier assignment."""
    out = []
    for line in text.splitlines():
        i = line.find("%")
        out.append(line if i < 0 else line[:i])
    return "\n".join(out)


def parse_mod_param(path: Path) -> dict:
    """{parameter: value-as-written} for every literal assignment in the file.

    An indexed target keeps its index in the name only when it is not the
    trivial "(1,:)" the template uses for a single-PFT layer, so hc_H(1,:) and
    hc_H are the same parameter and do not split into two rows.
    """
    text = strip_comments(path.read_text(errors="replace"))
    out = {}
    for m in ASSIGN.finditer(text):
        name, idx, rhs = m.group(1), (m.group(2) or "").strip(), m.group(3).strip()
        if not LITERAL.match(rhs):
            continue
        if idx and idx not in ("1,:", "1,:,:", "1", ":"):
            name = f"{name}({idx})"
        val = re.sub(r"\s+", " ", rhs).strip()
        # Later assignments win: MOD_PARAM sometimes sets a default then
        # overrides it, and the override is what the run used.
        out[name] = val
    return out


def as_float(v: str):
    try:
        return float(v)
    except (TypeError, ValueError):
        pass
    s = str(v).strip().strip("[]").strip()
    if s.lower() in ("nan", "-nan"):
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return None


def vec(v: str):
    s = str(v).strip().strip("[]")
    parts = [p for p in re.split(r"[,\s]+", s) if p]
    try:
        return np.array([float(p) for p in parts], float)
    except ValueError:
        return None


def read_forcing(path: Path) -> dict:
    """Lat/Lon/Zbas and the timing constants, as the run actually received them."""
    want = ("Lat", "Lon", "Zbas", "DeltaGMT", "t_bef", "t_aft")
    out = {}
    try:
        import scipy.io as sio
        m = sio.loadmat(path, squeeze_me=True,
                        variable_names=want)          # v7 and earlier
        for k in want:
            if k in m:
                out[f"mat_{k}"] = float(np.asarray(m[k]).ravel()[0])
        return out
    except NotImplementedError:
        pass                                          # -v7.3, fall through
    except Exception as e:                            # noqa: BLE001
        print(f"    ! {path.name}: {type(e).__name__}: {e}", file=sys.stderr)
        return out
    try:
        import h5py
        with h5py.File(path, "r") as f:
            for k in want:
                if k in f:
                    out[f"mat_{k}"] = float(np.array(f[k]).ravel()[0])
    except Exception as e:                            # noqa: BLE001
        print(f"    ! {path.name}: {type(e).__name__}: {e}", file=sys.stderr)
    return out


def rank_dir(p: Path) -> tuple:
    """Prefer the ERA5-Land fixed arm: Sl_H differs between arms and periods."""
    s = str(p).lower().replace("\\", "/")
    return (0 if "era5" in s else 1, 0 if "fixed" in s else 1, s)


def find_runs(model_run: Path, stations: set) -> dict:
    """station -> the MOD_PARAM to read, one per station."""
    best = {}
    for f in Path(model_run).glob("*/**/MOD_PARAM_*.m"):
        st = f.parent
        while st.parent != Path(model_run) and st.parent != st:
            st = st.parent
        if st.name not in stations:
            continue
        if st.name not in best or rank_dir(f) < rank_dir(best[st.name]):
            best[st.name] = f
    return best


def fetched_zr95(path: Path | None) -> dict:
    """station -> the UNCAPPED Schenk & Jackson D95, to detect where it was cut."""
    if not path or not Path(path).is_file():
        print(f"  ! fetched root-depth table not found ({path}); capping "
              f"cannot be identified", file=sys.stderr)
        return {}
    out = {}
    with Path(path).open(newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            sid = (r.get("StationID") or r.get("station") or "").strip()
            for key in ("ZR95_H_mm", "ZR95_mm", "zr95_mm", "ZR95"):
                if r.get(key):
                    try:
                        out[sid] = float(r[key])
                    except ValueError:
                        pass
                    break
    print(f"  fetched root depth: {len(out)} station(s)")
    return out


def tower_overlap(model_dir: Path, tower_dir: Path | None, stations) -> pd.DataFrame:
    """Usable model/tower annual pairs per station and variable.

    Uses figure_skill_maps' own readers so this subset is exactly the one the
    skill figures drew, not a second definition of "has tower data".
    """
    if not tower_dir or not Path(tower_dir).is_dir():
        print(f"  ! tower dir not found ({tower_dir}); the flux-comparison "
              f"subset cannot be identified", file=sys.stderr)
        return pd.DataFrame(columns=["station", "variable", "n"])
    try:
        import figure_skill_maps as SK
        sites = SK.read_sites()
        model = SK.read_model(Path(model_dir) / "era5_annual.csv", "annual")
        tower = SK.read_tower(Path(tower_dir), "annual", sites)
    except Exception as e:                            # noqa: BLE001
        print(f"  ! tower overlap unavailable: {type(e).__name__}: {e}",
              file=sys.stderr)
        return pd.DataFrame(columns=["station", "variable", "n"])

    rows = []
    for (sid, year, period, var), (f, d) in model.items():
        if sid not in stations:
            continue
        obs = tower.get((sid, year, period, var))
        if obs is None or not np.isfinite(obs) or not np.isfinite(f):
            continue
        rows.append({"station": sid, "variable": var})
    if not rows:
        return pd.DataFrame(columns=["station", "variable", "n"])
    d = (pd.DataFrame(rows).groupby(["station", "variable"])
           .size().reset_index(name="n"))
    print(f"  tower overlap: {d['station'].nunique()} station(s) with any "
          f"usable pair")
    return d


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--model-run", type=Path, default=None)
    ap.add_argument("--tower-dir", type=Path, default=None)
    ap.add_argument("--root-depth", type=Path, default=None,
                    help="root_depth_schenk_jackson.csv, for the uncapped D95")
    ap.add_argument("--min-n", type=int, default=3,
                    help="annual pairs a station needs before it counts as "
                         "entering the flux comparison")
    ap.add_argument("--out-prefix", default="")
    a = ap.parse_args(argv)

    import os
    mr = a.model_run or os.environ.get("MODEL_RUN")
    if not mr or not Path(mr).is_dir():
        print(f"ERROR: MODEL_RUN='{mr}' is not a directory", file=sys.stderr)
        return 1
    try:
        results = Path(a.results or resolve_out(".", create=False))
    except NoResultsDir as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    p = results / "station_metrics.csv"
    if not p.is_file():
        print(f"ERROR: {p} not found -- run station_metrics.py first",
              file=sys.stderr)
        return 1
    m = pd.read_csv(p, low_memory=False)
    m = m[(m["freq"] == "annual") & (m["subset"] == "all")]
    per = {ds: set(g["station"].unique()) for ds, g in m.groupby("dataset")}
    era5 = sorted(per.get("era5", ()))
    gcm = set.intersection(*per.values()) if per else set()
    print(f"fleets: {len(era5)} ERA5, {len(gcm)} with the full GCM set")

    runs = find_runs(Path(mr), set(era5))
    print(f"MOD_PARAM: {len(runs)}/{len(era5)} stations")
    missing = [s for s in era5 if s not in runs]
    if missing:
        print(f"  ! no MOD_PARAM for: {', '.join(missing)}", file=sys.stderr)

    fetched = fetched_zr95(a.root_depth or
                           (Path(os.environ.get("TC_INPUT_DATA", "")) /
                            "root_depth" / "root_depth_schenk_jackson.csv"))

    rows, long_rows = [], []
    for sid in era5:
        rec = {"station": sid, "in_gcm": sid in gcm}
        f = runs.get(sid)
        if f is not None:
            rec["mod_param"] = str(f)
            rec["arm_is_era5"] = "era5" in str(f).lower()
            par = parse_mod_param(f)
            for k, v in par.items():
                long_rows.append({"station": sid, "parameter": k, "value": v})
            # The handful the table needs directly, resolved to numbers here so
            # the local build does no MATLAB parsing at all.
            zs = vec(par.get("Zs", ""))
            if zs is not None and zs.size > 1:
                rec["soil_depth_mm"] = float(zs[-1])
                rec["ms"] = int(zs.size - 1)
                th = np.diff(zs)
                for key, name in (("sand", "Psan_Zs"), ("clay", "Pcla_Zs"),
                                  ("org", "Porg_Zs")):
                    vv = vec(par.get(name, ""))
                    if vv is not None and vv.size:
                        w = th if th.size == vv.size else None
                        rec[f"{key}_pct"] = float(
                            100.0 * (np.average(vv, weights=w) if w is not None
                                     else np.mean(vv)))
            for key, name in (("hc_m", "hc_H"), ("zatm_m", "zatm"),
                              ("zr95_mm", "ZR95_H"), ("sl_h", "Sl_H"),
                              ("ase", "aSE_H"), ("kbot", "Kbot")):
                if name in par:
                    rec[key] = as_float(par[name])
            if rec.get("sl_h"):
                rec["lma_g_m2"] = 1.0 / (rec["sl_h"] * 0.5)
            rec["kbot_free"] = (rec.get("kbot") is None or
                                (isinstance(rec.get("kbot"), float) and
                                 np.isnan(rec["kbot"])))
            # The ERA5 layout puts Meteo_*.mat inside the arm directory; the GCM
            # layout puts one per model, a level above the two arms. Look in
            # both, nearest first, or every GCM-only station reports no Zbas.
            met = (sorted(f.parent.glob("Meteo_*.mat")) or
                   sorted(f.parent.parent.glob("Meteo_*.mat")))
            if met:
                rec.update(read_forcing(met[0]))
                rec["forcing_mat"] = met[0].name
        rec["zr95_fetched_mm"] = fetched.get(sid)
        # Capped where the uncapped D95 exceeds the column: that is the exact
        # condition build_model_run applies, so it is reproduced, not inferred
        # from a difference that rounding could fake.
        if rec.get("zr95_fetched_mm") and rec.get("soil_depth_mm"):
            rec["zr95_capped"] = bool(rec["zr95_fetched_mm"] > rec["soil_depth_mm"])
        rows.append(rec)

    ov = tower_overlap(results, a.tower_dir, set(era5))
    used = set()
    if not ov.empty:
        used = set(ov[ov["n"] >= a.min_n]["station"].unique())
    for rec in rows:
        rec["tower_used"] = rec["station"] in used
        if not ov.empty:
            g = ov[ov["station"] == rec["station"]]
            for _, r in g.iterrows():
                rec[f"tower_n_{r['variable']}"] = int(r["n"])

    d = pd.DataFrame(rows)
    out1 = resolve_out(f"{a.out_prefix}station_inputs.csv")
    d.to_csv(out1, index=False)
    print(f"\n-> {out1}  ({len(d)} stations, {len(d.columns)} columns)")

    L = pd.DataFrame(long_rows)
    out2 = resolve_out(f"{a.out_prefix}mod_param_values.csv")
    L.to_csv(out2, index=False)
    print(f"-> {out2}  ({len(L)} rows, "
          f"{L['parameter'].nunique() if len(L) else 0} distinct parameters)")

    out3 = resolve_out(f"{a.out_prefix}tower_overlap.csv")
    ov.to_csv(out3, index=False)
    print(f"-> {out3}  ({len(ov)} rows)")

    n_cap = int(d.get("zr95_capped", pd.Series(dtype=bool)).fillna(False).sum())
    print(f"\n   ZR95 capped at the column depth: {n_cap} station(s)")
    if n_cap:
        print("     " + ", ".join(d[d["zr95_capped"].fillna(False)]["station"]))
    print(f"   entering the flux comparison (n >= {a.min_n}): {len(used)}")
    non_era5 = d[d.get("arm_is_era5") == False]["station"].tolist()  # noqa: E712
    if non_era5:
        print(f"   ! read from a non-ERA5 arm (LMA is a GCM-period mean): "
              f"{', '.join(non_era5)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
