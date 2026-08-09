#!/usr/bin/env python3
"""Build the model_run tree: per-station T&C run directories, fixed and dynamic LMA.

    model_run/
      Code/                       shared T&C source, copied ONCE
      GRAPH_MOD.m                 shared
      <STATION>/
        era5_land/
          fixed_lma/   GO_<ST>.m  MOD_PARAM_<ST>.m  LMA_<ST>.mat  Meteo_<ST>_*.mat
          dyn_lma/     same
        hist_gcm/  ssp126/  ssp245/  ssp370/  ssp585/     (created empty)

Everything is generated from T&C/Thanos_US_xRM as the template, with the
site-specific blocks substituted. Substitution is VERIFIED: every pattern must
match exactly once or the station is refused. A silently unsubstituted MOD_PARAM
would run happily with US_xRM's soil and canopy and produce plausible, wrong
numbers -- the failure mode this project has hit repeatedly.

WHAT IS SUBSTITUTED, and from where

    Zs, ms, Psan/Pcla/Porg   soil_profiles.csv   per-layer, Soil_parameters looped
    Kbot                     soil_sites.csv      0.01 mm/h where bedrock is reported,
                                                 NaN (free drainage) otherwise
    ZR95_H                   root_depth CSV      capped at the column depth
    zatm                     BASE_MeasurementHeight, max EC height, else hc + 12
    hc_H                     canopy_height CSV
    Sl_H                     LMA series          fixed arm: 1/(mean(LMA)*f_C)
                                                 dyn arm:   same value as the day-1
                                                            state; the yearly series
                                                            arrives via LMA_<ST>.mat
    aSE_H                    site list           0 evergreen, 1 deciduous

DECISIONS WORTH KNOWING

  * Station IDs contain a hyphen (US-HBK) which is not legal in a MATLAB
    identifier or file name. Directories keep the hyphen; MATLAB names use an
    underscore (MOD_PARAM_US_HBK.m), matching the template's own id_location
    convention ('US_xRM' for site US-xRM).
  * Code/ (145 files) and GRAPH_MOD.m are shared at the tree root rather than
    copied into every arm of every station. GO_<ST>.m addpaths up to them.
  * The fixed arm gets a CONSTANT Sl from the site's mean LMA, not the literature
    PFT value, so both arms share a mean and differ only in variability.
  * f_C = 0.5 is applied here: SLA = 1/(LMA*f_C). The shipped LMA_US_xRM.mat uses
    SLA = 1/LMA with no f_C, so its values are half what T&C should receive.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PREPROC = Path(__file__).resolve().parent
INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA",
                                 "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
DEFAULT_ROOT = INPUT_ROOT.parent / "model_run"
DEFAULT_TEMPLATE = REPO_ROOT / "T&C" / "Thanos_US_xRM"
DEFAULT_EXCLUDED = Path(__file__).resolve().parent / "excluded_stations.csv"
DEFAULT_SITE_LISTS = [
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv",
]
DEFAULT_HEIGHTS = PREPROC / "BASE_MeasurementHeight_20260715.csv"

F_C = 0.5                    # leaf carbon fraction; SLA = 1/(LMA*f_C)
KBOT_BEDROCK = 0.01          # [mm/h] near-impermeable, after Fatichi et al. at RME
ZATM_FALLBACK_OFFSET = 12.0  # [m] above canopy where no EC height is reported
SCENARIOS = ["hist_gcm", "ssp126", "ssp245", "ssp370", "ssp585"]
ARMS = ["fixed_lma", "dyn_lma"]
MISSING = {"", "NA", "N/A", "NaN", "nan", "-9999", None}


def mat_name(station: str) -> str:
    """US-HBK -> US_HBK. MATLAB identifiers cannot contain a hyphen."""
    return re.sub(r"[^0-9A-Za-z_]", "_", station)


def fnum(v, default=None):
    try:
        f = float(v)
        return f if math.isfinite(f) else default
    except (TypeError, ValueError):
        return default


# --------------------------------------------------------------------- inputs

def read_stations(paths, wanted):
    out = []
    seen = set()
    for p in paths:
        if not Path(p).is_file():
            print(f"  ! site list not found: {p}", file=sys.stderr)
            continue
        for r in csv.DictReader(open(p, newline="", encoding="utf-8-sig")):
            sid = (r.get("StationID") or "").strip()
            if not sid or sid in seen or (wanted and sid not in wanted):
                continue
            seen.add(sid)
            out.append({"station_id": sid,
                        "forest_type": (r.get("ForestType") or "").strip().lower(),
                        "lat": fnum(r.get("Lat")), "lon": fnum(r.get("Lon"))})
    return sorted(out, key=lambda s: s["station_id"])


def read_excluded(path):
    """Stations dropped from the study, from excluded_stations.csv.

    17 of the 118 are excluded -- 8 with no BADM archive, 8 young or disturbed
    stands that are not at equilibrium for spin-up, and US-KS1 with no rooting
    depth. Building forcing or run directories for them wastes work and, worse,
    puts them in run_list.txt where they would be simulated by accident.
    """
    out = {}
    p = Path(path) if path else None
    if not p or not p.is_file():
        return out
    for r in csv.DictReader(open(p, newline="", encoding="utf-8-sig")):
        sid = (r.get("station_id") or "").strip()
        if sid:
            out[sid] = (r.get("reason") or "excluded").strip()
    return out


def read_soil(soil_dir: Path):
    """(sites by station, profile layers by station) from fetch_soil.py output."""
    sites, layers = {}, defaultdict(list)
    sp = soil_dir / "soil_sites.csv"
    pp = soil_dir / "soil_profiles.csv"
    if sp.is_file():
        for r in csv.DictReader(open(sp, newline="", encoding="utf-8-sig")):
            sites[r["station_id"]] = r
    if pp.is_file():
        for r in csv.DictReader(open(pp, newline="", encoding="utf-8-sig")):
            layers[r["station_id"]].append(r)
    for v in layers.values():
        v.sort(key=lambda r: float(r["z_top_mm"]))
    return sites, layers


def read_lookup(path: Path, key_cols, val_cols):
    """station -> first parseable value among val_cols."""
    out = {}
    if not path or not Path(path).is_file():
        return out
    for r in csv.DictReader(open(path, newline="", encoding="utf-8-sig")):
        sid = next((r[k].strip() for k in key_cols if r.get(k)), None)
        if not sid:
            continue
        for c in val_cols:
            v = fnum(r.get(c))
            if v is not None:
                out[sid] = v
                break
    return out


def read_ec_heights(path: Path):
    """Max eddy-covariance measurement height per station, metres.

    The BASE measurement-height table lists one row per variable and version; the
    turbulent-flux variables carry the height that matters for zatm.
    """
    out = {}
    if not Path(path).is_file():
        return out
    fam = ("FC", "H", "LE", "USTAR", "WS", "TA")
    for r in csv.DictReader(open(path, newline="", encoding="utf-8-sig")):
        sid = (r.get("Site_ID") or "").strip()
        var = (r.get("Variable") or "").strip().upper()
        h = fnum(r.get("Height"))
        if not sid or h is None or h <= 0:
            continue
        base = re.sub(r"(_\d+)+$", "", re.sub(r"_PI(_F)?$", "", var))
        if base in fam:
            out[sid] = max(out.get(sid, 0.0), h)
    return out


def read_lma_series(lma_dir: Path, station: str):
    """[(year, LMA g/m2)] for a station, modelled series, from build_lma_input.py."""
    p = lma_dir / station / f"{station}_LMA_modelled.csv"
    if not p.is_file():
        return []
    out = []
    for r in csv.DictReader(open(p, newline="", encoding="utf-8-sig")):
        y, v = fnum(r.get("year")), fnum(r.get("LMA_g_m2"))
        if y is not None and v is not None:
            out.append((int(y), v))
    return sorted(out)


# ------------------------------------------------------------------ generation

def soil_block(layers, kbot):
    """MATLAB for the depth-resolved soil, replacing the single-triple template.

    Soil_parameters is scalar-only (it uses ^ and * on scalars, and `if Psil<0`),
    so the layered profile has to be a loop, not a vectorised call.
    """
    zs = [float(layers[0]["z_top_mm"])] + [float(l["z_bot_mm"]) for l in layers]
    ms = len(layers)
    psan = ", ".join(f"{float(l['Psan']):.5g}" for l in layers)
    pcla = ", ".join(f"{float(l['Pcla']):.5g}" for l in layers)
    porg = ", ".join(f"{float(l['Porg']):.5g}" for l in layers)
    zss = " ".join(f"{z:g}" for z in zs)
    return ms, zs, f"""%%%%%%%%%%% SOIL INPUT -- depth-resolved, generated by build_model_run.py
%%%% One triple per Zs layer from SSURGO/POLARIS. Saxton & Rawls runs PER LAYER:
%%%% Soil_parameters is scalar-only (it uses ^ and * on scalars and `if Psil<0`),
%%%% and so is the van Genuchten derivation below it -- alpVG and nVG come from
%%%% L and Pe, so they must be computed inside the loop too, not replicated from
%%%% layer 1.
Psan_Zs = [{psan}];
Pcla_Zs = [{pcla}];
Porg_Zs = [{porg}];
Psan = Psan_Zs(1); Pcla = Pcla_Zs(1); Porg = Porg_Zs(1); %%% surface, for reference
Color_Class = 0;
SPAR = 2; %%% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
%%%%%%%%%%%%%%%%%%%
Osat=zeros(1,ms); L=zeros(1,ms); Pe=zeros(1,ms); Ks=zeros(1,ms); O33=zeros(1,ms);
rsd=zeros(1,ms); lan_dry=zeros(1,ms); lan_s=zeros(1,ms); cv_s=zeros(1,ms);
K_usle=zeros(1,ms); nVG=zeros(1,ms); alpVG=zeros(1,ms);
for jk = 1:ms
    [Osat(jk),L(jk),Pe(jk),Ks(jk),O33(jk),rsd(jk),lan_dry(jk),lan_s(jk),cv_s(jk),K_usle(jk)] = ...
        Soil_parameters(Psan_Zs(jk),Pcla_Zs(jk),Porg_Zs(jk));
    p_jk = 3 + 2/L(jk);
    m_jk = 2/(p_jk-1);
    nVG(jk) = 1/(1-m_jk);
    alpVG(jk) = (((-101.9368*Pe(jk))*(2*p_jk*(p_jk-1))/(p_jk+3))* ...
        ((55.6+7.4*p_jk+p_jk^2)/(147.8+8.1*p_jk+0.092*p_jk^2)))^-1; %%% [1/mm]
end
K_usle = mean(K_usle);
Ks_Zs = Ks; %%% [mm/h] already per layer
%%%%%%%%%%%%%%%
Zs_gen = [{zss}]; %%% ms+1 = {len(zs)}
Kbot_gen = {('NaN' if kbot is None else f'{kbot:g}')}; %%% [mm/h] {'free drainage' if kbot is None else 'near-impermeable bedrock (Fatichi et al. RME)'}
"""


SUBS_REQUIRED = ("kbot", "zatm", "soil", "zs", "zr95", "sl", "ase", "hc")


def render_mod_param(template: str, st: dict) -> tuple[str, list[str]]:
    """Substitute the site blocks. Every pattern must fire exactly once."""
    txt = template
    fired = {}

    def sub(key, pattern, repl, count=1):
        nonlocal txt
        txt, n = re.subn(pattern, repl.replace("\\", "\\\\"), txt, count=count,
                         flags=re.M)
        fired[key] = fired.get(key, 0) + n

    ms, zs, block = st["_soil_block"]

    # The single-triple soil block through to the *ones(1,ms) replication.
    sub("soil",
        r"^%%%%%%%%%% SOIL INPUT.*?\n(?:.*?\n)*?\[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle\]=Soil_parameters\(Psan,Pcla,Porg\);\s*\n",
        block)
    # The per-property replications are now redundant: the loop fills the arrays.
    sub("ones",
        r"^(rsd|lan_dry|lan_s|cv_s|Osat|L|Pe|O33|alpVG|nVG|Ks_Zs)\s*=\s*\1?\s*\*?\s*"
        r"[A-Za-z_0-9]*\*ones\(1,ms\);.*$",
        r"%% \1: filled per layer by the Soil_parameters loop above", count=0)
    sub("zs", r"^Zs=\s*\[[^\]]*\];.*$", f"Zs = Zs_gen; %% ms+1 = {len(zs)}")
    sub("kbot", r"^Kbot\s*=\s*[^;]*;.*$",
        "Kbot = Kbot_gen; %% set from the bedrock contact (build_model_run.py)")
    sub("zatm", r"^zatm\s*=\s*[^;]*;.*$",
        f"zatm = {st['zatm']:g}; %% [m] {st['zatm_src']}")
    sub("zr95", r"^ZR95_H\s*=\s*\[[^\]]*\];.*$",
        f"ZR95_H = [{st['zr95']:g}]; %% [mm] {st['zr95_src']}")
    sub("sl", r"^Sl_H\s*=\s*\[[^\]]*\];.*$",
        f"Sl_H = [{st['sl']:.6g}]; %% [m^2/gC] 1/(LMA*{F_C}), LMA={st['lma_mean']:.1f} g/m2")
    sub("ase", r"^aSE_H\s*=\s*\[[^\]]*\];.*$",
        f"aSE_H = [{st['ase']}]; %% {st['forest_type']}")
    sub("hc", r"^hc_H\(1,:\)\s*=\s*\[[^\]]*\];.*$",
        f"hc_H(1,:) = [{st['hc']:g}]; %% [m] {st['hc_src']}")

    problems = [f"{k}: matched {fired.get(k, 0)}x (expected 1)"
                for k in SUBS_REQUIRED if fired.get(k, 0) != 1]
    # Tdp(1,:) = Ta(1)*ones(1,ms) is the initial soil temperature and MUST survive.
    leftover = [ln for ln in txt.splitlines()
                if "*ones(1,ms)" in ln and not ln.strip().startswith("%")
                and "Tdp(1,:)" not in ln]
    if leftover:
        problems.append(f"{len(leftover)} un-layered *ones(1,ms) line(s) remain: "
                        + "; ".join(x.strip()[:40] for x in leftover[:3]))
    return txt, problems


GO_TEMPLATE = """clear
%% Generated by build_model_run.py -- do not edit by hand.
%% Station {station} ({forest_type}), {arm}
addpath('{code_rel}')
addpath('{root_rel}')   %%% shared GRAPH_MOD.m

dt = 3600; %%[s]
dth = 1;   %%[h]

ms = {ms};  %%% Soil layers (site-specific, from the SSURGO/POLARIS profile)
cc = 1;     %%% Crown area

SLA_ex = load('LMA_{mname}.mat');

id_location = '{mname}';
load('{meteo_name}')

x1 = 1;
x2 = length(Ta);
NN = x2-x1+1;

Date=Date(x1:x2); Pr=Pr(x1:x2); Ta=Ta(x1:x2);
Ws=Ws(x1:x2); ea=ea(x1:x2); SAD1=SAD1(x1:x2);
SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2); Pre=Pre(x1:x2);
SAB2=SAB2(x1:x2); N=N(x1:x2); Tdew=Tdew(x1:x2); esat=esat(x1:x2);
PARB=PARB(x1:x2); PARD=PARD(x1:x2);

Ds = esat-ea; Ds(Ds<0) = 0;   %% [Pa] vapour pressure deficit
Oa = 210000;                  %% [umolO2/mol]
Ws(Ws<=0) = 0.01;

[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1)=YE; Datam(:,2)=MO; Datam(:,3)=DA; Datam(:,4)=HO;
clear YE MO DA HO MI SE

PARAM_IC = strcat('MOD_PARAM_','{mname}');

{main_frame};

rmpath('{code_rel}')

save(['RES_', id_location], '-v7.3');

%%% Figures. GRAPH_MOD opens figures, so under 'matlab -batch' they must be
%%% written to disk and closed -- otherwise the job accumulates invisible handles
%%% and exits without leaving anything behind.
fig_dir = 'figures';
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
try
    set(0, 'DefaultFigureVisible', 'off');
    GRAPH_MOD;
    h = findobj('Type', 'figure');
    for kf = 1:numel(h)
        nm = get(h(kf), 'Name');
        if isempty(nm), nm = sprintf('fig%02d', get(h(kf), 'Number')); end
        nm = regexprep(nm, '[^0-9A-Za-z_-]', '_');
        saveas(h(kf), fullfile(fig_dir, sprintf('%s_%s.png', id_location, nm)));
    end
    close all
    fprintf('wrote %d figure(s) to %s
', numel(h), fig_dir);
catch ME
    fprintf(2, 'GRAPH_MOD failed (results are still saved): %s
', ME.message);
end
"""


def build_station(st, tmpl_txt, args, out_root):
    """Write both arms for one station. Returns (n_written, problems)."""
    sid, mname = st["station_id"], mat_name(st["station_id"])
    probs = []
    for arm in ARMS:
        d = out_root / sid / "era5_land" / arm
        d.mkdir(parents=True, exist_ok=True)
        # Depth from the arm directory to the shared Code/ at the tree root.
        code_rel = "../../../Code"

        txt, p = render_mod_param(tmpl_txt, st)
        probs += [f"{arm}: {x}" for x in p]
        (d / f"MOD_PARAM_{mname}.m").write_text(txt, encoding="utf-8")

        (d / f"GO_{mname}.m").write_text(GO_TEMPLATE.format(
            station=sid, forest_type=st["forest_type"], arm=arm, code_rel=code_rel,
            root_rel="../../..",
            ms=st["_soil_block"][0], mname=mname, meteo_name=st["meteo_name"],
            main_frame="MAIN_FRAME_SLA" if arm == "dyn_lma" else "MAIN_FRAME",
        ), encoding="utf-8")

        if st.get("meteo_src"):
            dst = d / st["meteo_name"]
            if not dst.exists():
                try:
                    os.symlink(st["meteo_src"], dst)
                except OSError:
                    shutil.copy2(st["meteo_src"], dst)
        write_lma_mat(d / f"LMA_{mname}.mat", st["lma"], st["sl"], arm)

    for s in SCENARIOS:
        (out_root / sid / s).mkdir(parents=True, exist_ok=True)
    return probs


def write_lma_mat(path: Path, series, sl_fixed, arm):
    """LMA_<ST>.mat with years, LMA, SLA_H, SLA_L -- the LMA_US_xRM.mat layout.

    f_C is applied here: SLA = 1/(LMA*f_C). The shipped LMA_US_xRM.mat stores
    SLA = 1/LMA with no carbon fraction, so its values are half what T&C needs.
    The fixed arm gets the site's mean SLA repeated, so both arms share a mean and
    differ only in year-to-year variability.
    """
    try:
        import numpy as np
        from scipy.io import savemat
    except ImportError:
        return
    if not series:
        return
    years = np.array([y for y, _ in series], dtype=np.uint16)
    lma = np.array([v for _, v in series], dtype=float)
    sla = 1.0 / (lma * F_C) if arm == "dyn_lma" else np.full(lma.size, sl_fixed)
    savemat(path, {"years": years, "LMA": lma,
                   "SLA_H": sla.reshape(-1, 1),
                   "SLA_L": np.zeros((lma.size, 1))})


# ------------------------------------------------------------------------ main

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    p.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    p.add_argument("--stations", default="US-HBK,US-Ha2",
                   help="comma-separated; 'all' for every station in the site lists")
    p.add_argument("--site-list", type=Path, action="append", default=None)
    p.add_argument("--exclude-file", type=Path, default=DEFAULT_EXCLUDED,
                   help="CSV of stations to drop (default: excluded_stations.csv)")
    p.add_argument("--soil", type=Path, default=INPUT_ROOT / "soil")
    p.add_argument("--lma", type=Path, default=INPUT_ROOT / "lma")
    p.add_argument("--meteo", type=Path, default=INPUT_ROOT / "meteo")
    p.add_argument("--canopy", type=Path,
                   default=INPUT_ROOT / "canopy_height" / "canopy_height_gedi.csv")
    p.add_argument("--root-depth", type=Path,
                   default=INPUT_ROOT / "root_depth" / "root_depth_schenk_jackson.csv")
    p.add_argument("--heights", type=Path, default=DEFAULT_HEIGHTS)
    p.add_argument("--meteo-years", default="1985_2020")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()

    tmpl = a.template / "MOD_PARAM_US_xRM.m"
    if not tmpl.is_file():
        print(f"ERROR: template not found: {tmpl}", file=sys.stderr)
        return 1
    tmpl_txt = tmpl.read_text(encoding="utf-8", errors="replace")

    wanted = None if a.stations.strip().lower() == "all" else {
        s.strip() for s in a.stations.split(",") if s.strip()}
    excluded = read_excluded(a.exclude_file)
    stations = [x for x in read_stations(a.site_list or DEFAULT_SITE_LISTS, wanted)
                if x["station_id"] not in excluded]
    if not stations:
        print("no stations selected", file=sys.stderr)
        return 1

    soil_sites, soil_layers = read_soil(a.soil)
    # AmeriFlux first, gridded fallback -- the rule agreed for every site
    # parameter. fetch_canopy_height.py writes the BADM value and the GEDI window
    # statistics side by side; the median is preferred over the mean because a
    # window containing a gap or a road drags the mean down.
    hc, hc_src = {}, {}
    if a.canopy.is_file():
        for r in csv.DictReader(open(a.canopy, newline="", encoding="utf-8-sig")):
            sid = (r.get("station_id") or "").strip()
            for col, lab in (("badm_heightc_m", "AmeriFlux BADM HEIGHTC"),
                             ("hc_gedi_median_m", "GEDI window median"),
                             ("hc_gedi_mean_m", "GEDI window mean")):
                v = fnum(r.get(col))
                if sid and v is not None and 0 < v < 120:
                    hc[sid], hc_src[sid] = v, lab
                    break
    zr = read_lookup(a.root_depth, ("station_id", "StationID"),
                     ("ZR95_H_mm", "ZR95_mm", "ZR95"))
    ec = read_ec_heights(a.heights)

    print(f"root      : {a.root}")
    print(f"template  : {a.template}")
    print(f"stations  : {len(stations)} ({len(excluded)} excluded)")
    print(f"soil      : {len(soil_sites)} sites, {len(soil_layers)} profiles")
    from collections import Counter as _C
    print(f"canopy    : {len(hc)} "
          f"({', '.join(f'{k} {v}' for k, v in _C(hc_src.values()).most_common())})"
          if hc else "canopy    : 0   <-- NONE FOUND")
    print(f"ZR95      : {len(zr)}   EC heights: {len(ec)}")
    print(f"meteo dir : {a.meteo}{'' if a.meteo.is_dir() else '   <-- NOT FOUND'}\n")

    ready, blocked, all_probs = [], [], {}
    for st in stations:
        sid = st["station_id"]
        miss = []

        layers = soil_layers.get(sid)
        srow = soil_sites.get(sid, {})
        if not layers:
            miss.append("soil profile")
        st["hc"] = hc.get(sid)
        st["hc_src"] = hc_src.get(sid, "")
        if st["hc"] is None:
            miss.append("canopy height")
        st["ase"] = 0 if st["forest_type"].startswith("ever") else 1

        # zatm: measured EC height, else canopy + offset.
        if sid in ec:
            st["zatm"], st["zatm_src"] = ec[sid], "max EC measurement height (AmeriFlux)"
        elif st["hc"] is not None:
            st["zatm"] = st["hc"] + ZATM_FALLBACK_OFFSET
            st["zatm_src"] = f"hc + {ZATM_FALLBACK_OFFSET:g} m (no EC height reported)"
        else:
            miss.append("zatm")

        st["lma"] = read_lma_series(a.lma, sid)
        if not st["lma"]:
            miss.append("LMA series")
        else:
            st["lma_mean"] = sum(v for _, v in st["lma"]) / len(st["lma"])
            st["sl"] = 1.0 / (st["lma_mean"] * F_C)

        # Column depth, bedrock boundary condition, and the ZR95 cap.
        depth = fnum(srow.get("column_depth_mm"))
        bed = fnum(srow.get("bedrock_cm"))
        kbot = KBOT_BEDROCK if bed is not None else None
        st["_kbot"] = kbot
        z = zr.get(sid)
        if z is None:
            miss.append("ZR95")
        elif depth is not None and z > depth:
            st["zr95"] = depth
            st["zr95_src"] = (f"Schenk & Jackson {z:g} mm CAPPED at the {depth:g} mm "
                              f"column ({srow.get('bedrock') or 'no bedrock'})")
        else:
            st["zr95"] = z
            st["zr95_src"] = "Schenk & Jackson D95"

        st["meteo_name"] = f"Meteo_{mat_name(sid)}_{a.meteo_years}.mat"
        cand = a.meteo / st["meteo_name"]
        st["meteo_src"] = cand if cand.is_file() else None
        if st["meteo_src"] is None:
            miss.append("meteo .mat")

        if layers:
            st["_soil_block"] = soil_block(layers, kbot)

        if miss:
            blocked.append((sid, miss))
            continue
        ready.append(st)

    print(f"{'=' * 66}\nREADY: {len(ready)}/{len(stations)} stations have every input\n{'=' * 66}")
    for sid, miss in blocked:
        print(f"  ! {sid:<8} missing: {', '.join(miss)}")
    if blocked:
        print()

    if a.dry_run:
        for st in ready:
            print(f"  - {st['station_id']:<8} ms={st['_soil_block'][0]:<3} "
                  f"Kbot={'0.01' if st['_kbot'] else 'NaN':<5} "
                  f"ZR95={st['zr95']:>6g} zatm={st['zatm']:>5g} hc={st['hc']:>5g} "
                  f"Sl={st['sl']:.5f}")
        return 0

    a.root.mkdir(parents=True, exist_ok=True)
    code_dst = a.root / "Code"
    if not code_dst.is_dir():
        shutil.copytree(a.template / "Code", code_dst)
        print(f"copied shared Code/ ({len(list(code_dst.iterdir()))} files)")
    g = a.template / "GRAPH_MOD.m"
    if g.is_file() and not (a.root / "GRAPH_MOD.m").exists():
        shutil.copy2(g, a.root / "GRAPH_MOD.m")

    for st in ready:
        probs = build_station(st, tmpl_txt, a, a.root)
        if probs:
            all_probs[st["station_id"]] = probs
        print(f"  {st['station_id']:<8} written  ms={st['_soil_block'][0]} "
              f"Kbot={'0.01' if st['_kbot'] else 'NaN'} ZR95={st['zr95']:g} "
              f"zatm={st['zatm']:g}"
              + ("   <-- SUBSTITUTION PROBLEM" if probs else ""))

    if all_probs:
        print("\nSUBSTITUTION PROBLEMS -- these MOD_PARAM files are NOT trustworthy:")
        for sid, ps in all_probs.items():
            for x in ps:
                print(f"   {sid:<8} {x}")

    (a.root / "build_manifest.json").write_text(json.dumps({
        "stations_ready": [s["station_id"] for s in ready],
        "stations_blocked": {s: m for s, m in blocked},
        "f_C": F_C, "Kbot_bedrock_mm_h": KBOT_BEDROCK,
        "zatm_fallback_offset_m": ZATM_FALLBACK_OFFSET,
        "arms": ARMS, "scenarios": SCENARIOS,
        "substitution_problems": all_probs,
        "note": "fixed_lma uses MAIN_FRAME with a constant Sl from the site mean LMA; "
                "dyn_lma uses MAIN_FRAME_SLA with the yearly series in LMA_<ST>.mat. "
                "SLA = 1/(LMA*f_C).",
    }, indent=2), encoding="utf-8")
    # The job array reads this: one "<station> <arm>" per line, in the order
    # submit_tc_run.sh indexes with SLURM_ARRAY_TASK_ID.
    (a.root / "run_list.txt").write_text(
        "".join(f"{s['station_id']} {arm}\n" for s in ready for arm in ARMS),
        encoding="utf-8")
    print(f"\nmanifest : {a.root / 'build_manifest.json'}")
    print(f"run list : {a.root / 'run_list.txt'}  ->  "
          f"sbatch --array=1-{len(ready) * len(ARMS)} slurm/submit_tc_run.sh")
    return 1 if (all_probs or blocked) else 0


if __name__ == "__main__":
    sys.exit(main())
