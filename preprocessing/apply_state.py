#!/usr/bin/env python3
"""Write a harvested initial state into run directories that already exist.

    python apply_state.py --root $MODEL_RUN --from '*/historical/*/*_lma' \
        --ic initial_state.csv --ic-key 'historical/{gcm}/spinup'

    python apply_state.py --root $MODEL_RUN --from '*/ssp*/*/*_lma' \
        --ic initial_state.csv --ic-key 'historical/{gcm}/{arm}'

WHY THIS EXISTS RATHER THAN RE-RUNNING THE BUILDER

Between a pre-spin-up build and a post-spin-up one, exactly four lines change:

    LAI_H(1,:)   B_H(1,:,:)   PHE_S_H(1,:)   AgeL_H(1,:)

Everything else in the arm directory is already right and already there --
MOD_PARAM's soil/root/canopy substitution, Sl_H (the GCM's own 1985-2014 mean,
which does not depend on the initial state), GO's ms and its '../Meteo_*.mat'
load, the MAIN_FRAME/MAIN_FRAME_SLA choice, and LMA_<ST>.mat carrying the yearly
SLA series for the dynamic arm.

build_gcm_model_run.py would regenerate all of that identically, but to do so it
needs two things that are NOT part of a run: the era5_land MOD_PARAM it patches
as a template, and the PLSR projection CSVs it rebuilds the LMA series from.
Requiring those on the run cluster would break the property the tree was
restructured to have -- that model_run is self-contained and a plain rsync of it
is everything a run needs. Job 60692281 is what that costs: 92 stations blocked
on "no era5_land MOD_PARAM" for a rebuild that would have changed four lines.

So: patch what changes, in place, from model_run alone.

IDEMPOTENT. The IC patterns match the line whether or not it already carries a
"%% spun-up IC from ..." comment, so re-applying with a different key simply
replaces the state. That is what makes the ssp round work: the same directories
are re-pointed from the spin-up state to each arm's own historical end state.

NO FALLBACKS. A directory with no MOD_PARAM, a station with no row for the
requested key, or an IC line that does not match exactly once is an error naming
the directory. Nothing is defaulted, and a run left on the template's initial
pools while every log says it was restarted is precisely the failure this whole
exercise exists to remove.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_model_run import apply_ic, mat_name, read_ic_table   # noqa: E402
from gcm_variables import SCENARIOS                             # noqa: E402


def lma_gap(d: Path, station: str, scenario: str):
    """Years the dynamic arm will simulate that LMA_<ST>.mat does not carry.

    MAIN_FRAME_SLA looks each simulated year up in that file and STOPS if it is
    missing, so a short series is a mid-run failure, not a degraded result. Job
    60693242 lost five US-MtB arms at 2011 because its projection ends in 2010.
    Checking here costs a millisecond; discovering it 70 years into an 86-year
    ssp arm costs the better part of an hour, times however many arms share the
    problem.

    Only the dynamic arm reads the series -- the fixed arm takes its Sl_H from
    MOD_PARAM -- but the pair is useless with one half missing, so the caller
    refuses both.
    """
    if scenario not in SCENARIOS:
        return None                       # era5_land, or a name we do not price
    f = d / f"LMA_{mat_name(station)}.mat"
    if not f.is_file():
        return f"no {f.name}"
    try:
        from scipy.io import loadmat
        yrs = {int(y) for y in loadmat(f)["years"].ravel()}
    except Exception as e:                # noqa: BLE001 -- report, never guess
        return f"cannot read {f.name}: {e}"
    lo, hi = SCENARIOS[scenario]
    gap = [y for y in range(lo, hi + 1) if y not in yrs]
    if gap:
        return (f"LMA series covers {len(yrs)}/{hi-lo+1} years of {scenario}; "
                f"missing {gap[0]}{'..' + str(gap[-1]) if len(gap) > 1 else ''}")
    return None


def parse_arm(rel: Path):
    """<ST>/<scenario>/<GCM>/<arm>  or  <ST>/era5_land/<arm>  ->  the four parts."""
    p = rel.parts
    if len(p) == 4:
        return dict(station=p[0], scenario=p[1], gcm=p[2], arm=p[3])
    if len(p) == 3:
        return dict(station=p[0], scenario=p[1], gcm="", arm=p[2])
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True, help="model_run root")
    ap.add_argument("--from", dest="pattern", required=True,
                    help="glob under the root selecting ARM directories, e.g. "
                         "'*/historical/*/*_lma'")
    ap.add_argument("--ic", type=Path, default=None,
                    help="initial_state.csv (default <root>/initial_state.csv)")
    ap.add_argument("--ic-key", required=True,
                    help="which harvested state each arm takes, as it appears in "
                         "the 'key' column. {station}/{scenario}/{gcm}/{arm} are "
                         "substituted.")
    ap.add_argument("--run-list", default=None,
                    help="write the patched arms to <root>/<name> as a run list "
                         "for submit_gcm_tc_run.sh")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    if not a.root.is_dir():
        print(f"ERROR: model_run root not found: {a.root}", file=sys.stderr)
        return 1
    ic_path = a.ic or (a.root / "initial_state.csv")
    table = read_ic_table(ic_path)

    dirs = sorted(p for p in a.root.glob(a.pattern) if p.is_dir())
    if not dirs:
        print(f"ERROR: '{a.pattern}' matched no directory under {a.root}. Nothing "
              f"was patched, which is not the same as nothing needing it.",
              file=sys.stderr)
        return 1
    print(f"model_run    : {a.root}\npattern      : {a.pattern}\n"
          f"initial state: {ic_path}  ({len(table)} row(s))\n"
          f"key          : {a.ic_key}\nmatched      : {len(dirs)}\n")

    # PASS 1 -- decide, write nothing. The pair rule below needs to know every
    # refusal before any file is touched.
    ok, bad = [], []
    for d in dirs:
        rel = d.relative_to(a.root)
        parts = parse_arm(rel)
        if parts is None:
            bad.append((rel, "path is not <ST>/<scen>/<GCM>/<arm> or <ST>/<scen>/<arm>"))
            continue
        mp = d / f"MOD_PARAM_{mat_name(parts['station'])}.m"
        if not mp.is_file():
            bad.append((rel, f"no {mp.name}"))
            continue
        if parts["arm"].startswith("dyn_lma"):
            why = lma_gap(d, parts["station"], parts["scenario"])
            if why:
                bad.append((rel, why))
                continue
        key = a.ic_key.format(**parts)
        rec = table.get((parts["station"], key))
        if rec is None:
            bad.append((rel, f"no harvested state '{key}' in {ic_path.name}"))
            continue
        ok.append((d, rel, parts, mp, rec))

    # THE PAIR RULE. The experiment measures dyn MINUS fixed, so a fixed arm
    # whose dynamic twin was refused measures nothing -- it would run, cost an
    # hour, and produce a result with no counterpart. US-MtB and US-SHC lost
    # their dynamic arms to short LMA series in jobs 60693242/60700026; without
    # this, their fixed arms would still have been built and launched.
    refused_dyn = {(rel.parent, rel.name[len("dyn_lma"):])
                   for rel, _ in bad if rel.name.startswith("dyn_lma")}
    keep = []
    for d, rel, parts, mp, rec in ok:
        if rel.name.startswith("fixed_lma"):
            twin = (rel.parent, rel.name[len("fixed_lma"):])
            if twin in refused_dyn:
                bad.append((rel, f"dynamic twin dyn_lma{twin[1]} was refused -- an "
                                 f"unpaired fixed arm measures nothing"))
                continue
        keep.append((d, rel, parts, mp, rec))

    # PASS 2 -- write.
    done, runs = 0, []
    for d, rel, parts, mp, rec in keep:
        try:
            txt = apply_ic(mp.read_text(encoding="utf-8"), rec, str(rel))
        except SystemExit as e:
            bad.append((rel, str(e)))
            continue
        if not a.dry_run:
            mp.write_text(txt, encoding="utf-8")
        done += 1
        runs.append(f"{parts['station']} {parts['scenario']} {parts['gcm']} "
                    f"{parts['arm']}" if parts["gcm"] else
                    f"{parts['station']} {parts['arm']}")

    if bad:
        print(f"REFUSED -- {len(bad)} directory(ies) could not be patched:")
        for rel, why in bad[:20]:
            print(f"  ! {str(rel):<44}{why}")
        if len(bad) > 20:
            print(f"  ... and {len(bad) - 20} more")
        print()

    if a.run_list and runs and not a.dry_run:
        (a.root / a.run_list).write_text("".join(r + "\n" for r in runs),
                                         encoding="utf-8")
        print(f"run list : {a.root / a.run_list}  ({len(runs)} arms)")
        print(f"next     : RUN_LIST={a.run_list} sbatch --array=1-{len(runs)}%NN "
              f"slurm/submit_gcm_tc_run.sh")
        print(f"           (on Amarel add -p main; MaxSubmitPU on 'main' is 500, "
              f"so chunk with OFFSET beyond that)")

    verb = "would be patched" if a.dry_run else "patched"
    print(f"\n{done}/{len(dirs)} arm(s) {verb}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
