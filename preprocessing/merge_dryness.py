#!/usr/bin/env python3
"""Concatenate the 16 dryness_part_*.csv into station_dryness.csv."""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from results_dir import resolve_out                              # noqa: E402

root = Path(resolve_out(".", create=False))
parts = sorted(root.glob("dryness_part_*.csv"))
if not parts:
    print("ERROR: no dryness_part_*.csv found", file=sys.stderr)
    raise SystemExit(1)
d = pd.concat([pd.read_csv(p, low_memory=False) for p in parts],
              ignore_index=True)
out = resolve_out("station_dryness.csv")
d.to_csv(out, index=False)
print(f"merged {len(parts)} parts -> {out}  ({len(d)} rows)")
print(d.groupby("dataset").agg(rows=("phi", "size"),
                               stations=("station", "nunique"),
                               gcms=("gcm", "nunique")).to_string())
med = d.groupby(["dataset", "station"], as_index=False)["phi"].median()
print("\nphi across stations:")
print(med.groupby("dataset")["phi"].describe()[["count", "min", "50%", "max"]]
         .to_string(float_format=lambda x: f"{x:.3f}"))
frac = float((med["phi"] > 1).mean())
print(f"\nfraction water-limited (phi > 1): {100*frac:.1f}%")
if frac in (0.0, 1.0):
    print("  WARNING: every station on one side of phi = 1 -- that is what a "
          "unit or longitude error looks like", file=sys.stderr)
