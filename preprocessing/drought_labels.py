"""Which SPEI value labels a given model time step, at each frequency.

The accumulation window has to match the step being labelled, otherwise a
"drought month" is being judged by a year of antecedent conditions and the
composite means something else entirely. So:

    monthly   -> SPEI-3 for THAT month
    seasonal  -> SPEI-3 for the LAST month of the season (DJF->Feb, MAM->May,
                 JJA->Aug, SON->Nov). A 3-month accumulation ending in the
                 season's final month covers exactly that season.
    annual    -> SPEI-12 for SEPTEMBER, the end of the water year falling in
                 that calendar year. A 12-month accumulation ending in September
                 spans Oct-Sep, so it carries the snowpack and the growing
                 season that actually drove that year's fluxes -- a December
                 SPEI-12 would instead be dominated by the autumn after them.

SEASONAL YEARS FOLLOW THE FLUX TABLE. analyze_period_effect files DJF under the
year of its January, so DJF 2002 is Dec 2001 + Jan/Feb 2002 and its label is
SPEI-3 at February 2002 -- the same convention on both sides, which is what
makes the join correct.

ACCUMULATION OFFSETS ARE REAL. An N-month index has no value for the first N-1
months, and the stacks are trimmed rather than padded: SPEI-3 starts 1980-03 and
SPEI-12 starts 1980-12. A missing key is therefore normal near the start of the
record and is reported as absent, never as zero.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from era5_predictors import DEFAULT_ERA5_ROOT, Era5Monthly      # noqa: E402

INDEX_FOR = {"monthly": "SPEI3_ts", "seasonal": "SPEI3_ts", "annual": "SPEI12_ts"}
SEASON_END_MONTH = {"DJF": 2, "MAM": 5, "JJA": 8, "SON": 11}
WATER_YEAR_END = 9          # September closes the Oct-Sep water year


class NoLabel(Exception):
    """No SPEI value exists for this step, and why."""


def key_for(freq: str, year: int, period) -> int:
    """The year*100+month key whose SPEI labels this (year, period) step."""
    if freq == "monthly":
        m = int(period)
        if not 1 <= m <= 12:
            raise NoLabel(f"month {period} out of range")
        return year * 100 + m
    if freq == "seasonal":
        p = str(period)
        if p not in SEASON_END_MONTH:
            raise NoLabel(f"unknown season {period!r}")
        return year * 100 + SEASON_END_MONTH[p]
    if freq == "annual":
        return year * 100 + WATER_YEAR_END
    raise NoLabel(f"unknown freq {freq!r}")


def station_spei(sites: dict, freq: str, root: Path = DEFAULT_ERA5_ROOT) -> dict:
    """{(station, year, period): spei} for every step the index can label.

    sites is {station_id: (lat, lon)}. Stations whose pixel is all-NaN are
    omitted; the caller reports them rather than treating them as non-drought.
    """
    index = INDEX_FOR.get(freq)
    if index is None:
        raise NoLabel(f"unknown freq {freq!r}")
    store = Era5Monthly(root)
    out: dict = {}
    try:
        for sid, (lat, lon) in sorted(sites.items()):
            ser = store.pixel_series(lat, lon)
            vals = np.asarray(ser["si"][index], dtype=float)
            keys = np.asarray(ser["si_time"][index], dtype=int)
            if not np.isfinite(vals).any():
                continue
            by_key = {int(k): float(v) for k, v in zip(keys, vals)
                      if np.isfinite(v)}
            for k, v in by_key.items():
                y, m = divmod(k, 100)
                if freq == "monthly":
                    out[(sid, y, m)] = v
                elif freq == "seasonal":
                    for p, em in SEASON_END_MONTH.items():
                        if m == em:
                            out[(sid, y, p)] = v
                elif m == WATER_YEAR_END:
                    out[(sid, y, "ANN")] = v
    finally:
        store.close()
    return out


def is_drought(spei: float, threshold: float = -1.0) -> bool:
    return np.isfinite(spei) and spei <= threshold
