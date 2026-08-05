"""ERA5-Land variable registry for the T&C dynamic-LMA forcing pipeline.

Single source of truth for which ERA5-Land variables we pull, what they are called
in the returned netCDF, their native units, their time convention, and which T&C
forcing field they eventually feed (see CLAUDE.md section 4).

Nothing here converts units. Files are stored in the native ERA5-Land format and
units; conversion happens later, in the forcing-builder step.
"""

from __future__ import annotations

# CDS collection IDs.
TIMESERIES_DATASET = "reanalysis-era5-land-timeseries"
GRIDDED_DATASET = "reanalysis-era5-land"

# Time conventions, spelled out once so the metadata sidecars stay self-describing.
INSTANT = "instantaneous"
ACCUM = "accumulated_from_00utc"

ACCUMULATION_NOTE = (
    "ERA5-Land accumulated fields reset at 00 UTC each day: the value stored at hour H "
    "is the accumulation from 00 UTC up to H, NOT the flux over the preceding hour. "
    "De-accumulate with a within-day difference before use, and verify against a known "
    "day the first time a new download batch is processed."
)

# The time-series collection does NOT return one merged file. It splits the response
# into one netCDF per variable GROUP, delivered together in a zip. These are the group
# names the CDS uses in the archive member names
# (reanalysis-era5-land-timeseries-sfc-<group><random>.nc), so the download step can map
# each member to a predictable output filename.
GROUP_2M = "2m-temperature"
GROUP_PRESSURE_PRECIP = "pressure-precipitation"
GROUP_RADIATION = "radiation-heat"
GROUP_WIND = "wind"

# Ordered so requests and metadata always list variables the same way.
VARIABLES = [
    {
        "cds_name": "2m_temperature",
        "short_name": "t2m",
        "long_name": "2 metre temperature",
        "units": "K",
        "time_convention": INSTANT,
        "timeseries_group": GROUP_2M,
        "tc_field": "Ta",
        "notes": "Air temperature. T&C expects degrees Celsius (subtract 273.15).",
    },
    {
        "cds_name": "surface_pressure",
        "short_name": "sp",
        "long_name": "Surface pressure",
        "units": "Pa",
        "time_convention": INSTANT,
        "timeseries_group": GROUP_PRESSURE_PRECIP,
        "tc_field": "Pre",
        "notes": (
            "Surface pressure at the ERA5-Land grid-cell elevation, which may differ from "
            "the tower elevation in complex terrain. Confirm the unit T&C expects for Pre "
            "(Pa vs mbar) against the reference site forcing before converting."
        ),
    },
    {
        "cds_name": "2m_dewpoint_temperature",
        "short_name": "d2m",
        "long_name": "2 metre dewpoint temperature",
        "units": "K",
        "time_convention": INSTANT,
        "timeseries_group": GROUP_2M,
        "tc_field": "Tdew",
        "notes": (
            "Dewpoint temperature. Feeds ea via Tetens; also an input to the radiation "
            "partition routine. T&C expects degrees Celsius."
        ),
    },
    {
        "cds_name": "total_precipitation",
        "short_name": "tp",
        "long_name": "Total precipitation",
        "units": "m",
        "time_convention": ACCUM,
        "timeseries_group": GROUP_PRESSURE_PRECIP,
        "tc_field": "Pr",
        "notes": (
            "De-accumulate, then multiply by 1000 for mm/h. Also drives the N=1 "
            "overcast rule in the radiation partition."
        ),
    },
    {
        "cds_name": "10m_u_component_of_wind",
        "short_name": "u10",
        "long_name": "10 metre U wind component",
        "units": "m s**-1",
        "time_convention": INSTANT,
        "timeseries_group": GROUP_WIND,
        "tc_field": "Ws",
        "notes": "Eastward component. Ws = sqrt(u10^2 + v10^2).",
    },
    {
        "cds_name": "10m_v_component_of_wind",
        "short_name": "v10",
        "long_name": "10 metre V wind component",
        "units": "m s**-1",
        "time_convention": INSTANT,
        "timeseries_group": GROUP_WIND,
        "tc_field": "Ws",
        "notes": (
            "Northward component. Ws = sqrt(u10^2 + v10^2) at 10 m, whereas the T&C run "
            "reference height is zatm = 31 m in MOD_PARAM -- decide whether to apply a "
            "log-profile adjustment in the forcing builder."
        ),
    },
    {
        "cds_name": "surface_solar_radiation_downwards",
        "short_name": "ssrd",
        "long_name": "Surface solar radiation downwards",
        "units": "J m**-2",
        "time_convention": ACCUM,
        "timeseries_group": GROUP_RADIATION,
        "tc_field": "Rsw",
        "notes": (
            "De-accumulate, then divide by 3600 for W/m2. Feeds "
            "C_Automatic_Radiation_Partition to produce SAB1/SAB2/SAD1/SAD2, PARB/PARD "
            "and the derived cloudiness N. Already hourly, so the ERA5-Land path skips "
            "the daily-to-hourly shortwave disaggregation used for the GCM path."
        ),
    },
]

CDS_VARIABLE_NAMES = [v["cds_name"] for v in VARIABLES]

# Groups the time-series response will be split across, given the variables above.
# One output netCDF per group, per station. Sorted longest-first so prefix matching
# against archive member names never stops at a shorter group that is also a prefix.
TIMESERIES_GROUPS = sorted({v["timeseries_group"] for v in VARIABLES})
_GROUPS_BY_MATCH_ORDER = sorted(TIMESERIES_GROUPS, key=len, reverse=True)

# Archive members look like: reanalysis-era5-land-timeseries-sfc-<group><random>.nc
ARCHIVE_MEMBER_PREFIX = "reanalysis-era5-land-timeseries-sfc-"


def group_of_member(member_name: str) -> str | None:
    """Map a CDS archive member filename to one of TIMESERIES_GROUPS, or None."""
    stem = member_name.rsplit("/", 1)[-1]
    if stem.endswith(".nc"):
        stem = stem[:-3]
    if stem.startswith(ARCHIVE_MEMBER_PREFIX):
        stem = stem[len(ARCHIVE_MEMBER_PREFIX):]
    for group in _GROUPS_BY_MATCH_ORDER:
        if stem.startswith(group):
            return group
    return None

# Fields the download carries implicitly rather than as a requested variable.
COORDINATE_NOTES = {
    "time": {
        "tc_field": "Date",
        "units": "UTC",
        "notes": (
            "Hourly timestamps in UTC. Set DeltaGMT = 0 and keep the true station "
            "longitude for the solar geometry. For ERA5-Land the radiation timestamp "
            "convention corresponds to t_bef = 1, t_aft = 0 (CLAUDE.md section 4); this "
            "is a per-product constant, set once and reused across all sites."
        ),
    }
}

# Variables T&C does NOT need from ERA5-Land, recorded so they are not re-added by mistake.
DELIBERATELY_OMITTED = {
    "surface_thermal_radiation_downwards": (
        "Longwave is not a T&C input -- it is computed internally by "
        "Incoming_Longwave(Ta, ea, N). Downloading strd would be dead weight."
    ),
    "total_cloud_cover": (
        "N is derived from the clearness index inside the radiation partition, not "
        "downloaded. No ERA5 cloud product is needed."
    ),
    "relative_humidity": (
        "Not an ERA5-Land field, and not needed: ea comes from Tdew via Tetens, and "
        "U = ea/esat is never referenced by the run."
    ),
}
