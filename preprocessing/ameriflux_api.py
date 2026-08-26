"""AmeriFlux Data Services endpoints and the BADM -> T&C parameter wish list.

Endpoint URLs and the data_download request/response shape follow the `amerifluxr` R
package (chuhousen/amerifluxr), which is the reference client for these services.

Nothing here is guessed about payload *contents*: the download step unpacks whatever the
archive holds and reports it, and the BADM inspector discovers the variables actually
present rather than assuming a fixed schema.
"""

from __future__ import annotations

BASE = "https://amfcdn.lbl.gov/api/v1"

ENDPOINTS = {
    # Public, no credentials needed.
    "sitemap": f"{BASE}/site_display/AmeriFlux",
    "site_ccby4": f"{BASE}/site_availability/AmeriFlux/BIF/CCBY4.0",
    "data_year": f"{BASE}/data_availability/AmeriFlux",
    # Which sites actually HAVE a FLUXNET product. Discovered 2026-08-26 while
    # diagnosing job 38863: site_availability/AmeriFlux lists the real product
    # names -- BASE-BADM, BIF, FLUXNET -- and FLUXNET exists under CCBY4.0 only,
    # for 407 sites against BASE-BADM's 514. 86 of our 118 stations are in it.
    "site_fluxnet": f"{BASE}/site_availability/AmeriFlux/FLUXNET/CCBY4.0",
    "variables": f"{BASE}/fp_var?limits=True",
    # Requires an AmeriFlux account (user_id + user_email in the POST body).
    "data_download": f"{BASE}/data_download",
}

# BASE-BADM bundles the BASE measurement files and that site's BADM metadata in one
# archive. It is the half-hourly (or hourly) tower record as the site team submits it.
DATA_PRODUCT = "BASE-BADM"

# The ONEFlux-processed product, matching amerifluxr's amf_download_fluxnet(). This is
# what carries PARTITIONED carbon: BASE has FC (net CO2 flux) and no GPP column at all,
# verified against AMF_US-Ha2_BASE_HH_16-5.csv -- zero GPP or RECO fields. GPP is a
# derived quantity and only exists in this product.
DATA_PRODUCT_FLUXNET = "FLUXNET"

# FLUXNET REQUIRES A VARIANT; BASE-BADM DOES NOT. This is what job 39563 was
# missing: all four batches came back with a well-formed response whose
# data_urls list was empty, even though every one of the 86 sites requested was
# confirmed to have a FLUXNET product. The download endpoint will not build an
# archive without knowing which variant to package.
#
# FULLSET IS THE DEFAULT BECAUSE SUBSET RETURNS NOTHING. Tested 2026-08-26
# against US-HBK with identical credentials: SUBSET gave
# "number_of_sites_downloaded": 0 while FULLSET returned an 89.1 MB archive.
# A BASE-BADM control on the same site and credentials succeeded, which is what
# narrowed it from "our request is malformed" to "this variant is not published
# for AmeriFlux-FLUXNET sites". SUBSET is left selectable in case that changes.
#
# FULLSET is ~89 MB per site, so roughly 7.7 GB across our 86 sites -- large,
# but it is the only variant that exists.
DATA_VARIANTS = ["FULLSET", "SUBSET"]
DEFAULT_DATA_VARIANT = "FULLSET"

# FLUXNET IS CCBY4.0 ONLY. Requesting it for a LEGACY site returns nothing, so the
# station list must be filtered rather than split. Checked 2026-08-26 against
# ENDPOINTS["site_ccby4"]: 109 of our 118 stations qualify; US-Blk, US-CZ2, US-CZ3,
# US-CZ4, US-LPH, US-MRf, US-NR2, US-SB3 and US-WBW are LEGACY and cannot supply it.
#
# Policy eligibility is necessary, NOT sufficient -- a CCBY4.0 site still only has a
# FLUXNET product if ONEFlux has processed it. There is no working availability
# endpoint to pre-check that (data_availability/AmeriFlux now 404s), so the request
# itself is the test and any site that comes back empty must be reported by name.

# Sites are shared under one of two policies; requesting the wrong one returns nothing
# for that site, so the station list is split by policy before requesting.
POLICY_CCBY4 = "CCBY4.0"
POLICY_LEGACY = "LEGACY"

# Values the service accepts for intended_use.
INTENDED_USE_CHOICES = [
    "Research - Multi-site synthesis",
    "Research - Remote sensing",
    "Research - Land model/Earth system model",
    "Research - Other",
    "Education (Teacher or Student)",
    "Other",
]
DEFAULT_INTENDED_USE = "Research - Land model/Earth system model"
DEFAULT_DESCRIPTION = (
    "Driving and validating the Tethys-Chloris (T&C) ecohydrological model at forested "
    "AmeriFlux sites across CONUS ecoregions, to test how remote-sensing-derived dynamic "
    "leaf mass per area (LMA) alters simulated water, energy and carbon fluxes relative "
    "to a fixed-LMA baseline."
)

# ----------------------------------------------------------------------------------
# What we are hoping to find in BADM, expressed as patterns rather than exact names.
#
# BADM variable naming is not fully documented in one place and varies by site, so the
# inspector MATCHES on these patterns and also reports everything it did not match. That
# way a renamed or unexpected variable shows up as "other available" instead of being
# silently reported as absent.
#
# `tc_use` ties each item back to the parameters in CLAUDE.md sections 5-7.
# `fallback` records what happens when a site does not report it.
# ----------------------------------------------------------------------------------
PARAMETER_TARGETS = [
    {
        "key": "lma",
        "patterns": ["LMA"],
        "tc_use": "Sl_H = 1/(LMA x f_C) — GROUND TRUTH for the study's central variable",
        "fallback": "PLSR remote-sensing LMA only, with no in-situ check",
    },
    {
        "key": "sapwood_area",
        "patterns": ["GRP_SA", "SA_MAX", "SAPWOOD"],
        "tc_use": "Axyl_H — sapwood area per unit ground area",
        "fallback": "PFT-prescribed (15 cm2/m2 at US_xRM)",
    },
    {
        "key": "terrain",
        "patterns": ["TERRAIN", "ASPECT", "SURFACE_HOMOGENEITY"],
        "tc_use": "Slo_top / SvF — plot-scale runs assume flat; a sloping site breaks that",
        "fallback": "assume flat (Slo_top = 0), as at US_xRM",
    },
    {
        "key": "snow_cover",
        "patterns": ["SNOW_COVER"],
        "tc_use": "validation target for the snow module (TminS/TmaxS, sublimation)",
        "fallback": "no observational check on simulated snow duration",
    },
    {
        "key": "nep",
        "patterns": ["GRP_NEP", "NEP"],
        "tc_use": "reported net ecosystem production — independent flux validation target",
        "fallback": "tower BASE fluxes only",
    },
    {
        "key": "canopy_height",
        "patterns": ["HEIGHTC"],
        "tc_use": "hc — canopy height (High vegetation layer)",
        "fallback": "global canopy-height product (Potapov 2021 / Simard 2011)",
    },
    {
        "key": "lai",
        "patterns": ["LAI"],
        "tc_use": "validation target for LAI = Sl*B(1); also seeds/checks spin-up",
        "fallback": "spin-up equilibrium only, no observational check",
    },
    {
        # MUST come before "biomass": matching takes the first hit, and the generic
        # BIOMASS pattern would otherwise swallow every ROOT_BIOMASS_* variable and
        # report root biomass as absent while inflating above-ground coverage.
        "key": "root_biomass",
        "patterns": ["ROOT_BIOMASS"],
        "tc_use": "seeds the fine-root carbon pool B_H(3); PROFILE_MAX bounds sampling depth",
        "fallback": "spin-up only",
    },
    {
        "key": "biomass",
        "patterns": ["AG_BIOMASS", "BIOMASS", "AGB"],
        "tc_use": "seed/validate initial carbon pools B_H(1:8)",
        "fallback": "spin-up from a plausible guess; NBCD/GEDI for a cross-check",
    },
    {
        "key": "soil_texture",
        "patterns": ["SOIL_TEX"],
        "tc_use": "Psan / Pcla — Saxton & Rawls inputs, per layer",
        "fallback": "POLARIS (CONUS 30 m) / SSURGO / SoilGrids",
    },
    {
        "key": "soil_chem",
        "patterns": ["SOIL_CHEM"],
        "tc_use": "Porg — organic fraction (OM = SOC x 1.72)",
        "fallback": "POLARIS / SoilGrids SOC",
    },
    {
        "key": "soil_depth",
        "patterns": ["SOIL_DEPTH", "DEPTH_TO_BEDROCK", "BEDROCK"],
        "tc_use": "soil column depth — do NOT inherit US_xRM's 1 m for deep-rooted forest",
        "fallback": "Pelletier 2016 / Shangguan 2017 / SoilGrids BDTICM",
    },
    {
        # NOT just "ROOT": that matched ROOT_BIOMASS_* and ROOT_PROD_* and reported
        # 12/110 coverage for a variable BADM does not contain at all. Verified across
        # 110 stations: no rooting-DEPTH variable exists in the BADM schema.
        "key": "root_depth",
        "patterns": ["ROOT_DEPTH", "ROOTING_DEPTH", "ZR95"],
        "tc_use": "ZR95 — rooting depth (must stay <= soil column depth)",
        "fallback": "NOT IN BADM — Schenk & Jackson 2002 / Fan 2017 / PFT lookup, always",
    },
    {
        "key": "species",
        "patterns": ["SPP", "SPECIES"],
        "tc_use": "sanity-check the PFT choice and the deciduous/evergreen split",
        "fallback": "IGBP class alone",
    },
    {
        "key": "disturbance",
        "patterns": ["DOM_DIST", "DM_", "DISTURBANCE", "MGMT"],
        "tc_use": "flag sites whose fluxes reflect harvest/fire rather than climate",
        "fallback": "none — unflagged disturbance can contaminate the LMA signal",
    },
    {
        "key": "elevation",
        "patterns": ["LOCATION_ELEV"],
        "tc_use": "Zbas — site elevation (orthometric) for the radiation partition",
        "fallback": "Copernicus GLO-30",
    },
    {
        "key": "igbp",
        "patterns": ["IGBP"],
        "tc_use": "PFT selection and the phenology switch aSE (0 evergreen, 1 deciduous)",
        "fallback": "already in the ecoregion pairing CSVs",
    },
    {
        "key": "utc_offset",
        "patterns": ["UTC_OFFSET"],
        "tc_use": "align tower observations with the UTC ERA5-Land forcing for validation",
        "fallback": "infer from longitude (error-prone near boundaries)",
    },
]


def classify_member(name: str) -> str:
    """Label a file unpacked from an AmeriFlux archive.

    Deliberately permissive: unrecognised members are kept and reported as "other"
    rather than discarded, so a change in packaging shows up in the manifest instead of
    silently losing data.
    """
    stem = name.rsplit("/", 1)[-1].upper()
    if "BIF" in stem or "BADM" in stem:
        return "badm"
    if "BASE" in stem and stem.endswith(".CSV"):
        return "base"
    if stem.endswith(".PDF") or "README" in stem or stem.endswith(".TXT"):
        return "doc"
    return "other"
