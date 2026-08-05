#!/usr/bin/env python3
"""Download AmeriFlux BASE measurements + BADM metadata for the study's stations.

Pulls the BASE-BADM product from AmeriFlux Data Services for the stations in the
ecoregion pairing CSVs, unpacks each archive into a per-station directory, and records
what was found. Site-level metadata from the public site_display endpoint is saved
alongside, so coordinates/IGBP/elevation are available even for stations whose flux data
we cannot download.

Credentials: an AmeriFlux account is required for the data download (registration is
free at https://ameriflux.lbl.gov). Supply them via environment variables:

    export AMF_USER_ID=your_ameriflux_username
    export AMF_USER_EMAIL=you@example.edu

Downloading also means accepting the AmeriFlux data use policy, which is an explicit,
recorded act -- pass --agree-policy (or set AMF_AGREE_POLICY=1).

Examples:
    python download_ameriflux.py --dry-run                    # plan only, no credentials
    python download_ameriflux.py --metadata-only              # public site info, no account
    python download_ameriflux.py --stations US-HBK,US-Ha2 --agree-policy --is-test
    python download_ameriflux.py --agree-policy               # full run, all stations
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from ameriflux_api import (
    DATA_PRODUCT,
    DEFAULT_DESCRIPTION,
    DEFAULT_INTENDED_USE,
    ENDPOINTS,
    INTENDED_USE_CHOICES,
    POLICY_CCBY4,
    POLICY_LEGACY,
    classify_member,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SITE_LISTS = [
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
    REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv",
]
INPUT_ROOT = Path(os.environ.get("TC_INPUT_DATA", "/vol_efthymios/NFS07/dd1136/T_and_C/input_data"))
DEFAULT_OUT = INPUT_ROOT / "ameriflux"

TIMEOUT = 120


def log(msg: str = "") -> None:
    print(msg, flush=True)


# --------------------------------------------------------------------------------------
# HTTP helpers -- urllib only, so the download step adds no dependency beyond the stdlib.
# The service hands back ftp:// URLs, which requests cannot fetch but urllib can.
# --------------------------------------------------------------------------------------


def get_json(url: str):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
        return json.loads(resp.read().decode("utf-8"))


def post_json(url: str, payload: dict):
    body = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
        return json.loads(resp.read().decode("utf-8"))


def fetch_file(url: str, dest: Path) -> None:
    """Download to a .part file then rename, so an interrupted run leaves no half file."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=TIMEOUT) as resp, part.open("wb") as fh:
        while chunk := resp.read(1 << 20):
            fh.write(chunk)
    part.replace(dest)


# --------------------------------------------------------------------------------------
# Station list
# --------------------------------------------------------------------------------------


def read_stations(paths: list[Path], wanted: set[str] | None) -> list[dict]:
    stations: dict[str, dict] = {}
    for path in paths:
        if not path.exists():
            raise SystemExit(f"site list not found: {path}")
        with path.open(newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                sid = (row.get("StationID") or "").strip()
                if not sid or (wanted is not None and sid not in wanted):
                    continue
                stations.setdefault(sid, {
                    "station_id": sid,
                    "station_name": (row.get("StationName") or "").strip(),
                    "forest_type": (row.get("ForestType") or "").strip(),
                    "igbp": (row.get("IGBP") or "").strip(),
                    "us_l3name": (row.get("US_L3NAME") or "").strip(),
                    "lat": row.get("Lat", ""),
                    "lon": row.get("Lon", ""),
                })
    return sorted(stations.values(), key=lambda s: s["station_id"])


# --------------------------------------------------------------------------------------
# Public metadata
# --------------------------------------------------------------------------------------


def fetch_site_metadata(station_ids: set[str], out_dir: Path) -> dict:
    """Site info from the public endpoint. No account needed."""
    log("fetching public site metadata ...")
    records = get_json(ENDPOINTS["sitemap"])
    by_id = {r.get("SITE_ID"): r for r in records if isinstance(r, dict)}
    selected = {sid: by_id[sid] for sid in sorted(station_ids) if sid in by_id}

    missing = sorted(station_ids - set(selected))
    if missing:
        log(f"  ! not present in the AmeriFlux site registry: {missing}")

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "site_metadata.json").write_text(
        json.dumps(selected, indent=2), encoding="utf-8")

    # Flat CSV for quick inspection alongside the JSON.
    rows = []
    for sid, rec in selected.items():
        loc = rec.get("GRP_LOCATION") or {}
        clim = rec.get("GRP_CLIM_AVG") or {}
        rows.append({
            "SITE_ID": sid,
            "SITE_NAME": rec.get("SITE_NAME", ""),
            "IGBP": rec.get("IGBP", ""),
            "STATE": rec.get("STATE", ""),
            "LOCATION_LAT": loc.get("LOCATION_LAT", ""),
            "LOCATION_LONG": loc.get("LOCATION_LONG", ""),
            "LOCATION_ELEV": loc.get("LOCATION_ELEV", ""),
            "CLIMATE_KOEPPEN": clim.get("CLIMATE_KOEPPEN", ""),
            "MAT": clim.get("MAT", ""),
            "MAP": clim.get("MAP", ""),
            "TOWER_BEGAN": rec.get("TOWER_BEGAN", ""),
            "TOWER_END": rec.get("TOWER_END", ""),
            "URL_AMERIFLUX": rec.get("URL_AMERIFLUX", ""),
        })
    if rows:
        with (out_dir / "site_metadata.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    log(f"  + site metadata for {len(selected)}/{len(station_ids)} station(s)")
    return selected


def fetch_policy_map(station_ids: set[str]) -> dict[str, str]:
    """Which sites are shared under CC-BY-4.0; the rest fall under the legacy policy."""
    try:
        entries = get_json(ENDPOINTS["site_ccby4"])
    except Exception as exc:
        log(f"  ! could not read the CC-BY-4.0 site list ({exc}); assuming LEGACY for all")
        return {sid: POLICY_LEGACY for sid in station_ids}

    ccby4 = set()
    for entry in entries:
        if isinstance(entry, (list, tuple)) and entry:
            ccby4.add(str(entry[0]))
        elif isinstance(entry, dict):
            ccby4.add(str(entry.get("SITE_ID", "")))
        elif isinstance(entry, str):
            ccby4.add(entry)
    return {sid: (POLICY_CCBY4 if sid in ccby4 else POLICY_LEGACY) for sid in station_ids}


# --------------------------------------------------------------------------------------
# Data download
# --------------------------------------------------------------------------------------


def extract_urls(response) -> list[str]:
    """Pull download URLs out of the response without assuming one exact shape."""
    urls: list[str] = []
    payload = response.get("data_urls") if isinstance(response, dict) else None
    if payload is None and isinstance(response, dict):
        for key in ("data_url", "urls", "download_urls"):
            if key in response:
                payload = response[key]
                break
    if payload is None:
        return urls
    if isinstance(payload, (str, dict)):
        payload = [payload]
    for item in payload:
        if isinstance(item, str):
            urls.append(item)
        elif isinstance(item, dict):
            for key in ("url", "URL", "link"):
                if item.get(key):
                    urls.append(str(item[key]))
                    break
    return urls


def station_of(filename: str, station_ids: set[str]) -> str | None:
    """AmeriFlux archives are named AMF_<SITE_ID>_<PRODUCT>_... ."""
    parts = filename.split("_")
    for part in parts:
        if part in station_ids:
            return part
    return None


def unpack(archive: Path, dest: Path) -> list[dict]:
    """Unpack an archive, keeping every member and labelling what each one is."""
    found: list[dict] = []
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as zf:
        for info in zf.infolist():
            if info.is_dir():
                continue
            name = Path(info.filename).name
            target = dest / name
            with zf.open(info) as src, target.open("wb") as out:
                while chunk := src.read(1 << 20):
                    out.write(chunk)
            found.append({
                "file": name,
                "kind": classify_member(info.filename),
                "bytes": target.stat().st_size,
            })
    return found


def download_for_policy(policy: str, sids: list[str], out_dir: Path, args) -> list[Path]:
    payload = {
        "user_id": args.user_id,
        "user_email": args.user_email,
        "data_product": DATA_PRODUCT,
        "data_policy": policy,
        "site_ids": sids,
        "intended_use": args.intended_use,
        "description": args.description,
        "is_test": "true" if args.is_test else "",
    }
    log(f"  requesting {len(sids)} site(s) under {policy} ...")
    response = post_json(ENDPOINTS["data_download"], payload)

    urls = extract_urls(response)
    if not urls:
        log(f"  ! no download URLs returned for {policy}. Response keys: "
            f"{sorted(response) if isinstance(response, dict) else type(response).__name__}")
        (out_dir / f"_response_{policy.replace('.', '')}.json").write_text(
            json.dumps(response, indent=2), encoding="utf-8")
        return []

    archive_dir = out_dir / "_archives"
    archive_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for url in urls:
        name = url.split("?")[0].rstrip("/").split("/")[-1]
        dest = archive_dir / name
        if dest.exists() and dest.stat().st_size > 0 and not args.overwrite:
            log(f"  = {name} already downloaded, skipping")
            paths.append(dest)
            continue
        log(f"  > {name}")
        try:
            fetch_file(url, dest)
        except (urllib.error.URLError, OSError) as exc:
            log(f"  X failed to fetch {name}: {exc}")
            continue
        log(f"  + {name} ({dest.stat().st_size / 1e6:.1f} MB)")
        paths.append(dest)
    return paths


# --------------------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Download AmeriFlux BASE measurements + BADM metadata.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--site-list", type=Path, action="append", default=None,
                   help="station CSV (repeatable); defaults to the deciduous + evergreen lists")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--stations", default=None,
                   help="comma-separated StationIDs — use this for test runs")
    p.add_argument("--user-id", default=os.environ.get("AMF_USER_ID"))
    p.add_argument("--user-email", default=os.environ.get("AMF_USER_EMAIL"))
    p.add_argument("--intended-use", default=DEFAULT_INTENDED_USE, choices=INTENDED_USE_CHOICES)
    p.add_argument("--description", default=DEFAULT_DESCRIPTION)
    p.add_argument("--agree-policy", action="store_true",
                   default=os.environ.get("AMF_AGREE_POLICY") == "1",
                   help="accept the AmeriFlux data use policy (required to download)")
    p.add_argument("--is-test", action="store_true",
                   help="mark as a test download so site teams are not emailed; "
                        "use this while shaking out the pipeline")
    p.add_argument("--metadata-only", action="store_true",
                   help="fetch only the public site metadata (no account needed)")
    p.add_argument("--batch-size", type=int, default=25,
                   help="site IDs per download request; a failed batch costs only its "
                        "own sites, and completed archives are skipped on re-run")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args(argv)

    site_lists = args.site_list or DEFAULT_SITE_LISTS
    wanted = {s.strip() for s in args.stations.split(",") if s.strip()} if args.stations else None
    stations = read_stations(site_lists, wanted)
    if not stations:
        raise SystemExit("no stations selected")
    ids = {s["station_id"] for s in stations}
    if wanted:
        for miss in sorted(wanted - ids):
            log(f"  ! {miss} is not in the site lists")

    log(f"stations       : {len(stations)}")
    log(f"product        : {DATA_PRODUCT}")
    log(f"output         : {args.out}")
    log(f"intended use   : {args.intended_use}")
    log(f"test download  : {bool(args.is_test)}")
    log("")

    if args.dry_run:
        log("[dry-run] would fetch public site metadata for: " + ", ".join(sorted(ids)))
        if not args.metadata_only:
            log("[dry-run] would POST to " + ENDPOINTS["data_download"])
            log("[dry-run] payload (credentials omitted):")
            log(json.dumps({
                "data_product": DATA_PRODUCT,
                "data_policy": "<per-site: CCBY4.0 or LEGACY>",
                "site_ids": sorted(ids),
                "intended_use": args.intended_use,
                "description": args.description,
                "is_test": "true" if args.is_test else "",
            }, indent=2))
        return 0

    args.out.mkdir(parents=True, exist_ok=True)
    meta = fetch_site_metadata(ids, args.out)

    if args.metadata_only:
        log("\nmetadata-only run complete.")
        return 0

    if not args.user_id or not args.user_email:
        raise SystemExit(
            "AmeriFlux credentials missing. Register free at https://ameriflux.lbl.gov "
            "then set:\n    export AMF_USER_ID=<username>\n    export AMF_USER_EMAIL=<email>"
        )
    if not args.agree_policy:
        raise SystemExit(
            "Downloading requires accepting the AmeriFlux data use policy "
            "(https://ameriflux.lbl.gov/data/data-policy/). Re-run with --agree-policy "
            "once you have read it."
        )

    policies = fetch_policy_map(ids)
    by_policy: dict[str, list[str]] = {}
    for sid, pol in policies.items():
        by_policy.setdefault(pol, []).append(sid)
    for pol, sids in sorted(by_policy.items()):
        log(f"  {pol}: {len(sids)} site(s)")
    log("")

    # Requests are chunked rather than sent as one 118-site POST: the service does not
    # document a site_ids limit, and a batch that fails takes only its own sites with it.
    # Archives already on disk are skipped, so re-running resumes.
    archives: list[Path] = []
    failed_batches: list[tuple[str, list[str]]] = []
    for pol, sids in sorted(by_policy.items()):
        ordered = sorted(sids)
        batches = [ordered[i:i + args.batch_size]
                   for i in range(0, len(ordered), args.batch_size)]
        for n, batch in enumerate(batches, 1):
            if len(batches) > 1:
                log(f"  [{pol} batch {n}/{len(batches)}]")
            try:
                archives.extend(download_for_policy(pol, batch, args.out, args))
            except Exception as exc:
                log(f"  X {pol} batch {n} failed: {exc}")
                failed_batches.append((pol, batch))

    log("")
    rows = []
    unpacked_ids = set()
    for archive in archives:
        sid = station_of(archive.name, ids)
        dest = args.out / (sid if sid else "_network")
        try:
            members = unpack(archive, dest)
        except zipfile.BadZipFile as exc:
            log(f"  X {archive.name} is not a readable zip ({exc})")
            continue
        if sid:
            unpacked_ids.add(sid)
        kinds = sorted({m["kind"] for m in members})
        log(f"  + {archive.name} -> {dest.name}/ ({len(members)} files: {', '.join(kinds)})")
        for m in members:
            rows.append({"station_id": sid or "", "archive": archive.name, **m})

    if rows:
        with (args.out / "files_manifest.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    (args.out / "download_provenance.json").write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "data_product": DATA_PRODUCT,
        "intended_use": args.intended_use,
        "description": args.description,
        "is_test": bool(args.is_test),
        "requested_stations": sorted(ids),
        "stations_with_data": sorted(unpacked_ids),
        "site_metadata_found": sorted(meta),
        "policies": policies,
        "citation": (
            "Cite the individual site DOIs and the AmeriFlux data policy: "
            "https://ameriflux.lbl.gov/data/data-policy/"
        ),
    }, indent=2), encoding="utf-8")

    no_data = sorted(ids - unpacked_ids)
    log(f"\nstations with downloaded data: {len(unpacked_ids)}/{len(ids)}")
    if no_data:
        log(f"no data returned for: {no_data}")
        log("  (a site may have no shared BASE data, or a different data policy)")
    log(f"manifest: {args.out / 'files_manifest.csv'}")

    if failed_batches:
        log(f"\n{len(failed_batches)} batch(es) FAILED:")
        for pol, batch in failed_batches:
            log(f"  {pol}: {', '.join(batch)}")
        log("Re-run the same command -- archives already downloaded are skipped.")
        return 1

    log("\nNext: python inspect_ameriflux_badm.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
