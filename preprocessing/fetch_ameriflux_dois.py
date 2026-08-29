#!/usr/bin/env python3
"""Fetch each station's AmeriFlux data DOI(s) and cache them in the repo.

Every AmeriFlux site mints a DOI for its BASE product, and separately for the
ONEFlux-processed FLUXNET product where one exists. Both must be cited when both
are used, so both are collected; a site with only BASE simply has an empty
FLUXNET column.

WHY SCRAPE THE SITE PAGE. The public API (site_display/AmeriFlux, the sitemap
this project already uses) carries name, coordinates, IGBP and tower years but
NO DOI -- checked field by field against a response. The DOI appears only on the
site's own page, labelled "AmeriFlux BASE:" and "AmeriFlux FLUXNET:", each
followed by a doi.org link. Those two labels are what this matches, so a page
that changes layout yields nothing for that site and says so, rather than
silently pairing a station with whichever DOI happened to appear first.

THE PREFIX IS NOT ASSUMED. All AmeriFlux DOIs seen so far are 10.17190/AMF/...,
but the pattern here captures whatever follows "https://doi.org/" after each
label, so a site issued under a different prefix is recorded correctly instead of
being dropped by an over-tight regex.

CACHED IN THE REPO on purpose. It is ~90 HTTP requests, it changes only when a
site releases a new version, and the citation list must be reproducible months
later when a page may have moved. --refresh re-fetches; without it, stations
already in the CSV are skipped.

    python fetch_ameriflux_dois.py --stations US-Ha2,US-HBK
    python fetch_ameriflux_dois.py --site-lists --refresh
"""
from __future__ import annotations

import argparse
import csv
import html
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

PREPROC = Path(__file__).resolve().parent
REPO_ROOT = PREPROC.parent
SITE_LISTS = [REPO_ROOT / "T&C" / "dynamic_lma_test" / "deciduous_ameriflux.csv",
              REPO_ROOT / "T&C" / "dynamic_lma_test" / "evergreen_ameriflux.csv"]
DEFAULT_OUT = PREPROC / "ameriflux_dois.csv"
SITEINFO = "https://ameriflux.lbl.gov/sites/siteinfo/{sid}"

# The two labels as they appear in the rendered page, each followed by the link.
# Tags are stripped before matching so markup between label and href does not
# matter; [^<]{0,400} bounds the gap so a label cannot capture a DOI belonging to
# a later section.
LABELS = {"base": r"AmeriFlux\s+BASE\s*:", "fluxnet": r"AmeriFlux\s+FLUXNET\s*:"}
DOI_AFTER = r"[\s\S]{0,400}?https?://doi\.org/(10\.\d{4,9}/[^\s\"'<>)]+)"
FIELDS = ["station", "doi_base", "doi_fluxnet", "citation_base", "citation_fluxnet"]


def strip_tags(t: str) -> str:
    return html.unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", t)))


def scrape(sid: str, timeout: float = 60.0, tries: int = 3) -> dict:
    """{'doi_base':…, 'doi_fluxnet':…, 'citation_*':…} for one station."""
    url = SITEINFO.format(sid=sid)
    last = None
    for attempt in range(tries):
        try:
            rq = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            raw = urllib.request.urlopen(rq, timeout=timeout).read()
            break
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last = e
            time.sleep(2 * (attempt + 1))
    else:
        raise RuntimeError(f"{sid}: {type(last).__name__}: {last}")

    text = strip_tags(raw.decode("utf-8", "replace"))
    out = {"station": sid, "doi_base": "", "doi_fluxnet": "",
           "citation_base": "", "citation_fluxnet": ""}
    for key, label in LABELS.items():
        m = re.search(label + DOI_AFTER, text)
        if not m:
            continue
        out[f"doi_{key}"] = m.group(1).rstrip(".,;")
        # The formatted citation the site asks you to use follows the DOI as
        # "Citation: <authors> (<year>), …". Captured verbatim so the reference
        # list is the site's own wording, not a reconstruction from metadata.
        c = re.search(r"Citation:\s*(.{20,400}?)\s*(?:https?://doi\.org/|To cite|"
                      r"AmeriFlux (?:BASE|FLUXNET)\s*:|Find global)",
                      text[m.end():])
        if c:
            out[f"citation_{key}"] = c.group(1).strip()
    return out


def read_site_lists() -> list:
    ids = []
    for p in SITE_LISTS:
        if not p.exists():
            continue
        with p.open(newline="", encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                sid = (r.get("StationID") or "").strip()
                if sid:
                    ids.append(sid)
    return sorted(set(ids))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stations", default=None, help="comma-separated station IDs")
    ap.add_argument("--station-file", type=Path, default=None,
                    help="file with one station ID per line")
    ap.add_argument("--site-lists", action="store_true",
                    help="use every station in the deciduous/evergreen lists")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--refresh", action="store_true",
                    help="re-fetch stations already present in the cache")
    ap.add_argument("--sleep", type=float, default=0.4,
                    help="seconds between requests; be polite to the server")
    a = ap.parse_args(argv)

    if a.stations:
        want = [s.strip() for s in a.stations.split(",") if s.strip()]
    elif a.station_file:
        want = [l.strip() for l in a.station_file.read_text().splitlines() if l.strip()]
    elif a.site_lists:
        want = read_site_lists()
    else:
        print("ERROR: give --stations, --station-file or --site-lists",
              file=sys.stderr)
        return 1

    have = {}
    if a.out.is_file():
        with a.out.open(newline="", encoding="utf-8") as fh:
            have = {r["station"]: r for r in csv.DictReader(fh)}
    todo = want if a.refresh else [s for s in want if s not in have]
    print(f"stations: {len(want)}   cached: {len(want) - len(todo)}   "
          f"to fetch: {len(todo)}")

    failed = []
    for i, sid in enumerate(todo, 1):
        try:
            rec = scrape(sid)
        except RuntimeError as e:
            failed.append(str(e))
            print(f"  [{i:3}/{len(todo)}] {sid:9} FAILED")
            continue
        have[sid] = rec
        n = sum(bool(rec[f"doi_{k}"]) for k in ("base", "fluxnet"))
        print(f"  [{i:3}/{len(todo)}] {sid:9} {n} DOI(s)  "
              f"{rec['doi_base'] or '-':<22} {rec['doi_fluxnet'] or '-'}")
        time.sleep(a.sleep)

    with a.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for sid in sorted(have):
            w.writerow({k: have[sid].get(k, "") for k in FIELDS})
    print(f"\n-> {a.out}  ({len(have)} stations)")

    got = [s for s in want if have.get(s, {}).get("doi_base")]
    print(f"   BASE DOI:    {len(got)}/{len(want)}")
    print(f"   FLUXNET DOI: {sum(1 for s in want if have.get(s, {}).get('doi_fluxnet'))}"
          f"/{len(want)}")
    missing = [s for s in want if not have.get(s, {}).get("doi_base")]
    if missing:
        print(f"   NO BASE DOI: {', '.join(missing)}", file=sys.stderr)
    if failed:
        print(f"\n{len(failed)} request(s) failed:", file=sys.stderr)
        for f in failed:
            print(f"  {f}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
