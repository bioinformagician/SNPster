#!/usr/bin/env python3
"""
pgs_fetch_papers.py

Audit and download legally/publicly accessible papers referenced by a
PGS Catalog-style metadata workbook.

Primary strategy
----------------
1. Read the "Publications" sheet from the workbook.
2. Deduplicate works by DOI/PMID.
3. Batch-resolve DOI records in OpenAlex (up to 100 DOIs/request).
4. Record OA status, license, version, repository/preprint locations, and
   whether OpenAlex has a cached full-text PDF.
5. Optionally search OpenAlex by title for an open preprint when the journal
   version has no public copy.
6. Download PDFs from the supported OpenAlex content API.
7. Optionally fall back to a direct public OA PDF URL.
8. Resume safely: valid PDFs already present in --output-dir are detected and
   skipped before OpenAlex metadata/search/download requests are made.

The script never bypasses authentication/paywalls and does not use Sci-Hub.

Notes on rights
---------------
- "explicit_open_license" means OpenAlex reports a clear open license
  (e.g. CC BY / CC BY-NC / CC0 / public domain).
- "public_free_to_read_copy" means OpenAlex reports an OA/free-to-read copy,
  but the reuse/redistribution license is unclear. A local research copy may
  be accessible, but this script does not claim reuse rights.
- OpenAlex PDFs retain their original copyright. Downloading a copy does not
  grant additional redistribution/reuse rights.

Examples
--------
# Fast local inventory only; no network:
python pgs_fetch_papers.py pgs_all_metadata.xlsx --inventory-only

# Audit live OA availability without downloading:
python pgs_fetch_papers.py pgs_all_metadata.xlsx \
    --audit-only \
    --openalex-api-key "$OPENALEX_API_KEY"

# Download public copies cached by OpenAlex:
python pgs_fetch_papers.py pgs_all_metadata.xlsx \
    --download \
    --openalex-api-key "$OPENALEX_API_KEY"

# Also try direct public OA PDF URLs when OpenAlex has no cached PDF:
python pgs_fetch_papers.py pgs_all_metadata.xlsx \
    --download --direct-oa \
    --openalex-api-key "$OPENALEX_API_KEY"

# Only download works with an explicit open license:
python pgs_fetch_papers.py pgs_all_metadata.xlsx \
    --download --explicit-license-only \
    --openalex-api-key "$OPENALEX_API_KEY"
"""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import json
import os
import re
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Iterable, Optional
from urllib.parse import quote

import pandas as pd
import requests


OPENALEX_API = "https://api.openalex.org"
OPENALEX_CONTENT = "https://content.openalex.org"

# OpenAlex commonly uses these normalized license strings.
OPEN_LICENSE_PREFIXES = (
    "cc0",
    "cc-by",
    "public-domain",
)

PREPRINT_HOST_HINTS = (
    "biorxiv",
    "medrxiv",
    "arxiv",
    "research square",
    "ssrn",
    "osf preprints",
    "preprints.org",
)

DEFAULT_USER_AGENT = (
    "PGS-paper-fetcher/1.0 "
    "(research metadata audit; contact via configured email if supplied)"
)


@dataclass
class Publication:
    pgp_ids: str
    first_author: str
    title: str
    journal: str
    publication_date: str
    doi: str
    pmid: str


@dataclass
class AuditRow:
    pgp_ids: str
    first_author: str
    title: str
    journal: str
    publication_date: str
    doi: str
    pmid: str

    matched: bool = False
    matched_by: str = ""
    openalex_id: str = ""
    openalex_title: str = ""
    openalex_type: str = ""
    oa_status: str = ""
    is_oa: bool = False

    best_oa_license: str = ""
    best_oa_version: str = ""
    best_oa_source: str = ""
    best_oa_source_type: str = ""
    best_oa_pdf_url: str = ""
    best_oa_landing_url: str = ""

    has_cached_pdf: bool = False
    cached_pdf_url: str = ""

    preprint_found: bool = False
    preprint_openalex_id: str = ""
    preprint_title: str = ""
    preprint_pdf_url: str = ""
    preprint_license: str = ""

    explicit_open_license: bool = False
    public_fetchable: bool = False
    access_class: str = ""
    chosen_method: str = ""
    chosen_url: str = ""
    chosen_version: str = ""
    chosen_license: str = ""

    downloaded: bool = False
    skipped_existing: bool = False
    local_path: str = ""
    error: str = ""


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return str(value).strip()


def normalize_doi(value: Any) -> str:
    s = clean_text(value)
    if not s:
        return ""
    s = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", s, flags=re.I)
    s = re.sub(r"^doi:\s*", "", s, flags=re.I)
    return s.strip()


def normalize_pmid(value: Any) -> str:
    s = clean_text(value)
    if not s:
        return ""
    try:
        f = float(s)
        if f.is_integer():
            return str(int(f))
    except ValueError:
        pass
    return s


def normalize_title(title: str) -> str:
    title = title.casefold()
    title = re.sub(r"<[^>]+>", " ", title)
    title = re.sub(r"[^\w\s]", " ", title)
    title = re.sub(r"\s+", " ", title).strip()
    return title


def title_similarity(a: str, b: str) -> float:
    a2 = normalize_title(a)
    b2 = normalize_title(b)
    if not a2 or not b2:
        return 0.0
    return difflib.SequenceMatcher(None, a2, b2).ratio()


def first_author_surname(value: str) -> str:
    s = normalize_title(value)
    if not s:
        return ""
    # PGS sheets often store "Surname Initials".
    return s.split()[0]


def candidate_first_author_surname(work: dict[str, Any]) -> str:
    auth = work.get("authorships") or []
    if not auth:
        return ""
    name = clean_text((auth[0].get("author") or {}).get("display_name"))
    if not name:
        return ""
    return normalize_title(name).split()[-1]


def is_explicit_open_license(license_value: str) -> bool:
    s = clean_text(license_value).casefold()
    return any(s.startswith(prefix) for prefix in OPEN_LICENSE_PREFIXES)


def is_preprint_location(location: dict[str, Any]) -> bool:
    version = clean_text(location.get("version")).casefold()
    if version == "submittedversion":
        return True
    source = location.get("source") or {}
    source_name = clean_text(source.get("display_name")).casefold()
    return any(hint in source_name for hint in PREPRINT_HOST_HINTS)


def best_location(work: dict[str, Any]) -> dict[str, Any]:
    return work.get("best_oa_location") or {}


def all_oa_locations(work: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        loc for loc in (work.get("locations") or [])
        if bool(loc.get("is_oa"))
    ]


def location_score(loc: dict[str, Any]) -> tuple:
    """
    Prefer:
    1) explicit open license
    2) repository/preprint copy
    3) published > accepted > submitted version
    4) direct PDF
    """
    source = loc.get("source") or {}
    source_type = clean_text(source.get("type")).casefold()
    version = clean_text(loc.get("version"))
    version_rank = {
        "publishedVersion": 3,
        "acceptedVersion": 2,
        "submittedVersion": 1,
    }.get(version, 0)
    return (
        int(is_explicit_open_license(clean_text(loc.get("license")))),
        int(source_type == "repository" or is_preprint_location(loc)),
        version_rank,
        int(bool(loc.get("pdf_url"))),
    )


def select_direct_oa_pdf(work: dict[str, Any]) -> Optional[dict[str, Any]]:
    candidates = [
        loc for loc in all_oa_locations(work)
        if clean_text(loc.get("pdf_url"))
    ]
    if not candidates:
        return None
    return sorted(candidates, key=location_score, reverse=True)[0]


def load_publications(xlsx_path: Path) -> tuple[list[Publication], dict[str, int]]:
    df = pd.read_excel(xlsx_path, sheet_name="Publications")

    required = {
        "PGS Publication/Study (PGP) ID",
        "First Author",
        "Title",
        "Journal Name",
        "Publication Date",
        "digital object identifier (doi)",
        "PubMed ID (PMID)",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            "Publications sheet is missing required columns: "
            + ", ".join(sorted(missing))
        )

    work_map: dict[str, dict[str, Any]] = {}

    for _, row in df.iterrows():
        pgp = clean_text(row["PGS Publication/Study (PGP) ID"])
        doi = normalize_doi(row["digital object identifier (doi)"])
        pmid = normalize_pmid(row["PubMed ID (PMID)"])
        title = clean_text(row["Title"])

        # DOI is the strongest dedupe key; PMID is the fallback.
        if doi:
            key = "doi:" + doi.casefold()
        elif pmid:
            key = "pmid:" + pmid
        else:
            key = "title:" + normalize_title(title)

        if key not in work_map:
            work_map[key] = {
                "pgp_ids": [],
                "first_author": clean_text(row["First Author"]),
                "title": title,
                "journal": clean_text(row["Journal Name"]),
                "publication_date": clean_text(row["Publication Date"]),
                "doi": doi,
                "pmid": pmid,
            }

        work_map[key]["pgp_ids"].append(pgp)

    works = [
        Publication(
            pgp_ids=";".join(sorted(set(v["pgp_ids"]))),
            first_author=v["first_author"],
            title=v["title"],
            journal=v["journal"],
            publication_date=v["publication_date"],
            doi=v["doi"],
            pmid=v["pmid"],
        )
        for v in work_map.values()
    ]

    stats = {
        "metadata_rows": len(df),
        "unique_works": len(works),
        "rows_with_doi": int(df["digital object identifier (doi)"].notna().sum()),
        "rows_with_pmid": int(df["PubMed ID (PMID)"].notna().sum()),
        "works_with_doi": sum(bool(w.doi) for w in works),
        "works_with_pmid": sum(bool(w.pmid) for w in works),
        "works_without_any_id": sum(not w.doi and not w.pmid for w in works),
    }
    return works, stats


class HttpClient:
    def __init__(
        self,
        api_key: str = "",
        email: str = "",
        sleep_seconds: float = 0.12,
    ) -> None:
        self.api_key = api_key
        self.email = email
        self.sleep_seconds = sleep_seconds
        self.session = requests.Session()
        ua = DEFAULT_USER_AGENT
        if email:
            ua += f" ({email})"
        self.session.headers.update(
            {
                "User-Agent": ua,
                "Accept": "application/json",
            }
        )

    def get_json(
        self,
        url: str,
        params: Optional[dict[str, Any]] = None,
        max_attempts: int = 6,
    ) -> dict[str, Any]:
        params = dict(params or {})
        if self.api_key and "api_key" not in params:
            params["api_key"] = self.api_key

        delay = 1.0
        last_error: Optional[Exception] = None
        for attempt in range(max_attempts):
            try:
                r = self.session.get(url, params=params, timeout=45)
                if r.status_code == 429 or 500 <= r.status_code < 600:
                    retry_after = r.headers.get("Retry-After")
                    if retry_after:
                        try:
                            delay = max(delay, float(retry_after))
                        except ValueError:
                            pass
                    time.sleep(delay)
                    delay = min(delay * 2, 30)
                    continue
                r.raise_for_status()
                time.sleep(self.sleep_seconds)
                return r.json()
            except (requests.RequestException, ValueError) as exc:
                last_error = exc
                if attempt + 1 == max_attempts:
                    break
                time.sleep(delay)
                delay = min(delay * 2, 30)

        raise RuntimeError(f"GET failed after retries: {url}: {last_error}")

    def download_pdf(
        self,
        url: str,
        dest: Path,
        params: Optional[dict[str, Any]] = None,
        max_attempts: int = 5,
    ) -> None:
        params = dict(params or {})
        delay = 1.0

        for attempt in range(max_attempts):
            try:
                with self.session.get(
                    url,
                    params=params,
                    timeout=90,
                    stream=True,
                    allow_redirects=True,
                    headers={"Accept": "application/pdf,*/*"},
                ) as r:
                    if r.status_code == 429 or 500 <= r.status_code < 600:
                        time.sleep(delay)
                        delay = min(delay * 2, 30)
                        continue
                    r.raise_for_status()

                    tmp = dest.with_suffix(dest.suffix + ".part")
                    first = b""
                    with open(tmp, "wb") as fh:
                        for chunk in r.iter_content(chunk_size=1024 * 1024):
                            if not chunk:
                                continue
                            if not first:
                                first = chunk[:8]
                            fh.write(chunk)

                    # Some OA links return HTML landing pages despite ending in .pdf.
                    if not first.startswith(b"%PDF"):
                        tmp.unlink(missing_ok=True)
                        ctype = r.headers.get("Content-Type", "")
                        raise RuntimeError(
                            f"Response was not a PDF (Content-Type={ctype!r})"
                        )

                    tmp.replace(dest)
                    time.sleep(self.sleep_seconds)
                    return
            except requests.RequestException as exc:
                if attempt + 1 == max_attempts:
                    raise RuntimeError(f"PDF download failed: {url}: {exc}")
                time.sleep(delay)
                delay = min(delay * 2, 30)

        raise RuntimeError(f"PDF download failed after retries: {url}")


def chunks(items: list[str], size: int) -> Iterable[list[str]]:
    for i in range(0, len(items), size):
        yield items[i : i + size]


def openalex_batch_lookup(
    dois: list[str],
    client: HttpClient,
) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}

    select_fields = ",".join(
        [
            "id",
            "doi",
            "display_name",
            "publication_year",
            "type",
            "ids",
            "authorships",
            "open_access",
            "best_oa_location",
            "locations",
            "has_content",
            "content_urls",
        ]
    )

    for batch in chunks(dois, 100):
        payload = client.get_json(
            f"{OPENALEX_API}/works",
            params={
                "filter": "doi:" + "|".join(batch),
                "per_page": 100,
                "select": select_fields,
            },
        )
        for work in payload.get("results") or []:
            doi_url = clean_text(work.get("doi"))
            doi = normalize_doi(doi_url)
            if doi:
                result[doi.casefold()] = work

    return result


def openalex_title_search(
    publication: Publication,
    client: HttpClient,
    preprint_only: bool = False,
) -> Optional[dict[str, Any]]:
    params = {
        "search": publication.title,
        "per_page": 8,
        "select": ",".join(
            [
                "id",
                "doi",
                "display_name",
                "publication_year",
                "type",
                "ids",
                "authorships",
                "open_access",
                "best_oa_location",
                "locations",
                "has_content",
                "content_urls",
            ]
        ),
    }
    payload = client.get_json(f"{OPENALEX_API}/works", params=params)

    target_author = first_author_surname(publication.first_author)
    candidates = []

    for work in payload.get("results") or []:
        sim = title_similarity(publication.title, clean_text(work.get("display_name")))
        if sim < 0.92:
            continue

        if target_author:
            cand_author = candidate_first_author_surname(work)
            if cand_author and target_author != cand_author:
                # Strong title matches are still allowed; weaker ones must agree on author.
                if sim < 0.975:
                    continue

        if preprint_only:
            work_type = clean_text(work.get("type")).casefold()
            preprintish = work_type == "preprint" or any(
                is_preprint_location(loc)
                for loc in (work.get("locations") or [])
            )
            if not preprintish:
                continue

        candidates.append((sim, work))

    if not candidates:
        return None

    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]


def openalex_work_id(work: dict[str, Any]) -> str:
    oid = clean_text(work.get("id"))
    if oid.startswith("https://openalex.org/"):
        return oid.rsplit("/", 1)[-1]
    return oid


def source_details(loc: dict[str, Any]) -> tuple[str, str]:
    source = loc.get("source") or {}
    return (
        clean_text(source.get("display_name")),
        clean_text(source.get("type")),
    )


def populate_from_work(
    row: AuditRow,
    work: dict[str, Any],
    matched_by: str,
) -> None:
    row.matched = True
    row.matched_by = matched_by
    row.openalex_id = openalex_work_id(work)
    row.openalex_title = clean_text(work.get("display_name"))
    row.openalex_type = clean_text(work.get("type"))

    oa = work.get("open_access") or {}
    row.is_oa = bool(oa.get("is_oa"))
    row.oa_status = clean_text(oa.get("oa_status"))

    best = best_location(work)
    source_name, source_type = source_details(best)

    row.best_oa_license = clean_text(best.get("license"))
    row.best_oa_version = clean_text(best.get("version"))
    row.best_oa_source = source_name
    row.best_oa_source_type = source_type
    row.best_oa_pdf_url = clean_text(best.get("pdf_url"))
    row.best_oa_landing_url = clean_text(best.get("landing_page_url"))

    has_content = work.get("has_content") or {}
    content_urls = work.get("content_urls") or {}
    row.has_cached_pdf = bool(has_content.get("pdf"))
    row.cached_pdf_url = clean_text(content_urls.get("pdf"))

    licenses = [
        clean_text(loc.get("license"))
        for loc in all_oa_locations(work)
    ]
    row.explicit_open_license = any(
        is_explicit_open_license(lic) for lic in licenses
    )

    direct = select_direct_oa_pdf(work)
    if direct and not row.best_oa_pdf_url:
        row.best_oa_pdf_url = clean_text(direct.get("pdf_url"))
        d_source, d_type = source_details(direct)
        row.best_oa_source = row.best_oa_source or d_source
        row.best_oa_source_type = row.best_oa_source_type or d_type
        row.best_oa_license = row.best_oa_license or clean_text(direct.get("license"))
        row.best_oa_version = row.best_oa_version or clean_text(direct.get("version"))

    row.public_fetchable = bool(
        row.has_cached_pdf or row.best_oa_pdf_url
    )

    if row.explicit_open_license:
        row.access_class = "explicit_open_license"
    elif row.is_oa:
        row.access_class = "public_free_to_read_copy"
    else:
        row.access_class = "closed_or_no_public_copy"


def maybe_attach_preprint(
    row: AuditRow,
    publication: Publication,
    client: HttpClient,
) -> Optional[dict[str, Any]]:
    candidate = openalex_title_search(publication, client, preprint_only=True)
    if not candidate:
        return None

    oa = candidate.get("open_access") or {}
    if not oa.get("is_oa"):
        return None

    loc = select_direct_oa_pdf(candidate) or best_location(candidate)
    has_content = candidate.get("has_content") or {}
    content_urls = candidate.get("content_urls") or {}

    row.preprint_found = True
    row.preprint_openalex_id = openalex_work_id(candidate)
    row.preprint_title = clean_text(candidate.get("display_name"))
    row.preprint_pdf_url = (
        clean_text(content_urls.get("pdf"))
        if bool(has_content.get("pdf"))
        else clean_text(loc.get("pdf_url"))
    )
    row.preprint_license = clean_text(loc.get("license"))

    if is_explicit_open_license(row.preprint_license):
        row.explicit_open_license = True
        row.access_class = "explicit_open_license"
    else:
        row.access_class = "public_free_to_read_copy"

    if row.preprint_pdf_url:
        row.public_fetchable = True

    return candidate


def choose_download_source(
    row: AuditRow,
    primary_work: Optional[dict[str, Any]],
    preprint_work: Optional[dict[str, Any]],
    direct_oa: bool,
    explicit_license_only: bool,
) -> None:
    if explicit_license_only and not row.explicit_open_license:
        return

    # 1) Preferred: OpenAlex cached PDF for the primary work.
    if primary_work:
        has_content = primary_work.get("has_content") or {}
        content_urls = primary_work.get("content_urls") or {}
        if bool(has_content.get("pdf")) and clean_text(content_urls.get("pdf")):
            row.chosen_method = "openalex_content_pdf"
            row.chosen_url = clean_text(content_urls.get("pdf"))
            row.chosen_version = row.best_oa_version
            row.chosen_license = row.best_oa_license
            return

    # 2) OpenAlex cached PDF for a matched preprint.
    if preprint_work:
        has_content = preprint_work.get("has_content") or {}
        content_urls = preprint_work.get("content_urls") or {}
        if bool(has_content.get("pdf")) and clean_text(content_urls.get("pdf")):
            best = best_location(preprint_work)
            row.chosen_method = "openalex_preprint_content_pdf"
            row.chosen_url = clean_text(content_urls.get("pdf"))
            row.chosen_version = clean_text(best.get("version")) or "submittedVersion"
            row.chosen_license = clean_text(best.get("license"))
            return

    if not direct_oa:
        return

    # 3) Direct OA location for primary work.
    if primary_work:
        loc = select_direct_oa_pdf(primary_work)
        if loc:
            license_value = clean_text(loc.get("license"))
            if (not explicit_license_only) or is_explicit_open_license(license_value):
                row.chosen_method = "direct_oa_pdf"
                row.chosen_url = clean_text(loc.get("pdf_url"))
                row.chosen_version = clean_text(loc.get("version"))
                row.chosen_license = license_value
                return

    # 4) Direct OA preprint location.
    if preprint_work:
        loc = select_direct_oa_pdf(preprint_work)
        if loc:
            license_value = clean_text(loc.get("license"))
            if (not explicit_license_only) or is_explicit_open_license(license_value):
                row.chosen_method = "direct_preprint_pdf"
                row.chosen_url = clean_text(loc.get("pdf_url"))
                row.chosen_version = clean_text(loc.get("version")) or "submittedVersion"
                row.chosen_license = license_value


def publication_file_token(publication: Publication) -> str:
    """Stable suffix used to recognize the same paper across reruns."""
    if publication.doi:
        return hashlib.sha1(publication.doi.encode("utf-8")).hexdigest()[:8]
    if publication.pmid:
        return publication.pmid
    return hashlib.sha1(publication.title.encode("utf-8")).hexdigest()[:8]


def safe_filename(publication: Publication, suffix: str = ".pdf") -> str:
    year = ""
    m = re.search(r"\b(19|20)\d{2}\b", publication.publication_date)
    if m:
        year = m.group(0)

    author = re.sub(r"[^A-Za-z0-9_-]+", "_", publication.first_author)[:40]
    pgp = publication.pgp_ids.split(";")[0] if publication.pgp_ids else "PGP"
    file_token = publication_file_token(publication)

    bits = [b for b in [pgp, year, author, file_token] if b]
    return "__".join(bits) + suffix


def is_valid_pdf(path: Path) -> bool:
    """Return True only for a non-trivial file that actually contains a PDF header."""
    try:
        if not path.is_file() or path.stat().st_size <= 1024:
            return False
        with path.open("rb") as fh:
            return b"%PDF-" in fh.read(1024)
    except OSError:
        return False


def find_existing_pdf(publication: Publication, out_dir: Path) -> Optional[Path]:
    """
    Find a previously downloaded PDF before doing any network work.

    We first check the exact current filename.  We then match by the stable
    DOI/PMID-derived suffix so a file is still recognized if metadata such as
    author/year/PGP naming changed between workbook versions.
    """
    expected = out_dir / safe_filename(publication)
    candidates = [expected]

    token = publication_file_token(publication)
    candidates.extend(sorted(out_dir.glob(f"*__{token}.pdf")))

    seen: set[Path] = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        if is_valid_pdf(candidate):
            return candidate
    return None


def write_audit(rows: list[AuditRow], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([asdict(r) for r in rows]).to_csv(out_csv, index=False)


def print_inventory(stats: dict[str, int]) -> None:
    print("\nMetadata inventory")
    print("------------------")
    print(f"Publication rows:              {stats['metadata_rows']}")
    print(f"Unique works after dedupe:     {stats['unique_works']}")
    print(f"Rows with DOI:                 {stats['rows_with_doi']}")
    print(f"Rows with PMID:                {stats['rows_with_pmid']}")
    print(f"Unique works with DOI:         {stats['works_with_doi']}")
    print(f"Unique works with PMID:        {stats['works_with_pmid']}")
    print(f"Unique works without DOI/PMID: {stats['works_without_any_id']}")


def print_audit_summary(rows: list[AuditRow]) -> None:
    n = len(rows)
    matched = sum(r.matched for r in rows)
    is_oa = sum(r.is_oa for r in rows)
    explicit = sum(r.explicit_open_license for r in rows)
    cached = sum(r.has_cached_pdf for r in rows)
    public_fetchable = sum(r.public_fetchable for r in rows)
    preprints = sum(r.preprint_found for r in rows)
    chosen = sum(bool(r.chosen_url) for r in rows)
    downloaded = sum(r.downloaded for r in rows)
    skipped_existing = sum(r.skipped_existing for r in rows)
    errors = sum(bool(r.error) for r in rows)

    print("\nLive access audit")
    print("-----------------")
    print(f"Unique works:                       {n}")
    print(f"Matched in OpenAlex:                {matched}")
    print(f"OpenAlex says OA/free-to-read:      {is_oa}")
    print(f"Explicit open license found:        {explicit}")
    print(f"OpenAlex cached PDF available:      {cached}")
    print(f"Public PDF route found:             {public_fetchable}")
    print(f"Matched open preprint found:        {preprints}")
    print(f"Selected for download:              {chosen}")
    print(f"Downloaded/already present:         {downloaded}")
    print(f"Existing PDFs skipped pre-network:  {skipped_existing}")
    print(f"Rows with errors:                   {errors}")

    print(
        "\nInterpretation:\n"
        "  explicit open license = strongest reuse-rights signal.\n"
        "  OA/free-to-read without an explicit license may be readable/downloadable,\n"
        "  but downstream redistribution/reuse rights should not be assumed."
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Audit/download public full text for PGS metadata publications."
    )
    p.add_argument("workbook", type=Path)
    p.add_argument("--output-dir", type=Path, default=Path("pgs_papers"))
    p.add_argument(
        "--audit-csv",
        type=Path,
        default=None,
        help="Default: <output-dir>/paper_access_audit.csv",
    )
    p.add_argument(
        "--openalex-api-key",
        default=os.getenv("OPENALEX_API_KEY", ""),
        help="OpenAlex API key; defaults to OPENALEX_API_KEY.",
    )
    p.add_argument(
        "--email",
        default=os.getenv("CONTACT_EMAIL", ""),
        help="Optional contact email included in User-Agent.",
    )
    p.add_argument(
        "--inventory-only",
        action="store_true",
        help="Read workbook and print identifier coverage; no network.",
    )
    p.add_argument(
        "--audit-only",
        action="store_true",
        help="Query availability but do not download.",
    )
    p.add_argument(
        "--download",
        action="store_true",
        help="Download selected public PDFs.",
    )
    p.add_argument(
        "--direct-oa",
        action="store_true",
        help=(
            "If OpenAlex has no cached PDF, also try a direct public OA PDF URL. "
            "Off by default because host-specific automated-access terms vary."
        ),
    )
    p.add_argument(
        "--explicit-license-only",
        action="store_true",
        help="Only choose PDFs whose OA metadata contains an explicit open license.",
    )
    p.add_argument(
        "--no-preprint-search",
        action="store_true",
        help="Do not title-search for a public preprint when the primary work is closed.",
    )
    p.add_argument(
        "--sleep",
        type=float,
        default=0.12,
        help="Delay between API requests/downloads in seconds.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if not args.workbook.exists():
        raise FileNotFoundError(args.workbook)

    publications, stats = load_publications(args.workbook)
    print_inventory(stats)

    if args.inventory_only:
        return 0

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    audit_csv = args.audit_csv or (out_dir / "paper_access_audit.csv")
    metadata_cache = out_dir / "openalex_metadata_cache.json"

    # Resume support: when downloading, scan the output directory before any
    # OpenAlex requests. Existing valid PDFs are excluded from DOI resolution,
    # title/preprint searches, and download attempts.
    existing_pdfs: dict[str, Path] = {}
    if args.download:
        for pub in publications:
            existing = find_existing_pdf(pub, out_dir)
            if existing is not None:
                existing_pdfs[publication_file_token(pub)] = existing

        print(
            f"\nFound {len(existing_pdfs)} existing valid PDF(s) in {out_dir}. "
            "They will be skipped before any OpenAlex network calls."
        )

    pending_publications = [
        pub
        for pub in publications
        if publication_file_token(pub) not in existing_pdfs
    ]

    client = HttpClient(
        api_key=args.openalex_api_key,
        email=args.email,
        sleep_seconds=args.sleep,
    )

    # Resolve only papers that are still missing locally.
    all_dois = [p.doi for p in pending_publications if p.doi]
    print(
        f"\nResolving {len(all_dois)} DOI-bearing missing works in OpenAlex "
        f"({len(pending_publications)} total works still pending)..."
    )
    doi_lookup = openalex_batch_lookup(all_dois, client) if all_dois else {}

    rows: list[AuditRow] = []

    # Preserve metadata collected during previous partial runs instead of
    # overwriting the cache with only this run's newly processed records.
    work_cache: dict[str, dict[str, Any]] = {}
    if metadata_cache.exists():
        try:
            cached = json.loads(metadata_cache.read_text(encoding="utf-8"))
            if isinstance(cached, dict):
                work_cache.update(cached)
        except (OSError, json.JSONDecodeError):
            pass

    for idx, pub in enumerate(publications, start=1):
        row = AuditRow(
            pgp_ids=pub.pgp_ids,
            first_author=pub.first_author,
            title=pub.title,
            journal=pub.journal,
            publication_date=pub.publication_date,
            doi=pub.doi,
            pmid=pub.pmid,
        )

        existing_dest = existing_pdfs.get(publication_file_token(pub))
        if args.download and existing_dest is not None:
            row.downloaded = True
            row.skipped_existing = True
            row.chosen_method = "existing_file"
            row.local_path = str(existing_dest)
            rows.append(row)

            # Keep checkpoint/progress behaviour identical to normal records.
            if idx % 25 == 0 or idx == len(publications):
                write_audit(rows, audit_csv)
                metadata_cache.write_text(
                    json.dumps(work_cache, indent=2, ensure_ascii=False),
                    encoding="utf-8",
                )
                print(
                    f"Processed {idx}/{len(publications)} "
                    f"(skipped existing PDF: {existing_dest.name})",
                    flush=True,
                )
            continue

        primary_work: Optional[dict[str, Any]] = None
        preprint_work: Optional[dict[str, Any]] = None

        try:
            if pub.doi:
                primary_work = doi_lookup.get(pub.doi.casefold())
                if primary_work:
                    populate_from_work(row, primary_work, "doi")

            # Fallback for DOI-less or unmatched records.
            if primary_work is None:
                primary_work = openalex_title_search(pub, client, preprint_only=False)
                if primary_work:
                    populate_from_work(row, primary_work, "title")

            # If the matched journal work is not publicly fetchable, try to find
            # a separate preprint record with a near-identical title.
            if (
                not args.no_preprint_search
                and not row.public_fetchable
                and pub.title
            ):
                preprint_work = maybe_attach_preprint(row, pub, client)

            choose_download_source(
                row=row,
                primary_work=primary_work,
                preprint_work=preprint_work,
                direct_oa=args.direct_oa,
                explicit_license_only=args.explicit_license_only,
            )

            if primary_work:
                work_cache[row.openalex_id or pub.doi or pub.pmid] = primary_work
            if preprint_work:
                work_cache[
                    "preprint:" + (row.preprint_openalex_id or pub.title)
                ] = preprint_work

            if args.download and row.chosen_url:
                dest = out_dir / safe_filename(pub)

                if is_valid_pdf(dest):
                    row.downloaded = True
                    row.skipped_existing = True
                    row.local_path = str(dest)
                else:
                    if row.chosen_method.startswith("openalex_"):
                        if not args.openalex_api_key:
                            raise RuntimeError(
                                "OpenAlex content downloads require an API key. "
                                "Set OPENALEX_API_KEY or pass --openalex-api-key."
                            )
                        work_id = (
                            row.preprint_openalex_id
                            if "preprint" in row.chosen_method
                            else row.openalex_id
                        )
                        url = f"{OPENALEX_CONTENT}/works/{work_id}.pdf"
                        client.download_pdf(
                            url,
                            dest,
                            params={"api_key": args.openalex_api_key},
                        )
                    else:
                        client.download_pdf(row.chosen_url, dest)

                    row.downloaded = True
                    row.local_path = str(dest)

        except Exception as exc:
            row.error = f"{type(exc).__name__}: {exc}"

        rows.append(row)

        # Resume-friendly checkpoint.
        if idx % 25 == 0 or idx == len(publications):
            write_audit(rows, audit_csv)
            metadata_cache.write_text(
                json.dumps(work_cache, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            print(f"Processed {idx}/{len(publications)}", flush=True)

    print_audit_summary(rows)
    print(f"\nAudit CSV: {audit_csv}")
    print(f"Metadata cache: {metadata_cache}")
    if args.download:
        print(f"PDF directory: {out_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
