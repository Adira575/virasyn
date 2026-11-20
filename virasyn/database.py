from __future__ import annotations

import csv
import io
import os
import threading
import time
from typing import Any

import requests
import streamlit as st

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA = "https://data.rcsb.org/rest/v1/core/entry"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
EUROPE_PMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
OPENALEX_WORKS = "https://api.openalex.org/works"
CROSSREF_WORKS = "https://api.crossref.org/works"
ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
INTERPRO_SEARCH = "https://www.ebi.ac.uk/interpro/api/entry/interpro/"
USER_AGENT = "VIRASYN/0.4 (viral research discovery application)"


def _secret(name: str, default: str = "") -> str:
    value = os.getenv(name, "")
    if value:
        return value
    try:
        return str(st.secrets.get(name, default))
    except Exception:
        return default


def _ncbi_params() -> dict[str, str]:
    params = {"tool": "VIRASYN"}
    api_key = _secret("VIRASYN_NCBI_API_KEY")
    if api_key:
        params["api_key"] = api_key
    email = _secret("VIRASYN_CONTACT_EMAIL")
    if email:
        params["email"] = email
    return params


def _get_json(url: str, params: dict[str, Any] | None = None) -> dict:
    response = requests.get(url, params=params, timeout=35, headers={"User-Agent": USER_AGENT})
    response.raise_for_status()
    return response.json()


_NCBI_LOCK = threading.Lock()
_LAST_NCBI_REQUEST = 0.0


def _get_ncbi_json(url: str, params: dict[str, Any] | None = None) -> dict:
    """Serialize E-utilities calls to stay below the unauthenticated request rate."""
    global _LAST_NCBI_REQUEST
    with _NCBI_LOCK:
        elapsed = time.monotonic() - _LAST_NCBI_REQUEST
        if elapsed < 0.36:
            time.sleep(0.36 - elapsed)
        response = requests.get(url, params=params, timeout=35, headers={"User-Agent": USER_AGENT})
        _LAST_NCBI_REQUEST = time.monotonic()
    response.raise_for_status()
    return response.json()


def _entrez_search(db: str, term: str, offset: int, limit: int,
                   date_from: str = "", date_to: str = "",
                   sort: str = "Best match") -> tuple[int, list[str]]:
    params = {
        **_ncbi_params(), "db": db, "term": term, "retmode": "json",
        "retstart": offset, "retmax": limit,
    }
    if db == "pubmed":
        params["sort"] = "date" if sort == "Newest first" else "relevance"
        if date_from:
            params["mindate"] = date_from.replace("-", "/")
        if date_to:
            params["maxdate"] = date_to.replace("-", "/")
        if date_from or date_to:
            params["datetype"] = "pdat"
    payload = _get_ncbi_json(f"{NCBI_BASE}/esearch.fcgi", params)
    result = payload.get("esearchresult", {})
    if result.get("error"):
        raise RuntimeError("; ".join(result["error"]))
    return int(result.get("count", 0)), result.get("idlist", [])


def _entrez_summary(db: str, ids: list[str]) -> list[dict]:
    if not ids:
        return []
    payload = _get_ncbi_json(
        f"{NCBI_BASE}/esummary.fcgi",
        {**_ncbi_params(), "db": db, "id": ",".join(ids), "retmode": "json"},
    )
    result = payload.get("result", {})
    rows = []
    for uid in result.get("uids", ids):
        row = result.get(uid, {})
        if not row:
            continue
        if db == "pubmed":
            article_ids = row.get("articleids", [])
            doi = next((x.get("value", "") for x in article_ids if x.get("idtype") == "doi"), "")
            authors = row.get("authors", [])
            rows.append({
                "id": str(uid), "title": row.get("title", "Untitled record"),
                "date": row.get("pubdate", ""),
                "journal": row.get("fulljournalname", row.get("source", "")),
                "authors": ", ".join(a.get("name", "") for a in authors[:4]),
                "doi": doi, "url": f"https://pubmed.ncbi.nlm.nih.gov/{uid}/",
                "database": "PubMed",
            })
        elif db == "taxonomy":
            taxid = str(row.get("taxid", uid))
            rows.append({
                "id": taxid, "title": row.get("scientificname", "Taxonomy record"),
                "organism": row.get("scientificname", ""), "rank": row.get("rank", ""),
                "lineage": row.get("lineage", ""),
                "url": f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}",
                "database": "NCBI Taxonomy",
            })
        else:
            accession = row.get("accessionversion") or row.get("caption") or str(uid)
            rows.append({
                "id": str(uid), "accession": accession,
                "title": row.get("title", "Untitled sequence record"),
                "organism": row.get("organism", ""),
                "length": row.get("slen", row.get("length", "")),
                "date": row.get("updatedate", row.get("createdate", "")),
                "url": f"https://www.ncbi.nlm.nih.gov/{'nuccore' if db == 'nuccore' else 'protein'}/{uid}",
                "database": "NCBI Nucleotide" if db == "nuccore" else "NCBI Protein",
            })
    return rows


@st.cache_data(ttl=1800, show_spinner=False)
def search_pubmed(term: str, offset: int = 0, limit: int = 20,
                  date_from: str = "", date_to: str = "",
                  sort: str = "Best match") -> dict:
    count, ids = _entrez_search("pubmed", term, offset, limit, date_from, date_to, sort)
    return {"count": count, "results": _entrez_summary("pubmed", ids), "source": "NCBI PubMed"}


@st.cache_data(ttl=1800, show_spinner=False)
def search_ncbi_sequences(db: str, term: str, offset: int = 0, limit: int = 20) -> dict:
    if db not in {"nuccore", "protein"}:
        raise ValueError("Only nuccore and protein databases are supported.")
    count, ids = _entrez_search(db, term, offset, limit)
    return {
        "count": count, "results": _entrez_summary(db, ids),
        "source": "NCBI Nucleotide" if db == "nuccore" else "NCBI Protein",
    }


@st.cache_data(ttl=1800, show_spinner=False)
def search_ncbi_taxonomy(term: str, offset: int = 0, limit: int = 20) -> dict:
    count, ids = _entrez_search("taxonomy", term, offset, limit)
    return {"count": count, "results": _entrez_summary("taxonomy", ids), "source": "NCBI Taxonomy"}


@st.cache_data(ttl=1800, show_spinner=False)
def search_rcsb(term: str, offset: int = 0, limit: int = 20) -> dict:
    query = {
        "query": {"type": "terminal", "service": "full_text", "parameters": {"value": term}},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": offset, "rows": limit},
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }
    response = requests.post(
        RCSB_SEARCH, json=query, timeout=35,
        headers={"User-Agent": USER_AGENT, "Content-Type": "application/json"},
    )
    response.raise_for_status()
    payload = response.json()
    identifiers = [x.get("identifier", "") for x in payload.get("result_set", [])]
    results = []
    for identifier in identifiers:
        try:
            data = _get_json(f"{RCSB_DATA}/{identifier}")
            struct = data.get("struct", {})
            release = data.get("rcsb_accession_info", {}).get("initial_release_date", "")
            exptl = data.get("exptl", [])
            method = ", ".join(sorted({x.get("method", "") for x in exptl if x.get("method")}))
            results.append({
                "id": identifier, "title": struct.get("title", "PDB structure"),
                "date": release, "method": method,
                "url": f"https://www.rcsb.org/structure/{identifier}",
                "database": "RCSB Protein Data Bank",
            })
        except requests.RequestException:
            results.append({
                "id": identifier, "title": "PDB structure record", "date": "", "method": "",
                "url": f"https://www.rcsb.org/structure/{identifier}",
                "database": "RCSB Protein Data Bank",
            })
    return {"count": payload.get("total_count", len(results)), "results": results,
            "source": "RCSB Protein Data Bank"}


@st.cache_data(ttl=1800, show_spinner=False)
def search_uniprot(term: str, offset: int = 0, limit: int = 20) -> dict:
    params = {
        "query": term, "format": "tsv", "size": limit, "offset": offset,
        "fields": "accession,id,protein_name,organism_name,length",
    }
    response = requests.get(UNIPROT_SEARCH, params=params, timeout=35, headers={"User-Agent": USER_AGENT})
    response.raise_for_status()
    total = int(response.headers.get("x-total-results", "0"))
    reader = csv.DictReader(io.StringIO(response.text), delimiter="\t")
    results = []
    for row in reader:
        accession = row.get("Entry", "")
        results.append({
            "id": accession, "accession": accession,
            "title": row.get("Protein names", "Protein entry"),
            "organism": row.get("Organism", ""), "length": row.get("Length", ""),
            "entry_name": row.get("Entry Name", ""),
            "url": f"https://www.uniprot.org/uniprotkb/{accession}/entry",
            "database": "UniProtKB",
        })
    return {"count": total, "results": results, "source": "UniProtKB"}


@st.cache_data(ttl=1800, show_spinner=False)
def search_europe_pmc(term: str, page: int = 1, limit: int = 20) -> dict:
    payload = _get_json(EUROPE_PMC_SEARCH, {
        "query": term, "format": "json", "resultType": "lite",
        "page": page, "pageSize": limit,
    })
    items = payload.get("resultList", {}).get("result", [])
    results = []
    for item in items:
        pmcid, pmid, doi = item.get("pmcid", ""), item.get("pmid", ""), item.get("doi", "")
        url = f"https://europepmc.org/article/MED/{pmid}" if pmid else (
            f"https://europepmc.org/articles/{pmcid}" if pmcid else f"https://doi.org/{doi}"
        )
        results.append({
            "id": str(pmcid or pmid or doi),
            "title": item.get("title", "Untitled publication"),
            "date": item.get("firstPublicationDate", ""),
            "journal": item.get("journalTitle", ""),
            "authors": item.get("authorString", ""),
            "doi": doi, "pmid": pmid, "pmcid": pmcid,
            "is_open_access": item.get("isOpenAccess", "N") == "Y",
            "url": url, "database": "Europe PMC",
        })
    return {"count": int(payload.get("hitCount", 0)), "results": results, "source": "Europe PMC"}


@st.cache_data(ttl=1800, show_spinner=False)
def search_openalex(term: str, page: int = 1, limit: int = 20) -> dict:
    payload = _get_json(OPENALEX_WORKS, {
        "search": term, "page": page, "per-page": limit,
        "select": "id,doi,title,publication_year,authorships,primary_location,cited_by_count,open_access",
    })
    results = []
    for item in payload.get("results", []):
        authors = ", ".join(
            a.get("author", {}).get("display_name", "")
            for a in item.get("authorships", [])[:4] if a.get("author")
        )
        doi = item.get("doi") or ""
        source = ((item.get("primary_location") or {}).get("source") or {}).get("display_name", "")
        results.append({
            "id": item.get("id", ""), "title": item.get("title", "Untitled work"),
            "date": str(item.get("publication_year", "")), "authors": authors,
            "journal": source, "doi": doi,
            "cited_by_count": item.get("cited_by_count", 0),
            "is_open_access": (item.get("open_access") or {}).get("is_oa", False),
            "url": doi or item.get("id", ""), "database": "OpenAlex",
        })
    return {"count": (payload.get("meta") or {}).get("count", 0),
            "results": results, "source": "OpenAlex"}


@st.cache_data(ttl=1800, show_spinner=False)
def search_crossref(term: str, offset: int = 0, limit: int = 20) -> dict:
    params = {
        "query.bibliographic": term, "rows": limit, "offset": offset,
        "select": "DOI,title,published,container-title,author,URL,type",
    }
    mailto = _secret("VIRASYN_CONTACT_EMAIL")
    if mailto:
        params["mailto"] = mailto
    payload = _get_json(CROSSREF_WORKS, params)
    message = payload.get("message", {})
    results = []
    for item in message.get("items", []):
        doi = item.get("DOI", "")
        date_parts = (item.get("published") or {}).get("date-parts") or [[""]]
        year = str(date_parts[0][0]) if date_parts[0] else ""
        authors = ", ".join(
            " ".join(x for x in [a.get("given", ""), a.get("family", "")] if x)
            for a in item.get("author", [])[:4]
        )
        results.append({
            "id": doi, "title": (item.get("title") or ["Untitled work"])[0],
            "date": year, "authors": authors,
            "journal": (item.get("container-title") or [""])[0], "doi": doi,
            "url": item.get("URL") or (f"https://doi.org/{doi}" if doi else ""),
            "database": "Crossref",
        })
    return {"count": int(message.get("total-results", 0)), "results": results, "source": "Crossref"}


@st.cache_data(ttl=1800, show_spinner=False)
def search_ena(term: str, offset: int = 0, limit: int = 20) -> dict:
    clean_term = term.replace('"', "").strip()
    query = f'description="*{clean_term}*"'
    response = requests.get(
        ENA_SEARCH,
        params={
            "result": "sequence", "query": query, "format": "json",
            "limit": limit, "offset": offset,
            "fields": "accession,description,scientific_name,length",
        },
        timeout=35, headers={"User-Agent": USER_AGENT},
    )
    response.raise_for_status()
    results = []
    for item in response.json():
        accession = item.get("accession", "")
        results.append({
            "id": accession, "accession": accession,
            "title": item.get("description", "ENA sequence record"),
            "organism": item.get("scientific_name", ""), "length": item.get("length", ""),
            "url": f"https://www.ebi.ac.uk/ena/browser/view/{accession}",
            "database": "European Nucleotide Archive",
        })
    return {"count": len(results), "results": results, "source": "European Nucleotide Archive"}


@st.cache_data(ttl=1800, show_spinner=False)
def search_interpro(term: str, page: int = 1, limit: int = 20) -> dict:
    payload = _get_json(INTERPRO_SEARCH, {"search": term, "page": page, "page_size": limit})
    results = []
    for item in payload.get("results", []):
        metadata = item.get("metadata", {})
        accession = metadata.get("accession", item.get("accession", ""))
        name_obj = metadata.get("name") or {}
        name = name_obj.get("name", "") if isinstance(name_obj, dict) else str(name_obj)
        desc_obj = metadata.get("description") or []
        description = "; ".join(
            x.get("text", "") for x in desc_obj if isinstance(x, dict) and x.get("text")
        )
        results.append({
            "id": accession, "accession": accession,
            "title": name or accession or "InterPro entry",
            "description": description, "type": metadata.get("type", ""),
            "url": f"https://www.ebi.ac.uk/interpro/entry/InterPro/{accession}/",
            "database": "InterPro",
        })
    return {"count": (payload.get("meta") or {}).get("count", len(results)),
            "results": results, "source": "InterPro"}


@st.cache_data(ttl=1800, show_spinner=False)
def search_all_sources(term: str, page: int = 1, limit: int = 5) -> dict:
    """Federated search across public literature, sequence, protein and structure sources."""
    offset = (page - 1) * limit
    searches = [
        ("PubMed", lambda: search_pubmed(term, offset, limit)),
        ("Europe PMC", lambda: search_europe_pmc(term, page, limit)),
        ("OpenAlex", lambda: search_openalex(term, page, limit)),
        ("Crossref", lambda: search_crossref(term, offset, limit)),
        ("NCBI Nucleotide", lambda: search_ncbi_sequences("nuccore", term, offset, limit)),
        ("NCBI Protein", lambda: search_ncbi_sequences("protein", term, offset, limit)),
        ("NCBI Taxonomy", lambda: search_ncbi_taxonomy(term, offset, limit)),
        ("ENA", lambda: search_ena(term, offset, limit)),
        ("RCSB PDB", lambda: search_rcsb(term, offset, limit)),
        ("UniProtKB", lambda: search_uniprot(term, offset, limit)),
        ("InterPro", lambda: search_interpro(term, page, limit)),
    ]
    groups, errors = {}, {}
    for name, search in searches:
        try:
            groups[name] = search()
        except (requests.RequestException, ValueError, RuntimeError, KeyError) as exc:
            errors[name] = str(exc)
    return {"groups": groups, "errors": errors}


def fetch_sequence_fasta(db: str, uid: str) -> str:
    if db not in {"nuccore", "protein"}:
        raise ValueError("FASTA download supports NCBI nucleotide and protein records.")
    response = requests.get(
        f"{NCBI_BASE}/efetch.fcgi",
        params={**_ncbi_params(), "db": db, "id": uid, "rettype": "fasta", "retmode": "text"},
        timeout=35, headers={"User-Agent": USER_AGENT},
    )
    response.raise_for_status()
    text = response.text.strip()
    if not text.startswith(">"):
        raise RuntimeError("NCBI did not return a FASTA record for this identifier.")
    return text + "\n"
