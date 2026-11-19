from __future__ import annotations

import csv
import io
import os
from typing import Any

import requests
import streamlit as st

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA = "https://data.rcsb.org/rest/v1/core/entry"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
USER_AGENT = "VIRASYN/0.3 (academic research web application)"


def _ncbi_params() -> dict[str, str]:
    params = {"tool": "VIRASYN"}
    api_key = os.getenv("VIRASYN_NCBI_API_KEY")
    if not api_key:
        try:
            api_key = st.secrets.get("VIRASYN_NCBI_API_KEY", "")
        except Exception:
            api_key = ""
    if api_key:
        params["api_key"] = api_key
    return params


def _get_json(url: str, params: dict[str, Any] | None = None) -> dict:
    response = requests.get(url, params=params, timeout=30, headers={"User-Agent": USER_AGENT})
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
    payload = _get_json(f"{NCBI_BASE}/esearch.fcgi", params)
    result = payload.get("esearchresult", {})
    if result.get("error"):
        raise RuntimeError("; ".join(result["error"]))
    return int(result.get("count", 0)), result.get("idlist", [])


def _entrez_summary(db: str, ids: list[str]) -> list[dict]:
    if not ids:
        return []
    payload = _get_json(
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
        "count": count,
        "results": _entrez_summary(db, ids),
        "source": "NCBI Nucleotide" if db == "nuccore" else "NCBI Protein",
    }


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
        RCSB_SEARCH, json=query, timeout=30,
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
    response = requests.get(UNIPROT_SEARCH, params=params, timeout=30,
                            headers={"User-Agent": USER_AGENT})
    response.raise_for_status()
    total = int(response.headers.get("x-total-results", "0"))
    reader = csv.DictReader(io.StringIO(response.text), delimiter="\t")
    results = []
    for row in reader:
        accession = row.get("Entry", "")
        results.append({
            "id": accession, "accession": accession,
            "title": row.get("Protein names", "Protein entry"),
            "organism": row.get("Organism", ""),
            "length": row.get("Length", ""),
            "entry_name": row.get("Entry Name", ""),
            "url": f"https://www.uniprot.org/uniprotkb/{accession}/entry",
            "database": "UniProtKB",
        })
    return {"count": total, "results": results, "source": "UniProtKB"}


@st.cache_data(ttl=1800, show_spinner=False)
def search_all_sources(term: str, offset: int = 0, limit: int = 5) -> dict:
    """Small federated search across independent upstream services."""
    searches = [
        ("PubMed", lambda: search_pubmed(term, offset, limit)),
        ("NCBI Nucleotide", lambda: search_ncbi_sequences("nuccore", term, offset, limit)),
        ("NCBI Protein", lambda: search_ncbi_sequences("protein", term, offset, limit)),
        ("RCSB PDB", lambda: search_rcsb(term, offset, limit)),
        ("UniProtKB", lambda: search_uniprot(term, offset, limit)),
    ]
    groups, errors = {}, {}
    for name, search in searches:
        try:
            groups[name] = search()
        except (requests.RequestException, ValueError, RuntimeError) as exc:
            errors[name] = str(exc)
    return {"groups": groups, "errors": errors}


def fetch_sequence_fasta(db: str, uid: str) -> str:
    if db not in {"nuccore", "protein"}:
        raise ValueError("FASTA download supports NCBI nucleotide and protein records.")
    response = requests.get(
        f"{NCBI_BASE}/efetch.fcgi",
        params={**_ncbi_params(), "db": db, "id": uid, "rettype": "fasta", "retmode": "text"},
        timeout=30, headers={"User-Agent": USER_AGENT},
    )
    response.raise_for_status()
    text = response.text.strip()
    if not text.startswith(">"):
        raise RuntimeError("NCBI did not return a FASTA record for this identifier.")
    return text + "\n"
