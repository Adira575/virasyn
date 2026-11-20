"""Safe client for retrieving published models from the AlphaFold Protein Structure Database."""
from __future__ import annotations

import re
from urllib.parse import urlparse

import requests

API_BASE = "https://alphafold.ebi.ac.uk/api/prediction"
TIMEOUT = 30
ACCESSION_RE = re.compile(r"^[A-Z0-9]{6,10}$")


def lookup_alphafold_prediction(accession: str) -> list[dict]:
    """Return AlphaFold Database model metadata for a UniProt accession."""
    accession = accession.strip().upper()
    if not ACCESSION_RE.fullmatch(accession):
        raise ValueError("Use a valid UniProt accession (6–10 letters or digits).")

    try:
        response = requests.get(f"{API_BASE}/{accession}", timeout=TIMEOUT,
                                headers={"Accept": "application/json", "User-Agent": "VIRASYN/1.0"})
        if response.status_code == 404:
            return []
        response.raise_for_status()
        payload = response.json()
    except requests.RequestException as exc:
        raise RuntimeError(f"AlphaFold Database is temporarily unreachable: {exc}") from exc
    except ValueError as exc:
        raise RuntimeError("AlphaFold Database returned an unreadable response.") from exc

    if not isinstance(payload, list):
        raise RuntimeError("AlphaFold Database returned an unexpected response format.")
    return [item for item in payload if isinstance(item, dict)]


def _safe_structure_url(model: dict, coordinate_format: str) -> str:
    key = "pdbUrl" if coordinate_format.upper() == "PDB" else "cifUrl"
    url = model.get(key)
    parsed = urlparse(url or "")
    if parsed.scheme != "https" or parsed.hostname not in {
        "alphafold.ebi.ac.uk", "ftp.ebi.ac.uk"
    }:
        raise ValueError(f"The AlphaFold record does not contain a trusted {coordinate_format.upper()} download URL.")
    return url


def download_alphafold_structure(model: dict, coordinate_format: str = "PDB") -> tuple[str, str]:
    """Download a model coordinate file and return its text and safe filename."""
    fmt = coordinate_format.upper()
    if fmt not in {"PDB", "MMCIF"}:
        raise ValueError("Coordinate format must be PDB or mmCIF.")
    url = _safe_structure_url(model, fmt)
    try:
        response = requests.get(url, timeout=TIMEOUT,
                                headers={"User-Agent": "VIRASYN/1.0"})
        response.raise_for_status()
    except requests.RequestException as exc:
        raise RuntimeError(f"Could not download AlphaFold coordinates: {exc}") from exc

    text = response.text.strip()
    if not text or (fmt == "PDB" and not (text.startswith("HEADER") or "ATOM" in text)):
        raise RuntimeError("The downloaded file does not look like a valid structure coordinate file.")
    accession = re.sub(r"[^A-Z0-9_-]", "_", str(model.get("uniprotAccession", "protein")))
    return text + "\n", f"alphafold_{accession}.{ 'pdb' if fmt == 'PDB' else 'cif' }"
