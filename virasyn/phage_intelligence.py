"""Evidence-integrated bacteriophage RBP, receptor and bacterial-host intelligence.

This module deliberately reports evidence indices and ranked hypotheses, not calibrated
probabilities. A BLAST match or receptor annotation is not proof of adsorption or infection.
"""
from __future__ import annotations

from io import StringIO
import pandas as pd

HOST_COLUMNS = ["phage_id", "host_taxid", "host_name", "evidence_type", "evidence_level", "source", "citation"]
RECEPTOR_COLUMNS = ["subject_id", "receptor", "receptor_class", "host_taxid", "host_name", "evidence_type", "citation"]

_LEVEL_WEIGHT = {
    "experimental": 1.0,
    "curated": 0.75,
    "computational": 0.4,
    "inferred": 0.2,
}
_TYPE_WEIGHT = {
    "productive_infection": 1.0,
    "host_range": 0.95,
    "adsorption": 0.85,
    "receptor_binding": 0.8,
    "crispr_spacer": 0.45,
    "sequence_similarity": 0.25,
    "other": 0.2,
}


def _read_csv(text: str, required: list[str], label: str) -> pd.DataFrame:
    try:
        frame = pd.read_csv(StringIO(text), dtype=str).fillna("")
    except Exception as exc:
        raise ValueError(f"{label} CSV could not be parsed: {exc}") from exc
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} CSV is missing required columns: {', '.join(missing)}")
    if frame.empty:
        raise ValueError(f"{label} CSV contains no data rows.")
    return frame


def score_host_evidence(host_evidence_csv: str) -> tuple[pd.DataFrame, dict]:
    """Rank bacterial hosts from supplied phage-host evidence records.

    The evidence index is a transparent weighted sum, not a probability or a
    model trained to generalize to unseen phages.
    """
    frame = _read_csv(host_evidence_csv, HOST_COLUMNS, "Host evidence")
    frame["evidence_level_norm"] = frame["evidence_level"].str.strip().str.lower()
    frame["evidence_type_norm"] = frame["evidence_type"].str.strip().str.lower()
    unknown_levels = sorted(set(frame["evidence_level_norm"]) - set(_LEVEL_WEIGHT))
    if unknown_levels:
        raise ValueError("Unsupported evidence_level values: " + ", ".join(unknown_levels))
    frame["level_weight"] = frame["evidence_level_norm"].map(_LEVEL_WEIGHT)
    frame["type_weight"] = frame["evidence_type_norm"].map(lambda x: _TYPE_WEIGHT.get(x, _TYPE_WEIGHT["other"]))
    frame["evidence_weight"] = frame["level_weight"] * frame["type_weight"]
    frame["host_taxid"] = frame["host_taxid"].astype(str).str.strip()
    frame["host_name"] = frame["host_name"].astype(str).str.strip()
    if (frame["host_taxid"] == "").any() or (frame["host_name"] == "").any():
        raise ValueError("Every host evidence row needs host_taxid and host_name.")
    ranked = frame.groupby(["host_taxid", "host_name"], as_index=False).agg(
        evidence_index=("evidence_weight", "sum"),
        evidence_records=("evidence_weight", "size"),
        experimental_records=("evidence_level_norm", lambda x: int((x == "experimental").sum())),
        evidence_types=("evidence_type", lambda x: "; ".join(sorted(set(v for v in x if v)))),
        sources=("source", lambda x: "; ".join(sorted(set(v for v in x if v)))),
        citations=("citation", lambda x: "; ".join(sorted(set(v for v in x if v)))),
    )
    ranked = ranked.sort_values(["experimental_records", "evidence_index", "evidence_records"], ascending=False)
    ranked["evidence_tier"] = ranked.apply(
        lambda row: "Experimental support present" if row["experimental_records"] else
        ("Multiple non-experimental records" if row["evidence_records"] > 1 else "Single non-experimental record"),
        axis=1,
    )
    return ranked.reset_index(drop=True), {
        "records": len(frame),
        "unique_hosts": int(ranked["host_taxid"].nunique()),
        "experimental_records": int((frame["evidence_level_norm"] == "experimental").sum()),
        "index_is_probability": False,
    }


def map_rbp_receptor_evidence(candidate_hits: pd.DataFrame, receptor_csv: str) -> pd.DataFrame:
    """Attach curated receptor annotations to RBP reference hits by exact subject ID."""
    annotations = _read_csv(receptor_csv, RECEPTOR_COLUMNS, "RBP-receptor annotation")
    if annotations["subject_id"].duplicated().any():
        raise ValueError("RBP-receptor annotation contains duplicate subject_id values.")
    if candidate_hits is None or candidate_hits.empty:
        return pd.DataFrame(columns=list(candidate_hits.columns if candidate_hits is not None else []) + RECEPTOR_COLUMNS[1:])
    if "sseqid" not in candidate_hits.columns:
        raise ValueError("Candidate hit table must contain BLAST reference identifier column sseqid.")
    merged = candidate_hits.merge(annotations, left_on="sseqid", right_on="subject_id", how="left")
    return merged


def build_integrated_report(
    candidate_hits: pd.DataFrame,
    host_evidence_csv: str | None = None,
    receptor_annotation_csv: str | None = None,
) -> dict:
    """Create transparent ranked host and receptor evidence tables."""
    host_table = pd.DataFrame()
    host_summary = {"records": 0, "unique_hosts": 0, "experimental_records": 0, "index_is_probability": False}
    receptor_table = pd.DataFrame()
    if host_evidence_csv and host_evidence_csv.strip():
        host_table, host_summary = score_host_evidence(host_evidence_csv)
    if receptor_annotation_csv and receptor_annotation_csv.strip():
        receptor_table = map_rbp_receptor_evidence(candidate_hits, receptor_annotation_csv)
    return {
        "host_evidence": host_table,
        "host_summary": host_summary,
        "rbp_receptor_evidence": receptor_table,
        "interpretation": (
            "Ranked evidence hypotheses only. Scores are not probabilities. RBP homology, receptor annotations, "
            "host association, adsorption, receptor binding and productive infection are distinct claims. "
            "Independent validation is required."
        ),
    }


def csv_template(kind: str) -> str:
    if kind == "host":
        return ",".join(HOST_COLUMNS) + "\n"
    if kind == "receptor":
        return ",".join(RECEPTOR_COLUMNS) + "\n"
    raise ValueError("Template kind must be 'host' or 'receptor'.")
