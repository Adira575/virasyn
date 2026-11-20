"""PhARIS-inspired, reference-database-driven RBP homology workflow.

This module is an independent implementation. It does not bundle or reproduce
the original PhARIS reference sequences or notebook.
"""
from __future__ import annotations

from io import StringIO
from pathlib import Path
import shutil
import subprocess
import tempfile

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from virasyn.core import predict_six_frame_orfs

BLAST_FIELDS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend",
                "qlen", "sstart", "send", "slen", "evalue", "bitscore"]


def _read_nucleotide_fasta(text: str):
    records = list(SeqIO.parse(StringIO(text), "fasta"))
    if not records:
        raise ValueError("The genome FASTA contains no valid records.")
    for record in records:
        sequence = str(record.seq).upper().replace("U", "T")
        if not sequence or set(sequence) - set("ACGTUNRYKMSWBDHV"):
            raise ValueError(f"Record {record.id!r} is not a valid nucleotide sequence.")
        record.seq = record.seq.__class__(sequence)
    return records


def _read_reference_fasta(text: str):
    records = list(SeqIO.parse(StringIO(text), "fasta"))
    if not records:
        raise ValueError("The RBP reference FASTA contains no valid protein records.")
    alphabet = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")
    for record in records:
        seq = str(record.seq).upper().rstrip("*")
        if not seq or set(seq) - alphabet:
            raise ValueError(f"Reference {record.id!r} is not a valid protein sequence.")
        record.seq = record.seq.__class__(seq)
    return records


def analyse_pharis_genome(
    genome_fasta: str,
    reference_fasta: str,
    metadata_csv: str | None = None,
    min_orf_aa: int = 100,
    min_identity: float = 80.0,
    min_query_coverage: float = 80.0,
    genetic_code: int = 11,
) -> dict:
    """Predict candidate S. aureus-associated RBP homologues by local BLASTP.

    A match is sequence-similarity evidence to the supplied reference library;
    it is not proof of adsorption, receptor binding, infection, or host range.
    """
    if min_orf_aa < 30:
        raise ValueError("Minimum ORF length must be at least 30 amino acids.")
    if not 0 <= min_identity <= 100 or not 0 <= min_query_coverage <= 100:
        raise ValueError("Identity and coverage thresholds must be between 0 and 100.")

    genomes = _read_nucleotide_fasta(genome_fasta)
    references = _read_reference_fasta(reference_fasta)
    orfs = []
    for record in genomes:
        orfs.extend(predict_six_frame_orfs(
            str(record.seq), sequence_id=record.id, min_aa=min_orf_aa,
            genetic_code=genetic_code, alternative_starts=True, max_orfs=20000,
        ))
    if not orfs:
        raise ValueError(
            "No start-to-stop ORFs met the selected minimum length. Try a lower minimum, "
            "or confirm the sequence and genetic code."
        )

    with tempfile.TemporaryDirectory(prefix="virasyn_pharis_") as temp_dir:
        work = Path(temp_dir)
        query_path = work / "candidate_orfs.faa"
        reference_path = work / "rbp_reference.faa"
        db_prefix = work / "rbp_reference_db"
        output_path = work / "blast_hits.tsv"
        SeqIO.write(
            [SeqIO.SeqRecord(
                seq=Seq(item["protein_sequence"]),
                id=item["orf_id"], description=f"{item['sequence_id']} {item['strand']} frame {item['frame']}",
            ) for item in orfs], query_path, "fasta"
        )
        SeqIO.write(references, reference_path, "fasta")
        makeblastdb, blastp = shutil.which("makeblastdb"), shutil.which("blastp")
        if not makeblastdb or not blastp:
            raise RuntimeError(
                "NCBI BLAST+ is required. Use the VIRASYN Docker image or install makeblastdb and blastp on the server."
            )
        try:
            subprocess.run(
                [makeblastdb, "-in", str(reference_path), "-dbtype", "prot", "-parse_seqids", "-out", str(db_prefix)],
                check=True, capture_output=True, text=True, timeout=180,
            )
            outfmt = "6 " + " ".join(BLAST_FIELDS)
            subprocess.run(
                [blastp, "-query", str(query_path), "-db", str(db_prefix),
                 "-outfmt", outfmt, "-out", str(output_path), "-evalue", "1e-5",
                 "-max_target_seqs", "100", "-num_threads", "2"],
                check=True, capture_output=True, text=True, timeout=900,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError("BLAST+ exceeded the analysis time limit. Try a smaller genome or reference library.") from exc
        except subprocess.CalledProcessError as exc:
            detail = (exc.stderr or exc.stdout or str(exc))[-2000:]
            raise RuntimeError(f"BLAST+ failed: {detail}") from exc
        if output_path.exists() and output_path.stat().st_size:
            hits = pd.read_csv(output_path, sep="\\t", names=BLAST_FIELDS)
        else:
            hits = pd.DataFrame(columns=BLAST_FIELDS)

    orf_table = pd.DataFrame([{
        key: value for key, value in item.items() if key not in {"protein_sequence", "nucleotide_sequence"}
    } for item in orfs])
    if hits.empty:
        hits["query_coverage_pct"] = pd.Series(dtype=float)
        hits["subject_coverage_pct"] = pd.Series(dtype=float)
    else:
        hits["query_coverage_pct"] = (hits["length"] / hits["qlen"] * 100).clip(upper=100).round(2)
        hits["subject_coverage_pct"] = (hits["length"] / hits["slen"] * 100).clip(upper=100).round(2)
        hits = hits.merge(
            orf_table[["orf_id", "sequence_id", "strand", "frame", "start_1based", "end_1based", "length_aa"]],
            left_on="qseqid", right_on="orf_id", how="left",
        )
        hits = hits.sort_values(["bitscore", "evalue"], ascending=[False, True])

    if metadata_csv and not hits.empty:
        metadata = pd.read_csv(StringIO(metadata_csv))
        if "subject_id" not in metadata.columns:
            raise ValueError("Reference metadata CSV must contain a subject_id column matching FASTA record IDs.")
        if metadata["subject_id"].duplicated().any():
            raise ValueError("Reference metadata contains duplicate subject_id values.")
        hits = hits.merge(metadata, left_on="sseqid", right_on="subject_id", how="left")

    hits["passes_thresholds"] = (
        (hits["pident"] >= min_identity) &
        (hits["query_coverage_pct"] >= min_query_coverage)
    ) if not hits.empty else pd.Series(dtype=bool)
    ranked = hits[hits["passes_thresholds"]].copy() if not hits.empty else hits.copy()
    if not ranked.empty:
        ranked = ranked.drop_duplicates(subset=["qseqid"], keep="first")

    return {
        "orfs": orf_table,
        "all_hits": hits,
        "candidate_hits": ranked,
        "summary": {
            "genome_records": len(genomes),
            "reference_records": len(references),
            "orfs_screened": len(orfs),
            "blast_hits": len(hits),
            "candidate_orfs": int(ranked["qseqid"].nunique()) if not ranked.empty else 0,
            "min_orf_aa": min_orf_aa,
            "min_identity": min_identity,
            "min_query_coverage": min_query_coverage,
            "genetic_code": genetic_code,
        },
    }


def reference_fasta_manifest(reference_fasta: str) -> pd.DataFrame:
    """Return an auditable inventory of reference FASTA identifiers and lengths."""
    records = _read_reference_fasta(reference_fasta)
    return pd.DataFrame([{
        "subject_id": record.id, "description": record.description,
        "length_aa": len(record.seq),
    } for record in records])
