from collections import Counter
from datetime import datetime, timezone
from io import StringIO
import json
from pathlib import Path
import shutil
import subprocess

import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

NUCLEOTIDE_ALPHABET = set("ACGTUNRYKMSWBDHV")
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")


def parse_fasta(text: str) -> list[dict]:
    """Parse FASTA records and infer nucleotide/protein type."""
    parsed = list(SeqIO.parse(StringIO(text), "fasta"))
    if not parsed:
        raise ValueError("No FASTA records found. Check for a valid >header.")
    records = []
    for record in parsed:
        sequence = "".join(str(record.seq).upper().replace("*", "").split())
        if not sequence:
            continue
        symbols = set(sequence)
        if symbols <= NUCLEOTIDE_ALPHABET:
            kind = "nucleotide"
        elif symbols <= PROTEIN_ALPHABET:
            kind = "protein"
        else:
            invalid = sorted(symbols - NUCLEOTIDE_ALPHABET - PROTEIN_ALPHABET)
            raise ValueError(f"Sequence {record.id!r} has unsupported symbols: {invalid}")
        records.append({
            "id": record.id,
            "description": record.description,
            "sequence": sequence,
            "type": kind,
        })
    if not records:
        raise ValueError("FASTA contained no non-empty sequences.")
    return records


def sequence_summary(sequence: str) -> dict:
    sequence = sequence.upper()
    is_nucleotide = set(sequence) <= NUCLEOTIDE_ALPHABET
    return {
        "length": len(sequence),
        "sequence_type": "nucleotide" if is_nucleotide else "protein",
        "gc_percent": round(100 * (sequence.count("G") + sequence.count("C")) / len(sequence), 2)
        if is_nucleotide and sequence else None,
        "most_common_symbols": Counter(sequence).most_common(10),
    }


def parse_pdb(text: str, filename: str = "uploaded.pdb") -> dict:
    """Extract basic protein chain metadata from a PDB coordinate file."""
    structure = PDBParser(QUIET=True).get_structure("uploaded", StringIO(text))
    chains, residue_count = [], 0
    first_model = next(structure.get_models(), None)
    if first_model is not None:
        for chain in first_model:
            residues = [res for res in chain if is_aa(res, standard=False)]
            if residues:
                chains.append(chain.id)
                residue_count += len(residues)
    return {
        "file": filename,
        "model_count": len(list(structure.get_models())),
        "protein_chains": ", ".join(chains),
        "amino_acid_residue_count": residue_count,
    }


def run_blastp(query_path: Path, reference_fasta: Path, workdir: Path) -> pd.DataFrame:
    """Run BLASTP against a temporary database made from the supplied FASTA."""
    blastp, makeblastdb = shutil.which("blastp"), shutil.which("makeblastdb")
    if not blastp or not makeblastdb:
        raise RuntimeError("NCBI BLAST+ was not found. Install blastp and makeblastdb or use Docker.")
    prefix, output = workdir / "virasyn_ref", workdir / "blastp_hits.tsv"
    fields = "qseqid sseqid pident length qlen slen evalue bitscore"
    subprocess.run(
        [makeblastdb, "-in", str(reference_fasta), "-dbtype", "prot", "-out", str(prefix)],
        check=True, capture_output=True, text=True,
    )
    subprocess.run(
        [blastp, "-query", str(query_path), "-db", str(prefix),
         "-outfmt", f"6 {fields}", "-out", str(output), "-max_target_seqs", "50"],
        check=True, capture_output=True, text=True,
    )
    columns = fields.split()
    if not output.exists() or output.stat().st_size == 0:
        return pd.DataFrame(columns=columns + ["query_coverage_percent"])
    results = pd.read_csv(output, sep="\t", names=columns)
    results["query_coverage_percent"] = (100 * results["length"] / results["qlen"]).round(2)
    return results.sort_values(["bitscore", "evalue"], ascending=[False, True])


def make_report(virus_group: str, sequences: list[dict], structures: list[dict]) -> str:
    payload = {
        "platform": "VIRASYN",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "virus_group_context": virus_group,
        "sequences": sequences,
        "structures": structures,
        "limitations": [
            "No validated host-range or receptor-binding model is included.",
            "Sequence similarity is not proof of biological function or infection.",
            "Different virus groups require separate biologically appropriate models.",
            "PDB processing currently extracts metadata only.",
        ],
    }
    return json.dumps(payload, indent=2)
