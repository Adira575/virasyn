from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from io import StringIO
import json
from pathlib import Path
import shutil
import subprocess

import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

NUCLEOTIDE_ALPHABET = set("ACGTUNRYKMSWBDHV")
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta(text: str, filename: str = "") -> list[dict]:
    """Parse FASTA and classify using file extension when informative, otherwise alphabet."""
    parsed = list(SeqIO.parse(StringIO(text), "fasta"))
    if not parsed:
        raise ValueError("No FASTA records found. Check for a valid >header.")
    suffix = Path(filename).suffix.lower()
    hint = "protein" if suffix in {".faa", ".pep", ".aa"} else (
        "nucleotide" if suffix in {".fna", ".ffn"} else None
    )
    records = []
    for record in parsed:
        sequence = "".join(str(record.seq).upper().replace("*", "").split())
        if not sequence:
            continue
        symbols = set(sequence)
        if hint:
            kind = hint
        elif symbols <= NUCLEOTIDE_ALPHABET:
            kind = "nucleotide"
        elif symbols <= PROTEIN_ALPHABET:
            kind = "protein"
        else:
            invalid = sorted(symbols - NUCLEOTIDE_ALPHABET - PROTEIN_ALPHABET)
            raise ValueError(f"Sequence {record.id!r} contains unsupported symbols: {invalid}")
        allowed = NUCLEOTIDE_ALPHABET if kind == "nucleotide" else PROTEIN_ALPHABET
        invalid = sorted(symbols - allowed)
        if invalid:
            raise ValueError(f"Sequence {record.id!r} treated as {kind} but has symbols {invalid}. Check its file extension.")
        records.append({"id": record.id, "description": record.description, "sequence": sequence, "type": kind})
    if not records:
        raise ValueError("FASTA contained no non-empty sequences.")
    return records


def sequence_summary(sequence: str) -> dict:
    """Return general-purpose sequence composition metrics."""
    sequence = sequence.upper()
    is_nucleotide = set(sequence) <= NUCLEOTIDE_ALPHABET
    return {
        "length": len(sequence),
        "sequence_type": "nucleotide" if is_nucleotide else "protein",
        "gc_percent": round(100 * (sequence.count("G") + sequence.count("C")) / len(sequence), 2)
        if is_nucleotide and sequence else None,
        "ambiguous_symbols": sorted(set(sequence) - (NUCLEOTIDE_ALPHABET if is_nucleotide else PROTEIN_ALPHABET)),
        "most_common_symbols": Counter(sequence).most_common(10),
    }


def protein_properties(sequence: str) -> dict:
    """Compute sequence-derived properties for standard amino acids."""
    sequence = sequence.upper().replace("*", "")
    standard = "".join(letter for letter in sequence if letter in STANDARD_AA)
    result = {"length_aa": len(sequence), "analysed_standard_residues": len(standard)}
    if len(standard) < 5:
        result["note"] = "Too few standard residues for stable property estimates."
        return result
    analysis = ProteinAnalysis(standard)
    result.update({
        "molecular_weight_da": round(analysis.molecular_weight(), 2),
        "isoelectric_point": round(analysis.isoelectric_point(), 2),
        "instability_index": round(analysis.instability_index(), 2),
        "gravy": round(analysis.gravy(), 3),
        "aromaticity": round(analysis.aromaticity(), 3),
    })
    if len(standard) != len(sequence):
        result["note"] = "Non-standard residues were omitted from physicochemical estimates."
    return result


def predict_six_frame_orfs(
    nucleotide_sequence: str,
    sequence_id: str = "genome",
    min_aa: int = 45,
    genetic_code: int = 11,
    alternative_starts: bool = True,
    max_orfs: int = 20000,
) -> list[dict]:
    """
    Screen all six reading frames for start-to-stop ORFs.

    Coordinates are 1-based inclusive on the submitted sequence. Nested starts
    are retained. This is a screening method, not a complete gene caller for
    spliced or otherwise complex viral genomes.
    """
    dna = nucleotide_sequence.upper().replace("U", "T")
    if not dna or set(dna) - NUCLEOTIDE_ALPHABET:
        raise ValueError("ORF screening requires a nucleotide sequence using IUPAC DNA/RNA symbols.")
    if min_aa < 1:
        raise ValueError("Minimum protein length must be at least 1 amino acid.")
    try:
        table = CodonTable.unambiguous_dna_by_id[int(genetic_code)]
    except KeyError as exc:
        raise ValueError(f"Unsupported NCBI genetic code table: {genetic_code}") from exc

    starts = set(table.start_codons) if alternative_starts else {"ATG"}
    stops = set(table.stop_codons)
    genome_length = len(dna)
    results = []

    for strand, oriented in (("+", dna), ("-", str(Seq(dna).reverse_complement()))):
        for frame in range(3):
            active_starts = []
            for pos in range(frame, len(oriented) - 2, 3):
                codon = oriented[pos:pos + 3]
                if any(base not in "ACGT" for base in codon):
                    active_starts = []
                    continue
                if codon in starts:
                    active_starts.append(pos)
                if codon in stops:
                    for start_pos in active_starts:
                        coding_nt = oriented[start_pos:pos]
                        protein = str(Seq(coding_nt).translate(table=int(genetic_code)))
                        if protein and coding_nt[:3] in starts:
                            protein = "M" + protein[1:]
                        if len(protein) < min_aa:
                            continue
                        if strand == "+":
                            start_1based, end_1based = start_pos + 1, pos + 3
                        else:
                            start_1based = genome_length - (pos + 3) + 1
                            end_1based = genome_length - start_pos
                        results.append({
                            "sequence_id": sequence_id,
                            "strand": strand,
                            "frame": frame + 1 if strand == "+" else -(frame + 1),
                            "start_1based": start_1based,
                            "end_1based": end_1based,
                            "length_nt_including_stop": pos + 3 - start_pos,
                            "length_aa": len(protein),
                            "start_codon": oriented[start_pos:start_pos + 3],
                            "stop_codon": codon,
                            "protein_sequence": protein,
                            "nucleotide_sequence": oriented[start_pos:pos + 3],
                        })
                        if len(results) >= max_orfs:
                            break
                    active_starts = []
                    if len(results) >= max_orfs:
                        break
            if len(results) >= max_orfs:
                break
        if len(results) >= max_orfs:
            break

    results.sort(key=lambda item: (item["start_1based"], item["end_1based"], item["strand"]))
    for index, record in enumerate(results, start=1):
        record["orf_id"] = f"{sequence_id}_ORF{index:04d}"
    return results


def parse_structure(text: str, filename: str = "uploaded.pdb") -> dict:
    """Extract protein-chain metadata from PDB or mmCIF coordinate files."""
    suffix = Path(filename).suffix.lower()
    parser = MMCIFParser(QUIET=True) if suffix in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    structure = parser.get_structure("uploaded", StringIO(text))
    chains, residue_count, atom_count = [], 0, 0
    models = list(structure.get_models())
    if models:
        for chain in models[0]:
            residues = [res for res in chain if is_aa(res, standard=False)]
            if residues:
                chains.append(chain.id)
                residue_count += len(residues)
                atom_count += sum(len(list(res.get_atoms())) for res in residues)
    return {
        "file": filename,
        "format": suffix.lstrip(".").upper(),
        "model_count": len(models),
        "protein_chains": ", ".join(chains),
        "amino_acid_residue_count": residue_count,
        "protein_atom_count": atom_count,
    }


def render_structure_html(text: str) -> str:
    """Build an embeddable 3Dmol viewer for an uploaded PDB coordinate file."""
    try:
        import py3Dmol
    except ImportError as exc:
        raise RuntimeError("Install py3Dmol to enable the interactive structure viewer.") from exc
    view = py3Dmol.view(width=900, height=520)
    view.addModel(text, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.addStyle({"hetflag": True}, {"stick": {"radius": 0.18}})
    view.zoomTo()
    return view._make_html()


def run_blastp(query_path: Path, reference_fasta: Path, workdir: Path) -> pd.DataFrame:
    """Run BLASTP against a temporary database built from a supplied FASTA."""
    blastp, makeblastdb = shutil.which("blastp"), shutil.which("makeblastdb")
    if not blastp or not makeblastdb:
        raise RuntimeError("NCBI BLAST+ was not found. Install blastp and makeblastdb or use Docker.")
    prefix, output = workdir / "virasyn_protein_ref", workdir / "blastp_hits.tsv"
    fields = "qseqid sseqid pident length qlen slen qcovhsp evalue bitscore stitle"
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
        return pd.DataFrame(columns=columns)
    results = pd.read_csv(output, sep="\t", names=columns)
    return results.sort_values(["bitscore", "evalue"], ascending=[False, True])


def run_hmmscan(query_path: Path, hmm_path: Path, workdir: Path) -> pd.DataFrame:
    """Scan proteins against a user-supplied HMM profile library with HMMER."""
    hmmscan, hmmpress = shutil.which("hmmscan"), shutil.which("hmmpress")
    if not hmmscan or not hmmpress:
        raise RuntimeError("HMMER hmmscan and hmmpress are required for profile scanning.")
    database = workdir / "profiles.hmm"
    database.write_bytes(hmm_path.read_bytes())
    subprocess.run([hmmpress, "-f", str(database)], check=True, capture_output=True, text=True)
    domtbl = workdir / "hmmscan.domtbl"
    subprocess.run(
        [hmmscan, "--noali", "--domtblout", str(domtbl), str(database), str(query_path)],
        check=True, capture_output=True, text=True,
    )
    rows = []
    if domtbl.exists():
        for line in domtbl.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split(maxsplit=22)
            if len(fields) < 22:
                continue
            rows.append({
                "profile": fields[0],
                "target_accession": fields[1],
                "protein_id": fields[3],
                "full_sequence_evalue": float(fields[6]),
                "full_sequence_score": float(fields[7]),
                "domain_independent_evalue": float(fields[12]),
                "domain_score": float(fields[13]),
                "alignment_start": int(fields[17]),
                "alignment_end": int(fields[18]),
                "description": fields[22] if len(fields) > 22 else "",
            })
    return pd.DataFrame(rows)


def make_report(
    virus_group: str,
    sequence_summaries: list[dict],
    orf_summaries: list[dict],
    protein_summaries: list[dict],
    structure_summaries: list[dict],
) -> str:
    payload = {
        "platform": "VIRASYN",
        "report_schema_version": "0.2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "virus_group_context": virus_group,
        "sequences": sequence_summaries,
        "orf_candidates": orf_summaries,
        "protein_properties": protein_summaries,
        "structures": structure_summaries,
        "limitations": [
            "Six-frame ORF screening is preliminary and is not a splice-aware or virus-specific gene caller.",
            "Sequence similarity is not proof of biological function, host range, receptor binding or infection.",
            "Different virus groups require biologically appropriate reference datasets and separately validated models.",
            "Structure intake and visualisation do not establish protein function.",
        ],
    }
    return json.dumps(payload, indent=2)
