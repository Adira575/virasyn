import pytest
from Bio.Seq import Seq

from virasyn.core import (
    parse_fasta,
    predict_six_frame_orfs,
    protein_properties,
    sequence_summary,
)


def test_protein_fasta_extension_overrides_ambiguous_alphabet():
    records = parse_fasta(">protein_1 example\nACGTACGT\n", "sample.faa")
    assert records[0]["type"] == "protein"


def test_nucleotide_fasta():
    records = parse_fasta(">genome_1\nATGCGCAT\n", "sample.fna")
    assert records[0]["type"] == "nucleotide"


def test_invalid_fasta():
    with pytest.raises(ValueError):
        parse_fasta("not a fasta file")


def test_gc_summary():
    assert sequence_summary("ATGC")["gc_percent"] == 50.0


def test_six_frame_orf_finds_forward_orf():
    genome = "ATG" + "GCT" * 30 + "TAA"
    hits = predict_six_frame_orfs(
        genome, sequence_id="toy", min_aa=30, genetic_code=11
    )
    forward = [hit for hit in hits if hit["strand"] == "+" and hit["frame"] == 1]
    assert forward
    assert forward[0]["start_1based"] == 1
    assert forward[0]["end_1based"] == len(genome)
    assert forward[0]["length_aa"] == 31
    assert forward[0]["protein_sequence"].startswith("M")


def test_six_frame_orf_finds_reverse_orf():
    coding = "ATG" + "GCT" * 30 + "TAA"
    reverse_strand_genome = str(Seq(coding).reverse_complement())
    hits = predict_six_frame_orfs(
        reverse_strand_genome, sequence_id="reverse_toy", min_aa=30, genetic_code=11
    )
    reverse = [hit for hit in hits if hit["strand"] == "-" and hit["frame"] == -1]
    assert reverse
    assert reverse[0]["start_1based"] == 1
    assert reverse[0]["end_1based"] == len(reverse_strand_genome)


def test_protein_properties():
    props = protein_properties("MKWVTFISLLFLFSSAYS")
    assert props["length_aa"] == 18
    assert "molecular_weight_da" in props
