import pytest
from virasyn.core import parse_fasta, sequence_summary


def test_protein_fasta():
    records = parse_fasta(">protein_1 example\nMKWVTFISLL\n")
    assert len(records) == 1
    assert records[0]["id"] == "protein_1"
    assert records[0]["type"] == "protein"


def test_nucleotide_fasta():
    records = parse_fasta(">genome_1\nATGCGCAT\n")
    assert records[0]["type"] == "nucleotide"


def test_invalid_fasta():
    with pytest.raises(ValueError):
        parse_fasta("not a fasta file")


def test_gc_summary():
    assert sequence_summary("ATGC")["gc_percent"] == 50.0
