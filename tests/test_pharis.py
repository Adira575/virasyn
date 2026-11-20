from unittest.mock import Mock

import pandas as pd
import pytest

from virasyn import pharis


def test_nucleotide_fasta_validation():
    records = pharis._read_nucleotide_fasta(">genome1\nATGAAATAG\n")
    assert records[0].id == "genome1"
    with pytest.raises(ValueError):
        pharis._read_nucleotide_fasta(">bad\nMALFORMED\n")


def test_reference_manifest():
    manifest = pharis.reference_fasta_manifest(">RBP1 example\n" + "M" * 100 + "\n")
    assert manifest.loc[0, "subject_id"] == "RBP1"
    assert manifest.loc[0, "length_aa"] == 100


def test_reference_fasta_rejects_invalid_symbols():
    with pytest.raises(ValueError):
        pharis._read_reference_fasta(">bad\n" + "M" * 20 + "!" + "\n")


def test_analysis_returns_thresholded_hits(monkeypatch):
    fake_orf = {
        "orf_id": "genome1_ORF0001", "sequence_id": "genome1", "strand": "+",
        "frame": 1, "start_1based": 1, "end_1based": 303,
        "length_nt_including_stop": 303, "length_aa": 100,
        "start_codon": "ATG", "stop_codon": "TAA",
        "protein_sequence": "M" * 100, "nucleotide_sequence": "ATG" + "AAA" * 99 + "TAA",
    }
    monkeypatch.setattr(pharis, "predict_six_frame_orfs", lambda *args, **kwargs: [fake_orf])
    monkeypatch.setattr(pharis.shutil, "which", lambda name: f"/usr/bin/{name}")

    def fake_run(command, **kwargs):
        if command[0].endswith("blastp"):
            output_path = command[command.index("-out") + 1]
            pd.DataFrame([[
                "genome1_ORF0001", "RBP1", 95, 100, 1, 100, 100,
                1, 100, 100, "1e-30", 200
            ]]).to_csv(output_path, sep="\t", header=False, index=False)
        return Mock(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(pharis.subprocess, "run", fake_run)
    result = pharis.analyse_pharis_genome(
        ">genome1\n" + "A" * 303 + "\n",
        ">RBP1 reference\n" + "M" * 100 + "\n",
        "subject_id,cluster\nRBP1,cluster_A\n",
        min_orf_aa=90, min_identity=80, min_query_coverage=80,
    )
    assert result["summary"]["candidate_orfs"] == 1
    assert result["candidate_hits"].iloc[0]["cluster"] == "cluster_A"
