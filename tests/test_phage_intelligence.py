import pandas as pd
import pytest

from virasyn.phage_intelligence import (
    build_integrated_report, csv_template, score_host_evidence,
)

def test_host_evidence_ranks_experimental_support_first():
    data = (
        "phage_id,host_taxid,host_name,evidence_type,evidence_level,source,citation\n"
        "P1,1280,Staphylococcus aureus,productive_infection,experimental,study,doi:1\n"
        "P1,562,Escherichia coli,sequence_similarity,computational,tool,doi:2\n"
    )
    ranked, summary = score_host_evidence(data)
    assert ranked.iloc[0]["host_taxid"] == "1280"
    assert ranked.iloc[0]["experimental_records"] == 1
    assert summary["index_is_probability"] is False

def test_unknown_evidence_level_is_rejected():
    data = (
        "phage_id,host_taxid,host_name,evidence_type,evidence_level,source,citation\n"
        "P1,1280,Staphylococcus aureus,host_range,guess,db,ref\n"
    )
    with pytest.raises(ValueError, match="Unsupported evidence_level"):
        score_host_evidence(data)

def test_receptor_annotation_maps_by_exact_reference_id():
    hits = pd.DataFrame([{"sseqid": "RBP_1", "pident": 91.0}])
    receptor = (
        "subject_id,receptor,receptor_class,host_taxid,host_name,evidence_type,citation\n"
        "RBP_1,WTA,cell wall polymer,1280,Staphylococcus aureus,receptor_binding,doi:1\n"
    )
    result = build_integrated_report(hits, receptor_annotation_csv=receptor)
    assert result["rbp_receptor_evidence"].iloc[0]["receptor"] == "WTA"

def test_templates_have_expected_headers():
    assert "host_taxid" in csv_template("host")
    assert "subject_id" in csv_template("receptor")
