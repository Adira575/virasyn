import numpy as np
from io import BytesIO
import zipfile
import pandas as pd
import pytest

from virasyn.qpcr import (
    analyse_relative_expression,
    calculate_standard_curves,
    normalise_qpcr_wells,
    read_qpcr_table,
)


def _sample_frame():
    return pd.DataFrame([
        {"Sample": "C1", "Condition": "Control", "Target": "ACTB", "Cq": 20.0},
        {"Sample": "C1", "Condition": "Control", "Target": "GeneA", "Cq": 25.0},
        {"Sample": "C1", "Condition": "Control", "Target": "GeneA", "Cq": 25.2},
        {"Sample": "C2", "Condition": "Control", "Target": "ACTB", "Cq": 20.1},
        {"Sample": "C2", "Condition": "Control", "Target": "GeneA", "Cq": 25.1},
        {"Sample": "T1", "Condition": "Treatment", "Target": "ACTB", "Cq": 20.0},
        {"Sample": "T1", "Condition": "Treatment", "Target": "GeneA", "Cq": 22.0},
        {"Sample": "T2", "Condition": "Treatment", "Target": "ACTB", "Cq": 20.2},
        {"Sample": "T2", "Condition": "Treatment", "Target": "GeneA", "Cq": 22.1},
        {"Sample": "NTC_1", "Condition": "NTC", "Target": "GeneA", "Cq": 35.0},
    ])


def test_qpcr_import_keeps_undetermined_as_missing_not_zero():
    frame = pd.DataFrame([
        {"Sample": "S1", "Target": "GeneA", "Cq": "Undetermined"},
        {"Sample": "S2", "Target": "GeneA", "Cq": "28.4"},
    ])
    wells = normalise_qpcr_wells(frame, "Sample", "Target", "Cq")
    assert pd.isna(wells.loc[0, "Cq"])
    assert wells.loc[0, "detected"] is False or not bool(wells.loc[0, "detected"])
    assert wells.loc[1, "Cq"] == pytest.approx(28.4)
    assert wells.loc[0, "source_row"] == 2
    assert wells.loc[0, "raw_Cq"] == "Undetermined"


def test_efficiency_corrected_relative_expression_and_ntc_exclusion():
    wells = normalise_qpcr_wells(
        _sample_frame(), "Sample", "Target", "Cq", "Condition"
    )
    result = analyse_relative_expression(
        wells,
        reference_genes=["ACTB"],
        calibrator_condition="Control",
        efficiency_percent={"ACTB": 100.0, "GeneA": 100.0},
    )
    summary = result["group_summary"]
    treatment_gene = summary[
        summary["target"].eq("GeneA") & summary["condition"].eq("Treatment")
    ].iloc[0]
    assert treatment_gene["fold_change_vs_calibrator"] == pytest.approx(8.574187700290347, rel=0.001)
    assert result["ntc_qc"].iloc[0]["NTC_wells"] == 1
    assert result["ntc_qc"].iloc[0]["NTC_detected_wells"] == 1
    assert not result["sample_level"]["sample_id"].str.startswith("NTC").any()


def test_standard_curve_estimates_efficiency_from_slope():
    frame = pd.DataFrame({
        "Target": ["GeneA"] * 5,
        "Cq": [20.0, 23.3219, 26.6438, 29.9657, 33.2876],
        "Concentration": [1, 0.1, 0.01, 0.001, 0.0001],
    })
    result = calculate_standard_curves(frame, "Target", "Cq", "Concentration")
    row = result.iloc[0]
    assert row["efficiency_percent"] == pytest.approx(100.0, abs=1.0)
    assert row["R_squared"] > 0.999


def test_standard_curve_requires_three_distinct_concentrations():
    frame = pd.DataFrame({
        "Target": ["GeneA"] * 2,
        "Cq": [20.0, 21.0],
        "Concentration": [1, 1],
    })
    result = calculate_standard_curves(frame, "Target", "Cq", "Concentration")
    assert pd.isna(result.iloc[0]["efficiency_percent"])
    assert "at least 3" in result.iloc[0]["status"]


def test_sample_id_cannot_span_multiple_conditions():
    frame = pd.DataFrame([
        {"Sample": "S1", "Condition": "Control", "Target": "ACTB", "Cq": 20.0},
        {"Sample": "S1", "Condition": "Treatment", "Target": "GeneA", "Cq": 25.0},
    ])
    wells = normalise_qpcr_wells(frame, "Sample", "Target", "Cq", "Condition")
    with pytest.raises(ValueError, match="one condition"):
        analyse_relative_expression(
            wells, ["ACTB"], "Control", {"ACTB": 100.0, "GeneA": 100.0}
        )


def test_csv_loader_reads_instrument_export():
    content = b"Sample,Target,Cq\nS1,ACTB,20.1\nS1,GeneA,25.4\n"
    frame = read_qpcr_table(content, "plate.csv")
    assert list(frame.columns) == ["Sample", "Target", "Cq"]
    assert len(frame) == 2

def test_single_technical_well_is_flagged_as_unreplicated():
    frame = pd.DataFrame([
        {"Sample": "C1", "Condition": "Control", "Target": "ACTB", "Cq": 20.0},
        {"Sample": "C1", "Condition": "Control", "Target": "GeneA", "Cq": 25.0},
    ])
    wells = normalise_qpcr_wells(frame, "Sample", "Target", "Cq", "Condition")
    result = analyse_relative_expression(
        wells, ["ACTB"], "Control", {"ACTB": 100.0, "GeneA": 100.0}
    )
    gene_qc = result["technical_qc"].query("target == 'GeneA'").iloc[0]
    assert "Single detected well" in gene_qc["technical_qc"]


RDML_XML = b"""<?xml version="1.0" encoding="UTF-8"?>
<rdml version="1.3">
  <sample id="S1">
    <name>Control_1</name>
    <description>Control</description>
    <type>unkn</type>
  </sample>
  <sample id="N1">
    <name>NTC_1</name>
    <type>ntc</type>
  </sample>
  <target id="T1"><name>ACTB</name></target>
  <experiment id="Experiment1">
    <run id="Run1">
      <react id="A1">
        <sample id="S1"/>
        <data>
          <tar id="T1"/>
          <cq>20.5</cq>
          <adp><fpoints>
            <point><cyc>1</cyc><fluor>0.10</fluor></point>
            <point><cyc>2</cyc><fluor>0.14</fluor></point>
          </fpoints></adp>
          <mdp><fpoints>
            <point><tmp>70.0</tmp><fluor>0.20</fluor></point>
            <point><tmp>71.0</tmp><fluor>0.15</fluor></point>
          </fpoints></mdp>
        </data>
      </react>
      <react id="A2">
        <sample id="N1"/>
        <data><tar id="T1"/><cq>-1</cq></data>
      </react>
    </run>
  </experiment>
</rdml>
"""


def test_rdml_xml_maps_samples_targets_cq_and_curves():
    frame = read_qpcr_table(RDML_XML, "run.rdml")
    control = frame[frame["Sample ID"].eq("S1")].iloc[0]
    ntc = frame[frame["Sample ID"].eq("N1")].iloc[0]
    assert control["Sample"] == "Control_1"
    assert control["Sample description"] == "Control"
    assert control["Target"] == "ACTB"
    assert control["Cq"] == "20.5"
    assert len(pd.read_json(control["Amplification trace (JSON)"])) == 2
    assert len(pd.read_json(control["Melt curve trace (JSON)"])) == 2
    assert ntc["Cq"] == "Undetermined"
    assert ntc["Sample role"] == "NTC"


def test_rdml_zip_package_is_supported():
    archive_bytes = BytesIO()
    with zipfile.ZipFile(archive_bytes, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("rdml_data.xml", RDML_XML)
    frame = read_qpcr_table(archive_bytes.getvalue(), "instrument_run.rdml")
    assert len(frame) == 2
    assert set(frame["Target"]) == {"ACTB"}
