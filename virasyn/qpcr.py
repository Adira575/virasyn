from __future__ import annotations

from io import BytesIO
from pathlib import Path
import json
import xml.etree.ElementTree as ET
import zipfile
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


NON_DETECT_TOKENS = {
    "", "undetermined", "undetected", "not detected", "no ct", "no cq",
    "na", "n/a", "nan", "none", "-", "--",
}


def _rdml_local_name(tag: str) -> str:
    return str(tag).rsplit("}", 1)[-1]


def _rdml_children(element: ET.Element, name: str) -> list[ET.Element]:
    return [child for child in list(element) if _rdml_local_name(child.tag) == name]


def _rdml_text(element: ET.Element, name: str) -> str:
    nodes = _rdml_children(element, name)
    if not nodes:
        return ""
    return "".join(nodes[0].itertext()).strip()


def _rdml_trace(data_element: ET.Element, trace_name: str) -> str:
    trace_nodes = _rdml_children(data_element, trace_name)
    if not trace_nodes:
        return ""
    root = trace_nodes[0]
    points: list[dict[str, float]] = []
    x_name, y_name = ("cyc", "fluor") if trace_name == "adp" else ("tmp", "fluor")
    for node in root.iter():
        if _rdml_local_name(node.tag) != "point":
            continue
        x_value = node.attrib.get(x_name) or _rdml_text(node, x_name)
        y_value = node.attrib.get(y_name) or _rdml_text(node, y_name)
        try:
            points.append({x_name: float(x_value), "fluor": float(y_value)})
        except (TypeError, ValueError):
            continue
    return json.dumps(points, separators=(",", ":")) if points else ""


def _read_rdml(data: bytes) -> pd.DataFrame:
    """Read common RDML 1.x XML exports, including XML packaged inside a ZIP."""
    xml_payload: bytes | None = None
    if zipfile.is_zipfile(BytesIO(data)):
        with zipfile.ZipFile(BytesIO(data)) as archive:
            candidates = [
                name for name in archive.namelist()
                if name.lower().endswith(".xml") and not name.startswith("__MACOSX/")
            ]
            for name in candidates:
                candidate = archive.read(name)
                if b"<rdml" in candidate[:4096].lower():
                    xml_payload = candidate
                    break
            if xml_payload is None:
                raise ValueError("The RDML archive did not contain an XML file with an <rdml> root.")
    else:
        xml_payload = data

    try:
        root = ET.fromstring(xml_payload)
    except ET.ParseError as exc:
        raise ValueError(f"Could not parse RDML XML: {exc}") from exc
    if _rdml_local_name(root.tag) != "rdml":
        raise ValueError("The uploaded XML does not have an RDML root element.")

    sample_info: dict[str, dict[str, str]] = {}
    target_names: dict[str, str] = {}
    for node in root.iter():
        local = _rdml_local_name(node.tag)
        node_id = node.attrib.get("id", "")
        if local == "sample" and node_id and _rdml_children(node, "name"):
            raw_role = ""
            for role_node in _rdml_children(node, "type"):
                raw_role = (
                    role_node.attrib.get("type")
                    or "".join(role_node.itertext()).strip()
                    or raw_role
                )
            role_aliases = {
                "ntc": "NTC",
                "std": "Standard",
                "pos": "Positive control",
                "unkn": "Unknown sample",
                "opt": "Optical calibrator",
            }
            quantities: dict[str, str] = {}
            for quantity_node in _rdml_children(node, "quantity"):
                target_id = quantity_node.attrib.get("targetId", "")
                value = _rdml_text(quantity_node, "value")
                if not value:
                    value = "".join(quantity_node.itertext()).strip()
                if value:
                    quantities[target_id] = value
            sample_info[node_id] = {
                "name": _rdml_text(node, "name") or node_id,
                "description": _rdml_text(node, "description"),
                "role": role_aliases.get(raw_role.lower(), raw_role),
                "quantities": json.dumps(quantities),
            }
        elif local == "target" and node_id and _rdml_children(node, "name"):
            target_names[node_id] = _rdml_text(node, "name") or node_id

    rows: list[dict[str, Any]] = []
    for experiment in [node for node in root.iter() if _rdml_local_name(node.tag) == "experiment"]:
        experiment_id = experiment.attrib.get("id", "")
        for run in _rdml_children(experiment, "run"):
            run_id = run.attrib.get("id", "")
            for reaction in _rdml_children(run, "react"):
                reaction_id = reaction.attrib.get("id", "")
                sample_refs = _rdml_children(reaction, "sample")
                sample_id = sample_refs[0].attrib.get("id", "") if sample_refs else ""
                sample = sample_info.get(sample_id, {
                    "name": sample_id or "Unnamed sample",
                    "description": "",
                    "role": "",
                    "quantities": "{}",
                })
                try:
                    sample_quantities = json.loads(sample.get("quantities", "{}"))
                except (TypeError, json.JSONDecodeError):
                    sample_quantities = {}

                for data_element in _rdml_children(reaction, "data"):
                    target_refs = _rdml_children(data_element, "tar")
                    target_id = target_refs[0].attrib.get("id", "") if target_refs else ""
                    target = target_names.get(target_id, target_id or "Unknown target")
                    cq_value = _rdml_text(data_element, "cq")
                    if cq_value in {"-1", "-1.0"}:
                        cq_value = "Undetermined"
                    known_quantity = sample_quantities.get(target_id, sample_quantities.get("", ""))
                    rows.append({
                        "Sample": sample.get("name", sample_id or "Unnamed sample"),
                        "Sample ID": sample_id,
                        "Sample description": sample.get("description", ""),
                        "Sample role": sample.get("role", ""),
                        "Target": target,
                        "Target ID": target_id,
                        "Cq": cq_value,
                        "Known concentration": known_quantity,
                        "Well": reaction_id,
                        "Run": run_id,
                        "Experiment": experiment_id,
                        "Amplification trace (JSON)": _rdml_trace(data_element, "adp"),
                        "Melt curve trace (JSON)": _rdml_trace(data_element, "mdp"),
                    })

    if not rows:
        raise ValueError("No reaction/target Cq records were found in this RDML file.")
    return pd.DataFrame(rows)


def read_qpcr_table(
    data: bytes,
    filename: str,
    sheet_name: str | int | None = None,
) -> pd.DataFrame:
    """Read a tabular instrument export. Raw fluorescence traces are not Cq tables."""
    suffix = Path(filename).suffix.lower()
    if suffix == ".rdml":
        frame = _read_rdml(data)
    elif suffix in {".xlsx", ".xlsm"}:
        frame = pd.read_excel(BytesIO(data), sheet_name=sheet_name or 0, engine="openpyxl")
    elif suffix == ".tsv":
        frame = pd.read_csv(BytesIO(data), sep="\t")
    elif suffix in {".csv", ".txt"}:
        frame = pd.read_csv(BytesIO(data), sep=None, engine="python")
    else:
        raise ValueError("Supported qPCR exports are RDML, CSV, TSV, TXT, XLSX and XLSM.")
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        raise ValueError("The selected file or worksheet contains no tabular data.")
    frame.columns = [str(column).strip() for column in frame.columns]
    if len(set(frame.columns)) != len(frame.columns):
        raise ValueError("Column names must be unique. Rename duplicate headers in the export.")
    return frame


def detect_column(columns: list[str], candidates: tuple[str, ...]) -> str | None:
    """Choose a likely column by common qPCR header names."""
    lookup = {str(column).strip().lower(): str(column) for column in columns}
    for candidate in candidates:
        if candidate in lookup:
            return lookup[candidate]
    for candidate in candidates:
        for lower, original in lookup.items():
            if candidate in lower:
                return original
    return None


def _numeric_cq(values: pd.Series) -> pd.Series:
    text = values.astype("string").str.strip()
    text_lower = text.str.lower()
    numeric = pd.to_numeric(text.mask(text_lower.isin(NON_DETECT_TOKENS)), errors="coerce")
    # Zero, negative and non-numeric Cq entries are invalid/missing, never treated as zero copies.
    return numeric.where(numeric > 0, np.nan)


def normalise_qpcr_wells(
    frame: pd.DataFrame,
    sample_col: str,
    target_col: str,
    cq_col: str,
    condition_col: str | None = None,
    role_col: str | None = None,
) -> pd.DataFrame:
    """Map instrument-specific headings to a standard well-level table."""
    for column in (sample_col, target_col, cq_col):
        if column not in frame.columns:
            raise ValueError(f"Required column not found: {column}")
    if condition_col and condition_col not in frame.columns:
        raise ValueError(f"Condition column not found: {condition_col}")
    if role_col and role_col not in frame.columns:
        raise ValueError(f"Sample-role column not found: {role_col}")

    result = pd.DataFrame({
        "source_row": np.arange(2, len(frame) + 2),
        "raw_Cq": frame[cq_col].astype("string"),
        "sample_id": frame[sample_col].astype("string").str.strip(),
        "target": frame[target_col].astype("string").str.strip(),
        "Cq": _numeric_cq(frame[cq_col]),
    })
    if condition_col:
        result["condition"] = frame[condition_col].astype("string").str.strip()
    else:
        result["condition"] = "All samples"

    if role_col:
        result["sample_role"] = frame[role_col].astype("string").str.strip()
    else:
        result["sample_role"] = ""

    result = result[
        result["sample_id"].notna()
        & result["target"].notna()
        & result["sample_id"].ne("")
        & result["target"].ne("")
    ].copy()
    if result.empty:
        raise ValueError("No rows with both a sample identifier and target/gene name were found.")

    role_text = result["sample_role"].fillna("").str.lower()
    sample_text = result["sample_id"].fillna("").str.lower()
    result["is_ntc"] = (
        role_text.str.contains(r"ntc|no[- ]template", regex=True)
        | sample_text.str.contains(r"^ntc(?:$|[_\-\s0-9])|no[- ]template control", regex=True)
    )
    result["detected"] = result["Cq"].notna()
    result["well_status"] = np.select(
        [result["is_ntc"] & result["detected"], result["is_ntc"], result["detected"]],
        ["NTC amplified — investigate contamination", "NTC", "Detected"],
        default="Undetermined / no Cq",
    )
    return result.reset_index(drop=True)


def calculate_standard_curves(
    frame: pd.DataFrame,
    target_col: str,
    cq_col: str,
    concentration_col: str,
) -> pd.DataFrame:
    """Fit Cq against log10(concentration) for each assay and estimate efficiency."""
    for column in (target_col, cq_col, concentration_col):
        if column not in frame.columns:
            raise ValueError(f"Standard-curve column not found: {column}")
    rows: list[dict[str, Any]] = []
    work = frame[[target_col, cq_col, concentration_col]].copy()
    work.columns = ["target", "Cq", "concentration"]
    work["target"] = work["target"].astype("string").str.strip()
    work["Cq"] = _numeric_cq(work["Cq"])
    work["concentration"] = pd.to_numeric(work["concentration"], errors="coerce")
    work = work[
        work["target"].notna()
        & work["concentration"].gt(0)
        & work["Cq"].notna()
    ]

    for target, group in work.groupby("target"):
        x = np.log10(group["concentration"].to_numpy(dtype=float))
        y = group["Cq"].to_numpy(dtype=float)
        unique_levels = int(group["concentration"].nunique())
        row: dict[str, Any] = {
            "target": str(target),
            "n_wells": int(len(group)),
            "dilution_levels": unique_levels,
            "slope": np.nan,
            "intercept": np.nan,
            "R_squared": np.nan,
            "efficiency_percent": np.nan,
            "status": "Insufficient dilution levels (need at least 3)",
        }
        if unique_levels >= 3 and len(group) >= 3:
            slope, intercept = np.polyfit(x, y, 1)
            predicted = slope * x + intercept
            ss_res = float(np.sum((y - predicted) ** 2))
            ss_tot = float(np.sum((y - np.mean(y)) ** 2))
            r_squared = 1.0 if ss_tot == 0 and ss_res == 0 else (
                1 - ss_res / ss_tot if ss_tot > 0 else np.nan
            )
            efficiency = (10 ** (-1 / slope) - 1) * 100 if slope < 0 else np.nan
            status = "OK"
            if not np.isfinite(efficiency):
                status = "Invalid slope; efficiency not calculated"
            elif efficiency < 90 or efficiency > 110:
                status = "Review efficiency outside common 90–110% range"
            if np.isfinite(r_squared) and r_squared < 0.98:
                status = "Review standard-curve linearity (R² < 0.98)"
            row.update({
                "slope": float(slope),
                "intercept": float(intercept),
                "R_squared": float(r_squared),
                "efficiency_percent": float(efficiency) if np.isfinite(efficiency) else np.nan,
                "status": status,
            })
        rows.append(row)
    return pd.DataFrame(rows)


def _benjamini_hochberg(p_values: list[float]) -> list[float]:
    adjusted = [np.nan] * len(p_values)
    valid = [(i, p) for i, p in enumerate(p_values) if np.isfinite(p)]
    if not valid:
        return adjusted
    ordered = sorted(valid, key=lambda item: item[1])
    total = len(ordered)
    running = 1.0
    for rank_index in range(total - 1, -1, -1):
        original_index, p_value = ordered[rank_index]
        rank = rank_index + 1
        running = min(running, p_value * total / rank)
        adjusted[original_index] = min(running, 1.0)
    return adjusted


def analyse_relative_expression(
    wells: pd.DataFrame,
    reference_genes: list[str],
    calibrator_condition: str,
    efficiency_percent: dict[str, float],
    max_technical_sd: float = 1.0,
) -> dict[str, pd.DataFrame]:
    """
    Average technical wells within biological sample/target, create efficiency-corrected
    quantities, normalize to the geometric mean of reference-gene quantities, then compare
    group means. Welch tests use per-biological-sample log2 normalized quantities.
    """
    required = {"sample_id", "target", "condition", "Cq", "is_ntc"}
    missing = required.difference(wells.columns)
    if missing:
        raise ValueError("Normalised well table is missing: " + ", ".join(sorted(missing)))
    if not reference_genes:
        raise ValueError("Select at least one reference gene for normalization.")

    all_targets = set(wells["target"].dropna().astype(str))
    unknown_refs = set(reference_genes).difference(all_targets)
    if unknown_refs:
        raise ValueError("Reference gene(s) not found: " + ", ".join(sorted(unknown_refs)))
    condition_values = set(wells.loc[~wells["is_ntc"], "condition"].dropna().astype(str))
    if calibrator_condition not in condition_values:
        raise ValueError("The selected calibrator/control group has no non-NTC samples.")

    efficiencies: dict[str, float] = {}
    for target in all_targets:
        value = float(efficiency_percent.get(target, 100.0))
        if value <= -100 or value > 500:
            raise ValueError(f"Invalid amplification efficiency for {target}: {value}")
        efficiencies[target] = value

    sample_wells = wells.loc[~wells["is_ntc"]].copy()
    condition_counts = sample_wells.groupby("sample_id")["condition"].nunique()
    if condition_counts.gt(1).any():
        bad_samples = condition_counts[condition_counts.gt(1)].index.astype(str).tolist()
        raise ValueError(
            "Each biological sample ID must belong to one condition. Check: "
            + ", ".join(bad_samples[:10])
        )
    grouped = sample_wells.groupby(
        ["sample_id", "condition", "target"], dropna=False
    )["Cq"]
    technical = grouped.agg(
        n_technical_wells="count",
        mean_Cq="mean",
        sd_Cq="std",
        min_Cq="min",
        max_Cq="max",
    ).reset_index()
    technical["Cq_range"] = technical["max_Cq"] - technical["min_Cq"]
    technical["technical_qc"] = np.select(
        [
            technical["n_technical_wells"].eq(0).fillna(False).to_numpy(dtype=bool),
            technical["n_technical_wells"].eq(1).fillna(False).to_numpy(dtype=bool),
            technical["sd_Cq"].gt(max_technical_sd).fillna(False).to_numpy(dtype=bool),
        ],
        [
            "No detected wells",
            "Single detected well; technical variation cannot be estimated",
            f"Review: technical SD > {max_technical_sd:g} Cq",
        ],
        default="OK / no automatic exclusion",
    )

    technical["efficiency_percent"] = technical["target"].map(efficiencies)
    technical["amplification_factor"] = 1 + technical["efficiency_percent"] / 100
    technical["quantity_relative"] = np.power(
        technical["amplification_factor"], -technical["mean_Cq"]
    )
    technical.loc[technical["mean_Cq"].isna(), "quantity_relative"] = np.nan

    ref_rows = technical[
        technical["target"].isin(reference_genes)
        & technical["quantity_relative"].gt(0)
    ]
    reference_summary = ref_rows.groupby(["sample_id", "condition"]).agg(
        reference_quantity=("quantity_relative", lambda values: float(np.exp(np.log(values).mean()))),
        reference_genes_detected=("target", "nunique"),
    ).reset_index()

    sample_level = technical.merge(
        reference_summary, on=["sample_id", "condition"], how="left"
    )
    sample_level["normalized_quantity"] = (
        sample_level["quantity_relative"] / sample_level["reference_quantity"]
    )
    sample_level["log2_normalized_quantity"] = np.log2(
        sample_level["normalized_quantity"].where(sample_level["normalized_quantity"] > 0)
    )
    sample_level["normalization_status"] = "Normalized"
    sample_level.loc[
        sample_level["reference_quantity"].isna(),
        "normalization_status",
    ] = "Not normalized: no detected reference-gene Cq"
    partial_reference = (
        sample_level["reference_quantity"].notna()
        & sample_level["reference_genes_detected"].lt(len(reference_genes)).fillna(False)
    )
    sample_level.loc[
        partial_reference,
        "normalization_status",
    ] = (
        "Partial reference set: "
        + sample_level.loc[partial_reference, "reference_genes_detected"].astype(str)
        + f"/{len(reference_genes)} detected"
    )

    valid = sample_level.dropna(subset=["normalized_quantity"]).copy()
    group_summary = valid.groupby(["target", "condition"]).agg(
        n_biological_samples=("sample_id", "nunique"),
        mean_normalized_quantity=("normalized_quantity", "mean"),
        sd_normalized_quantity=("normalized_quantity", "std"),
        mean_log2_normalized_quantity=("log2_normalized_quantity", "mean"),
    ).reset_index()

    baseline = (
        group_summary[group_summary["condition"].eq(calibrator_condition)]
        [["target", "mean_normalized_quantity"]]
        .rename(columns={"mean_normalized_quantity": "calibrator_mean_quantity"})
    )
    group_summary = group_summary.merge(baseline, on="target", how="left")
    group_summary["fold_change_vs_calibrator"] = (
        group_summary["mean_normalized_quantity"] / group_summary["calibrator_mean_quantity"]
    )
    group_summary["log2_fold_change"] = np.log2(
        group_summary["fold_change_vs_calibrator"].where(
            group_summary["fold_change_vs_calibrator"] > 0
        )
    )

    p_values: list[float] = []
    for _, row in group_summary.iterrows():
        target = row["target"]
        condition = row["condition"]
        if condition == calibrator_condition:
            p_values.append(np.nan)
            continue
        control_values = valid.loc[
            valid["target"].eq(target)
            & valid["condition"].eq(calibrator_condition),
            "log2_normalized_quantity",
        ].dropna().to_numpy(dtype=float)
        treatment_values = valid.loc[
            valid["target"].eq(target)
            & valid["condition"].eq(condition),
            "log2_normalized_quantity",
        ].dropna().to_numpy(dtype=float)
        if len(control_values) >= 2 and len(treatment_values) >= 2:
            p_values.append(float(ttest_ind(
                treatment_values, control_values, equal_var=False, nan_policy="omit"
            ).pvalue))
        else:
            p_values.append(np.nan)
    group_summary["p_value_welch_log2"] = p_values
    group_summary["q_value_BH"] = _benjamini_hochberg(p_values)

    ntc = wells.loc[wells["is_ntc"]].copy()
    ntc_qc = pd.DataFrame([{
        "NTC_wells": int(len(ntc)),
        "NTC_detected_wells": int(ntc["Cq"].notna().sum()),
        "status": (
            "WARNING: one or more NTC wells have a Cq value"
            if ntc["Cq"].notna().any()
            else "No NTC amplification detected in supplied table"
        ),
    }])

    return {
        "sample_level": sample_level,
        "group_summary": group_summary,
        "technical_qc": technical,
        "ntc_qc": ntc_qc,
    }
