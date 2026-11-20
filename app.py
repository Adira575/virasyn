from __future__ import annotations

from pathlib import Path
from io import BytesIO
import json
import tempfile

import pandas as pd
import streamlit as st

from virasyn.core import (
    make_report, parse_fasta, parse_structure, predict_six_frame_orfs,
    protein_properties, render_structure_html, run_blastp, run_hmmscan,
    sequence_summary,
)
from virasyn.database import (
    fetch_sequence_fasta, search_all_sources, search_crossref, search_ena,
    search_europe_pmc, search_interpro, search_ncbi_sequences,
    search_ncbi_taxonomy, search_openalex, search_pubmed, search_rcsb,
    search_uniprot,
)
from virasyn.qpcr import (
    analyse_relative_expression, calculate_standard_curves, detect_column,
    normalise_qpcr_wells, read_qpcr_table,
)

st.set_page_config(
    page_title="VIRASYN | Viral Systems Intelligence",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
.block-container {max-width: 1480px; padding-top: 1.3rem; padding-bottom: 3rem;}
.hero {padding: 1.35rem 1.6rem; border-radius: 18px; background: linear-gradient(120deg,#0c1b2b,#123c50 68%,#126b72); color:#f4fbff; margin-bottom:1.1rem;}
.hero h1 {font-size:2.55rem; letter-spacing:.11em; margin:0;}
.hero p {font-size:1rem; opacity:.9; margin:.35rem 0 0;}
div[data-testid="stMetric"] {background:rgba(120,150,170,.08); border:1px solid rgba(120,150,170,.18); padding:12px; border-radius:12px;}
.search-box {padding:1rem; border:1px solid rgba(120,150,170,.25); border-radius:14px; background:rgba(120,150,170,.04);}
.record-title {font-size:1.04rem; font-weight:650;}
.source-tag {font-size:.78rem; letter-spacing:.04em; text-transform:uppercase; color:#17818c;}
</style>
<div class="hero">
  <h1>VIRASYN</h1>
  <p>Viral Systems Intelligence · Discover literature, genomes, proteins and structures.</p>
</div>
""", unsafe_allow_html=True)

st.caption("An integrated viral research discovery and analysis workspace.")

database_tab, analysis_tab, qpcr_tab, workbench_tab, alphafold_tab, pharis_tab, phage_tab, resources_tab, about_tab = st.tabs(
    [
        "🔎 Research database",
        "🧬 Sequence & structure analysis",
        "🧪 PCR / qPCR analysis",
        "🧰 Bioinformatics workbench",
        "🧬 AlphaFold",
        "🦠 PhARIS RBP intelligence",
        "🧫 Phage intelligence",
        "📚 Database directory",
        "ℹ️ About VIRASYN",
    ]
)
st.caption("Build: virasyn-phage-intelligence-2026-11-20-r6")

def render_record(record: dict, kind: str):
    with st.container(border=True):
        title = record.get("title") or record.get("accession") or record.get("id") or "Untitled record"
        st.write(title)
        meta = []
        for key in ("authors", "journal", "date", "organism", "length", "method", "doi",
                    "rank", "lineage", "description", "type", "cited_by_count", "is_open_access"):
            value = record.get(key)
            if value:
                meta.append(f"**{key.replace('_', ' ').title()}:** {value}")
        if meta:
            st.markdown(" · ".join(meta))
        st.caption(f"{record.get('database', kind)} · ID: {record.get('accession', record.get('id', ''))}")
        st.markdown(f"[Open original record ↗]({record.get('url', '#')})")
        if kind == "NCBI Nucleotide":
            if st.button("Download nucleotide FASTA", key=f"fna_{record.get('id')}"):
                try:
                    fasta = fetch_sequence_fasta("nuccore", record["id"])
                    st.download_button("Save FASTA", fasta, f"{record.get('accession', record['id'])}.fna",
                                       "text/plain", key=f"dl_fna_{record.get('id')}")
                except Exception as exc:
                    st.error(str(exc))
        elif kind == "NCBI Protein":
            if st.button("Download protein FASTA", key=f"faa_{record.get('id')}"):
                try:
                    fasta = fetch_sequence_fasta("protein", record["id"])
                    st.download_button("Save FASTA", fasta, f"{record.get('accession', record['id'])}.faa",
                                       "text/plain", key=f"dl_faa_{record.get('id')}")
                except Exception as exc:
                    st.error(str(exc))


with database_tab:
    st.subheader("Search viral research literature and molecular data")
    st.write(
        "Search publications, viral nucleotide records, protein records and "
        "three-dimensional structures from connected public scientific databases."
    )
    with st.form("research_search_form"):
        query = st.text_input(
            "Search terms",
            placeholder='e.g. bacteriophage receptor binding protein OR influenza hemagglutinin',
            help="Use scientific names, accession IDs, protein names, genes, receptors or PubMed-style search terms.",
        )
        source = st.selectbox(
            "Search database",
            ["All sources", "PubMed articles", "Europe PMC articles", "OpenAlex scholarly works",
             "Crossref records", "NCBI nucleotide genomes", "NCBI protein sequences",
             "NCBI taxonomy", "ENA nucleotide sequences", "RCSB PDB structures",
             "UniProtKB proteins", "InterPro entries"],
        )
        c1, c2, c3, c4 = st.columns([1, 1, 1, 1])
        with c1:
            page_size = st.selectbox("Results per page", [10, 20, 50], index=1)
        with c2:
            page_number = st.number_input("Page", min_value=1, value=1, step=1)
        with c3:
            date_from = st.text_input("Published from (PubMed only)", placeholder="YYYY-MM-DD")
        with c4:
            date_to = st.text_input("Published to (PubMed only)", placeholder="YYYY-MM-DD")
        sort = st.selectbox("Sort publications", ["Best match", "Newest first"])
        offset = int((page_number - 1) * page_size)
        submitted = st.form_submit_button("Search research databases", type="primary", use_container_width=True)

    if submitted:
        if not query.strip():
            st.warning("Enter a search term first.")
        else:
            try:
                with st.spinner("Searching scientific databases..."):
                    if source == "All sources":
                        bundle = search_all_sources(query.strip(), int(page_number), 5)
                        st.session_state["virasyn_search_bundle"] = bundle
                        st.session_state["virasyn_search_query"] = query.strip()
                        st.session_state["virasyn_search_source"] = source
                    elif source == "PubMed articles":
                        result = search_pubmed(query.strip(), offset, page_size, date_from, date_to, sort)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    elif source == "Europe PMC articles":
                        result = search_europe_pmc(query.strip(), page_number, page_size)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    elif source == "OpenAlex scholarly works":
                        result = search_openalex(query.strip(), page_number, page_size)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    elif source == "Crossref records":
                        result = search_crossref(query.strip(), offset, page_size)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    elif source == "NCBI nucleotide genomes":
                        result = search_ncbi_sequences("nuccore", query.strip(), offset, page_size)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    elif source == "NCBI protein sequences":
                        result = search_ncbi_sequences("protein", query.strip(), offset, page_size)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    elif source == "NCBI taxonomy":
                        result = search_ncbi_taxonomy(query.strip(), offset, page_size)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    elif source == "ENA nucleotide sequences":
                        result = search_ena(query.strip(), offset, page_size)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    elif source == "InterPro entries":
                        result = search_interpro(query.strip(), page_number, page_size)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    elif source == "RCSB PDB structures":
                        result = search_rcsb(query.strip(), offset, page_size)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    else:
                        result = search_uniprot(query.strip(), offset, page_size)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
            except Exception as exc:
                st.error(f"Database search failed: {exc}")

    current_source = st.session_state.get("virasyn_search_source")
    if current_source:
        st.divider()
        st.markdown(f"### Results for: {st.session_state.get('virasyn_search_query', '')}")
        if current_source == "All sources":
            bundle = st.session_state.get("virasyn_search_bundle", {})
            if bundle.get("errors"):
                with st.expander("Some sources could not be reached"):
                    for name, error in bundle["errors"].items():
                        st.warning(f"{name}: {error}")
            total_display = 0
            for name, payload in bundle.get("groups", {}).items():
                records = payload.get("results", [])
                total_display += len(records)
                with st.expander(f"{name} · {payload.get('count', len(records)):,} records", expanded=True):
                    if not records:
                        st.info("No results returned from this source.")
                    for record in records:
                        render_record(record, name)
            st.caption(f"Showing up to five results per source · {total_display} records displayed.")
        else:
            result = st.session_state.get("virasyn_search_result", {})
            records = result.get("results", [])
            st.metric("Matching records", f"{result.get('count', 0):,}")
            st.caption(f"Source: {result.get('source', current_source)} · Results are retrieved live and cached briefly.")
            if not records:
                st.info("No records found. Try a broader query or different search terms.")
            for record in records:
                kind = (
                    "NCBI Nucleotide" if current_source == "NCBI nucleotide genomes" else
                    "NCBI Protein" if current_source == "NCBI protein sequences" else
                    "NCBI Taxonomy" if current_source == "NCBI taxonomy" else
                    "Europe PMC" if current_source == "Europe PMC articles" else
                    "OpenAlex" if current_source == "OpenAlex scholarly works" else
                    "Crossref" if current_source == "Crossref records" else
                    "European Nucleotide Archive" if current_source == "ENA nucleotide sequences" else
                    "InterPro" if current_source == "InterPro entries" else
                    current_source
                )
                render_record(record, kind)


with qpcr_tab:
    st.subheader("Real-time PCR (qPCR / RT-qPCR) analysis")
    st.write(
        "Import a well-level Cq/Ct export, map the instrument-specific columns, inspect technical "
        "replicates and calculate normalized relative expression. Technical wells are grouped within "
        "each biological sample; they are not counted as independent biological replicates."
    )
    st.info(
        "Supported inputs: RDML, CSV, TSV, TXT and Excel XLSX/XLSM files containing reported Cq/Ct values. "
        "RDML files may also include amplification and melt-curve traces, which VIRASYN can display. "
        "Cq is not yet recalculated from raw fluorescence traces, and proprietary binary instrument projects "
        "still need export to a supported format."
    )
    qpcr_file = st.file_uploader(
        "Upload qPCR plate/result export",
        type=["rdml", "csv", "tsv", "txt", "xlsx", "xlsm"],
        key="qpcr_input_file",
    )

    if qpcr_file is not None:
        try:
            qpcr_sheet = None
            if Path(qpcr_file.name).suffix.lower() in {".xlsx", ".xlsm"}:
                workbook = pd.ExcelFile(BytesIO(qpcr_file.getvalue()), engine="openpyxl")
                qpcr_sheet = st.selectbox(
                    "Worksheet", workbook.sheet_names, key="qpcr_sheet_name"
                )
            qpcr_raw = read_qpcr_table(
                qpcr_file.getvalue(), qpcr_file.name, sheet_name=qpcr_sheet
            )
            st.caption(f"Loaded {len(qpcr_raw):,} rows and {len(qpcr_raw.columns)} columns.")
            with st.expander("Preview original instrument export", expanded=True):
                st.dataframe(qpcr_raw.head(15), use_container_width=True, hide_index=True)

            if "Amplification trace (JSON)" in qpcr_raw.columns or "Melt curve trace (JSON)" in qpcr_raw.columns:
                trace_rows = []
                for trace_index, trace_row in qpcr_raw.iterrows():
                    amp_text = trace_row.get("Amplification trace (JSON)", "")
                    melt_text = trace_row.get("Melt curve trace (JSON)", "")
                    if (isinstance(amp_text, str) and amp_text) or (isinstance(melt_text, str) and melt_text):
                        trace_label = (
                            f"{trace_row.get('Sample', trace_row.get('Sample ID', 'Sample'))} · "
                            f"{trace_row.get('Target', 'Target')} · {trace_row.get('Well', trace_index)}"
                        )
                        trace_rows.append((trace_label, trace_index))
                if trace_rows:
                    with st.expander("RDML amplification and melt-curve trace viewer", expanded=False):
                        trace_labels = [item[0] for item in trace_rows]
                        selected_trace_label = st.selectbox(
                            "Reaction / target",
                            trace_labels,
                            key="qpcr_rdml_trace_select",
                        )
                        selected_trace_index = dict(trace_rows)[selected_trace_label]
                        selected_trace = qpcr_raw.loc[selected_trace_index]
                        t1, t2 = st.columns(2)
                        amp_text = selected_trace.get("Amplification trace (JSON)", "")
                        melt_text = selected_trace.get("Melt curve trace (JSON)", "")
                        with t1:
                            st.caption("Amplification: cycle vs fluorescence")
                            if isinstance(amp_text, str) and amp_text:
                                amp_frame = pd.DataFrame(json.loads(amp_text))
                                if not amp_frame.empty and {"cyc", "fluor"}.issubset(amp_frame.columns):
                                    st.line_chart(amp_frame.set_index("cyc")["fluor"], use_container_width=True)
                                else:
                                    st.info("No amplification-curve points were included in this RDML reaction.")
                            else:
                                st.info("No amplification-curve points were included in this RDML reaction.")
                        with t2:
                            st.caption("Melt curve: temperature vs fluorescence")
                            if isinstance(melt_text, str) and melt_text:
                                melt_frame = pd.DataFrame(json.loads(melt_text))
                                if not melt_frame.empty and {"tmp", "fluor"}.issubset(melt_frame.columns):
                                    st.line_chart(melt_frame.set_index("tmp")["fluor"], use_container_width=True)
                                else:
                                    st.info("No melt-curve points were included in this RDML reaction.")
                            else:
                                st.info("No melt-curve points were included in this RDML reaction.")
                        st.caption("Stored trace visualization only; VIRASYN does not yet call or recalculate Cq from these raw curves.")

            qpcr_columns = list(qpcr_raw.columns)

            def _qpcr_index(candidate):
                return qpcr_columns.index(candidate) if candidate in qpcr_columns else 0

            c1, c2, c3 = st.columns(3)
            with c1:
                qpcr_sample_col = st.selectbox(
                    "Biological sample ID",
                    qpcr_columns,
                    index=_qpcr_index(detect_column(qpcr_columns, (
                        "sample", "sample id", "sample name", "specimen"
                    ))),
                    help="Use the same ID for technical wells from one biological sample. Different biological samples need different IDs.",
                    key="qpcr_sample_col",
                )
            with c2:
                qpcr_target_col = st.selectbox(
                    "Target / gene / assay",
                    qpcr_columns,
                    index=_qpcr_index(detect_column(qpcr_columns, (
                        "target", "gene", "assay", "detector"
                    ))),
                    key="qpcr_target_col",
                )
            with c3:
                qpcr_cq_col = st.selectbox(
                    "Cq / Ct value",
                    qpcr_columns,
                    index=_qpcr_index(detect_column(qpcr_columns, (
                        "cq", "ct", "quantification cycle", "cycle threshold"
                    ))),
                    key="qpcr_cq_col",
                )

            optional_columns = ["— None —"] + qpcr_columns
            c4, c5 = st.columns(2)
            with c4:
                qpcr_condition_choice = st.selectbox(
                    "Experimental condition / group",
                    optional_columns,
                    index=(
                        optional_columns.index(detect_column(qpcr_columns, (
                            "condition", "group", "treatment", "experimental group"
                        )))
                        if detect_column(qpcr_columns, (
                            "condition", "group", "treatment", "experimental group"
                        )) in optional_columns else 0
                    ),
                    key="qpcr_condition_col",
                )
            with c5:
                qpcr_role_choice = st.selectbox(
                    "Sample role (optional, for NTC identification)",
                    optional_columns,
                    index=(
                        optional_columns.index(detect_column(qpcr_columns, (
                            "sample type", "role", "control type", "well type"
                        )))
                        if detect_column(qpcr_columns, (
                            "sample type", "role", "control type", "well type"
                        )) in optional_columns else 0
                    ),
                    key="qpcr_role_col",
                )

            qpcr_condition_col = None if qpcr_condition_choice == "— None —" else qpcr_condition_choice
            qpcr_role_col = None if qpcr_role_choice == "— None —" else qpcr_role_choice
            try:
                qpcr_wells = normalise_qpcr_wells(
                    qpcr_raw,
                    sample_col=qpcr_sample_col,
                    target_col=qpcr_target_col,
                    cq_col=qpcr_cq_col,
                    condition_col=qpcr_condition_col,
                    role_col=qpcr_role_col,
                )
                qpcr_targets = sorted(qpcr_wells["target"].dropna().astype(str).unique().tolist())
                qpcr_conditions = sorted(
                    qpcr_wells.loc[~qpcr_wells["is_ntc"], "condition"].dropna().astype(str).unique().tolist()
                )
                detected_count = int(qpcr_wells.loc[~qpcr_wells["is_ntc"], "Cq"].notna().sum())
                ntc_detected = int(qpcr_wells.loc[qpcr_wells["is_ntc"], "Cq"].notna().sum())

                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Wells / rows", f"{len(qpcr_wells):,}")
                m2.metric("Targets", len(qpcr_targets))
                m3.metric("Detected sample wells", detected_count)
                m4.metric("NTC detections", ntc_detected)

                if not qpcr_condition_col:
                    st.warning(
                        "No condition column selected. VIRASYN will report normalized sample-level "
                        "quantities, but a treatment-versus-control fold change cannot be interpreted."
                    )
                st.caption(
                    "Undetermined, blank and non-numeric Cq entries remain missing values; they are never converted to zero."
                )

                qpcr_refs_default = [
                    gene for gene in qpcr_targets
                    if gene.lower() in {"gapdh", "actb", "b2m", "hprt1", "hprt", "rplp0", "tbp", "18s", "18s rrna"}
                ]
                qpcr_refs = st.multiselect(
                    "Reference gene(s) for normalization",
                    qpcr_targets,
                    default=qpcr_refs_default,
                    help="Choose one or more experimentally justified, stable reference genes. Do not assume a housekeeping gene is stable without evidence.",
                    key="qpcr_reference_genes",
                )
                qpcr_calibrator = st.selectbox(
                    "Control / calibrator condition",
                    qpcr_conditions or ["All samples"],
                    key="qpcr_calibrator",
                )
                qpcr_max_sd = st.slider(
                    "Technical replicate SD review threshold (Cq)",
                    min_value=0.25, max_value=3.0, value=1.0, step=0.25,
                    help="This flags spread for review; it does not automatically remove wells.",
                    key="qpcr_max_sd",
                )

                st.markdown("#### Optional assay standard curves")
                st.caption(
                    "Upload a separate standard dilution-series table with target, positive concentration "
                    "and Cq/Ct columns. At least three distinct concentrations are needed for each assay."
                )
                qpcr_standard_file = st.file_uploader(
                    "Standard-curve export (optional)",
                    type=["rdml", "csv", "tsv", "txt", "xlsx", "xlsm"],
                    key="qpcr_standard_file",
                )
                qpcr_standard_results = pd.DataFrame()
                qpcr_standard_map = {}
                if qpcr_standard_file is not None:
                    try:
                        qpcr_standard_raw = read_qpcr_table(
                            qpcr_standard_file.getvalue(), qpcr_standard_file.name
                        )
                        st.dataframe(qpcr_standard_raw.head(10), use_container_width=True, hide_index=True)
                        sc_columns = list(qpcr_standard_raw.columns)
                        sc_target = st.selectbox(
                            "Standard curve target column", sc_columns,
                            index=(
                                sc_columns.index(detect_column(sc_columns, ("target", "gene", "assay")))
                                if detect_column(sc_columns, ("target", "gene", "assay")) in sc_columns else 0
                            ),
                            key="qpcr_sc_target",
                        )
                        sc_cq = st.selectbox(
                            "Standard curve Cq/Ct column", sc_columns,
                            index=sc_columns.index(detect_column(sc_columns, ("cq", "ct")))
                            if detect_column(sc_columns, ("cq", "ct")) in sc_columns else 0,
                            key="qpcr_sc_cq",
                        )
                        sc_concentration = st.selectbox(
                            "Known standard concentration column", sc_columns,
                            index=sc_columns.index(detect_column(sc_columns, (
                                "concentration", "copies", "quantity"
                            )))
                            if detect_column(sc_columns, (
                                "concentration", "copies", "quantity"
                            )) in sc_columns else 0,
                            key="qpcr_sc_concentration",
                        )
                        qpcr_standard_results = calculate_standard_curves(
                            qpcr_standard_raw, sc_target, sc_cq, sc_concentration
                        )
                        st.dataframe(qpcr_standard_results, use_container_width=True, hide_index=True)
                        qpcr_standard_map = {
                            str(row["target"]): float(row["efficiency_percent"])
                            for _, row in qpcr_standard_results.iterrows()
                            if pd.notna(row["efficiency_percent"])
                        }
                    except Exception as standard_error:
                        st.error(f"Could not process standard-curve file: {standard_error}")

                efficiency_rows = [
                    {
                        "target": target,
                        "efficiency_percent": qpcr_standard_map.get(target, 100.0),
                        "source": "Standard curve" if target in qpcr_standard_map else "Assumed 100% — edit if known",
                    }
                    for target in qpcr_targets
                ]
                efficiency_frame = pd.DataFrame(efficiency_rows)
                st.markdown("#### Amplification efficiency by assay")
                efficiency_edit = st.data_editor(
                    efficiency_frame,
                    hide_index=True,
                    use_container_width=True,
                    disabled=["target", "source"],
                    column_config={
                        "efficiency_percent": st.column_config.NumberColumn(
                            "Efficiency (%)", min_value=0.0, max_value=300.0, step=0.5,
                            help="Use the value measured from a validated standard curve where available."
                        )
                    },
                    key="qpcr_efficiency_editor",
                )
                if not qpcr_standard_map:
                    st.warning(
                        "No usable standard-curve efficiencies were imported. The table starts at an "
                        "explicit 100% efficiency assumption for every target. Edit these values using "
                        "your assay validation data; results are exploratory until efficiencies are established."
                    )
                else:
                    st.caption(
                        "Standard-curve efficiencies were prefilled. Review slope, R² and efficiency status "
                        "before accepting the values."
                    )

                st.markdown("#### Run analysis")
                if st.button(
                    "Analyse qPCR data",
                    type="primary",
                    use_container_width=True,
                    disabled=not qpcr_refs,
                    key="qpcr_analyse_button",
                ):
                    try:
                        efficiency_map = {
                            str(row["target"]): float(row["efficiency_percent"])
                            for _, row in efficiency_edit.iterrows()
                        }
                        qpcr_result = analyse_relative_expression(
                            qpcr_wells,
                            reference_genes=qpcr_refs,
                            calibrator_condition=qpcr_calibrator,
                            efficiency_percent=efficiency_map,
                            max_technical_sd=qpcr_max_sd,
                        )
                        sample_level = qpcr_result["sample_level"]
                        group_summary = qpcr_result["group_summary"]
                        technical_qc = qpcr_result["technical_qc"]
                        ntc_qc = qpcr_result["ntc_qc"]

                        st.success("qPCR calculations completed. Review the quality flags before interpreting the biology.")
                        out_of_range_efficiencies = {
                            target: value for target, value in efficiency_map.items()
                            if value < 90 or value > 110
                        }
                        if out_of_range_efficiencies:
                            st.warning(
                                "Review assay efficiencies outside the common 90–110% range: "
                                + ", ".join(f"{target}={value:g}%" for target, value in out_of_range_efficiencies.items())
                            )
                        if not group_summary.empty and group_summary["n_biological_samples"].min() < 2:
                            st.warning(
                                "Some groups have fewer than two biological samples. Inferential p-values are unavailable or underpowered."
                            )
                        if not technical_qc.empty and technical_qc["technical_qc"].str.startswith("Review").any():
                            st.warning("Some technical replicate groups exceed the selected SD threshold. They are flagged, not automatically excluded.")
                        partial_ref_rows = sample_level["normalization_status"].str.startswith("Partial reference set").sum()
                        if partial_ref_rows:
                            st.warning(
                                f"{partial_ref_rows} sample-target rows used fewer than all selected reference genes. "
                                "Review the reference-gene detection table before interpreting these samples."
                            )
                        if ntc_detected:
                            st.error("One or more NTC wells have a Cq value. Investigate possible contamination before interpreting results.")
                        elif int(ntc_qc.iloc[0]["NTC_wells"]) == 0:
                            st.info("No NTC rows were detected in the supplied table. Confirm that appropriate negative controls were run and reported.")

                        r1, r2, r3 = st.columns(3)
                        r1.metric("Biological samples", sample_level["sample_id"].nunique())
                        r2.metric("Normalized sample-target rows", int(sample_level["normalized_quantity"].notna().sum()))
                        r3.metric("Targets analysed", group_summary["target"].nunique() if not group_summary.empty else 0)

                        st.markdown("##### Group-level relative expression")
                        st.dataframe(group_summary, use_container_width=True, hide_index=True)
                        fold_plot = group_summary[
                            group_summary["condition"].ne(qpcr_calibrator)
                        ].pivot_table(
                            index="target",
                            columns="condition",
                            values="fold_change_vs_calibrator",
                            aggfunc="first",
                        )
                        if not fold_plot.empty:
                            st.markdown("##### Fold change versus calibrator")
                            st.bar_chart(fold_plot, use_container_width=True)

                        with st.expander("Per-biological-sample normalized quantities", expanded=True):
                            st.dataframe(sample_level, use_container_width=True, hide_index=True)
                        with st.expander("Technical replicate quality-control table"):
                            st.dataframe(technical_qc, use_container_width=True, hide_index=True)
                        with st.expander("NTC quality-control summary"):
                            st.dataframe(ntc_qc, use_container_width=True, hide_index=True)
                        if not qpcr_standard_results.empty:
                            with st.expander("Standard-curve fit details"):
                                st.dataframe(qpcr_standard_results, use_container_width=True, hide_index=True)

                        report_buffer = BytesIO()
                        with pd.ExcelWriter(report_buffer, engine="openpyxl") as writer:
                            qpcr_raw.to_excel(writer, sheet_name="OriginalInput", index=False)
                            qpcr_wells.to_excel(writer, sheet_name="WellData", index=False)
                            sample_level.to_excel(writer, sheet_name="SampleNormalised", index=False)
                            group_summary.to_excel(writer, sheet_name="GroupSummary", index=False)
                            technical_qc.to_excel(writer, sheet_name="TechnicalQC", index=False)
                            ntc_qc.to_excel(writer, sheet_name="NTC_QC", index=False)
                            efficiency_edit.to_excel(writer, sheet_name="Efficiencies", index=False)
                            if not qpcr_standard_results.empty:
                                qpcr_standard_results.to_excel(writer, sheet_name="StandardCurves", index=False)
                        st.download_button(
                            "Download complete qPCR analysis workbook",
                            report_buffer.getvalue(),
                            "virasyn_qpcr_analysis.xlsx",
                            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            use_container_width=True,
                        )
                        st.download_button(
                            "Download normalized sample-level data (CSV)",
                            sample_level.to_csv(index=False),
                            "virasyn_qpcr_sample_normalised.csv",
                            "text/csv",
                        )
                        st.download_button(
                            "Download group summary (CSV)",
                            group_summary.to_csv(index=False),
                            "virasyn_qpcr_group_summary.csv",
                            "text/csv",
                        )
                    except Exception as analysis_error:
                        st.error(f"qPCR analysis could not be completed: {analysis_error}")
            except Exception as import_error:
                st.error(f"Could not read or map this qPCR export: {import_error}")
        except Exception as file_error:
            st.error(f"Could not load qPCR file: {file_error}")
    else:
        st.markdown(
            """
            **Expected data layout (example):**

            | Sample | Condition | Target | Cq |
            |---|---|---|---:|
            | Control_1 | Control | ACTB | 19.8 |
            | Control_1 | Control | GeneA | 25.4 |
            | Treatment_1 | Treatment | ACTB | 20.1 |
            | Treatment_1 | Treatment | GeneA | 22.7 |

            Upload your instrument's CSV, TSV or Excel Cq/Ct export to begin.
            """
        )


with workbench_tab:
    st.subheader("VIRASYN Bioinformatics Workbench")
    st.write(
        "A consolidated catalogue of research workflows spanning sequence, genome, transcriptome, "
        "protein, structure, docking, comparative genomics and virology. The status column distinguishes "
        "functions available in the current release from partial integrations and planned development."
    )

    module_catalogue = pd.DataFrame([
        {
            "Module": "Genome analysis",
            "Capabilities": "FASTA import, sequence summaries, GC content, six-frame ORF screening and downloadable ORFs. FASTQ quality control, FastQC integration and automated genome annotation are planned.",
            "Current status": "Partial",
        },
        {
            "Module": "Sequence analysis",
            "Capabilities": "BLASTP against user-provided protein references and federated sequence-database discovery. Multiple sequence alignment and phylogenetic analysis are planned.",
            "Current status": "Partial",
        },
        {
            "Module": "Protein analysis",
            "Capabilities": "Protein physicochemical properties and optional HMMER profile scanning with user-provided HMM libraries. Automated motif/domain annotation and richer functional annotation are planned.",
            "Current status": "Partial",
        },
        {
            "Module": "Structural bioinformatics",
            "Capabilities": "PDB/mmCIF metadata parsing and PDB 3D visualization. Structural alignment, TM-align integration and RMSD comparison are planned.",
            "Current status": "Partial",
        },
        {
            "Module": "Molecular docking",
            "Capabilities": "Protein and ligand preparation, docking workflows, pose scoring and protein–ligand interaction analysis.",
            "Current status": "Planned",
        },
        {
            "Module": "Transcriptomics",
            "Capabilities": "RNA-seq quality control, differential expression, PCA, volcano plots and heatmaps.",
            "Current status": "Planned",
        },
        {
            "Module": "Comparative genomics",
            "Capabilities": "Genome-to-genome comparison, gene conservation, orthology and synteny analysis, with automated annotation links.",
            "Current status": "Planned",
        },
        {
            "Module": "Virology",
            "Capabilities": "Viral literature, genome, protein and taxonomy discovery plus preliminary sequence characterization. Validated candidate attachment-protein workflows and virus–host evidence integration remain planned.",
            "Current status": "Partial",
        },
        {
            "Module": "Data visualization",
            "Capabilities": "Analysis tables, qPCR fold-change charts and interactive PDB viewing. Publication-ready figure templates, PCA plots, volcano plots and heatmaps are planned.",
            "Current status": "Partial",
        },
        {
            "Module": "PCR / qPCR analysis",
            "Capabilities": "RDML and tabular Cq/Ct import, technical-replicate QC, standard curves, efficiency-aware relative expression, NTC flags and Excel/CSV reports.",
            "Current status": "Available",
        },
        {
            "Module": "PhARIS-inspired RBP intelligence",
            "Capabilities": "Six-frame ORF screening, local BLASTP against an authorized RBP reference FASTA, identity/query-coverage thresholds, optional cluster metadata and downloadable candidate reports.",
            "Current status": "Available · reference library supplied by user",
        },
        {
            "Module": "Bacteria-focused Phage Intelligence",
            "Capabilities": "Integrated RBP homology discovery, curated bacterial host evidence ranking, exact-ID receptor annotation mapping, structure-file inventory and reproducible exports. Calibrated universal binding prediction and structural docking are not yet implemented.",
            "Current status": "Partial · evidence integration available",
        },
    ])
    st.dataframe(
        module_catalogue,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Module": st.column_config.TextColumn("Module", width="medium"),
            "Capabilities": st.column_config.TextColumn("Capabilities", width="large"),
            "Current status": st.column_config.TextColumn("Current status", width="small"),
        },
    )
    st.caption(
        "Available = implemented in the current release. Partial = some component is available, "
        "but the full workflow is not yet implemented. Planned = not currently available in VIRASYN."
    )

    st.divider()
    st.subheader("Next five development priorities")
    st.markdown("""
1. **RNA-seq analysis** — FastQC, differential expression, PCA, volcano plots and heatmaps.
2. **Multiple sequence alignment and phylogenetics** — aligned sequence workflows, model selection and evolutionary trees.
3. **Protein structure comparison** — structural superposition, TM-align and RMSD reporting.
4. **Molecular docking** — protein/ligand preparation, docking workflows and protein–ligand interaction analysis.
5. **Comparative genomics and automated genome annotation** — genome comparison, gene conservation, orthology and synteny.
""")
    st.info(
        "These five items are the planned development sequence, not features that are already active. "
        "New workflows will be added with reproducible parameters, inspectable intermediate results, "
        "clear software provenance and downloadable outputs."
    )



with alphafold_tab:
    st.subheader("AlphaFold protein structure lookup")
    st.write(
        "Retrieve published AlphaFold Database predictions by UniProt accession, inspect "
        "model confidence metadata, visualize the predicted structure and download the coordinates."
    )
    st.info(
        "This module retrieves existing AlphaFold Database models. It does not run a new "
        "AlphaFold prediction. For a new sequence prediction, use the official AlphaFold Server."
    )
    from virasyn.alphafold import lookup_alphafold_prediction, download_alphafold_structure

    with st.form("alphafold_lookup_form"):
        alphafold_accession = st.text_input(
            "UniProt accession",
            placeholder="e.g. P0DTC2",
            help="Enter a UniProt accession for a protein with a model in the AlphaFold Protein Structure Database.",
        )
        alphafold_format = st.selectbox("Coordinate format", ["PDB", "mmCIF"])
        alphafold_submit = st.form_submit_button("Find AlphaFold structure", type="primary")

    if alphafold_submit:
        if not alphafold_accession.strip():
            st.warning("Enter a UniProt accession.")
        else:
            try:
                with st.spinner("Querying the AlphaFold Protein Structure Database..."):
                    predictions = lookup_alphafold_prediction(alphafold_accession.strip())
                if not predictions:
                    st.warning("No AlphaFold Database model was found for this accession. Check the accession or try the UniProt entry.")
                    st.session_state.pop("virasyn_alphafold_model", None)
                else:
                    st.session_state["virasyn_alphafold_model"] = predictions[0]
                    st.success(f"Found {len(predictions)} AlphaFold model record(s). Showing the first model.")
            except (ValueError, RuntimeError) as exc:
                st.error(f"AlphaFold lookup failed: {exc}")

    model = st.session_state.get("virasyn_alphafold_model")
    if model:
        st.markdown(f"**UniProt accession:** [{model.get('uniprotAccession', alphafold_accession if 'alphafold_accession' in locals() else '')}](https://www.uniprot.org/uniprotkb/{model.get('uniprotAccession', '')}/entry)")
        model_id = model.get("modelEntityId", model.get("entryId", "Open model"))
        st.markdown(f"**AlphaFold model:** [{model_id}](https://alphafold.com/entry/{model_id})")
        m1.metric("Mean pLDDT", f"{float(model.get('globalMetricValue', model.get('confidenceScore'))):.1f}" if model.get('globalMetricValue', model.get('confidenceScore')) is not None else "Not supplied")
        m2.metric("Model version", str(model.get("latestVersion", "Not supplied")))
        m3.metric("UniProt", model.get("uniprotAccession", "—"))
        st.caption("pLDDT is a local confidence estimate (0–100), not experimental validation or a measure of biological activity.")
        try:
            with st.spinner(f"Downloading {alphafold_format} coordinates..."):
                structure_text, filename = download_alphafold_structure(model, alphafold_format)
            st.download_button(
                f"Download {alphafold_format} structure",
                structure_text,
                filename,
                "chemical/x-pdb" if alphafold_format == "PDB" else "chemical/x-mmcif",
            )
            if alphafold_format == "PDB":
                try:
                    import streamlit.components.v1 as components
                    components.html(render_structure_html(structure_text), height=560, scrolling=False)
                except Exception as viewer_error:
                    st.warning(f"Structure viewer could not render this model: {viewer_error}")
            else:
                st.caption("The mmCIF coordinates are available to download. Switch to PDB format to use the built-in 3D viewer.")
        except (ValueError, RuntimeError) as exc:
            st.error(f"Could not retrieve model coordinates: {exc}")

    st.markdown("[Open AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/) · [Open AlphaFold Server for new sequence predictions](https://alphafoldserver.com/)")




with pharis_tab:
    st.subheader("PhARIS-inspired RBP intelligence")
    st.write(
        "Screen phage genomes for candidate receptor-binding proteins (RBPs), then compare "
        "translated ORFs against an authorized, versioned RBP reference FASTA using local BLASTP."
    )
    st.info(
        "This is an independent, reference-driven workflow inspired by the published PhARIS approach. "
        "It does not include the original PhARIS reference sequences. Upload a reference library you are "
        "authorized to use. Sequence similarity is not proof of receptor binding, adsorption or productive infection."
    )

    from virasyn.pharis import analyse_pharis_genome, reference_fasta_manifest

    ph_genome = st.file_uploader(
        "Upload phage genome FASTA", type=["fa", "fasta", "fna", "fas"], key="pharis_genome"
    )
    ph_reference = st.file_uploader(
        "Upload authorized RBP reference protein FASTA", type=["fa", "fasta", "faa", "fas"], key="pharis_reference"
    )
    ph_metadata = st.file_uploader(
        "Optional reference annotation CSV", type=["csv"], key="pharis_metadata",
        help="Optional columns: subject_id (must match FASTA record IDs), cluster, RBP_name, strain, receptor, evidence, citation.",
    )
    st.download_button(
        "Download reference metadata template",
        "subject_id,cluster,RBP_name,strain,receptor,evidence,citation,source_version\\n",
        "virasyn_rbp_reference_metadata_template.csv",
        "text/csv",
    )
    if ph_reference is not None:
        try:
            manifest = reference_fasta_manifest(ph_reference.getvalue().decode("utf-8"))
            st.caption(f"Reference library loaded: {len(manifest)} protein sequences.")
            with st.expander("Reference library inventory"):
                st.dataframe(manifest, use_container_width=True, hide_index=True)
        except (UnicodeDecodeError, ValueError) as exc:
            st.error(f"Reference FASTA could not be read: {exc}")

    p1, p2, p3 = st.columns(3)
    with p1:
        ph_min_aa = st.number_input("Minimum ORF length (aa)", min_value=30, max_value=2000, value=100, step=10)
    with p2:
        ph_min_identity = st.slider("Minimum sequence identity (%)", min_value=30, max_value=100, value=80)
    with p3:
        ph_min_coverage = st.slider("Minimum query coverage (%)", min_value=30, max_value=100, value=80)
    ph_code = st.selectbox("NCBI genetic code", [11, 4, 1], index=0, help="11 = bacterial, archaeal and plant plastid; select the code appropriate to the phage.")

    if st.button("Run RBP discovery workflow", type="primary", use_container_width=True, disabled=(ph_genome is None or ph_reference is None)):
        try:
            with st.spinner("Translating six reading frames and running BLASTP against the reference library..."):
                ph_result = analyse_pharis_genome(
                    ph_genome.getvalue().decode("utf-8"),
                    ph_reference.getvalue().decode("utf-8"),
                    ph_metadata.getvalue().decode("utf-8") if ph_metadata is not None else None,
                    min_orf_aa=int(ph_min_aa),
                    min_identity=float(ph_min_identity),
                    min_query_coverage=float(ph_min_coverage),
                    genetic_code=int(ph_code),
                )
            st.session_state["virasyn_pharis_result"] = ph_result
            st.success("RBP homology screening completed.")
        except (UnicodeDecodeError, ValueError, RuntimeError) as exc:
            st.error(f"RBP analysis could not be completed: {exc}")

    ph_result = st.session_state.get("virasyn_pharis_result")
    if ph_result:
        summary = ph_result["summary"]
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Genome records", summary["genome_records"])
        c2.metric("ORFs screened", summary["orfs_screened"])
        c3.metric("BLAST hits", summary["blast_hits"])
        c4.metric("Candidate ORFs", summary["candidate_orfs"])
        st.caption(
            f"Parameters: minimum {summary['min_orf_aa']} aa · identity ≥ {summary['min_identity']}% · "
            f"query coverage ≥ {summary['min_query_coverage']}% · genetic code {summary['genetic_code']}"
        )
        st.markdown("**Top reference-supported candidate matches**")
        candidates = ph_result["candidate_hits"]
        if candidates.empty:
            st.info("No matches passed the selected identity and query-coverage thresholds. This does not exclude an RBP; the reference library may be incomplete or divergent.")
        else:
            display_cols = [c for c in [
                "qseqid", "sequence_id", "strand", "frame", "start_1based", "end_1based",
                "sseqid", "pident", "query_coverage_pct", "subject_coverage_pct",
                "evalue", "bitscore", "cluster", "RBP_name", "strain", "receptor", "evidence", "citation"
            ] if c in candidates.columns]
            st.dataframe(candidates[display_cols], use_container_width=True, hide_index=True)
        st.download_button(
            "Download all BLASTP hits (CSV)",
            ph_result["all_hits"].to_csv(index=False),
            "virasyn_pharis_all_hits.csv",
            "text/csv",
        )
        st.download_button(
            "Download screened ORFs (CSV)",
            ph_result["orfs"].to_csv(index=False),
            "virasyn_pharis_orfs.csv",
            "text/csv",
        )
        st.warning(
            "Interpretation: a BLASTP match indicates homology to a sequence in the uploaded reference set. "
            "It does not independently establish binding specificity or host range. Confirm candidates with "
            "the reference study, receptor assays and appropriate experimental validation."
        )
    st.markdown(
        "[Original PhARIS repository](https://github.com/JKrusche1/PhARIS) · "
        "[PhARIS study](https://doi.org/10.1016/j.celrep.2025.115369)"
    )



with phage_tab:
    st.subheader("VIRASYN Phage Intelligence · Bacteria")
    st.write(
        "An integrated research workflow connecting candidate RBP discovery, curated receptor annotations, "
        "bacterial host-interaction evidence and structure intake. Use a broad, versioned RBP reference library "
        "covering the bacterial taxa relevant to your study."
    )
    st.warning(
        "This release produces reference-supported RBP candidates and evidence-ranked host/receptor hypotheses. "
        "It does not contain a universal, trained RBP–receptor binding model. Evidence indices are not probabilities; "
        "sequence homology does not prove receptor binding, adsorption or productive infection."
    )
    from virasyn.pharis import analyse_pharis_genome, reference_fasta_manifest
    from virasyn.phage_intelligence import build_integrated_report, csv_template

    with st.expander("Workflow and evidence model", expanded=False):
        st.markdown("""
**Analysis stages**
1. Translate uploaded bacteriophage genome records in six frames and screen ORFs against a user-authorized RBP reference library using BLASTP.
2. Attach reference metadata such as RBP family, cluster, bacterial host and receptor.
3. Integrate independently curated phage–bacterium evidence (experimental infection/host range, adsorption, receptor binding, CRISPR and computational evidence).
4. Map RBP reference identifiers to receptor annotations and rank the evidence transparently.
5. Inspect optional RBP/receptor coordinate files and export analysis tables for downstream structural work.

**Evidence vocabulary**
- evidence_level: experimental, curated, computational or inferred.
- evidence_type: productive_infection, host_range, adsorption, receptor_binding, crispr_spacer, sequence_similarity or another explicitly labelled type.
- The host evidence index is a transparent weighted evidence summary, not a calibrated probability or an ML prediction.
""")

    phage_genome = st.file_uploader(
        "Bacteriophage genome FASTA (one or multiple records)",
        type=["fa", "fasta", "fna", "fas"], key="phageintel_genome"
    )
    rbp_reference = st.file_uploader(
        "Broad-spectrum RBP reference protein FASTA",
        type=["fa", "fasta", "faa", "fas"], key="phageintel_rbp_reference",
        help="Use a curated, permission-cleared library with coverage across your target bacterial groups."
    )
    rbp_metadata = st.file_uploader(
        "Optional RBP reference metadata CSV", type=["csv"], key="phageintel_rbp_metadata",
        help="Must include subject_id matching reference FASTA IDs. Suggested: cluster, RBP_name, strain, receptor, evidence, citation, source_version."
    )
    host_evidence_file = st.file_uploader(
        "Optional phage–bacterial host evidence CSV", type=["csv"], key="phageintel_host_evidence"
    )
    receptor_file = st.file_uploader(
        "Optional RBP–receptor annotation CSV", type=["csv"], key="phageintel_receptor_map"
    )

    t1, t2 = st.columns(2)
    with t1:
        st.download_button("Download host-evidence CSV template", csv_template("host"),
                           "virasyn_host_evidence_template.csv", "text/csv")
    with t2:
        st.download_button("Download RBP–receptor CSV template", csv_template("receptor"),
                           "virasyn_rbp_receptor_template.csv", "text/csv")

    if rbp_reference is not None:
        try:
            ref_manifest = reference_fasta_manifest(rbp_reference.getvalue().decode("utf-8"))
            st.caption(f"RBP reference inventory: {len(ref_manifest)} protein sequences.")
            with st.expander("Inspect reference inventory"):
                st.dataframe(ref_manifest, use_container_width=True, hide_index=True)
        except (UnicodeDecodeError, ValueError) as exc:
            st.error(f"RBP reference library is invalid: {exc}")

    a1, a2, a3 = st.columns(3)
    with a1:
        intel_min_aa = st.number_input("Minimum translated ORF (aa)", min_value=30, max_value=2000, value=80, step=10, key="intel_min_aa")
    with a2:
        intel_identity = st.slider("Minimum RBP identity (%)", 30, 100, 35, key="intel_identity")
    with a3:
        intel_coverage = st.slider("Minimum query coverage (%)", 30, 100, 40, key="intel_coverage")
    intel_code = st.selectbox("Translation table", [11, 4, 1], index=0, key="intel_genetic_code")

    structure_uploads = st.file_uploader(
        "Optional RBP/receptor structure files for inventory (PDB or mmCIF)",
        type=["pdb", "ent", "cif", "mmcif"], accept_multiple_files=True, key="phageintel_structures"
    )

    if st.button("Run integrated phage intelligence", type="primary", use_container_width=True,
                 disabled=(phage_genome is None or rbp_reference is None), key="run_phage_intelligence"):
        try:
            with st.spinner("Screening ORFs against the RBP reference library and integrating evidence..."):
                discovery = analyse_pharis_genome(
                    phage_genome.getvalue().decode("utf-8"),
                    rbp_reference.getvalue().decode("utf-8"),
                    rbp_metadata.getvalue().decode("utf-8") if rbp_metadata else None,
                    min_orf_aa=int(intel_min_aa), min_identity=float(intel_identity),
                    min_query_coverage=float(intel_coverage), genetic_code=int(intel_code),
                )
                integrated = build_integrated_report(
                    discovery["candidate_hits"],
                    host_evidence_file.getvalue().decode("utf-8") if host_evidence_file else None,
                    receptor_file.getvalue().decode("utf-8") if receptor_file else None,
                )
                structure_inventory = []
                for structure_file in (structure_uploads or []):
                    try:
                        structure_inventory.append(parse_structure(
                            structure_file.getvalue().decode("utf-8", errors="replace"), structure_file.name
                        ))
                    except Exception as structure_error:
                        structure_inventory.append({"file": structure_file.name, "error": str(structure_error)})
                st.session_state["virasyn_phage_intelligence"] = {
                    "discovery": discovery, "integrated": integrated,
                    "structure_inventory": structure_inventory,
                }
            st.success("Integrated analysis completed.")
        except (UnicodeDecodeError, ValueError, RuntimeError) as exc:
            st.error(f"Phage intelligence workflow failed: {exc}")

    result = st.session_state.get("virasyn_phage_intelligence")
    if result:
        discovery, integrated = result["discovery"], result["integrated"]
        summary = discovery["summary"]
        metrics = st.columns(4)
        metrics[0].metric("Genome records", summary["genome_records"])
        metrics[1].metric("ORFs screened", summary["orfs_screened"])
        metrics[2].metric("RBP reference hits", summary["blast_hits"])
        metrics[3].metric("Candidate ORFs", summary["candidate_orfs"])
        st.caption(
            f"Discovery thresholds: ORF ≥ {summary['min_orf_aa']} aa · identity ≥ {summary['min_identity']}% · "
            f"query coverage ≥ {summary['min_query_coverage']}% · translation table {summary['genetic_code']}"
        )

        st.markdown("### 1 · RBP discovery and classification")
        if discovery["candidate_hits"].empty:
            st.info("No RBP homologues passed the selected thresholds. This does not exclude divergent RBPs or candidates absent from the reference library.")
        else:
            st.dataframe(discovery["candidate_hits"], use_container_width=True, hide_index=True)
        st.download_button("Export RBP candidate table", discovery["candidate_hits"].to_csv(index=False),
                           "virasyn_rbp_candidates.csv", "text/csv")
        st.download_button("Export all sequence-similarity hits", discovery["all_hits"].to_csv(index=False),
                           "virasyn_rbp_all_hits.csv", "text/csv")

        st.markdown("### 2 · Bacterial host evidence ranking")
        host_table = integrated["host_evidence"]
        if host_table.empty:
            st.info("No host evidence table was supplied. Upload curated phage–host evidence to rank known bacterial associations.")
        else:
            st.dataframe(host_table, use_container_width=True, hide_index=True)
            st.download_button("Export host evidence ranking", host_table.to_csv(index=False),
                               "virasyn_bacterial_host_evidence.csv", "text/csv")
            st.caption(
                f"{integrated['host_summary']['unique_hosts']} bacterial taxa · "
                f"{integrated['host_summary']['experimental_records']} experimental evidence records. "
                "Evidence index is not a probability."
            )

        st.markdown("### 3 · RBP–receptor specificity evidence")
        receptor_table = integrated["rbp_receptor_evidence"]
        if receptor_table.empty:
            st.info("No receptor annotations were mapped. Supply an annotation CSV whose subject_id values exactly match RBP reference FASTA identifiers.")
        else:
            st.dataframe(receptor_table, use_container_width=True, hide_index=True)
            st.download_button("Export RBP–receptor evidence", receptor_table.to_csv(index=False),
                               "virasyn_rbp_receptor_evidence.csv", "text/csv")

        st.markdown("### 4 · Structural data inventory")
        structures = result["structure_inventory"]
        if structures:
            st.dataframe(pd.DataFrame(structures), use_container_width=True, hide_index=True)
        else:
            st.caption("No optional structure files were supplied. Use the AlphaFold tab for available UniProt models; structural similarity and docking are planned extensions.")

        st.markdown("### 5 · Interpretation and reproducibility")
        st.warning(integrated["interpretation"])
        report = {
            "platform": "VIRASYN", "module": "Bacteria-focused phage intelligence",
            "discovery_summary": summary, "host_evidence_summary": integrated["host_summary"],
            "candidate_rbps": discovery["candidate_hits"].to_dict(orient="records"),
            "host_evidence": host_table.to_dict(orient="records"),
            "rbp_receptor_evidence": receptor_table.to_dict(orient="records"),
            "structure_inventory": structures,
            "limitations": [integrated["interpretation"], "No calibrated universal RBP–receptor binding model is included in this release."],
        }
        st.download_button("Download integrated JSON research report", json.dumps(report, indent=2, default=str),
                           "virasyn_phage_intelligence_report.json", "application/json")


with resources_tab:
    st.subheader("A wider directory of public virology and molecular-biology resources")
    st.write(
        "VIRASYN directly searches the providers marked as integrated. Other specialist resources are listed "
        "as external portals and can be opened in a new tab; they are not yet queried by VIRASYN."
    )
    resource_groups = {
        "Literature and scholarly discovery": [
            ("PubMed", "https://pubmed.ncbi.nlm.nih.gov/", "Integrated · biomedical literature"),
            ("Europe PMC", "https://europepmc.org/", "Integrated · life-sciences literature and full-text links"),
            ("OpenAlex", "https://openalex.org/", "Integrated · scholarly works and citation metadata"),
            ("Crossref", "https://search.crossref.org/", "Integrated · DOI and publication metadata"),
        ],
        "Viral genomes, sequence archives and taxonomy": [
            ("NCBI Virus", "https://www.ncbi.nlm.nih.gov/labs/virus/", "External portal · curated virus sequence discovery"),
            ("NCBI Nucleotide", "https://www.ncbi.nlm.nih.gov/nuccore", "Integrated · nucleotide records"),
            ("NCBI Taxonomy", "https://www.ncbi.nlm.nih.gov/taxonomy", "Integrated · organism taxonomy"),
            ("European Nucleotide Archive", "https://www.ebi.ac.uk/ena/browser/home", "Integrated · sequence archive"),
            ("ICTV Virus Taxonomy", "https://ictv.global/taxonomy", "External portal · official virus taxonomy"),
            ("Nextstrain", "https://nextstrain.org/", "External portal · pathogen genomic epidemiology"),
        ],
        "Proteins, domains and functional annotation": [
            ("UniProtKB", "https://www.uniprot.org/", "Integrated · protein knowledgebase"),
            ("InterPro", "https://www.ebi.ac.uk/interpro/", "Integrated · protein families and domains"),
            ("Pfam", "https://www.ebi.ac.uk/interpro/entry/pfam/", "External portal · protein families"),
            ("AlphaFold Protein Structure Database", "https://alphafold.ebi.ac.uk/", "Integrated · retrieve published predicted protein structures by UniProt accession"),
            ("NCBI Conserved Domain Database", "https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml", "External portal · conserved domains"),
        ],
        "Macromolecular structures": [
            ("RCSB Protein Data Bank", "https://www.rcsb.org/", "Integrated · experimental and computed structures"),
            ("PDBe", "https://www.ebi.ac.uk/pdbe/", "External portal · Europe-based structure archive access"),
        ],
        "Phage and virus-specialist resources": [
            ("PHROGs", "https://phrogs.lmge.uca.fr/", "External portal · phage protein orthologous groups"),
            ("PhagesDB", "https://phagesdb.org/", "External portal · bacteriophage genomes and metadata"),
            ("INPHARED", "https://github.com/RyanCook94/inphared", "External portal · curated phage genome dataset"),
            ("PhARIS", "https://github.com/JKrusche1/PhARIS", "External research tool · S. aureus phage RBP identification; VIRASYN provides an independent compatible workflow"),
            ("ViralZone", "https://viralzone.expasy.org/", "External portal · virus biology knowledge resource"),
            ("Virus-Host DB", "https://www.genome.jp/virushostdb/", "External portal · virus-host associations"),
            ("BV-BRC", "https://www.bv-brc.org/", "External portal · bacterial and viral research resources"),
            ("VEuPathDB", "https://veupathdb.org/", "External portal · eukaryotic pathogen bioinformatics"),
        ],
    }
    for group_name, entries in resource_groups.items():
        st.markdown(f"### {group_name}")
        cols = st.columns(2)
        for index, (name, url, status) in enumerate(entries):
            with cols[index % 2]:
                with st.container(border=True):
                    st.markdown(f"**[{name}]({url})**")
                    st.caption(status)
    st.info(
        "Database coverage and licensing differ. VIRASYN does not copy all external database contents into its own store. "
        "It searches selected public APIs and links users back to source records. Additional adapters can be added modularly."
    )

with analysis_tab:
    st.subheader("Sequence & structure analysis")
    st.caption("Upload FASTA sequences or PDB/mmCIF structures to begin. Optional reference files enable BLASTP and HMMER searches.")
    st.warning(
        "Analysis tools are exploratory. A sequence hit, predicted ORF or structural match "
        "does not establish receptor binding, host range or productive infection."
    )
    with st.sidebar:
        st.header("Analysis setup")
        virus_group = st.selectbox(
            "Virus group",
            ["Bacteriophage", "Animal virus", "Plant virus", "Archaeal virus", "Other / exploratory"],
            key="analysis_virus_group",
        )
        genetic_code = st.selectbox(
            "Genetic code for ORF screening",
            [11, 1, 2, 4],
            index=0,
            format_func=lambda x: {
                11: "11 · Bacterial / archaeal / plastid",
                1: "1 · Standard",
                2: "2 · Vertebrate mitochondrial",
                4: "4 · Mold / protozoan mitochondrial",
            }[x],
            key="analysis_genetic_code",
        )
        min_aa = st.slider("Minimum ORF length (amino acids)", 20, 300, 45, 5, key="analysis_min_aa")
        alternative_starts = st.checkbox("Include genetic-code alternative start codons", True)
        max_upload_mb = st.select_slider("Maximum upload size per file (MB)", options=[10, 25, 50, 100], value=50)

    st.subheader("Upload genomes, proteins and structures")
    uploaded = st.file_uploader(
        "Genome/protein FASTA and PDB/mmCIF structures",
        type=["fa", "fasta", "fna", "faa", "ffn", "pep", "pdb", "ent", "cif", "mmcif"],
        accept_multiple_files=True,
        key="analysis_primary_uploads",
    )
    reference_fasta = st.file_uploader(
        "Optional BLASTP reference proteins (FASTA)",
        type=["fa", "fasta", "faa"], key="analysis_blast_reference",
    )
    hmm_library = st.file_uploader(
        "Optional HMMER profile library (.hmm)", type=["hmm"], key="analysis_hmm_library",
    )
    run = st.button("Run VIRASYN analysis", type="primary", disabled=not uploaded, use_container_width=True)

    if run:
        too_large = [f.name for f in uploaded if f.size > max_upload_mb * 1024 * 1024]
        if too_large:
            st.error("Files exceed the selected upload limit: " + ", ".join(too_large))
            st.stop()

        sequence_rows, nucleotide_records, protein_records, structure_records = [], [], [], []
        pdb_for_view, errors = [], []
        for file in uploaded:
            suffix = Path(file.name).suffix.lower()
            text = file.getvalue().decode("utf-8", errors="replace")
            if suffix in {".pdb", ".ent", ".cif", ".mmcif"}:
                try:
                    structure_records.append(parse_structure(text, file.name))
                    if suffix in {".pdb", ".ent"}:
                        pdb_for_view.append((file.name, text))
                except Exception as exc:
                    errors.append(f"{file.name}: {exc}")
                continue
            try:
                for record in parse_fasta(text, file.name):
                    sequence_rows.append({
                        "file": file.name, "sequence_id": record["id"],
                        "sequence_type": record["type"], "length": len(record["sequence"]),
                    })
                    (nucleotide_records if record["type"] == "nucleotide" else protein_records).append(record)
            except Exception as exc:
                errors.append(f"{file.name}: {exc}")

        orf_records = []
        for genome in nucleotide_records:
            try:
                orf_records.extend(predict_six_frame_orfs(
                    genome["sequence"], sequence_id=genome["id"],
                    min_aa=min_aa, genetic_code=genetic_code,
                    alternative_starts=alternative_starts,
                ))
            except Exception as exc:
                errors.append(f"ORF screening for {genome['id']}: {exc}")

        analysis_proteins = list(protein_records) + [
            {"id": x["orf_id"], "sequence": x["protein_sequence"], "type": "protein"}
            for x in orf_records
        ]
        for error in errors:
            st.error(error)

        metrics = st.columns(4)
        metrics[0].metric("Sequences", len(sequence_rows))
        metrics[1].metric("Genome records", len(nucleotide_records))
        metrics[2].metric("ORF candidates", len(orf_records))
        metrics[3].metric("Protein records", len(analysis_proteins))

        if virus_group in {"Animal virus", "Plant virus"} and nucleotide_records:
            st.info("Six-frame ORF screening is preliminary and does not model splicing, RNA editing, overlapping transcripts or polyprotein processing.")
        if sequence_rows:
            st.subheader("Sequence inventory")
            st.dataframe(pd.DataFrame(sequence_rows), use_container_width=True, hide_index=True)
        if nucleotide_records:
            st.subheader("Genome composition")
            st.dataframe(pd.DataFrame([
                {"sequence_id": r["id"], **sequence_summary(r["sequence"])} for r in nucleotide_records
            ]), use_container_width=True, hide_index=True)
        if orf_records:
            st.subheader("Six-frame ORF candidates")
            orf_table = pd.DataFrame([
                {k: v for k, v in r.items() if k not in {"protein_sequence", "nucleotide_sequence"}}
                for r in orf_records
            ])
            st.dataframe(orf_table, use_container_width=True, hide_index=True)
            st.download_button(
                "Download predicted proteins (FASTA)",
                "".join(f">{r['orf_id']} {r['sequence_id']} strand={r['strand']} frame={r['frame']}\n{r['protein_sequence']}\n" for r in orf_records),
                "virasyn_predicted_orfs.faa", "text/plain",
            )
            st.download_button("Download ORF table (CSV)", orf_table.to_csv(index=False), "virasyn_orfs.csv", "text/csv")

        protein_rows = [{"protein_id": p["id"], **protein_properties(p["sequence"])} for p in analysis_proteins]
        if protein_rows:
            st.subheader("Protein physicochemical analysis")
            st.dataframe(pd.DataFrame(protein_rows), use_container_width=True, hide_index=True)

        if reference_fasta and analysis_proteins:
            with tempfile.TemporaryDirectory() as tmp:
                work = Path(tmp)
                q, ref = work / "query.faa", work / "reference.faa"
                q.write_text("".join(f">{p['id']}\n{p['sequence']}\n" for p in analysis_proteins), encoding="utf-8")
                ref.write_bytes(reference_fasta.getvalue())
                try:
                    hits = run_blastp(q, ref, work)
                    st.subheader("BLASTP reference comparison")
                    st.dataframe(hits, use_container_width=True, hide_index=True)
                    st.download_button("Download BLASTP hits", hits.to_csv(index=False), "virasyn_blastp_hits.csv", "text/csv")
                except Exception as exc:
                    st.warning(f"BLASTP could not run: {exc}")

        if hmm_library and analysis_proteins:
            with tempfile.TemporaryDirectory() as tmp:
                work = Path(tmp)
                q, hmm = work / "query.faa", work / "profiles.hmm"
                q.write_text("".join(f">{p['id']}\n{p['sequence']}\n" for p in analysis_proteins), encoding="utf-8")
                hmm.write_bytes(hmm_library.getvalue())
                try:
                    hits = run_hmmscan(q, hmm, work)
                    st.subheader("HMMER profile scan")
                    st.dataframe(hits, use_container_width=True, hide_index=True)
                    st.download_button("Download HMMER hits", hits.to_csv(index=False), "virasyn_hmmer_hits.csv", "text/csv")
                except Exception as exc:
                    st.warning(f"HMMER could not run: {exc}")

        if structure_records:
            st.subheader("Structure inventory")
            st.dataframe(pd.DataFrame(structure_records), use_container_width=True, hide_index=True)
            for filename, pdb_text in pdb_for_view:
                with st.expander(f"3D structure viewer · {filename}"):
                    try:
                        import streamlit.components.v1 as components
                        components.html(render_structure_html(pdb_text), height=540, scrolling=False)
                    except Exception as exc:
                        st.caption(f"3D viewer unavailable: {exc}")

        safe_orfs = [{k: v for k, v in r.items() if k not in {"protein_sequence", "nucleotide_sequence"}} for r in orf_records]
        report = make_report(virus_group, sequence_rows, safe_orfs, protein_rows, structure_records)
        st.download_button("Download complete analysis report (JSON)", report, "virasyn_report.json", "application/json")


with about_tab:
    st.subheader("A unified discovery and analysis workspace for virology")
    st.markdown("""
VIRASYN brings literature discovery and molecular data retrieval together with sequence and structure analysis.

**Directly integrated search providers**
- PubMed and Europe PMC — biomedical and life-sciences literature.
- OpenAlex and Crossref — scholarly work and DOI metadata.
- NCBI Nucleotide, Protein and Taxonomy.
- European Nucleotide Archive (ENA).
- RCSB Protein Data Bank.
- UniProtKB and InterPro.
- AlphaFold Protein Structure Database — published protein structure predictions, confidence metadata and coordinate downloads.
- PhARIS-inspired RBP intelligence — independent local six-frame ORF and BLASTP workflow using an authorized user-supplied reference library.
- Bacteria-focused Phage Intelligence — integrated RBP, receptor-annotation and curated bacterial host-evidence analysis.

VIRASYN retrieves records from upstream APIs when you search. It does not own, replace or mirror those databases. The Database Directory also links to specialist portals that are not yet directly queried by the app.
""")
    st.subheader("Scientific interpretation")
    st.markdown("""
VIRASYN is a bacteria-focused phage research platform. Its current Phage Intelligence module integrates reference-supported RBP homology with supplied host and receptor evidence; it is not yet a universal or calibrated RBP–receptor binding predictor. Sequence similarity, attachment, receptor binding, host range and productive infection are distinct levels of evidence.
""")
    st.subheader("Development roadmap")
    st.markdown("""
The PhARIS-inspired module is an independent implementation and does not bundle the original PhARIS reference sequences. The new bacteria-focused Phage Intelligence module adds evidence-ranked host and receptor hypotheses using researcher-supplied, provenance-tracked CSVs. Use only reference datasets you are authorized to use. Homology matches and evidence indices do not establish adsorption, receptor binding, infection or therapeutic suitability.

The next major development goal is to move VIRASYN towards a comprehensive viral research
knowledgebase through these planned capabilities:

- **Individual virus profile pages** that bring taxonomy, genome, protein, host and literature evidence together.
- **Saved research collections** for organizing records, searches and analysis outputs.
- **Automated genome annotation** with transparent methods, provenance and reviewable outputs.
- **Scientifically validated virus–host prediction models** developed for appropriate virus groups and evaluated
  against independent, well-curated evidence.

**Next five bioinformatics workbench priorities**

1. RNA-seq analysis — FastQC, differential expression, PCA, volcano plots and heatmaps.
2. Multiple sequence alignment and phylogenetics.
3. Protein structure comparison using TM-align and RMSD.
4. Molecular docking and protein–ligand interaction analysis.
5. Comparative genomics and automated genome annotation.

These are planned development areas, not capabilities of the current release. Each new analysis
workflow will need suitable reference data, reproducible methods and validation. Virus–host prediction
will require virus-group-appropriate datasets and independent validation rather than assuming one
universal model can cover bacteriophages, animal viruses and plant viruses.
""")
    st.caption("Cite the primary databases, records and software methods used in every analysis.")
