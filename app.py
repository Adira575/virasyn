from __future__ import annotations

from pathlib import Path
from io import BytesIO
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

database_tab, analysis_tab, qpcr_tab, resources_tab, about_tab = st.tabs(
    ["🔎 Research database", "🧬 Sequence & structure analysis", "🧪 PCR / qPCR analysis",
     "📚 Database directory", "ℹ️ About VIRASYN"]
)

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
        "Supported inputs: CSV, TSV, TXT and Excel XLSX/XLSM files containing reported Cq/Ct values. "
        "This release does not yet derive Cq from raw cycle-by-cycle fluorescence traces or parse "
        "proprietary binary instrument projects. Export a Cq/Ct table from your instrument software first."
    )
    qpcr_file = st.file_uploader(
        "Upload qPCR plate/result export",
        type=["csv", "tsv", "txt", "xlsx", "xlsm"],
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
                    type=["csv", "tsv", "txt", "xlsx", "xlsm"],
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
            ("AlphaFold Protein Structure Database", "https://alphafold.ebi.ac.uk/", "External portal · predicted protein structures"),
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

with about_tab:
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

VIRASYN retrieves records from upstream APIs when you search. It does not own, replace or mirror those databases. The Database Directory also links to specialist portals that are not yet directly queried by the app.
""")
    st.subheader("Scientific interpretation")
    st.markdown("""
VIRASYN is not yet a universal virus-host prediction system. Bacteriophage RBPs, animal-virus attachment/entry proteins and plant-virus host interactions require separate training data and methods. Sequence similarity, attachment, receptor binding, host range and productive infection are distinct levels of evidence.
""")
    st.subheader("Development roadmap")
    st.markdown("""
The next major development goal is to move VIRASYN towards a comprehensive viral research
knowledgebase through four planned capabilities:

- **Individual virus profile pages** that bring taxonomy, genome, protein, host and literature evidence together.
- **Saved research collections** for organizing records, searches and analysis outputs.
- **Automated genome annotation** with transparent methods, provenance and reviewable outputs.
- **Scientifically validated virus–host prediction models** developed for appropriate virus groups and evaluated
  against independent, well-curated evidence.

These are planned development areas, not capabilities of the current release. Virus–host prediction
will require virus-group-appropriate datasets and independent validation rather than assuming one
universal model can cover bacteriophages, animal viruses and plant viruses.
""")
    st.caption("Cite the primary databases, records and software methods used in every analysis.")
