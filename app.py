from __future__ import annotations

from pathlib import Path
import tempfile

import pandas as pd
import streamlit as st

from virasyn.core import (
    make_report, parse_fasta, parse_structure, predict_six_frame_orfs,
    protein_properties, render_structure_html, run_blastp, run_hmmscan,
    sequence_summary,
)
from virasyn.database import (
    fetch_sequence_fasta, search_all_sources, search_ncbi_sequences,
    search_pubmed, search_rcsb, search_uniprot,
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

database_tab, analysis_tab, about_tab = st.tabs(
    ["🔎 Research database", "🧬 Sequence & structure analysis", "ℹ️ About VIRASYN"]
)

def render_record(record: dict, kind: str):
    with st.container(border=True):
        title = record.get("title") or record.get("accession") or record.get("id") or "Untitled record"
        st.markdown(f'<div class="record-title">{title}</div>', unsafe_allow_html=True)
        meta = []
        for key in ("authors", "journal", "date", "organism", "length", "method", "doi"):
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
            ["All sources", "PubMed articles", "NCBI nucleotide genomes", "NCBI protein sequences",
             "RCSB PDB structures", "UniProtKB proteins"],
        )
        c1, c2, c3 = st.columns([1, 1, 1])
        with c1:
            page_size = st.selectbox("Results per page", [10, 20, 50], index=1)
        with c2:
            date_from = st.text_input("Published from (PubMed only)", placeholder="YYYY-MM-DD")
        with c3:
            date_to = st.text_input("Published to (PubMed only)", placeholder="YYYY-MM-DD")
        sort = st.selectbox("Sort publications", ["Best match", "Newest first"])
        submitted = st.form_submit_button("Search research databases", type="primary", use_container_width=True)

    if submitted:
        if not query.strip():
            st.warning("Enter a search term first.")
        else:
            try:
                with st.spinner("Searching scientific databases..."):
                    if source == "All sources":
                        bundle = search_all_sources(query.strip(), 0, min(5, page_size))
                        st.session_state["virasyn_search_bundle"] = bundle
                        st.session_state["virasyn_search_query"] = query.strip()
                        st.session_state["virasyn_search_source"] = source
                    elif source == "PubMed articles":
                        result = search_pubmed(query.strip(), 0, page_size, date_from, date_to, sort)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    elif source == "NCBI nucleotide genomes":
                        result = search_ncbi_sequences("nuccore", query.strip(), 0, page_size)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    elif source == "NCBI protein sequences":
                        result = search_ncbi_sequences("protein", query.strip(), 0, page_size)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    elif source == "RCSB PDB structures":
                        result = search_rcsb(query.strip(), 0, page_size)
                        st.session_state["virasyn_search_result"] = result
                        st.session_state["virasyn_search_source"] = source
                        st.session_state["virasyn_search_query"] = query.strip()
                    else:
                        result = search_uniprot(query.strip(), 0, page_size)
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
                    current_source
                )
                render_record(record, kind)


with analysis_tab:
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

**Live external data sources**
- PubMed — biomedical literature.
- NCBI Nucleotide — nucleotide sequence records.
- NCBI Protein — protein sequence records.
- RCSB Protein Data Bank — macromolecular structures.
- UniProtKB — protein knowledgebase.

VIRASYN retrieves records from these services when you search. It does not claim to own, replace or reproduce their underlying databases. Record pages link back to the original providers.
""")
    st.subheader("Scientific interpretation")
    st.markdown("""
VIRASYN is not yet a universal virus-host prediction system. Bacteriophage RBPs, animal-virus attachment/entry proteins and plant-virus host interactions require separate training data and methods. Sequence similarity, attachment, receptor binding, host range and productive infection are distinct levels of evidence.
""")
    st.caption("Cite the primary databases, records and software methods used in every analysis.")
