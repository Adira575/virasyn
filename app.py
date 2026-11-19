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

st.set_page_config(
    page_title="VIRASYN | Viral Systems Intelligence",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
.block-container {max-width: 1480px; padding-top: 1.6rem; padding-bottom: 3rem;}
.hero {padding: 1.4rem 1.6rem; border-radius: 18px; background: linear-gradient(120deg,#0c1b2b,#123c50 68%,#126b72); color:#f4fbff; margin-bottom:1.2rem;}
.hero h1 {font-size:2.6rem; letter-spacing:.10em; margin:0;}
.hero p {font-size:1rem; opacity:.88; margin:.35rem 0 0;}
div[data-testid="stMetric"] {background:rgba(120,150,170,.08); border:1px solid rgba(120,150,170,.18); padding:12px; border-radius:12px;}
</style>
<div class="hero">
  <h1>VIRASYN</h1>
  <p>Viral Systems Intelligence · Decode. Discover. Understand.</p>
</div>
""", unsafe_allow_html=True)

st.warning(
    "Research prototype. This release performs sequence screening, protein analysis, "
    "reference searches and structure inspection. It does not claim a universal or "
    "experimentally validated virus-host prediction model."
)

with st.sidebar:
    st.header("Analysis setup")
    virus_group = st.selectbox(
        "Virus group",
        ["Bacteriophage", "Animal virus", "Plant virus", "Archaeal virus", "Other / exploratory"],
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
    )
    min_aa = st.slider("Minimum ORF length (amino acids)", 20, 300, 45, 5)
    alternative_starts = st.checkbox("Include genetic-code alternative start codons", True)
    max_upload_mb = st.select_slider("Maximum upload size per file (MB)", options=[10, 25, 50, 100], value=50)
    st.divider()
    st.caption("VIRASYN v0.2 · Analysis workspace")

workspace_tab, guide_tab, methods_tab = st.tabs(
    ["Analysis workspace", "How to use", "Methods & limitations"]
)

with workspace_tab:
    st.subheader("1 · Upload viral genomes, proteins and structures")
    st.write("Submit multiple files in one analysis.")
    uploaded = st.file_uploader(
        "Genome/protein FASTA and PDB/mmCIF structures",
        type=["fa", "fasta", "fna", "faa", "ffn", "pep", "pdb", "ent", "cif", "mmcif"],
        accept_multiple_files=True,
        key="primary_uploads",
    )
    reference_fasta = st.file_uploader(
        "Optional BLASTP reference proteins (FASTA)",
        type=["fa", "fasta", "faa"],
        key="blast_reference",
        help="Use a curated reference set. FASTA deflines can carry metadata such as host=..., role=..., receptor=....",
    )
    hmm_library = st.file_uploader(
        "Optional HMMER profile library (.hmm)",
        type=["hmm"],
        key="hmm_library",
        help="VIRASYN runs hmmpress before hmmscan on a supplied HMM profile library.",
    )
    run = st.button("Run VIRASYN analysis", type="primary", disabled=not uploaded, use_container_width=True)

    if run:
        too_large = [f.name for f in uploaded if f.size > max_upload_mb * 1024 * 1024]
        if too_large:
            st.error("These files exceed the selected upload limit: " + ", ".join(too_large))
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
                        "file": file.name,
                        "sequence_id": record["id"],
                        "sequence_type": record["type"],
                        "length": len(record["sequence"]),
                    })
                    if record["type"] == "nucleotide":
                        nucleotide_records.append(record)
                    else:
                        protein_records.append(record)
            except Exception as exc:
                errors.append(f"{file.name}: {exc}")

        orf_records = []
        for genome in nucleotide_records:
            try:
                candidates = predict_six_frame_orfs(
                    genome["sequence"], sequence_id=genome["id"],
                    min_aa=min_aa, genetic_code=genetic_code,
                    alternative_starts=alternative_starts,
                )
                orf_records.extend(candidates)
            except Exception as exc:
                errors.append(f"ORF screening for {genome['id']}: {exc}")

        analysis_proteins = list(protein_records)
        for orf in orf_records:
            analysis_proteins.append({
                "id": orf["orf_id"],
                "description": f"predicted_ORF source={orf['sequence_id']} strand={orf['strand']} frame={orf['frame']}",
                "sequence": orf["protein_sequence"],
                "type": "protein",
            })

        for error in errors:
            st.error(error)

        st.success("Analysis completed for the files that could be processed.")
        metrics = st.columns(4)
        metrics[0].metric("Sequences", len(sequence_rows))
        metrics[1].metric("Genome records", len(nucleotide_records))
        metrics[2].metric("ORF candidates", len(orf_records))
        metrics[3].metric("Protein records", len(analysis_proteins))

        if virus_group in {"Animal virus", "Plant virus"} and nucleotide_records:
            st.info(
                "Six-frame ORF screening is preliminary. It does not model eukaryotic "
                "splicing, RNA editing, overlapping transcripts, polyprotein processing "
                "or all virus-specific expression strategies. Use curated protein FASTA "
                "or a validated virus-specific gene caller for interpretation."
            )

        if sequence_rows:
            st.subheader("Sequence inventory")
            st.dataframe(pd.DataFrame(sequence_rows), use_container_width=True, hide_index=True)

        if nucleotide_records:
            st.subheader("Genome composition")
            genome_table = pd.DataFrame([
                {"sequence_id": r["id"], **sequence_summary(r["sequence"])}
                for r in nucleotide_records
            ])
            st.dataframe(genome_table, use_container_width=True, hide_index=True)

        if orf_records:
            st.subheader("Six-frame ORF candidates")
            orf_table = pd.DataFrame([
                {k: v for k, v in record.items() if k not in {"protein_sequence", "nucleotide_sequence"}}
                for record in orf_records
            ])
            st.dataframe(orf_table, use_container_width=True, hide_index=True)
            st.download_button(
                "Download predicted proteins (FASTA)",
                "".join(
                    f">{r['orf_id']} {r['sequence_id']} strand={r['strand']} frame={r['frame']}\n{r['protein_sequence']}\n"
                    for r in orf_records
                ),
                "virasyn_predicted_orfs.faa",
                "text/plain",
            )
            st.download_button(
                "Download ORF table (CSV)", orf_table.to_csv(index=False),
                "virasyn_orfs.csv", "text/csv",
            )

        protein_rows = [
            {"protein_id": protein["id"], **protein_properties(protein["sequence"])}
            for protein in analysis_proteins
        ]
        if protein_rows:
            st.subheader("Protein physicochemical analysis")
            protein_table = pd.DataFrame(protein_rows)
            st.dataframe(protein_table, use_container_width=True, hide_index=True)
            st.download_button(
                "Download protein properties (CSV)", protein_table.to_csv(index=False),
                "virasyn_protein_properties.csv", "text/csv",
            )

        if reference_fasta and analysis_proteins:
            st.subheader("BLASTP reference comparison")
            with tempfile.TemporaryDirectory() as tmp:
                work = Path(tmp)
                query_path, reference_path = work / "query.faa", work / "reference.faa"
                query_path.write_text(
                    "".join(f">{p['id']}\n{p['sequence']}\n" for p in analysis_proteins),
                    encoding="utf-8",
                )
                reference_path.write_bytes(reference_fasta.getvalue())
                try:
                    blast_results = run_blastp(query_path, reference_path, work)
                    st.dataframe(blast_results, use_container_width=True, hide_index=True)
                    st.download_button(
                        "Download BLASTP hits (CSV)", blast_results.to_csv(index=False),
                        "virasyn_blastp_hits.csv", "text/csv",
                    )
                except Exception as exc:
                    st.warning(f"BLASTP could not run: {exc}")

        if hmm_library and analysis_proteins:
            st.subheader("HMMER profile scan")
            with tempfile.TemporaryDirectory() as tmp:
                work = Path(tmp)
                query_path, hmm_path = work / "query.faa", work / "profiles.hmm"
                query_path.write_text(
                    "".join(f">{p['id']}\n{p['sequence']}\n" for p in analysis_proteins),
                    encoding="utf-8",
                )
                hmm_path.write_bytes(hmm_library.getvalue())
                try:
                    hmm_results = run_hmmscan(query_path, hmm_path, work)
                    st.dataframe(hmm_results, use_container_width=True, hide_index=True)
                    st.download_button(
                        "Download HMMER domains (CSV)", hmm_results.to_csv(index=False),
                        "virasyn_hmmer_hits.csv", "text/csv",
                    )
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

        safe_orfs = [
            {k: v for k, v in r.items() if k not in {"protein_sequence", "nucleotide_sequence"}}
            for r in orf_records
        ]
        report = make_report(virus_group, sequence_rows, safe_orfs, protein_rows, structure_records)
        st.download_button(
            "Download complete analysis report (JSON)", report,
            "virasyn_report.json", "application/json", use_container_width=True,
        )

with guide_tab:
    st.subheader("Quick start")
    st.markdown("""
1. Choose the closest virus group in the sidebar.
2. Upload genome FASTA, protein FASTA and/or PDB/mmCIF files.
3. Adjust the ORF length threshold and genetic code for nucleotide genomes.
4. Optionally supply a curated protein reference FASTA for BLASTP or an HMM profile library for HMMER.
5. Click Run VIRASYN analysis and download the results.

Use .fna for nucleotide genomes and .faa for protein sequences whenever possible. For ambiguous .fa or .fasta files, sequence type is inferred from the alphabet.
""")
    st.subheader("Reference FASTA metadata convention")
    st.code(">RefProtein01 host=Escherichia_coli role=tail_fiber receptor=OmpC virus_group=bacteriophage\nMKK...", language="text")
    st.caption("The defline is retained in BLASTP results. Metadata is not treated as validated truth; curate and cite your reference set.")

with methods_tab:
    st.subheader("Implemented in this release")
    st.markdown("""
- Six-frame start-to-stop ORF screening with selectable NCBI genetic code.
- Protein physicochemical summaries: molecular weight, estimated pI, instability index, GRAVY and aromaticity.
- BLASTP and HMMER searches against reference data supplied by the researcher.
- PDB/mmCIF chain and residue inventory, plus an interactive viewer for PDB files.
- CSV, FASTA and JSON exports.
""")
    st.subheader("Scientific limitations")
    st.markdown("""
- This is not a universal virus-host prediction model.
- Six-frame ORF screening is not splice-aware or virus-specific gene prediction.
- Similarity does not prove function, host range, receptor binding or infection.
- Structure intake and visualisation do not establish protein function.
- Different virus groups need separate datasets, biological assumptions and validation.
- Attachment, receptor binding, entry and productive infection are distinct endpoints.
""")

st.divider()
st.caption("VIRASYN · Independent research software · Cite original methods and databases used in each analysis.")
