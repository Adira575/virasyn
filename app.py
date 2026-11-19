from pathlib import Path
import tempfile

import pandas as pd
import streamlit as st

from virasyn.core import make_report, parse_fasta, parse_pdb, run_blastp, sequence_summary

st.set_page_config(page_title="VIRASYN | Viral Systems Intelligence", page_icon="🧬", layout="wide")
st.title("🧬 VIRASYN")
st.caption("Viral Systems Intelligence · Decode. Discover. Understand.")
st.info(
    "Research software prototype. VIRASYN currently provides sequence intake, "
    "basic PDB metadata and optional similarity search. It does not yet provide "
    "validated host-range, receptor-binding or infectivity predictions."
)

with st.sidebar:
    st.header("Analysis context")
    virus_group = st.selectbox(
        "Virus group",
        ["Bacteriophage", "Animal virus", "Plant virus", "Archaeal virus", "Other / exploratory"],
    )
    st.caption("This labels the analysis context; it does not select a trained model.")

st.header("1 · Upload sequences and structures")
uploads = st.file_uploader(
    "Upload FASTA (.fa, .fasta, .fna, .faa) and/or PDB (.pdb, .ent) files",
    type=["fa", "fasta", "fna", "faa", "fas", "pdb", "ent"],
    accept_multiple_files=True,
)
sequence_rows, protein_records, structure_rows = [], [], []

if uploads:
    for upload in uploads:
        name = upload.name
        text = upload.getvalue().decode("utf-8", errors="replace")
        if Path(name).suffix.lower() in {".pdb", ".ent"}:
            try:
                structure_rows.append(parse_pdb(text, name))
            except Exception as exc:
                st.error(f"Could not parse {name}: {exc}")
            continue
        try:
            for record in parse_fasta(text):
                sequence_rows.append({
                    "file": name,
                    "sequence_id": record["id"],
                    "sequence_type": record["type"],
                    "length": len(record["sequence"]),
                })
                if record["type"] == "protein":
                    protein_records.append(record)
        except Exception as exc:
            st.error(f"Could not parse {name}: {exc}")

    if sequence_rows:
        st.header("2 · Sequence overview")
        st.dataframe(pd.DataFrame(sequence_rows), use_container_width=True, hide_index=True)
        st.metric("Sequences parsed", len(sequence_rows))
        with st.expander("Protein sequence summaries"):
            for record in protein_records:
                st.markdown(f"**{record['id']}**")
                st.json(sequence_summary(record["sequence"]))

    if structure_rows:
        st.header("3 · PDB structure overview")
        st.dataframe(pd.DataFrame(structure_rows), use_container_width=True, hide_index=True)

    if protein_records:
        st.header("4 · Optional BLASTP reference comparison")
        st.write("Upload a curated protein FASTA reference set. Results depend on reference quality and biological scope.")
        ref_upload = st.file_uploader("Reference protein FASTA", type=["fa", "fasta", "faa"], key="ref")
        if ref_upload and st.button("Run BLASTP", type="primary"):
            with tempfile.TemporaryDirectory() as tmp:
                work = Path(tmp)
                query = work / "query.faa"
                reference = work / "reference.faa"
                query.write_text(
                    "".join(f">{r['id']}\n{r['sequence']}\n" for r in protein_records),
                    encoding="utf-8",
                )
                reference.write_bytes(ref_upload.getvalue())
                try:
                    hits = run_blastp(query, reference, work)
                    st.dataframe(hits, use_container_width=True, hide_index=True)
                    st.download_button(
                        "Download BLASTP results (CSV)",
                        hits.to_csv(index=False),
                        "virasyn_blastp_hits.csv",
                        "text/csv",
                    )
                except Exception as exc:
                    st.error(str(exc))

    if st.button("Prepare analysis summary"):
        summary = make_report(virus_group, sequence_rows, structure_rows)
        st.download_button(
            "Download summary (JSON)",
            summary,
            "virasyn_analysis_summary.json",
            "application/json",
        )

st.divider()
st.header("Analysis modules")
left, right = st.columns(2)
with left:
    st.subheader("Genome & protein analysis")
    st.write("FASTA validation, sequence statistics, similarity search and future annotation modules.")
    st.subheader("Virus-specific interaction analysis")
    st.write("Separate host and receptor workflows for bacteriophages, animal viruses and plant viruses.")
with right:
    st.subheader("Structure analysis")
    st.write("PDB intake is available; structure quality and alignment are planned.")
    st.subheader("Reproducible reporting")
    st.write("Export machine-readable summaries and analysis results.")
st.caption("Similarity is evidence for relatedness, not proof of biological function or infectivity.")
