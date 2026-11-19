# VIRASYN — Viral Systems Intelligence

**Decode. Discover. Understand.**

VIRASYN is a browser-based research discovery and analysis workspace. Its interface is inspired by literature databases such as PubMed: search, browse results, open original records and move from discovery into sequence and structure analysis.

## Live research database search

- **PubMed** — article title, authors, journal, date and DOI via NCBI E-utilities.
- **NCBI Nucleotide** — nucleotide sequence record metadata and FASTA retrieval.
- **NCBI Protein** — protein sequence record metadata and FASTA retrieval.
- **RCSB PDB** — structure search and links to original structures.
- **UniProtKB** — protein entry search and links to UniProt records.
- **All sources** — federated search across the providers above.

Results are retrieved live from each provider and cached briefly. VIRASYN links back to original records; it does not mirror or replace PubMed, GenBank, PDB or UniProt.

## Run locally

Requirements: Python 3.11 or newer. Optional BLAST+ and HMMER are needed for local reference searches.

~~~bash
git clone https://github.com/Adira575/virasyn.git
cd virasyn
python -m venv .venv
# Windows PowerShell
.venv\Scripts\Activate.ps1
# macOS/Linux
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
~~~

Open http://localhost:8501.

## Deploy

The application is configured for Streamlit Community Cloud. Select repository `Adira575/virasyn`, branch `main`, main file `app.py`.

## Analysis features

- Multi-file nucleotide/protein FASTA upload.
- Six-frame preliminary ORF screening with selectable NCBI genetic code.
- Protein physicochemical properties.
- Optional BLASTP against user-provided references.
- Optional HMMER profile scanning.
- PDB/mmCIF metadata and PDB 3D viewer.
- FASTA, CSV and JSON exports.

## Scientific limitations

VIRASYN is not a validated universal virus-host predictor. Six-frame ORF screening is not a complete virus-specific gene caller and does not model all eukaryotic splicing, RNA editing, transcript processing or polyprotein strategies. Similarity does not prove function. Attachment, receptor binding, host range and productive infection are distinct biological endpoints. Each virus group requires suitable reference data and independent validation.

## External data and attribution

Follow the terms, citation guidance and rate limits of each upstream service. NCBI E-utilities requests are cached to reduce repeated traffic. More than three requests per second requires an NCBI API key under NCBI's published usage guidance. Set `VIRASYN_NCBI_API_KEY` as an environment variable or Streamlit secret if you have one.

- NCBI E-utilities: https://eutilities.github.io/site/
- NCBI Datasets: https://www.ncbi.nlm.nih.gov/datasets/
- RCSB PDB APIs: https://www.rcsb.org/docs/programmatic-access/web-apis-overview
- UniProt API: https://www.uniprot.org/api-documentation

## Relationship to PhARIS

VIRASYN is an independent project inspired by the computational problem addressed by PhARIS; it is not an official PhARIS release. Check the upstream licence and obtain permission before reusing PhARIS code, reference data, logo or branding. Cite the original work when comparing methods.
