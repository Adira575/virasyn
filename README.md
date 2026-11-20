# VIRASYN — Viral Systems Intelligence

**Decode. Discover. Understand.**

VIRASYN is a browser-based viral research discovery and analysis workspace. It combines federated search across public literature and molecular-data services with local sequence/structure workflows and a floating research assistant.

## Research database integrations

The federated search can query:
- **PubMed** — biomedical publications.
- **Europe PMC** — life-sciences articles, preprints and links to full text.
- **OpenAlex** — scholarly works and citation metadata.
- **Crossref** — DOI and publication metadata.
- **NCBI Nucleotide** — nucleotide sequence records.
- **NCBI Protein** — protein sequence records.
- **NCBI Taxonomy** — organism taxonomy.
- **European Nucleotide Archive (ENA)** — nucleotide archive records.
- **RCSB Protein Data Bank** — macromolecular structures.
- **UniProtKB** — protein knowledgebase.
- **InterPro** — protein families, domains and functional signatures.

The Database Directory also links to specialist resources including ICTV taxonomy, NCBI Virus, PHROGs, PhagesDB, INPHARED, ViralZone, Virus-Host DB, BV-BRC, Nextstrain, PDBe, Pfam, AlphaFold DB and VEuPathDB. These directory links are external portals; they are not all queried directly by the app.

VIRASYN retrieves results from upstream providers and links to the source records. It does not mirror or replace those databases. Data coverage, licensing, citation requirements and query syntax differ by provider.

## Floating VIRASYN AI assistant

A compact bottom-right helper is available across the site. It can explain the interface, search sources, ORF results, protein properties, BLASTP/HMMER outputs and scientific limitations.

### Enable generative answers

The assistant has a basic quick-help mode without an API key. To enable generative answers, add a server-side OpenAI API key in **Streamlit Community Cloud → App settings → Secrets**:

~~~toml
OPENAI_API_KEY = "your_api_key"
VIRASYN_AI_MODEL = "gpt-6.2-astra"
~~~

The model setting is optional; use a model that is available to your API account. API usage is billed separately from a ChatGPT subscription. Never put the API key in browser-side code or commit it to GitHub.

The assistant is a research aid, not a clinical decision system. Verify scientific claims against original records and literature. It should not be treated as evidence that a predicted protein binds a receptor or that a virus infects a host.

## Analysis features

- Multi-file nucleotide/protein FASTA upload.
- Six-frame preliminary ORF screening with selectable NCBI genetic code.
- Protein physicochemical properties.
- Optional BLASTP against researcher-provided protein references.
- Optional HMMER profile scanning against researcher-provided HMM libraries.
- PDB/mmCIF metadata and PDB 3D viewer.
- FASTA, CSV and JSON exports.

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

The application is configured for Streamlit Community Cloud. Select repository `Adira575/virasyn`, branch `main`, main file `app.py`. Add optional API keys under App settings → Secrets, not as committed files.

## Scientific limitations

VIRASYN is not a validated universal virus-host predictor. Six-frame ORF screening is not a complete virus-specific gene caller and does not model all eukaryotic splicing, RNA editing, transcript processing or polyprotein strategies. Similarity does not prove function. Attachment, receptor binding, host range and productive infection are distinct biological endpoints. Each virus group requires suitable reference data and independent validation.

## External data and attribution

Follow the terms, citation guidance and rate limits of each upstream service. NCBI E-utilities requests are cached and rate-limited; add `VIRASYN_NCBI_API_KEY` and `VIRASYN_CONTACT_EMAIL` as Streamlit secrets where appropriate.

- NCBI E-utilities: https://eutilities.github.io/site/
- Europe PMC REST API: https://www.ebi.ac.uk/europepmc/webservices/rest/
- ENA Portal API: https://ena-docs.readthedocs.io/
- InterPro API: https://interpro-documentation.readthedocs.io/en/latest/api.html
- OpenAlex API: https://help.openalex.org/api/
- Crossref REST API: https://www.crossref.org/documentation/retrieve-metadata/rest-api/
- RCSB PDB APIs: https://www.rcsb.org/docs/programmatic-access/web-apis-overview
- ICTV taxonomy: https://ictv.global/taxonomy

## Relationship to PhARIS

VIRASYN is an independent project inspired by the computational problem addressed by PhARIS; it is not an official PhARIS release. Check the upstream licence and obtain permission before reusing PhARIS code, reference data, logo or branding. Cite the original work when comparing methods.
