# VIRASYN — Viral Systems Intelligence

Decode. Discover. Understand.

VIRASYN is an independent, modular research website for exploratory viral genome, protein and structure analysis. It is designed to support several virus groups while keeping their biological prediction tasks separate.

## Run the website locally

Requirements: Python 3.11 or newer. Optional system tools for reference searches: NCBI BLAST+ (blastp, makeblastdb) and HMMER (hmmscan, hmmpress).

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

Open the local URL printed by Streamlit, normally http://localhost:8501.

## Docker

~~~bash
docker build -t virasyn .
docker run --rm -p 8501:8501 virasyn
~~~

Then open http://localhost:8501. The Docker image installs BLAST+ and HMMER.

## Deploy a public website

For a simple public deployment, use Streamlit Community Cloud:
1. Sign in at https://share.streamlit.io/ using GitHub.
2. Choose repository Adira575/virasyn.
3. Choose branch main.
4. Set the main file path to app.py.
5. Deploy.

The repository includes packages.txt for BLAST+ and HMMER system packages. Deployment depends on hosting resource quotas and system-package support. A Docker-compatible host can use the supplied Dockerfile instead.

## Features in Release 0.2

- Multiple nucleotide/protein FASTA files.
- FASTA validation and sequence-type inference; use .fna for nucleotide and .faa for proteins where possible.
- Six-frame start-to-stop ORF screening, selectable NCBI genetic code, configurable minimum length and predicted protein FASTA export.
- Protein physicochemical summaries: molecular weight, estimated pI, instability index, GRAVY and aromaticity.
- Optional BLASTP comparison against a researcher-supplied protein FASTA.
- Optional HMMER profile scanning against a researcher-supplied HMM library.
- PDB and mmCIF metadata inspection; interactive PDB structure viewer.
- CSV, FASTA and JSON downloads.
- Automated unit tests and GitHub Actions.

## Scientific scope and limitations

VIRASYN is not a validated universal virus-host predictor. Six-frame ORF screening is preliminary and is not a complete virus-specific gene caller. It does not model eukaryotic splicing, RNA editing, overlapping transcripts, polyprotein processing or every viral expression strategy. For complex animal and plant virus genomes, use curated protein FASTA or a validated virus-specific gene-calling workflow.

BLASTP/HMMER results depend on the reference libraries selected by the researcher. Similarity is not proof of function. Attachment, receptor binding, host range and productive infection are distinct biological endpoints. Each virus group needs curated labels, a biologically appropriate model and independent validation.

## Tests

~~~bash
pip install -e ".[test]"
pytest -q
~~~

## Relationship to PhARIS

VIRASYN is an independent project inspired by the computational problem addressed by PhARIS; it is not an official PhARIS release. Before reusing PhARIS source code, reference sequences, logo or branding, inspect its licence and obtain permission where required. Cite the original paper and repository when comparing methods.

See ROADMAP.md for future development priorities.
