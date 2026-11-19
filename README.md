# VIRASYN — Viral Systems Intelligence

An independent, modular research platform for viral genome, protein and structure analysis.

## Scientific scope
VIRASYN is designed for multiple virus groups through separate analysis modules. Bacteriophage receptor-binding proteins, animal-virus attachment/entry proteins, and plant-virus host-interaction proteins are distinct biological tasks and require separate evidence, datasets and models.

## Current MVP
- Upload nucleotide or protein FASTA files.
- Infer sequence type and summarise sequences.
- Upload PDB structures and extract basic model, chain and residue metadata.
- Optionally run BLASTP against a user-supplied curated protein FASTA.
- Download JSON summaries and BLASTP CSV results.

## Run locally
```bash
python -m venv .venv
# Windows PowerShell
.venv\Scripts\Activate.ps1
# macOS/Linux
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

BLASTP requires NCBI BLAST+ executables `blastp` and `makeblastdb`. The Docker image installs BLAST+.

## Docker
```bash
docker build -t virasyn .
docker run --rm -p 8501:8501 virasyn
```

Open http://localhost:8501.

## Limitations
This starter does not provide validated host-range, receptor-binding, tropism, pathogenicity or infectivity predictions. Sequence similarity does not establish function. Genome ORF prediction and structural alignment are not yet implemented.

## Relationship to PhARIS
VIRASYN is an independent project inspired by the general computational problem addressed by PhARIS, not an official PhARIS release. Check upstream licences and seek permission before reusing source code, datasets, logos or branding. Cite the original work when comparing methods.

See [ROADMAP.md](ROADMAP.md).
