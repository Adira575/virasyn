# VIRASYN — Viral Systems Intelligence

**Decode. Discover. Understand.**

VIRASYN is a browser-based viral research discovery and analysis workspace. It combines federated search across public literature and molecular-data services with local sequence/structure workflows.

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

## PhARIS-inspired RBP intelligence

VIRASYN includes an independent, reference-driven workflow inspired by the published PhARIS method for screening bacteriophage genomes for candidate receptor-binding proteins associated with *Staphylococcus aureus*. It accepts a phage genome FASTA, an authorized RBP protein reference FASTA, and optional reference metadata CSV. The workflow screens six reading frames for start-to-stop ORFs, builds a temporary local BLASTP database, compares candidate proteins to the supplied reference library, calculates query and subject coverage, applies configurable identity/coverage thresholds, and exports candidate and all-hit tables.

The optional metadata CSV must include a `subject_id` column matching the reference FASTA record identifiers. Suggested additional columns are `cluster`, `RBP_name`, `strain`, `receptor`, `evidence`, and `citation`. The reference FASTA and its provenance/version should be maintained by the researcher. VIRASYN does not redistribute the original PhARIS reference sequences or notebook; obtain permission and follow the source terms before using those materials.

This workflow is a homology screen, not a validated binding-specificity predictor. Sequence similarity does not establish adsorption, receptor binding, productive infection, host range, or therapeutic suitability. Experimental validation is required. NCBI BLAST+ is needed; it is installed through `packages.txt` on Streamlit Community Cloud and in the Docker image.

Original project: https://github.com/JKrusche1/PhARIS  
Publication: https://doi.org/10.1016/j.celrep.2025.115369

## Bacteria-focused Phage Intelligence

See the [detailed architecture and validation plan](docs/PHAGE_INTELLIGENCE_ARCHITECTURE.md) for the data contracts, evidence model, current boundaries and model-development roadmap.

The Phage Intelligence workspace combines the existing reference-driven RBP discovery workflow with an evidence integration layer for bacterial host and receptor hypotheses.

### Current implemented capabilities
- Six-frame ORF screening from one or multiple bacteriophage nucleotide records.
- Local BLASTP comparison against a researcher-supplied, permission-cleared RBP protein reference FASTA.
- Configurable minimum ORF length, sequence identity, query coverage and translation table.
- Optional RBP metadata joined by exact FASTA subject identifier.
- Optional bacterial host evidence CSV with explicit evidence levels (experimental, curated, computational, inferred) and evidence types (productive infection, host range, adsorption, receptor binding, CRISPR spacer, sequence similarity).
- Transparent evidence-weighted host ranking; the index is not a probability and is not a trained predictor.
- Optional RBP-to-receptor annotation mapping by exact reference identifier.
- Optional PDB/mmCIF structure inventory and integrated JSON/CSV exports.
- Downloadable CSV templates for host evidence and receptor annotations.

### Required CSV schemas

Host evidence requires: `phage_id, host_taxid, host_name, evidence_type, evidence_level, source, citation`.

RBP–receptor annotations require: `subject_id, receptor, receptor_class, host_taxid, host_name, evidence_type, citation`. The `subject_id` must exactly match the protein identifier in the RBP reference FASTA.

### Scientific boundary
This is an evidence-integration and homology-screening release, not a universal RBP–receptor interaction model. The current host ranking summarizes supplied records; it does not infer unseen hosts from sequence alone. Structural file intake is an inventory step, not docking or a binding-energy calculation. Future work must add curated training labels, independent taxonomic holdouts, calibrated evaluation, structural interface methods and experimental validation before claiming generalized receptor specificity.

## AlphaFold protein structure lookup

The AlphaFold tab retrieves published predictions from the AlphaFold Protein Structure Database by UniProt accession. It displays available model metadata and the global confidence metric (pLDDT where supplied), supports PDB and mmCIF coordinate downloads, and renders PDB coordinates in the built-in 3D viewer. The app links to the official AlphaFold Server for submitting new sequence-prediction jobs; VIRASYN does not run AlphaFold inference itself. A predicted structure and its confidence scores are not experimental validation.

## PCR / qPCR analysis

VIRASYN includes a qPCR workbench for RDML, CSV, TSV, TXT, XLSX and XLSM exports containing instrument-reported Cq/Ct values. RDML XML (including zipped RDML packages) can also provide amplification and melt-curve traces for visualization. Users map sample, target, Cq/Ct, condition and optional sample-role columns.

Features include:
- Technical-well aggregation within biological sample × condition × target.
- Missing/undetermined Cq values are retained as missing, never converted to zero.
- NTC identification and positive-NTC warnings.
- Technical replicate mean, SD, range and review flags (no automatic outlier deletion).
- Optional standard-curve fitting by assay using Cq versus log10 concentration, with slope, R² and estimated amplification efficiency.
- Editable per-assay efficiency values when validated efficiency data are available.
- Efficiency-corrected relative quantities, normalized to the geometric mean of selected reference-gene quantities.
- Group fold changes versus a selected calibrator and Welch tests on per-biological-sample log2 normalized quantities, with Benjamini–Hochberg q-values.
- Downloadable Excel workbook and CSV outputs.

The module starts from Cq/Ct values already called by the instrument software. RDML amplification and melt-curve traces can be displayed, but the app does not yet recalculate Cq from raw cycle-by-cycle fluorescence traces or decode proprietary binary instrument projects. The 100% efficiency values shown when no standard curve is supplied are explicit assumptions, not measured assay performance. Review reference-gene stability, efficiency, controls and replicate flags before interpreting results. This is a research workflow, not a diagnostic test.

Methods are informed by the MIQE 2.0 recommendations for transparent qPCR data processing, efficiency-aware quantification, normalization and quality control: https://doi.org/10.1093/clinchem/hvaf043

## Bioinformatics workbench

The Bioinformatics Workbench page in the app presents this module catalogue. Status labels describe the current implementation rather than the long-term ambition.

| Module | Capabilities | Current status |
|---|---|---|
| Genome analysis | FASTA import, sequence summaries, GC content, six-frame ORF screening and downloadable ORFs. FASTQ quality control, FastQC integration and automated genome annotation are planned. | Partial |
| Sequence analysis | BLASTP against user-provided protein references and federated sequence-database discovery. Multiple sequence alignment and phylogenetic analysis are planned. | Partial |
| Protein analysis | Protein physicochemical properties and optional HMMER profile scanning with user-provided HMM libraries. Automated motif/domain annotation and richer functional annotation are planned. | Partial |
| Structural bioinformatics | PDB/mmCIF metadata parsing and PDB 3D visualization. Structural alignment, TM-align integration and RMSD comparison are planned. | Partial |
| Molecular docking | Protein and ligand preparation, docking workflows, pose scoring and protein–ligand interaction analysis. | Planned |
| Transcriptomics | RNA-seq quality control, differential expression, PCA, volcano plots and heatmaps. | Planned |
| Comparative genomics | Genome-to-genome comparison, gene conservation, orthology and synteny analysis, with automated annotation links. | Planned |
| Virology | Viral literature, genome, protein and taxonomy discovery plus preliminary sequence characterization. Validated candidate attachment-protein workflows and virus–host evidence integration remain planned. | Partial |
| Data visualization | Analysis tables, qPCR fold-change charts and interactive PDB viewing. Publication-ready figure templates, PCA plots, volcano plots and heatmaps are planned. | Partial |
| PCR / qPCR analysis | RDML and tabular Cq/Ct import, technical-replicate QC, standard curves, efficiency-aware relative expression, NTC flags and Excel/CSV reports. | Available |

### Next five development priorities

1. **RNA-seq analysis** — FastQC, differential expression, PCA, volcano plots and heatmaps.
2. **Multiple sequence alignment and phylogenetics.**
3. **Protein structure comparison** — TM-align and RMSD.
4. **Molecular docking** — protein/ligand preparation, docking workflows and protein–ligand interaction analysis.
5. **Comparative genomics and automated genome annotation.**

These are planned additions, not active features. Each workflow will need reproducible parameters, inspectable intermediate results, clear software provenance and appropriate validation.

## Development roadmap

The next major development goal is to move VIRASYN towards a comprehensive viral research
knowledgebase through four planned capabilities:

1. **Individual virus profile pages** that bring taxonomy, genome, protein, host and literature evidence together.
2. **Saved research collections** for organizing records, searches and analysis outputs.
3. **Automated genome annotation** with transparent methods, provenance and reviewable outputs.
4. **Scientifically validated virus–host prediction models** developed for appropriate virus groups and evaluated against independent, well-curated evidence.

These are planned development areas, not current-release features. Virus–host prediction will require
virus-group-appropriate datasets and independent validation; a single universal model should not be
assumed to cover bacteriophages, animal viruses and plant viruses.

## Navigation

The main workspace is organized into nine sections: Research database, Sequence & structure analysis, PCR / qPCR analysis, Bioinformatics workbench, AlphaFold, PhARIS RBP intelligence, Phage Intelligence, Database directory and About VIRASYN.

## Sequence and structure analysis

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

The application is configured for Streamlit Community Cloud. Select repository `Adira575/virasyn`, branch `main`, main file `app.py`. No AI API key is required. Optional NCBI API credentials and a contact email can be configured as server-side secrets.

## Scientific limitations

qPCR outputs are dependent on assay efficiency, validated reference genes, sample design, negative controls and instrument Cq calling. Technical replicates must not be counted as independent biological samples. The qPCR module flags high technical replicate spread and detected NTC wells but does not automatically exclude data or establish assay validity.


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
