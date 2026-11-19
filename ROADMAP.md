# VIRASYN development roadmap

## Release 0.2 — Interactive analysis MVP
- Multi-file FASTA upload and validation.
- Six-frame ORF screening with selectable NCBI genetic code.
- Protein physicochemical analysis.
- BLASTP and HMMER against user-supplied reference databases.
- PDB/mmCIF metadata and interactive PDB viewer.
- CSV, FASTA and JSON exports.
- Docker deployment and automated tests.

## Stage 2 — Genome and protein annotation
- Genome quality checks and virus-aware ORF callers.
- Automated functional annotation.
- Curated domain databases and database-version provenance.
- Viral protein family and architecture classification.

## Stage 3 — Virus-specific analysis modules
- Bacteriophages: tail fibres, tail spikes, receptor-binding proteins and bacterial host/receptor evidence.
- Animal viruses: attachment/entry proteins, host-cell receptor evidence and tropism datasets.
- Plant viruses: movement/coat proteins, host-factor interactions and vector-associated biology.
- Other virus groups only where suitable references and methods exist.

## Stage 4 — Advanced structure analysis
- Structure quality/confidence metadata and domain architecture.
- Structural similarity and domain comparisons.
- Interactive visualisation and reproducible export.

## Stage 5 — Host and interaction prediction
- Integrate existing tools only within their validated scope and licence.
- Separate host association, receptor recognition, entry and productive infection.
- Train new models only with curated labelled data.
- Use cluster-aware data splits to reduce homolog leakage and independent test sets.
- Report precision, recall, F1, MCC, PR-AUC, calibration and uncertainty.

## Stage 6 — Production readiness
- Background jobs, upload quotas and user/job isolation.
- Persistent job history and versioned reference databases.
- Deployment monitoring and reproducible analysis manifests.
