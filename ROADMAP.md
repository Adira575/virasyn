# VIRASYN development roadmap

## Stage 1 — Input and reproducible baseline (current)
- FASTA parsing and sequence-type detection.
- PDB metadata intake.
- Optional BLASTP against user-supplied protein references.
- Downloadable JSON report.

## Stage 2 — Genome and protein annotation
- Genome quality checks.
- ORF prediction with a method appropriate to genome type.
- Protein translation and annotation.
- Curated profile/domain databases with explicit version provenance.

## Stage 3 — Virus-specific modules
- Bacteriophages: tail fibre, tail spike and RBP candidates; bacterial host/receptor evidence.
- Animal viruses: attachment/entry proteins and host-cell receptor evidence.
- Plant viruses: coat/movement proteins and host-factor evidence.
- Other virus groups only where suitable reference data and validated methods exist.

## Stage 4 — Structure analysis
- PDB/mmCIF validation and quality metadata.
- Structure comparison using established software.
- Clear distinction between structural similarity and established biological function.

## Stage 5 — Host and interaction modelling
- Integrate existing tools only within their validated scope and licence.
- Distinguish host range, attachment, receptor binding and productive infection.
- Develop new models only with curated labelled datasets.
- Use cluster-aware splits to reduce homolog leakage and independent test sets.
- Evaluate precision, recall, F1, MCC, PR-AUC and calibration.

## Stage 6 — Production readiness
- Background jobs, file size limits and resource quotas.
- Temporary-file cleanup and job isolation.
- Database versioning, reproducible reports, tests and CI.
