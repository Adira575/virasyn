# VIRASYN Bacteria-Focused Phage Intelligence Architecture

## Product objective

Build a research platform for discovering bacteriophage receptor-binding proteins (RBPs), mapping candidate bacterial receptors, and ranking bacterial host-range hypotheses with traceable evidence.

The initial scope is bacteriophages and bacterial hosts. A host association is not equivalent to adsorption, receptor binding, or productive infection; these endpoints are stored and reported separately.

## Logical architecture

```text
Input and provenance
  ├── Phage genome FASTA (single or multi-record)
  ├── RBP protein reference FASTA + versioned metadata
  ├── Phage–bacterium evidence CSV
  ├── RBP–receptor annotation CSV
  └── Optional PDB/mmCIF coordinate files
          |
          v
Genome and protein processing
  ├── Six-frame ORF screening (preliminary gene calling)
  ├── Candidate protein translation
  └── Reference-library validation and inventory
          |
          v
RBP discovery
  ├── Local BLASTP homology search
  ├── Identity, query coverage, subject coverage and E-value
  └── RBP family / cluster metadata join
          |
          v
Evidence graph
  ├── RBP reference ID -> receptor annotation
  ├── Phage ID -> bacterial taxon (NCBI Taxonomy ID)
  ├── Evidence type and evidence level
  └── Source, citation, version and provenance
          |
          v
Interpretation and outputs
  ├── Per-phage bacterial host evidence ranking
  ├── RBP–receptor candidate table
  ├── Structural-file inventory
  └── CSV and JSON research reports
```

## Data contracts

### RBP reference FASTA

- FASTA record ID is the stable key used by BLASTP and metadata joins.
- Protein sequence should be curated and linked to an RBP family/cluster where known.
- Store database name, release/version, retrieval date, licence and citation outside the FASTA in a maintained manifest.

### RBP reference metadata CSV

Required:
- `subject_id`: exact FASTA record identifier.

Recommended:
- `cluster`, `RBP_name`, `strain`, `receptor`, `evidence`, `citation`, `source_version`, `retrieved_at`, `licence`.

### Host evidence CSV

Required:
- `phage_id`, `host_taxid`, `host_name`, `evidence_type`, `evidence_level`, `source`, `citation`.

Supported evidence levels:
- `experimental`
- `curated`
- `computational`
- `inferred`

Recommended evidence types:
- `productive_infection`
- `host_range`
- `adsorption`
- `receptor_binding`
- `crispr_spacer`
- `sequence_similarity`

### RBP–receptor annotation CSV

Required:
- `subject_id`, `receptor`, `receptor_class`, `host_taxid`, `host_name`, `evidence_type`, `citation`.

The join to RBP hits uses exact reference IDs, not fuzzy name matching.

## Current release boundary

Implemented:
- Six-frame ORF screening through the existing VIRASYN core.
- Local BLASTP against a user-provided RBP FASTA.
- Configurable discovery thresholds and translation table.
- Reference metadata, host evidence and receptor annotation imports.
- Per-phage host evidence aggregation with transparent evidence weighting.
- Optional PDB/mmCIF structure inventory.
- CSV templates and CSV/JSON exports.

Not yet implemented:
- Automated genome annotation with a phage-aware gene caller.
- Automated, broad-spectrum RBP family classifier trained on a curated non-redundant dataset.
- Protein language-model embeddings or a trained RBP–receptor interaction model.
- Receptor discovery from arbitrary bacterial genomes.
- Structural interface prediction, docking, binding-energy calculation or molecular dynamics.
- Phylogenetic and comparative-genomic visualization.
- Calibrated host-range prediction for unseen phages.

## Model-development plan

1. Curate a non-redundant RBP sequence set with experimentally supported labels and explicit negative/unknown labels.
2. Curate receptor annotations at molecule and strain level; distinguish receptor binding from whole-cell adsorption and productive infection.
3. Split evaluation by phage lineage, RBP cluster and bacterial taxon to reduce sequence leakage.
4. Establish baselines (BLAST/HMM profiles, k-mer models, protein embeddings) before training complex models.
5. Evaluate with precision-recall, AUROC where appropriate, calibration, taxonomic holdouts and per-clade performance.
6. Report abstentions and out-of-distribution warnings for poorly represented taxa.
7. Validate predictions against independent adsorption, receptor-knockout/complementation and infection assays.

## Scientific interpretation policy

- Never label a homology hit as experimentally confirmed RBP unless its source evidence supports that claim.
- Never label a host-evidence score as a probability.
- Keep receptor binding, adsorption, host range and productive infection as separate endpoints.
- Show source citations and provenance for every imported reference annotation.
- Missing evidence means unknown, not negative.
- AlphaFold confidence is structural-model confidence, not proof of receptor interaction.
