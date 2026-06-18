# Protein Annotator v0.3

Protein Annotator is a command-line workflow for protein functional annotation in non-model organisms.

It integrates multiple evidence sources:

```
DIAMOND best-hit search
eggNOG-mapper orthology-based annotation
Pfam / HMMER domain annotation
Source FASTA product annotation extraction
Multi-evidence annotation confidence scoring
Conservative functional name recommendation
Annotation conflict and ambiguity detection
```

This workflow is designed for proteomes from non-model organisms where gene annotation is incomplete, inconsistent, or difficult to interpret.

Protein Annotator does not only merge existing annotation results. In version 0.3, it introduces an interpretable multi-evidence annotation assessment framework to help users distinguish:

```
high-confidence gene-level annotations
family-level annotations
domain-level annotations
paralog-ambiguous annotations
annotations requiring manual review
weak or unsupported annotations
```

## Main features

Protein Annotator combines standard annotation tools with interpretable scoring and classification modules.

Core annotation workflow:

```
Protein FASTA
  -> DIAMOND best-hit search
  -> eggNOG-mapper annotation
  -> Pfam / HMMER domain annotation
  -> Integrated annotation table
  -> Multi-evidence confidence assessment
  -> High-confidence and review-needed annotation subsets
```

Version 0.3 adds the following assessment modules:

```
PACS
Protein Annotation Confidence Score.
A raw confidence score based on sequence similarity, orthology, domain evidence, functional annotation richness, evidence consistency and ambiguity penalties.

EGCS
Evidence Graph Consistency Score.
Different annotation sources are treated as evidence nodes, including source product, DIAMOND, eggNOG, Pfam/domain and functional terms. Their agreement is evaluated as an evidence graph.

PAI
Paralog Ambiguity Index.
Estimates whether available evidence supports an exact gene-level name or only a broader family-level annotation, especially for paralog-rich gene families.

DACS
Domain Architecture Consistency Status.
Evaluates whether detected protein domains support the predicted protein family or functional class.

CFNR
Conservative Functional Name Recommendation.
Provides a conservative recommended name and avoids over-specific paralog naming when the evidence only supports family-level annotation.

ACAD
Annotation Conflict and Ambiguity Detection.
Classifies annotations as consistent, domain-level consistent, family-level consistent, paralog ambiguous, potentially conflicting or unsupported.

WSCC
Weak-Supervised Confidence Calibration.
Uses informative source FASTA product annotations as weak labels to calibrate raw annotation confidence scores.
```

## Installation

Create a conda environment:

```
conda env create -f environment.yml
conda activate protein_anno
```

Install Protein Annotator:

```
git clone https://github.com/fungusL/protein-annotator-v0.3.git
cd protein-annotator-v0.3
python -m pip install -e .
```

Check installation:

```
protein-annotate --help
protein-annotate-simple --help
protein-add-source-annotation --help
```

## Prepare databases

Large databases are not included in this repository.

Prepare default databases:

```
protein-annotate prepare-db --db-dir /path/to/protein_annotator_databases
```

Load database paths:

```
source /path/to/protein_annotator_databases/protein_annotator_db_paths.sh
```

Check the environment:

```
protein-annotate doctor \
  --diamond-db "$PROTEIN_ANNOTATOR_DIAMOND_DB" \
  --pfam-hmm "$PROTEIN_ANNOTATOR_PFAM_HMM" \
  --eggnog-data-dir "$PROTEIN_ANNOTATOR_EGGNOG_DATA_DIR"
```

## Input

The main input is a protein FASTA file:

```
proteins.faa
```

The FASTA file can come from NCBI, genome annotation pipelines, transcriptome assemblies, BRAKER, AUGUSTUS, MAKER, TransDecoder or other protein prediction workflows.

If the FASTA header contains source product information, Protein Annotator can extract it automatically. For example:

```
>XP_035657343.1 semaphorin-5A-like, partial [Branchiostoma floridae]
```

Protein Annotator can extract:

```
source_protein_id
source_product
source_species
source_annotation_informative
source_header
```

## Recommended simple usage

For most users, run:

```
protein-annotate-simple \
  --protein-fasta proteins.faa \
  --outdir annotation_out \
  --prefix sample \
  --diamond-db "$PROTEIN_ANNOTATOR_DIAMOND_DB" \
  --pfam-hmm "$PROTEIN_ANNOTATOR_PFAM_HMM" \
  --eggnog-data-dir "$PROTEIN_ANNOTATOR_EGGNOG_DATA_DIR" \
  --threads 16
```

This command produces four main output files:

```
annotation_out/sample.annotation.tsv
annotation_out/sample.high_confidence.tsv
annotation_out/sample.review_needed.tsv
annotation_out/sample.summary.tsv
```

## Main outputs

### sample.annotation.tsv

This is the complete integrated annotation table.

It contains source FASTA product annotation, DIAMOND results, eggNOG-mapper results, Pfam/HMMER domain results and v0.3 multi-evidence assessment columns.

Important columns include:

```
gene_id
source_product
source_species
source_annotation_informative
diamond_best_hit
diamond_pident
diamond_evalue
diamond_bitscore
Preferred_name
Description
GOs
KEGG_ko
KEGG_Pathway
COG_category
PFAMs
pfam_domains_hmmscan
annotation_evidence_level
```

Version 0.3 also adds:

```
PACS_annotation_confidence_score
PACS_annotation_confidence_level

EGCS_evidence_graph_score
EGCS_evidence_graph_density
EGCS_supporting_edges
EGCS_conflicting_edges
EGCS_graph_status

PAI_paralog_ambiguity_index
PAI_paralog_resolution_level
PAI_detected_families
PAI_reason

DACS_domain_architecture_status

CFNR_recommended_name
CFNR_recommended_name_level
CFNR_recommended_name_reason

ACAD_annotation_consistency_status

WSCC_pseudo_label
WSCC_calibration_bin
WSCC_bin_positive_rate
WSCC_calibrated_confidence_score
WSCC_calibrated_confidence_level

manual_review_priority
```

### sample.high_confidence.tsv

This table contains high-confidence annotations.

In v0.3, high-confidence annotations are selected using multiple criteria:

```
DIAMOND + eggNOG + Pfam evidence are all present
WSCC calibrated confidence level is high or very high
Paralog ambiguity is not high
Manual review priority is not high
```

This file is useful when users want a conservative set of reliable gene annotations.

### sample.review_needed.tsv

This table contains annotations that should be manually checked.

Typical cases include:

```
potential annotation conflict
paralog ambiguity
family-level annotation only
source product more specific than supporting evidence
low-confidence but partially informative annotation
```

This file is useful for identifying genes that should not be blindly assigned to a specific gene-level name.

### sample.summary.tsv

This table summarizes annotation evidence levels and v0.3 assessment results.

It includes counts for:

```
annotation evidence levels
PACS confidence levels
EGCS graph status
PAI paralog resolution levels
DACS domain architecture status
ACAD consistency status
WSCC calibrated confidence levels
CFNR recommended name levels
manual review priority
```

## Interpretation of v0.3 scores

### PACS annotation confidence level

```
very_high_confidence
high_confidence
medium_confidence
low_confidence
weak_or_no_annotation
```

PACS is the raw interpretable confidence score before weak-supervised calibration.

### EGCS evidence graph status

```
high_graph_agreement
moderate_graph_agreement
weak_graph_agreement
low_graph_agreement
conflicting_evidence
insufficient_evidence
```

EGCS evaluates whether different annotation sources support similar functional interpretations.

### PAI paralog resolution level

```
clear_gene_level
clear_or_not_applicable
likely_gene_or_family_level
family_level_only
source_specific_evidence_family_level
evidence_specific_source_unspecific
unresolved
```

PAI is useful for paralog-rich families where exact gene-level naming may be unreliable.

### CFNR recommended name level

```
gene_like
family_like
orthology_like
functional_class
domain_containing
conflict_candidate
uncharacterized
```

CFNR provides conservative annotation recommendations. For example, if evidence supports a protein family but not a specific paralog, the recommended name will be family-level rather than over-specific.

### Manual review priority

```
low
medium
high
```

Genes with high or medium manual review priority should be checked before being used as key markers or functional conclusions.

## Keep intermediate files

By default, protein-annotate-simple removes intermediate files and keeps only the main output tables.

To keep intermediate files, add:

```
--keep-work
```

Example:

```
protein-annotate-simple \
  --protein-fasta proteins.faa \
  --outdir annotation_out \
  --prefix sample \
  --diamond-db "$PROTEIN_ANNOTATOR_DIAMOND_DB" \
  --pfam-hmm "$PROTEIN_ANNOTATOR_PFAM_HMM" \
  --eggnog-data-dir "$PROTEIN_ANNOTATOR_EGGNOG_DATA_DIR" \
  --threads 16 \
  --keep-work
```

This will keep:

```
annotation_out/_work/01_diamond/
annotation_out/_work/02_eggnog/
annotation_out/_work/03_pfam/
annotation_out/_work/04_merged/
annotation_out/_work/05_high_confidence/
```

## Detailed workflow usage

Advanced users can run the full workflow and keep all intermediate files:

```
protein-annotate run \
  --protein-fasta proteins.faa \
  --outdir full_annotation_out \
  --prefix sample \
  --diamond-db "$PROTEIN_ANNOTATOR_DIAMOND_DB" \
  --pfam-hmm "$PROTEIN_ANNOTATOR_PFAM_HMM" \
  --eggnog-data-dir "$PROTEIN_ANNOTATOR_EGGNOG_DATA_DIR" \
  --threads 16
```

The detailed output contains:

```
full_annotation_out/01_diamond/
full_annotation_out/02_eggnog/
full_annotation_out/03_pfam/
full_annotation_out/04_merged/
full_annotation_out/05_high_confidence/
```

## Add source FASTA annotation only

If users already have an annotation table and only want to add original FASTA header product names:

```
protein-add-source-annotation \
  --protein-fasta proteins.faa \
  --annotation-tsv sample.final_merged_annotation.tsv \
  --out-tsv sample.final_merged_annotation.with_source.tsv
```

This adds:

```
source_protein_id
source_product
source_species
source_annotation_informative
source_header
```

This step does not compare annotations. It only preserves original product names from the input FASTA.

## Example: amphioxus annotation

```
protein-annotate-simple \
  --protein-fasta Branchiostoma_floridae_GCF_000003815.2.protein.faa \
  --outdir amphioxus_annotation_out \
  --prefix amphioxus \
  --diamond-db "$PROTEIN_ANNOTATOR_DIAMOND_DB" \
  --pfam-hmm "$PROTEIN_ANNOTATOR_PFAM_HMM" \
  --eggnog-data-dir "$PROTEIN_ANNOTATOR_EGGNOG_DATA_DIR" \
  --threads 16
```

Main outputs:

```
amphioxus_annotation_out/amphioxus.annotation.tsv
amphioxus_annotation_out/amphioxus.high_confidence.tsv
amphioxus_annotation_out/amphioxus.review_needed.tsv
amphioxus_annotation_out/amphioxus.summary.tsv
```

## Citation

If you use this workflow, please cite the underlying tools and databases:

```
DIAMOND
eggNOG-mapper
eggNOG
Pfam
HMMER
```

If you use the v0.3 confidence assessment framework, please also cite this repository.

## Notes

This repository only contains workflow code.

Large databases, input proteomes and annotation outputs are not included.

The v0.3 scoring framework is designed to be interpretable and conservative. It is not intended to replace manual curation, phylogenetic analysis or experimental validation when exact paralog identity is critical.
