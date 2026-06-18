# Protein Annotator

Protein Annotator is a command-line workflow for protein functional annotation using multiple evidence sources.

It combines:

1. DIAMOND best-hit search
2. eggNOG-mapper orthology-based annotation
3. Pfam / HMMER domain annotation
4. Integrated annotation table generation
5. High-confidence annotation filtering
6. Optional source FASTA product annotation extraction

This workflow is designed for non-model organisms where gene annotation is incomplete or needs multi-evidence validation.

## Installation

Create a conda environment:

    conda env create -f environment.yml
    conda activate protein_anno

Install Protein Annotator:

    git clone https://github.com/YOUR_USERNAME/protein-annotator.git
    cd protein-annotator
    python -m pip install -e .

Check installation:

    protein-annotate --help
    protein-annotate-simple --help
    protein-add-source-annotation --help

## Prepare databases

Databases are not included in this repository because they are large.

Prepare default databases:

    protein-annotate prepare-db --db-dir /path/to/protein_annotator_databases

Load database paths:

    source /path/to/protein_annotator_databases/protein_annotator_db_paths.sh

Check the environment:

    protein-annotate doctor \
      --diamond-db "$PROTEIN_ANNOTATOR_DIAMOND_DB" \
      --pfam-hmm "$PROTEIN_ANNOTATOR_PFAM_HMM" \
      --eggnog-data-dir "$PROTEIN_ANNOTATOR_EGGNOG_DATA_DIR"

## Recommended simple usage

For most users, run:

    protein-annotate-simple \
      --protein-fasta proteins.faa \
      --outdir annotation_out \
      --prefix sample \
      --diamond-db "$PROTEIN_ANNOTATOR_DIAMOND_DB" \
      --pfam-hmm "$PROTEIN_ANNOTATOR_PFAM_HMM" \
      --eggnog-data-dir "$PROTEIN_ANNOTATOR_EGGNOG_DATA_DIR" \
      --threads 16

This command produces two main output files:

    annotation_out/sample.annotation.tsv
    annotation_out/sample.summary.tsv

The final annotation table contains:

- original source FASTA product annotation, if available
- DIAMOND best hit
- eggNOG annotation
- GO / KEGG / COG information
- Pfam domain annotation
- integrated evidence level

To keep intermediate files for debugging, add:

    --keep-work

## Detailed workflow usage

Advanced users can run the full workflow and keep all intermediate files:

    protein-annotate run \
      --protein-fasta proteins.faa \
      --outdir full_annotation_out \
      --prefix sample \
      --diamond-db "$PROTEIN_ANNOTATOR_DIAMOND_DB" \
      --pfam-hmm "$PROTEIN_ANNOTATOR_PFAM_HMM" \
      --eggnog-data-dir "$PROTEIN_ANNOTATOR_EGGNOG_DATA_DIR" \
      --threads 16

The detailed output contains:

    full_annotation_out/01_diamond/
    full_annotation_out/02_eggnog/
    full_annotation_out/03_pfam/
    full_annotation_out/04_merged/
    full_annotation_out/05_high_confidence/

## Add source FASTA annotation only

If users already have an annotation table and only want to add original FASTA header product names:

    protein-add-source-annotation \
      --protein-fasta proteins.faa \
      --annotation-tsv sample.final_merged_annotation.tsv \
      --out-tsv sample.final_merged_annotation.with_source.tsv

This adds:

- source_protein_id
- source_product
- source_species
- source_annotation_informative
- source_header

This step does not compare annotations. It only preserves original product names from the input FASTA.

## Example: amphioxus annotation

    protein-annotate-simple \
      --protein-fasta Branchiostoma_floridae_GCF_000003815.2.protein.faa \
      --outdir amphioxus_annotation_out \
      --prefix amphioxus \
      --diamond-db "$PROTEIN_ANNOTATOR_DIAMOND_DB" \
      --pfam-hmm "$PROTEIN_ANNOTATOR_PFAM_HMM" \
      --eggnog-data-dir "$PROTEIN_ANNOTATOR_EGGNOG_DATA_DIR" \
      --threads 16

Main outputs:

    amphioxus_annotation_out/amphioxus.annotation.tsv
    amphioxus_annotation_out/amphioxus.summary.tsv

## Output columns

The final annotation table includes source annotation and multi-evidence functional annotation columns, such as:

- gene_id
- source_product
- source_species
- source_annotation_informative
- diamond_best_hit
- diamond_pident
- diamond_evalue
- diamond_bitscore
- Preferred_name
- Description
- GOs
- KEGG_ko
- KEGG_Pathway
- COG_category
- PFAMs
- pfam_domains_hmmscan
- annotation_evidence_level

## Citation

If you use this workflow, please cite the underlying tools and databases:

- DIAMOND
- eggNOG-mapper
- eggNOG
- Pfam
- HMMER

## Notes

This repository only contains workflow code. Large databases and annotation outputs are not included.
