# Metadata Directory

This directory is reserved for the project metadata architecture.

## Source Of Truth

`analysis_questions.tsv` is the planned human-maintained source of truth for
analysis intent. Generated module input tables should be derived from it.

## Planned Tables

| File | Role |
|---|---|
| `samples.tsv` | scRNA sample design table |
| `sections.tsv` | Future spatial transcriptomics section table |
| `analysis_questions.tsv` | Tier 1 question schema |
| `comparisons.tsv` | Generated DEG/composition inputs |
| `communication_pairs.tsv` | Generated communication inputs |
| `trajectory_pairs.tsv` | Generated trajectory inputs |
| `scenic_targets.tsv` | Generated regulation inputs |
| `enrichment_targets.tsv` | Generated enrichment inputs |
| `gene_program_targets.tsv` | Generated 05 gene-program availability contract |
| `deconv_pairs.tsv` | Future generated ST deconvolution inputs |
| `spatial_pairs.tsv` | Future generated ST spatial inputs |

M1-M5 should replace the placeholder tables with generated, validated content.
