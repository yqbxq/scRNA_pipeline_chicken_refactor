# ST09 Spatial SVG Evidence

ST09 builds a spatially variable gene evidence layer for downstream interpretation. It does not replace communication, regulation, enrichment, trajectory, or velocity methods.

## Methods

- `09a_spatialde2_svg.R` runs the H5AD-first SpatialDE2 sidecar. If SpatialDE2 is unavailable, it writes `skipped_no_spatialde2` unless `SPATIALDE2_ALLOW_PROXY=yes`.
- `09b_sparkx_svg.R` is R-native for SPARK-X. It prefers spatial Seurat/RDS inputs and can use a H5AD adapter only for exploratory proxy output when `SPARKX_ALLOW_PROXY=yes`.
- `09c_svg_consensus_report.R` combines method outputs and writes global handoff tables.

Proxy outputs are always exploratory caution. They cannot satisfy I06/I07 PASS and cannot upgrade communication tiers.

## Core Outputs

ST09 consensus outputs live under:

`results/spatial/tables/09_svg/09c_consensus/`

Key files:

- `svg_consensus.tsv`
- `svg_evidence_tier.tsv`
- `svg_gene_sets.tsv`
- `svg_lr_support.tsv`
- `svg_regulation_support.tsv`
- `svg_enrichment_handoff.tsv`
- `svg_trajectory_support.tsv`
- `svg_question_gate_status.tsv`
- `svg_project_handoff_manifest.tsv`
- `svg_report.md`

## Gate Semantics

`I06_SVG` passes only when formal SVG evidence exists. `I07_SVG_by_stage` passes only when formal SVG evidence supports stage/condition comparison. Both gates are exploratory, never core.
