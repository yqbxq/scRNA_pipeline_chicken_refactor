# Communication Report Schema

`workflow/05single_script/07e_communication_eda.R` reads the 07d consensus
outputs and writes the review-facing communication report.

## Inputs

| Input | Meaning |
|---|---|
| `results/tables/communication/consensus/method_consensus.tsv` | Single authoritative ligand-receptor axis table from 07d. |
| `results/tables/communication/consensus/evidence_tiers.tsv` | Fallback compact table if the full consensus table is absent. |
| `metadata/communication_pairs.tsv` | Pair-level report expectations and evidence-tier defaults. |

## Outputs

| Output key | Path | Row semantics |
|---|---|---|
| `report_md` | `reports/eda/communication/report.md` | Markdown report with evidence disclaimer, three new panels, legacy sections, pair requirements, and filtered-axis preview. |
| `report_html` | `reports/eda/communication/report.html` | Simple HTML companion when `COMMUNICATION_REPORT_HTML=yes`. |
| `communication_eda_panels_tsv` | `results/tables/communication/consensus/communication_eda_panels.tsv` | One row per rendered report panel. |
| `communication_eda_filtered_tsv` | `results/tables/communication/consensus/communication_eda_filtered_axes.tsv` | Consensus axes after the configured report-tier filter. |

## Panels

| Panel | Purpose |
|---|---|
| Panel 1: Evidence Tier Distribution | Counts and plots `primary`, `exploratory`, `candidate`, and `blocked` axes. |
| Panel 2: Method Agreement | Summarizes CellChat, LIANA+, NicheNet/MultiNicheNet, and COMMOT hit patterns. |
| Panel 3: Downstream Target Chain | Reports NicheNet/MultiNicheNet downstream target support when available. |

The report keeps legacy sections for scDesign3 gate status, fallback summaries,
and missing gene-program diagnostics. These sections are summaries only and do
not override the 07d evidence tier.

## Filtering

`COMMUNICATION_REPORT_FILTER_TIER=auto` uses
`communication_pairs.tsv:evidence_tier_required` per pair, defaulting to
`exploratory` for legacy sheets. Explicit values `primary`, `exploratory`,
`candidate`, or `blocked` apply the same minimum tier to every axis.
