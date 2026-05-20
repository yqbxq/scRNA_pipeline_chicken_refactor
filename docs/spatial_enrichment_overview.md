# Spatial Module 06 Enrichment Overview

Spatial module 06a/06b consumes ST05 marker, pseudobulk, and spot-level DE manifests and runs region-level over-representation analysis.

| Script | Manifest | Purpose |
| --- | --- | --- |
| `06a_region_go_enrichment.R` | `spatial_06a_region_go` | GO BP/CC/MF enrichment for region or sub-region gene sets. |
| `06b_region_kegg_enrichment.R` | `spatial_06b_region_kegg` | KEGG pathway enrichment for the same ST05-derived gene sets. |
| `06c_region_enrichment_eda.R` | `spatial_06c_region_enrichment_eda` | Combined GO/KEGG status summary and long region-pathway score matrix. |

The scripts use the ST00 ortholog cache to map chicken symbols to human symbols, then use `org.Hs.eg.db` Entrez IDs for `clusterProfiler`. If the upstream ST05 manifests, ortholog map, Entrez mapping, or enrichment packages are unavailable, the module writes explicit skip statuses instead of failing the whole ST extension stage.

Default knobs are `SPATIAL_ENRICH_TOP_N=200`, `SPATIAL_ENRICH_MIN_GENES=10`, and `SPATIAL_ENRICH_QVALUE=0.1`. Outputs land under `results/spatial/tables/spatial_06a_region_go`, `results/spatial/tables/spatial_06b_region_kegg`, and `reports/eda/spatial_enrichment/`.

The `spatial_enrichment` EDA gate is set after `06c_region_enrichment_eda.R` so reviewers can inspect GO/KEGG coverage before neighborhood analysis.
