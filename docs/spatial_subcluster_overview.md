# Spatial Subcluster Overview

ST 04 implements the spatial "merge first, split later" drilldown. It reads `spatial_panorama_annotated.rds`, selects enabled `region_subset` rows from `config/spatial_object_layers.tsv`, rebuilds HVG/PCA/SNN clustering inside each selected region, annotates sub-regions, and projects the resulting `sub_region` labels back to `spatial_panorama_subannotated.rds`.

The template `region_gc` and `region_tc` rows are intentionally `enabled=no`. With no enabled rows, `04a_subcluster_build.R`, `04b_subcluster_annotate.R`, and `04c_subcluster_eda.R` write skipped/no-op manifests and let the pipeline continue to the review gate.

To enable a drilldown, set the row to `enabled=yes`, keep `layer_role=region_subset`, choose `selection_column=region`, and provide `selection_values` such as `GC_rich` or `TC_rich`. The selected subset must have at least `SPATIAL_SUBCLUSTER_MIN_SPOTS` spots, default `50`.

Sub-region labels are namespaced as `<parent_region>_<sub_label>`, for example `GC_rich_inner_ring`. Small SNN clusters are not dropped; labels below `SPATIAL_SUBCLUSTER_SMALL_CLUSTER_FRAC` get a `_merged_small` suffix so downstream grouping does not silently lose spots.

ST 04 always uses the `b1_seurat_snn` backend for these small region subsets. BayesSpace, SpaGCN, and STAGATE stay available for whole-panorama clustering but are intentionally ignored for subclustering.
