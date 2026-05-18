# Spatial Pipeline Overview

The main ST workflow currently runs:

| stage | manifest | purpose |
| --- | --- | --- |
| `01` | `spatial_01_build_objects` | Build per-section spatial Seurat objects. |
| `01a` / `01b` / `01c` | pre/post QC manifests | Review, filter, and re-check spatial spots. |
| `02` / `02a` / `02b` | normalization, integration, clustering manifests | Compare normalization, integration, and clustering choices. |
| `03` / `03a` | `spatial_03_region_annotation`, `spatial_03a_region_annotation_eda` | Assign whole-panorama regions and review diagnostics. |
| `04a` / `04b` / `04c` | `spatial_04a_subcluster_build`, `spatial_04b_subcluster_annotate`, `spatial_04c_subcluster_eda` | Optional region-subset subclustering and sub-region review. |
| `05*` | `spatial_05*` | Reserved for marker discovery, pseudobulk, composition, and spatial DE in the next ST module. |

The `spatial_region_annotation` gate is held twice: once after `03a` for region review, then reset and held again after `04c` for sub-region review. ST 04 does not add a new gate name.
