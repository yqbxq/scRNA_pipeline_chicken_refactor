# Spatial Region Annotation

Module 03 assigns cluster-level spatial regions with three evidence chains:

- cluster markers from `FindAllMarkers`
- curated marker-panel hits from active TSV files in `config/marker_panels/`
- region module scores calculated on the panorama object

The default panel template uses `GC_rich`, `TC_rich`, `stroma`, `vasculature`, and the reserved uncertain bucket. At runtime, `uncertain` is represented as `mixed_or_uncertain` in `panorama$region`.

Active marker panels must end in `.tsv`; `.tsv.template` files are copied for user guidance but are not loaded. Required columns are `layer_id`, `celltype`, `gene`, and `evidence_source`.

Manual cluster overrides can be written to `config/spatial_region_annotation_override.tsv` with:

```text
cluster_id	override_region	notes
```
