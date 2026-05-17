# Spatial Region Annotation

Module 03 assigns cluster-level spatial regions with two operational evidence signals:

- marker-panel hits among cluster markers from `FindAllMarkers`
- region module scores calculated from active TSV files in `config/marker_panels/`

The marker-panel signal combines exploratory cluster markers with curated panel membership, so it covers the planned marker and panel evidence without treating them as independent votes.

The default panel template uses `GC_rich`, `TC_rich`, `stroma`, `vasculature`, and the reserved uncertain bucket. At runtime, `uncertain` is represented as `mixed_or_uncertain` in `panorama$region`.

Active marker panels must end in `.tsv`; `.tsv.template` files are copied for user guidance but are not loaded. Required columns are `layer_id`, `celltype`, `gene`, and `evidence_source`.

Manual cluster overrides can be written to `config/spatial_region_annotation_override.tsv` with:

```text
cluster_id	override_region	notes
```
