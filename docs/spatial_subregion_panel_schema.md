# Spatial Subregion Panel Schema

Sub-region panels use the same active TSV schema as `spatial_region_panel.tsv`:

```text
layer_id	celltype	gene	evidence_source	tissue	panel_name	evidence_note	confidence_ceiling
```

Use one active file per enabled region subset, for example `config/marker_panels/region_gc.tsv` for the `region_gc` layer. The installed `.template` files are examples only; copy or generate a `.tsv` file before enabling the layer.

Required columns are `layer_id`, `celltype`, `gene`, and `evidence_source`. `celltype` is the sub-region label before prefixing, such as `inner_ring`; ST 04 writes the projected label as `GC_rich_inner_ring`.

Rows with `layer_id` matching the enabled layer are used. `layer_id` values `*`, `all`, and `ALL` are also accepted for shared rows.
