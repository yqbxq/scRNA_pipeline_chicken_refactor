# Marker Panel Directory

This directory is intentionally project-facing and biologically neutral.

- The workflow never hardcodes tissue-specific marker panels.
- Annotation meaning is injected only through TSV files placed here.
- If no panel TSV is present, annotation still runs and emits data-driven candidates with confidence `未定`.

## Files That Are Loaded

Only files ending with `.tsv` are treated as active marker panels.

- Safe examples or templates should use extensions such as `.tsv.template`.
- Multiple active TSV files can coexist.
- A panel file may target one layer or multiple layers.

## Required Columns

- `layer_id`
- `celltype`
- `gene`
- `evidence_source`

## Optional Columns

- `tissue`
- `panel_name`
- `evidence_note`
- `confidence_ceiling`

## Semantics

- `layer_id`: layer name from `object_layers.tsv`; use `*` for a global panel.
- `celltype`: user-defined annotation label. The workflow does not interpret its biology.
- `gene`: marker gene symbol or feature name expected in the object.
- `evidence_source`: free text such as `literature`, `curated`, `ortholog`, `candidate`.
- `tissue`: optional filter; use `*` or blank for all tissues.
- `confidence_ceiling`: optional cap on the final confidence. Allowed values: `确定`, `暂定`, `未定`.

## Example

See [example_panel.tsv.template](./example_panel.tsv.template). It uses fictional labels and fictional marker names on purpose.

## Notes

- One row represents one marker gene for one annotation label.
- The workflow first performs cluster-level marker discovery, then cross-checks against these panels, and finally runs Module Score as a validator.
- If a TSV is malformed, the workflow reports it and skips that file instead of silently inventing markers.
