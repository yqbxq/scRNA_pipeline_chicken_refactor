# Spatial Integration Signals

Module 03 evaluates `none`, `harmony`, and `cca` integration modes before clustering.

The default recommendation is `none`. In the initial SYF/F5 use case, section and stage can be confounded, so automatic batch correction may remove biological structure. Harmony and CCA reductions are still reported when their dependencies are available, and the reviewer can set `integration_mode` on the `panorama_st` row in `config/spatial_object_layers.tsv`.

`lisi_summary.tsv` reports section mixing for each reduction. When the `lisi` R package is unavailable, the pipeline writes an inverse-Simpson section-diversity proxy with `lisi_implementation=inverse_simpson_proxy` so the report remains reviewable.

`biological_consistency.tsv` checks whether high marker-module-score spots for each region concentrate in a coarse clustering on the same reduction. Treat this as a diagnostic, not a hard biological label.
