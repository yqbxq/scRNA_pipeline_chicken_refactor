# Spatial Replicate Gate

Module 05 uses a replicate gate only for formal sample-level statistics:

| Path | Gate |
| --- | --- |
| Region marker discovery | No replicate gate; exploratory spot-level marker discovery. |
| Region pseudobulk DE | Requires at least `min_biological_replicates` samples in both comparison groups. |
| Region composition propeller | Requires the same biological replicate gate. |
| Spot-level spatial DE | No replicate gate; exploratory Wilcoxon fallback. |

When the gate fails, pseudobulk rows write `status=skipped_replicate_gate` and point `exploratory_fallback` at the matching `05_spatial_de.R` result path. Composition still writes descriptive proportions and stack-bar figures.

To activate formal paths in a future expanded cohort, add enough samples per condition and set `min_biological_replicates` appropriately in `metadata/comparisons.tsv`. The module code does not need to change.
