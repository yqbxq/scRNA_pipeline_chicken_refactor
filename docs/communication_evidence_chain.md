# Communication Evidence Chain

ST08 communication tiers are still determined by scRNA communication evidence, deconvolution support, colocalization, neighborhood consistency, and COMMOT/LR spatial support.

ST09 SVG adds gene-level spatial support through `svg_lr_support.tsv`:

- ligand SVG tier
- receptor SVG tier
- receiver target overlap summary
- `svg_spatial_support_level`
- `svg_support_reason`

SVG support can enrich `reason` and supporting evidence fields. It must not convert `spatial_hypothesis` into `spatial_primary` by itself.

Interpretation examples:

- `spatial_primary` plus strong SVG support: core communication with SVG gene-level spatial support.
- `spatial_hypothesis` plus strong SVG support: hypothesis with SVG-compatible spatial expression pattern.
