# Spatial Regulation Pipeline

Spatial regulation sidecars consume H5AD-first spatial input and may also record ST09 SVG gene set support.

`08a_spatial_decoupler.R` passes `svg_gene_sets.tsv` to `spatial_decoupler.py` when the ST09 consensus output exists. The Python sidecar records this source in its manifest and H5AD `.uns`, but SVG evidence does not make decoupleR or pySCENIC results pass automatically.

`svg_regulation_support.tsv` links regulators and targets to SVG tiers. The intended interpretation is:

- high-confidence or condition-specific SVG targets can support a regulon spatial-expression note;
- single-method SVG targets are supplemental;
- proxy-only SVG targets require caution.
