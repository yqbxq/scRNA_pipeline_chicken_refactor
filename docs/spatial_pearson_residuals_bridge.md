# Spatial Pearson Residual Bridge

`workflow/04python/normalize_pearson_residuals.py` is the Python sidecar used by `02_normalize_spatial.R` for `m4_pearson_residuals`.

The bridge is invoked through `system2()` rather than reticulate. The R stage writes a gene-by-spot Matrix Market counts matrix plus feature and barcode files. Python builds an AnnData object, runs `scanpy.experimental.pp.normalize_pearson_residuals`, and writes a gene-by-spot Matrix Market residual matrix back for Seurat.

The script also accepts `--input` and `--output` h5ad paths for direct bridge testing. If the Python environment or scanpy call fails, ST 02 records `failed_py_bridge` for m4 while keeping the other normalization methods and selected default available.
