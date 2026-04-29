import glob
import os
import re

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scvelo as scv
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D


VELOCITY_INPUT_DIR = os.environ["VELOCITY_INPUT_DIR"]
VELOCITY_LOOM_DIR = os.environ["VELOCITY_LOOM_DIR"]
VELOCITY_OUTPUT_DIR = os.environ["VELOCITY_OUTPUT_DIR"]
SCVELO_THREADS = int(os.environ.get("SCVELO_THREADS", "16"))

UMAP_CSV = os.path.join(VELOCITY_INPUT_DIR, "velocity_umap.csv")
META_CSV = os.path.join(VELOCITY_INPUT_DIR, "velocity_metadata.csv")
FIGURE_DIR = os.path.join(VELOCITY_OUTPUT_DIR, "figures")

CELLTYPE_COLOR_POOL = [
    "#79BB7D",
    "#C4A9D6",
    "#FDBE6F",
    "#2F6CB3",
    "#E30073",
    "#A54F37",
    "#69706A",
    "#0E9B84",
    "#D95F02",
    "#6B63B7",
]

CLUSTER_COLORS = {
    "1": "#74BA73",
    "2": "#C8C6E4",
    "3": "#FDBE6F",
    "4": "#2E63AF",
    "5": "#E30073",
    "6": "#A54F37",
    "7": "#69706A",
    "8": "#0E9B84",
    "9": "#D95F02",
    "10": "#6B63B7",
    "11": "#D90445",
    "12": "#1893D1",
    "13": "#FDBA12",
    "14": "#4C2A7A",
    "15": "#D6457A",
}

os.makedirs(VELOCITY_OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params("scvelo")

RAINBOW_CMAP = LinearSegmentedColormap.from_list(
    "blue_to_red_rainbow",
    [
        "#2C7BB6", "#00A6CA", "#00CCBC", "#90EB9D",
        "#FFFF8C", "#F9D057", "#F29E2E", "#E76818", "#D7191C",
    ],
)


def normalize_obs_name(obs_name, fallback_sample=None):
    obs_name = str(obs_name)
    if ":" in obs_name:
        sample_id, barcode = obs_name.split(":", 1)
    else:
        sample_id, barcode = fallback_sample, obs_name

    barcode = barcode.strip()
    if barcode.endswith("x"):
        barcode = f"{barcode[:-1]}-1"
    elif not re.search(r"-\\d+$", barcode):
        barcode = f"{barcode}-1"

    if sample_id is None:
        raise ValueError(f"Cannot determine sample id for barcode: {obs_name}")
    return f"{sample_id}:{barcode}"


def assign_discrete_colors(categories, color_pool=None):
    clean_categories = [str(x) for x in categories if pd.notna(x) and str(x)]
    if color_pool is None:
        color_pool = CELLTYPE_COLOR_POOL

    if len(clean_categories) > len(color_pool):
        color_pool = color_pool + [plt.cm.tab20(i) for i in range(len(clean_categories) - len(color_pool))]

    return {
        category: color_pool[idx]
        for idx, category in enumerate(clean_categories)
    }


def load_looms():
    loom_files = sorted(glob.glob(os.path.join(VELOCITY_LOOM_DIR, "*.loom")))
    if not loom_files:
        raise FileNotFoundError("未在 loom 目录中找到 .loom 文件")

    adata_list = []
    for loom_path in loom_files:
        sample_name = os.path.splitext(os.path.basename(loom_path))[0]
        adata_i = scv.read(loom_path, cache=True)
        adata_i.obs_names = [normalize_obs_name(barcode, fallback_sample=sample_name) for barcode in adata_i.obs_names]
        adata_i.obs["sample_id"] = sample_name
        adata_list.append(adata_i)

    return ad.concat(adata_list, join="outer", merge="same")


def integrate_metadata(adata):
    umap_df = pd.read_csv(UMAP_CSV, index_col=0)
    meta_df = pd.read_csv(META_CSV, index_col=0)

    common_cells = adata.obs_names.intersection(meta_df.index).intersection(umap_df.index)
    if len(common_cells) == 0:
        raise ValueError("loom 结果与导出的 metadata/UMAP 没有重叠细胞")

    adata = adata[common_cells].copy()
    adata.obsm["X_umap"] = umap_df.loc[common_cells, ["UMAP_1", "UMAP_2"]].to_numpy()

    for column in meta_df.columns:
        adata.obs[column] = meta_df.loc[common_cells, column].values

    if "cell_type" in adata.obs:
        observed_celltypes = [
            value for value in pd.unique(adata.obs["cell_type"].astype(str))
            if value and value.lower() != "nan"
        ]
        adata.obs["cell_type"] = pd.Categorical(
            adata.obs["cell_type"].astype(str),
            categories=observed_celltypes,
            ordered=True,
        )
        celltype_colors = assign_discrete_colors(observed_celltypes)
        adata.uns["cell_type_colors"] = [
            celltype_colors[cell_type]
            for cell_type in adata.obs["cell_type"].cat.categories
        ]

    cluster_col = None
    for candidate in ("seurat_clusters", "predicted_cluster"):
        if candidate in adata.obs:
            cluster_col = candidate
            break
    if cluster_col is not None:
        adata.obs["cluster_for_plot"] = pd.Categorical(
            adata.obs[cluster_col].astype(str),
            categories=list(CLUSTER_COLORS.keys()),
            ordered=True,
        )
        adata.uns["cluster_for_plot_colors"] = [
            CLUSTER_COLORS[cluster_id]
            for cluster_id in adata.obs["cluster_for_plot"].cat.categories
        ]

    return adata


def annotate_cluster_numbers(ax, adata):
    if "cluster_for_plot" not in adata.obs:
        return

    coords = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"], index=adata.obs_names)
    coords["cluster"] = adata.obs["cluster_for_plot"].astype(str).values
    centers = (
        coords.dropna()
        .groupby("cluster", as_index=False)[["UMAP_1", "UMAP_2"]]
        .median()
    )
    for _, row in centers.iterrows():
        ax.text(
            row["UMAP_1"],
            row["UMAP_2"],
            row["cluster"],
            fontsize=11,
            fontweight="bold",
            color="white",
            ha="center",
            va="center",
            bbox=dict(boxstyle="circle,pad=0.25", fc="black", ec="none", alpha=0.65),
        )


def plot_stream(adata, color_key, file_name, title, palette=None, annotate_clusters=False):
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color=color_key,
        palette=palette,
        legend_loc="right margin",
        title=title,
        frameon=False,
        ax=ax,
        show=False,
    )
    if color_key == "cell_type":
        categories = list(adata.obs["cell_type"].cat.categories)
        palette = palette or adata.uns.get("cell_type_colors", [])
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markerfacecolor=palette[idx],
                markeredgecolor="none",
                markersize=8,
                label=label,
            )
            for idx, label in enumerate(categories)
        ]
        ax.legend(
            handles=handles,
            labels=[h.get_label() for h in handles],
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            borderaxespad=0,
            handletextpad=0.5,
            labelspacing=0.7,
            prop={"family": "DejaVu Sans Mono", "size": 12},
        )
    if annotate_clusters:
        annotate_cluster_numbers(ax, adata)
    fig.savefig(os.path.join(FIGURE_DIR, file_name), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_feature(adata, color_key, file_name, title, cmap):
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    scv.pl.scatter(
        adata,
        basis="umap",
        color=color_key,
        color_map=cmap,
        size=30,
        title=title,
        frameon=False,
        ax=ax,
        show=False,
    )
    fig.savefig(os.path.join(FIGURE_DIR, file_name), dpi=300, bbox_inches="tight")
    plt.close(fig)


def compute_moments_compat(adata, n_pcs=30, n_neighbors=30):
    try:
        scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
        return
    except TypeError as exc:
        if "write_knn_indices" not in str(exc):
            raise

    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, method="gauss")
    scv.pp.moments(adata, n_pcs=None, n_neighbors=None)


def main():
    adata = load_looms()
    adata = integrate_metadata(adata)

    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    compute_moments_compat(adata, n_pcs=30, n_neighbors=30)
    scv.tl.recover_dynamics(adata, n_jobs=SCVELO_THREADS)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)
    scv.tl.latent_time(adata)
    scv.tl.velocity_confidence(adata)

    if "cluster_for_plot" in adata.obs:
        plot_stream(
            adata,
            "cluster_for_plot",
            "Figure_5A_velocity_stream_clusters.png",
            "RNA velocity by clusters",
            palette=[CLUSTER_COLORS[x] for x in CLUSTER_COLORS],
            annotate_clusters=True,
        )

    if "cell_type" in adata.obs:
        plot_stream(
            adata,
            "cell_type",
            "Figure_5B_velocity_stream_celltypes.png",
            "RNA velocity by cell types",
            palette=list(adata.uns.get("cell_type_colors", [])),
        )

    plot_feature(adata, "latent_time", "Figure_5C_latent_time.png", "latent time", RAINBOW_CMAP)
    plot_feature(
        adata,
        "velocity_length",
        "Figure_5D_velocity_length.png",
        "dynamical velocity length",
        RAINBOW_CMAP,
    )
    plot_feature(
        adata,
        "velocity_confidence",
        "Figure_5E_velocity_confidence.png",
        "dynamical velocity confidence",
        RAINBOW_CMAP,
    )

    adata.write(os.path.join(VELOCITY_OUTPUT_DIR, "scvelo_result.h5ad"))


if __name__ == "__main__":
    main()
