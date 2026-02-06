"""
DEG 분석 및 Cell Type 예측

Wilcoxon 기반 DEG + Enrichr 기반 cell type 예측
"""
import os

import scanpy as sc
import gseapy
import matplotlib.pyplot as plt

from config.defaults import DEG
from config.species import get_species_config
from plotting.utils import save_figure, calc_figsize, plot_gene_cluster_panels


def deg_analysis(
    adata,
    save_path: str,
    species: str,
    groupby: str = "leiden",
    top_n: int = None,
    padj_threshold: float = None,
) -> sc.AnnData:
    """
    DEG 분석 및 Cell Type 예측

    Parameters
    ----------
    adata : AnnData
    save_path : str
    species : str
    groupby : str
    top_n : int, optional
    padj_threshold : float, optional

    Returns
    -------
    AnnData
    """
    os.makedirs(save_path, exist_ok=True)

    if top_n is None:
        top_n = DEG["top_n"]
    if padj_threshold is None:
        padj_threshold = DEG["padj_threshold"]

    config = get_species_config(species)

    # 1. 성염색체 유전자 제거
    print("[DEG] Removing sex chromosome genes...")
    try:
        annot = sc.queries.biomart_annotations(
            config["biomart_org"],
            ["ensembl_gene_id", "external_gene_name", "chromosome_name"]
        ).set_index("external_gene_name")

        chrY = adata.var_names.intersection(annot.index[annot.chromosome_name == "Y"])
        chrX = adata.var_names.intersection(annot.index[annot.chromosome_name == "X"])
        sex_genes = chrY.union(chrX)
        adata = adata[:, [g for g in adata.var_names if g not in sex_genes]]
        print(f"[DEG] Removed {len(sex_genes)} sex chromosome genes.")
    except Exception as e:
        print(f"[DEG] Warning: BioMart error: {e}")
        print("[DEG] Skipping sex chromosome gene removal.")

    # 2. DEG 분석
    print("[DEG] Running DEG analysis...")
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=DEG["method"], key_added="wilcoxon")
    deg_df = sc.get.rank_genes_groups_df(adata, group=None, key="wilcoxon")
    deg_df = deg_df[deg_df["pvals_adj"] < padj_threshold]

    if deg_df.empty:
        print("[DEG] No significant DEGs found.")
        return adata

    top_genes = (
        deg_df.groupby("group")
        .apply(lambda x: x.nlargest(top_n, "logfoldchanges"))["names"]
        .unique()
        .tolist()
    )
    print(f"[DEG] Selected {len(top_genes)} top genes.")

    if groupby not in adata.obs:
        raise KeyError(f"{groupby} not found in adata.obs")
    n_clusters = adata.obs[groupby].nunique()

    # 3. Enrichr로 cell type 예측
    print("[DEG] Predicting cell types via Enrichr...")
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", key_added="gsea_wilcoxon")
    pred = {}
    organism = "human" if species == "human" else "mouse"

    for cl in adata.obs[groupby].astype("category").cat.categories:
        glist = (
            sc.get.rank_genes_groups_df(adata, group=cl, key="gsea_wilcoxon")["names"]
            .dropna()
            .astype(str)
            .str.strip()
            .tolist()
        )
        if len(glist) == 0:
            pred[cl] = "Unassigned"
            continue

        try:
            enr = gseapy.enrichr(
                gene_list=glist[:300],
                gene_sets=config["enrichr_libs"],
                organism=organism,
                background=adata.var_names,
                cutoff=1.0,
                no_plot=True,
                outdir=None,
                verbose=False,
            )
            res = getattr(enr, "results", None)
            if res is None or res.shape[0] == 0:
                pred[cl] = "Unassigned"
                continue
            res = res.sort_values("P-value")
            pred[cl] = res["Term"].iloc[0]
        except Exception:
            pred[cl] = "Unassigned"

    adata.obs["cell_type"] = adata.obs[groupby].map(pred)
    print("[DEG] Cell type prediction complete.")

    # 4. 시각화
    # 4-1. rank_genes_groups
    fig_w, fig_h = calc_figsize(
        n_genes=top_n,
        n_clusters=n_clusters,
        base_size=(10.0, 5),
        gene_width=0.8,
        cluster_height=0.6,
        min_size=(6.0, 4.0),
        max_size=(20.0, 16.0),
    )
    sc.pl.rank_genes_groups(adata, n_genes=top_n, key="wilcoxon", sharey=False, show=False, figsize=(fig_w, fig_h))
    fig = plt.gcf()
    save_figure(fig, save_path, "rank_genes_groups.png")

    # 4-2. dotplot / matrixplot / stacked_violin / heatmap
    plot_gene_cluster_panels(sc.pl.dotplot, adata, top_genes, groupby, save_path, "dotplot", n_clusters)
    plot_gene_cluster_panels(sc.pl.matrixplot, adata, top_genes, groupby, save_path, "matrixplot", n_clusters)
    plot_gene_cluster_panels(sc.pl.stacked_violin, adata, top_genes, groupby, save_path, "stacked_violin", n_clusters)
    plot_gene_cluster_panels(sc.pl.heatmap, adata, top_genes, groupby, save_path, "heatmap", n_clusters)

    # 4-3. UMAP with cell_type
    if "cell_type" in adata.obs:
        n_celltypes = adata.obs["cell_type"].nunique(dropna=True)
        base_side = 20.0
        extra_for_legend = min(5.0, 0.25 * max(n_celltypes - 8, 0))
        fig_w_umap = base_side + extra_for_legend
        fig_h_umap = base_side

        with plt.rc_context({"figure.figsize": (fig_w_umap, fig_h_umap)}):
            sc.pl.umap(adata, color="cell_type", show=False, legend_loc="right margin")
            fig = plt.gcf()
            ax = plt.gca()

            H = fig.get_figheight()
            W = fig.get_figwidth()
            height_frac = 0.9
            width_frac = height_frac * (H / W)
            ax.set_position([0.08, 0.07, width_frac, height_frac])

            save_figure(fig, save_path, "predicted_cell_types.png")

    print(f"[DEG] Complete. Results: {save_path}")
    return adata
