import matplotlib
matplotlib.use("Agg")
import os
import math
import matplotlib.pyplot as plt
import scanpy as sc
import gseapy
import pandas as pd

def _biomart_org(species):
    if species == "human":
        return "hsapiens"
    if species == "mouse":
        return "mmusculus"
    raise ValueError(f"Invalid species: {species}")

def _enrichr_lib(species):
    if species == "human":
        return [
            "CellMarker_Augmented_2021",
            "PanglaoDB_Augmented_2021",
            "Descartes_Cell_Types_and_Tissue_2021",
            "Tabula_Sapiens",
            "Azimuth_Cell_Types_2021",
        ]
    if species == "mouse":
        return [
            "Tabula_Muris",
            "Mouse_Gene_Atlas",
            "PanglaoDB_Augmented_2021",
        ]
    raise ValueError(f"Invalid species: {species}")

def _chunk_list(lst, chunk_size):
    """리스트를 chunk_size 단위로 잘라 generator로 반환."""
    for i in range(0, len(lst), chunk_size):
        yield lst[i : i + chunk_size]

def _calc_figsize_genes_clusters(
    n_genes,
    n_clusters,
    base_size=(4.0, 3.0),
    gene_height=0.25,
    cluster_width=0.4,
    min_size=(4.0, 3.0),
    max_size=(30.0, 20.0),
):
    """
    gene 수와 cluster 수에 비례해서 figsize 계산.
    너무 작거나 너무 큰 경우를 막기 위해 min/max 를 둡니다.
    """
    width = base_size[0] + cluster_width * max(n_clusters - 1, 0)
    height = base_size[1] + gene_height * max(n_genes - 1, 0)

    width = max(min_size[0], min(max_size[0], width))
    height = max(min_size[1], min(max_size[1], height))
    return (width, height)

def _plot_gene_cluster_panels(
    plot_func,
    adata,
    genes,
    groupby,
    save_path,
    base_name,
    n_clusters,
    genes_per_fig=40,
):
    """
    dotplot / matrixplot / stacked_violin / heatmap 을 공통 로직으로 그림.
    gene 이 많으면 여러 장으로 나누고, 각 장마다 figsize를 동적으로 지정.
    """
    n_genes_total = len(genes)
    if n_genes_total == 0:
        return

    n_chunks = math.ceil(n_genes_total / genes_per_fig)

    for idx, genes_chunk in enumerate(_chunk_list(genes, genes_per_fig), start=1):
        n_genes_chunk = len(genes_chunk)
        fig_w, fig_h = _calc_figsize_genes_clusters(
            n_genes=n_genes_chunk,
            n_clusters=n_clusters,
            base_size=(4.0, 3.0),
            gene_height=0.25,
            cluster_width=0.4,
            min_size=(4.0, 3.0),
            max_size=(30.0, 20.0),
        )

        with plt.rc_context({"figure.figsize": (fig_w, fig_h)}):
            plot_func(
                adata,
                var_names=genes_chunk,
                groupby=groupby,
                use_raw=False,
                show=False,
            )
            plt.tight_layout()

            if n_chunks == 1:
                fname = f"{base_name}.png"
            else:
                fname = f"{base_name}_{idx}.png"

            plt.savefig(os.path.join(save_path, fname), dpi=300)
            plt.close()

def deg_analysis_with_sex_gene_filtering(
    adata,
    save_path,
    species,
    groupby="leiden",
    top_n=20,
    padj_threshold=0.05,
):
    os.makedirs(save_path, exist_ok=True)

    # ------------------------------------------------------------------ #
    # 1. 성별 염색체 유전자 제거
    # ------------------------------------------------------------------ #
    print("Removing sex chromosome genes...")
    org = _biomart_org(species)
    annot = sc.queries.biomart_annotations(
        org, ["ensembl_gene_id", "external_gene_name", "chromosome_name"]
    ).set_index("external_gene_name")

    chrY = adata.var_names.intersection(annot.index[annot.chromosome_name == "Y"])
    chrX = adata.var_names.intersection(annot.index[annot.chromosome_name == "X"])
    sex_genes = chrY.union(chrX)
    adata = adata[:, [g for g in adata.var_names if g not in sex_genes]]
    print(f"Removed {len(sex_genes)} sex chromosome genes.")

    # ------------------------------------------------------------------ #
    # 2. DEG 분석
    # ------------------------------------------------------------------ #
    print("Running DEG analysis...")
    sc.tl.rank_genes_groups(
        adata, groupby=groupby, method="wilcoxon", key_added="wilcoxon"
    )
    deg_df = sc.get.rank_genes_groups_df(adata, group=None, key="wilcoxon")
    deg_df = deg_df[deg_df["pvals_adj"] < padj_threshold]
    if deg_df.empty:
        print("No significant DEGs found.")
        return adata

    top_genes = (
        deg_df.groupby("group")
        .apply(lambda x: x.nlargest(top_n, "logfoldchanges"))["names"]
        .unique()
        .tolist()
    )
    print(f"Selected {len(top_genes)} top genes.")

    # cluster 수 계산
    if groupby not in adata.obs:
        raise KeyError(f"{groupby} not found in adata.obs")
    n_clusters = adata.obs[groupby].nunique()

    # 3. Enrichr로 cell type 예측
    print("Predicting cell types via Enrichr...")
    gene_sets = _enrichr_lib(species)
    organism = "human" if species == "human" else "mouse"

    sc.tl.rank_genes_groups(
        adata, groupby=groupby, method="wilcoxon", key_added="gsea_wilcoxon"
    )
    pred = {}

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
                gene_sets=gene_sets,
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
    print("Cell type prediction complete.")


    # 4. 플롯들 (모두 동적 크기 / 필요시 분할 저장)
    # 4-1. rank_genes_groups: top_n, n_clusters 기반 크기
    fig_w, fig_h = _calc_figsize_genes_clusters(
        n_genes=top_n,
        n_clusters=n_clusters,
        base_size=(4.0, 3.0),
        gene_height=0.35,
        cluster_width=0.8,
        min_size=(6.0, 4.0),
        max_size=(30.0, 20.0),
    )
    with plt.rc_context({"figure.figsize": (fig_w, fig_h)}):
        sc.pl.rank_genes_groups(
            adata, n_genes=top_n, key="wilcoxon", sharey=False, show=False
        )
        plt.tight_layout()
        plt.savefig(os.path.join(save_path, "rank_genes_groups.png"), dpi=300)
        plt.close()

    # 4-2. dotplot / matrixplot / stacked_violin / heatmap
    # gene 수가 많으면 자동으로 여러 장으로 나눔
    _plot_gene_cluster_panels(
        sc.pl.dotplot,
        adata,
        top_genes,
        groupby,
        save_path,
        base_name="dotplot",
        n_clusters=n_clusters,
        genes_per_fig=40,
    )

    _plot_gene_cluster_panels(
        sc.pl.matrixplot,
        adata,
        top_genes,
        groupby,
        save_path,
        base_name="matrixplot",
        n_clusters=n_clusters,
        genes_per_fig=40,
    )

    _plot_gene_cluster_panels(
        sc.pl.stacked_violin,
        adata,
        top_genes,
        groupby,
        save_path,
        base_name="stacked_violin",
        n_clusters=n_clusters,
        genes_per_fig=40,
    )

    _plot_gene_cluster_panels(
        sc.pl.heatmap,
        adata,
        top_genes,
        groupby,
        save_path,
        base_name="heatmap",
        n_clusters=n_clusters,
        genes_per_fig=40,
    )

    # 4-3. UMAP: cell_type 개수에 따라 가로 크기를 확장해서 legend가 잘리지 않도록
    if "cell_type" in adata.obs:
        n_celltypes = adata.obs["cell_type"].nunique(dropna=True)
    else:
        n_celltypes = n_clusters

    umap_width = 6.0 + 0.4 * max(n_celltypes - 5, 0)
    umap_width = max(6.0, min(20.0, umap_width))
    umap_height = 6.0

    with plt.rc_context({"figure.figsize": (umap_width, umap_height)}):
        sc.pl.umap(adata, color="cell_type", show=False)
        plt.tight_layout()
        plt.savefig(
            os.path.join(save_path, "predicted_cell_types.png"),
            dpi=300,
        )
        plt.close()

    print("DEG + Cell Type Prediction 완료 (동적 이미지 크기 적용)")
    return adata