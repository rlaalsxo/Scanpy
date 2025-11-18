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
    """리스트를 chunk_size 단위로 나누는 유틸."""
    for i in range(0, len(lst), chunk_size):
        yield lst[i : i + chunk_size]


def _calc_figsize_genes_clusters(
    n_genes,
    n_clusters,
    base_size=(3.5, 3.5),   # (기본 가로, 기본 세로)
    gene_width=0.7,        # gene 하나 늘어날 때 가로 증가량
    cluster_height=0.30,    # cluster 하나 늘어날 때 세로 증가량
    min_size=(4.5, 4.0),
    max_size=(18.0, 18.0),
):
    # gene 수가 적을 때는 가로를 너무 넓히지 않도록 약하게
    if n_genes <= 8:
        gene_width_eff = gene_width * 0.6
    else:
        gene_width_eff = gene_width

    # cluster 수가 적을 때는 세로를 너무 길게 만들지 않도록 약하게
    if n_clusters <= 8:
        cluster_height_eff = cluster_height * 0.7
    else:
        cluster_height_eff = cluster_height

    # 가로: gene 수에 비례
    width = base_size[0] + gene_width_eff * max(n_genes - 1, 0)
    # 세로: cluster 수에 비례
    height = base_size[1] + cluster_height_eff * max(n_clusters - 1, 0)

    # 최소/최대 사이즈 제한
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
    dotplot / matrixplot / stacked_violin / heatmap 공통 처리.

    - genes_per_fig 개 단위로 여러 장으로 나눔
    - 각 장마다 gene/cluster 수에 맞게 figsize 계산
    - bbox_inches="tight" 로 저장해서 위쪽/양옆 여백 제거
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
        )

        # scanpy 함수에 figsize 직접 전달
        plot_func(
            adata,
            var_names=genes_chunk,
            groupby=groupby,
            use_raw=False,
            show=False,
            figsize=(fig_w, fig_h),
        )
        fig = plt.gcf()
        fig.tight_layout()

        if n_chunks == 1:
            fname = f"{base_name}.png"
        else:
            fname = f"{base_name}_{idx}.png"

        fig.savefig(
            os.path.join(save_path, fname),
            dpi=300,
            bbox_inches="tight",  # 여백 자동 제거
        )
        plt.close(fig)


def deg_analysis_with_sex_gene_filtering(
    adata,
    save_path,
    species,
    groupby="leiden",
    top_n=20,
    padj_threshold=0.05,
):
    os.makedirs(save_path, exist_ok=True)

    # -------------------------------------------------------------
    # 1. 성염색체 유전자 제거
    # -------------------------------------------------------------
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

    # -------------------------------------------------------------
    # 2. DEG 분석
    # -------------------------------------------------------------
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

    if groupby not in adata.obs:
        raise KeyError(f"{groupby} not found in adata.obs")
    n_clusters = adata.obs[groupby].nunique()

    # -------------------------------------------------------------
    # 3. Enrichr 로 cell type 예측
    # -------------------------------------------------------------
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

    # -------------------------------------------------------------
    # 4-1. rank_genes_groups: gene/cluster 기반으로 크기 조절 + 여백 제거
    # -------------------------------------------------------------
    fig_w, fig_h = _calc_figsize_genes_clusters(
        n_genes=top_n,
        n_clusters=n_clusters,
        base_size=(10.0, 5),
        gene_width=0.8,
        cluster_height=0.6,
        min_size=(6.0, 4.0),
        max_size=(20.0, 16.0),
    )
    sc.pl.rank_genes_groups(
        adata,
        n_genes=top_n,
        key="wilcoxon",
        sharey=False,
        show=False,
        figsize=(fig_w, fig_h),
    )
    fig = plt.gcf()
    fig.tight_layout()
    fig.savefig(
        os.path.join(save_path, "rank_genes_groups.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # -------------------------------------------------------------
    # 4-2. dotplot / matrixplot / stacked_violin / heatmap
    # -------------------------------------------------------------
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

    # -------------------------------------------------------------
    # 4-3. UMAP: cell_type legend + 정사각형 지도
    # -------------------------------------------------------------
    if "cell_type" in adata.obs:
        n_celltypes = adata.obs["cell_type"].nunique(dropna=True)
    else:
        n_celltypes = n_clusters

    # 지도 부분은 정사각형, 오른쪽에 legend 공간만 약간 추가
    base_side = 20.0
    extra_for_legend = min(5.0, 0.25 * max(n_celltypes - 8, 0))
    fig_w_umap = base_side + extra_for_legend
    fig_h_umap = base_side

    with plt.rc_context({"figure.figsize": (fig_w_umap, fig_h_umap)}):
        sc.pl.umap(
            adata,
            color="cell_type",
            show=False,
            legend_loc="right margin",
        )
        fig = plt.gcf()
        ax = plt.gca()

        # 축 position 을 조정해서 실제 UMAP 축은 정사각형이 되도록 맞춤
        H = fig.get_figheight()
        W = fig.get_figwidth()
        height_frac = 0.9
        width_frac = height_frac * (H / W)  # width * W == height * H → 정사각형

        left = 0.08
        bottom = 0.07
        ax.set_position([left, bottom, width_frac, height_frac])

        fig.tight_layout()
        fig.savefig(
            os.path.join(save_path, "predicted_cell_types.png"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig)

    print("DEG + Cell Type Prediction 완료 (동적 이미지 + 여백 정리)")
    return adata
