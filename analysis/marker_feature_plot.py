"""
Marker Gene Feature Plot

마커 유전자의 UMAP 발현 시각화
"""
import os
import math

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from config.defaults import MARKER_FEATURE, DEG
from core.neighbors import detect_cluster_key, ensure_umap
from plotting.utils import save_figure


def marker_feature_plot(
    adata,
    save_path: str,
    genes: list = None,
    top_n_genes: int = None,
    ncols: int = None,
) -> sc.AnnData:
    """
    Marker Gene Feature Plot

    Parameters
    ----------
    adata : AnnData
    save_path : str
    genes : list, optional
        시각화할 유전자 리스트 (None이면 DEG에서 추출)
    top_n_genes : int, optional
        클러스터당 상위 유전자 수
    ncols : int, optional
        UMAP grid 열 수

    Returns
    -------
    AnnData
    """
    os.makedirs(save_path, exist_ok=True)

    if top_n_genes is None:
        top_n_genes = MARKER_FEATURE["top_n_genes"]
    if ncols is None:
        ncols = MARKER_FEATURE["ncols"]

    groupby = detect_cluster_key(adata)
    ensure_umap(adata)

    # 1. 유전자 리스트 결정
    if genes is not None:
        # 사용자 지정 유전자 - adata에 존재하는 것만 필터
        valid_genes = [g for g in genes if g in adata.var_names]
        missing = [g for g in genes if g not in adata.var_names]
        if missing:
            print(f"[MarkerFeature] Warning: genes not found: {missing}")
        if not valid_genes:
            print("[MarkerFeature] No valid genes found.")
            return adata
        selected_genes = valid_genes
        print(f"[MarkerFeature] Using {len(selected_genes)} user-specified genes")
    else:
        # DEG 결과에서 추출
        deg_key = "wilcoxon"
        if deg_key not in adata.uns:
            print("[MarkerFeature] DEG results not found. Running rank_genes_groups...")
            if groupby is None:
                raise ValueError("[MarkerFeature] No cluster key found.")
            sc.tl.rank_genes_groups(adata, groupby=groupby, method=DEG["method"], key_added=deg_key)

        clusters = adata.obs[groupby].astype("category").cat.categories
        selected_genes = []
        for cl in clusters:
            cl_genes = (
                sc.get.rank_genes_groups_df(adata, group=cl, key=deg_key)
                .query("pvals_adj < @DEG['padj_threshold']")
                .head(top_n_genes)["names"]
                .tolist()
            )
            selected_genes.extend(cl_genes)

        # 중복 제거 (순서 유지)
        seen = set()
        unique_genes = []
        for g in selected_genes:
            if g not in seen:
                seen.add(g)
                unique_genes.append(g)
        selected_genes = unique_genes
        print(f"[MarkerFeature] Selected {len(selected_genes)} genes from DEG (top {top_n_genes}/cluster)")

    if not selected_genes:
        print("[MarkerFeature] No genes to plot.")
        return adata

    # 2. Feature UMAP grid
    print("[MarkerFeature] Generating feature UMAP plots...")
    max_genes_per_page = ncols * 5  # 5행씩
    n_pages = math.ceil(len(selected_genes) / max_genes_per_page)

    for page in range(n_pages):
        start = page * max_genes_per_page
        end = min(start + max_genes_per_page, len(selected_genes))
        page_genes = selected_genes[start:end]

        n_genes = len(page_genes)
        nrows = math.ceil(n_genes / ncols)

        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(ncols * 4, nrows * 3.5),
        )
        if nrows == 1 and ncols == 1:
            axes = np.array([[axes]])
        elif nrows == 1:
            axes = axes[np.newaxis, :]
        elif ncols == 1:
            axes = axes[:, np.newaxis]

        for i, gene in enumerate(page_genes):
            row, col = divmod(i, ncols)
            ax = axes[row, col]
            sc.pl.umap(adata, color=gene, show=False, ax=ax, use_raw=False, frameon=False, title=gene)

        # 빈 subplot 제거
        for i in range(n_genes, nrows * ncols):
            row, col = divmod(i, ncols)
            axes[row, col].set_visible(False)

        suffix = f"_{page + 1}" if n_pages > 1 else ""
        save_figure(fig, save_path, f"feature_umap{suffix}.png")

    # 3. Violin plot (클러스터별 발현 분포)
    if groupby is not None:
        print("[MarkerFeature] Generating violin plots...")
        max_genes_violin = 20
        n_violin_pages = math.ceil(len(selected_genes) / max_genes_violin)

        for page in range(n_violin_pages):
            start = page * max_genes_violin
            end = min(start + max_genes_violin, len(selected_genes))
            page_genes = selected_genes[start:end]

            n_clusters = adata.obs[groupby].nunique()
            fig_h = max(3, len(page_genes) * 0.8)
            fig_w = max(8, n_clusters * 0.4 + 4)

            fig, axes = plt.subplots(
                len(page_genes), 1,
                figsize=(fig_w, fig_h * len(page_genes) / max(len(page_genes), 4)),
                sharex=True,
            )
            if len(page_genes) == 1:
                axes = [axes]

            for i, gene in enumerate(page_genes):
                sc.pl.violin(
                    adata, gene, groupby=groupby,
                    show=False, ax=axes[i], use_raw=False,
                    rotation=45, size=0,
                )
                axes[i].set_title("")
                axes[i].set_ylabel(gene, fontsize=8, rotation=0, labelpad=50, ha="right")

            fig.suptitle("Marker Gene Expression by Cluster", fontsize=12)
            suffix = f"_{page + 1}" if n_violin_pages > 1 else ""
            save_figure(fig, save_path, f"marker_violin{suffix}.png")

    # 4. 유전자 리스트 저장
    pd.DataFrame({"gene": selected_genes}).to_csv(
        os.path.join(save_path, "marker_genes.csv"), index=False
    )

    print(f"[MarkerFeature] Complete. Results: {save_path}")
    return adata
