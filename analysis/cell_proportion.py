"""
Cell Proportion 분석

샘플/배치별 세포 비율 비교
"""
import os

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from config.defaults import CELL_PROPORTION
from core.neighbors import detect_cluster_key, ensure_umap
from plotting.utils import save_figure


def cell_proportion_analysis(
    adata,
    save_path: str,
    groupby: str = None,
    batch_key: str = "sample",
) -> sc.AnnData:
    """
    Cell Proportion 분석

    Parameters
    ----------
    adata : AnnData
    save_path : str
    groupby : str, optional
        세포 그룹 키 (None이면 cell_type 또는 cluster key 자동 탐지)
    batch_key : str
        샘플/배치 키

    Returns
    -------
    AnnData
    """
    os.makedirs(save_path, exist_ok=True)

    # 1. groupby 키 결정
    if groupby is None:
        if "cell_type" in adata.obs:
            groupby = "cell_type"
        else:
            groupby = detect_cluster_key(adata)
    if groupby is None or groupby not in adata.obs:
        raise ValueError(f"[CellProportion] groupby key not found: {groupby}")

    # batch_key 확인
    if batch_key not in adata.obs:
        print(f"[CellProportion] Warning: '{batch_key}' not in adata.obs. Using groupby only.")
        batch_key = None

    print(f"[CellProportion] groupby: {groupby}, batch_key: {batch_key}")

    # 2. Crosstab 계산
    if batch_key is not None:
        counts = pd.crosstab(adata.obs[batch_key], adata.obs[groupby])
        fractions = counts.div(counts.sum(axis=1), axis=0)

        # CSV 저장
        counts.to_csv(os.path.join(save_path, "cell_proportion_counts.csv"))
        fractions.to_csv(os.path.join(save_path, "cell_proportion_fractions.csv"))
        print(f"[CellProportion] Samples: {counts.shape[0]}, Groups: {counts.shape[1]}")

        # 3. 시각화
        # 3-1. Stacked barplot (샘플별 세포 비율)
        print("[CellProportion] Generating stacked barplot...")
        n_samples = len(counts.index)
        n_groups = len(counts.columns)
        fig_h = max(4, n_samples * 0.6)
        fig_w = max(10, n_groups * 0.3 + 6)

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        fractions.plot(
            kind="barh",
            stacked=True,
            ax=ax,
            colormap="tab20",
            edgecolor="white",
            linewidth=0.5,
        )
        ax.set_xlabel("Proportion")
        ax.set_ylabel(batch_key.capitalize())
        ax.set_title(f"Cell Proportion by {batch_key.capitalize()}")
        ax.legend(
            title=groupby,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            fontsize=7,
        )
        ax.set_xlim(0, 1)
        save_figure(fig, save_path, "cell_proportion_stacked_bar.png")

        # 3-2. Proportion heatmap
        print("[CellProportion] Generating proportion heatmap...")
        fig, ax = plt.subplots(figsize=(max(8, n_groups * 0.5), max(4, n_samples * 0.5)))
        im = ax.imshow(fractions.values, cmap="YlOrRd", aspect="auto", vmin=0, vmax=fractions.values.max())
        ax.set_xticks(range(n_groups))
        ax.set_xticklabels(fractions.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(n_samples))
        ax.set_yticklabels(fractions.index, fontsize=8)
        ax.set_xlabel(groupby)
        ax.set_ylabel(batch_key.capitalize())
        ax.set_title(f"Cell Proportion Heatmap")

        # 셀 내 비율 텍스트 표시
        for i in range(n_samples):
            for j in range(n_groups):
                val = fractions.values[i, j]
                if val > 0.01:
                    color = "white" if val > fractions.values.max() * 0.6 else "black"
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6, color=color)

        plt.colorbar(im, ax=ax, label="Proportion")
        save_figure(fig, save_path, "cell_proportion_heatmap.png")

    else:
        # batch_key 없으면 전체 비율만 계산
        value_counts = adata.obs[groupby].value_counts()
        fractions_series = value_counts / value_counts.sum()

        df = pd.DataFrame({"count": value_counts, "fraction": fractions_series})
        df.to_csv(os.path.join(save_path, "cell_proportion_counts.csv"))

        # barplot
        fig, ax = plt.subplots(figsize=(max(8, len(value_counts) * 0.5), 5))
        bars = ax.bar(range(len(value_counts)), fractions_series.values, color=plt.cm.tab20(np.arange(len(value_counts)) % 20))
        ax.set_xticks(range(len(value_counts)))
        ax.set_xticklabels(value_counts.index, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Proportion")
        ax.set_title("Cell Proportion")
        save_figure(fig, save_path, "cell_proportion_bar.png")

    # 3-3. UMAP 참고 시각화
    print("[CellProportion] Generating UMAP plots...")
    ensure_umap(adata)

    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color=groupby, show=False, ax=ax, legend_loc="right margin")
    save_figure(fig, save_path, "umap_celltype.png")

    if batch_key is not None:
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(adata, color=batch_key, show=False, ax=ax, legend_loc="right margin")
        save_figure(fig, save_path, "umap_sample.png")

    print(f"[CellProportion] Complete. Results: {save_path}")
    return adata
