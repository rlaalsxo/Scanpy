"""
Cell Cycle 분석

S/G2M phase 스코어링 및 시각화
"""
import os

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from config.species import get_cell_cycle_genes
from config.defaults import NEIGHBORS
from core.neighbors import compute_pca, compute_neighbors, detect_cluster_key
from plotting.utils import save_figure


def _plot_phase_composition(adata, groupby_key: str, save_path: str):
    """
    phase별 구성비 요약 및 stacked bar plot
    """
    if "phase" not in adata.obs:
        print(f"[CellCycle] `phase` not found in adata.obs. Skip composition.")
        return
    if groupby_key not in adata.obs:
        print(f"[CellCycle] `{groupby_key}` not found. Skip composition.")
        return

    # 교차표
    comp_counts = pd.crosstab(adata.obs[groupby_key], adata.obs["phase"])
    if comp_counts.empty:
        return

    comp_frac = comp_counts.div(comp_counts.sum(axis=1), axis=0)

    # 테이블 저장
    comp_counts.to_csv(os.path.join(save_path, f"cell_cycle_composition_{groupby_key}_counts.csv"))
    comp_frac.to_csv(os.path.join(save_path, f"cell_cycle_composition_{groupby_key}_fraction.csv"))

    # Phase 순서 정렬
    phase_order = ["G1", "S", "G2M"]
    cols = [p for p in phase_order if p in comp_frac.columns] + [
        p for p in comp_frac.columns if p not in phase_order
    ]
    comp_frac = comp_frac[cols]

    # Stacked bar plot
    x = np.arange(comp_frac.shape[0])
    fig_width = max(6.0, comp_frac.shape[0] * 0.6)
    fig, ax = plt.subplots(figsize=(fig_width, 4.0))

    bottom = np.zeros_like(x, dtype=float)
    for ph in comp_frac.columns:
        vals = comp_frac[ph].values
        ax.bar(x, vals, bottom=bottom, label=ph)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(comp_frac.index, rotation=45, ha="right")
    ax.set_ylabel("Fraction of cells")
    ax.set_xlabel(groupby_key)
    ax.set_title(f"Cell cycle phase composition by {groupby_key}")
    ax.legend(title="Phase", bbox_to_anchor=(1.05, 1.0), loc="upper left")

    save_figure(fig, save_path, f"cell_cycle_composition_{groupby_key}_stacked_bar.png")


def _compare_umap_cell_cycle_effect(adata, save_path: str, color: str = "phase"):
    """
    Cell cycle 효과 제거 전후 UMAP 비교
    """
    if "S_score" not in adata.obs or "G2M_score" not in adata.obs:
        print("[CellCycle] S_score/G2M_score not found. Skip UMAP comparison.")
        return adata

    has_umap_before = "X_umap" in adata.obsm

    # Cell cycle score 회귀 제거
    adata_cc = adata.copy()
    sc.pp.regress_out(adata_cc, ["S_score", "G2M_score"])
    compute_pca(adata_cc)
    compute_neighbors(adata_cc, n_neighbors=10, n_pcs=NEIGHBORS["n_pcs"])
    sc.tl.umap(adata_cc)

    # Regressed UMAP 저장
    adata.obsm["X_umap_cell_cycle_regressed"] = adata_cc.obsm["X_umap"].copy()

    # 시각화
    if has_umap_before:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        sc.pl.umap(adata, color=color, ax=axes[0], show=False, title="Before regressing cell cycle")
        sc.pl.umap(adata_cc, color=color, ax=axes[1], show=False, title="After regressing cell cycle")
        save_figure(fig, save_path, "umap_cell_cycle_before_after.png")
    else:
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.umap(adata_cc, color=color, show=False, ax=ax, title="UMAP after regressing cell cycle")
        save_figure(fig, save_path, "umap_after_cell_cycle_regress.png")

    return adata


def score_cell_cycle(
    adata,
    save_path: str,
    species: str,
) -> sc.AnnData:
    """
    Cell Cycle 스코어링

    Parameters
    ----------
    adata : AnnData
    save_path : str
    species : str

    Returns
    -------
    AnnData
    """
    os.makedirs(save_path, exist_ok=True)

    # 1. Cell cycle gene set 로딩
    print("[CellCycle] Loading cell cycle genes...")
    s_genes, g2m_genes = get_cell_cycle_genes(species)
    s_genes = [g for g in s_genes if g in adata.var_names]
    g2m_genes = [g for g in g2m_genes if g in adata.var_names]

    print(f"[CellCycle] S genes used: {len(s_genes)}")
    print(f"[CellCycle] G2M genes used: {len(g2m_genes)}")

    # 2. Cell cycle score 계산
    print("[CellCycle] Scoring cell cycle...")
    sc.pp.scale(adata, max_value=10)
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    # 3. Violin plot
    fig, ax = plt.subplots(figsize=(10, 6))
    sc.pl.violin(adata, ["S_score", "G2M_score"], jitter=0.4, groupby="sample", rotation=60, show=False, ax=ax)
    save_figure(fig, save_path, "violin_cell_cycle.png")

    # 4. Scatter plot
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.scatter(adata, x="S_score", y="G2M_score", color="phase", show=False, ax=ax)
    save_figure(fig, save_path, "scatter_cell_cycle.png")

    # 5. Phase composition by sample
    if "sample" in adata.obs:
        _plot_phase_composition(adata, "sample", save_path)

    # 6. Phase composition by cluster
    cluster_key = detect_cluster_key(adata)
    if cluster_key:
        _plot_phase_composition(adata, cluster_key, save_path)

    # 7. Cell cycle 효과 제거 비교
    _compare_umap_cell_cycle_effect(adata, save_path)

    print(f"[CellCycle] Complete. Results: {save_path}")
    return adata
