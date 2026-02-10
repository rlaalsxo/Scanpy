"""
Pseudotime Gene Expression 분석

Pseudotime에 따른 유전자 발현 변화 시각화
"""
import os

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from config.defaults import PSEUDOTIME_EXPRESSION, DEG
from core.neighbors import detect_cluster_key
from plotting.utils import save_figure


def pseudotime_expression(
    adata,
    save_path: str,
    n_bins: int = None,
    top_n_genes: int = None,
) -> sc.AnnData:
    """
    Pseudotime Gene Expression 분석

    Parameters
    ----------
    adata : AnnData
    save_path : str
    n_bins : int, optional
    top_n_genes : int, optional

    Returns
    -------
    AnnData
    """
    os.makedirs(save_path, exist_ok=True)

    if n_bins is None:
        n_bins = PSEUDOTIME_EXPRESSION["n_bins"]
    if top_n_genes is None:
        top_n_genes = PSEUDOTIME_EXPRESSION["top_n_genes"]

    # 1. Pseudotime 확인
    if "dpt_pseudotime" not in adata.obs:
        raise ValueError("[PseudotimeExpr] Pseudotime not found. Run Trajectory analysis (Step 5) first.")

    pseudotime = adata.obs["dpt_pseudotime"].values
    valid_mask = ~np.isinf(pseudotime) & ~np.isnan(pseudotime)
    print(f"[PseudotimeExpr] Valid cells: {valid_mask.sum()} / {len(pseudotime)}")

    # 2. Pseudotime과 상관관계가 높은 유전자 선택
    print("[PseudotimeExpr] Computing gene-pseudotime correlations...")
    pt_valid = pseudotime[valid_mask]

    # expression matrix
    X = adata.X[valid_mask]
    if hasattr(X, "toarray"):
        X = X.toarray()

    correlations = []
    for i in range(X.shape[1]):
        expr = X[:, i]
        if expr.std() == 0:
            correlations.append(0.0)
            continue
        rho, _ = spearmanr(pt_valid, expr)
        correlations.append(rho if not np.isnan(rho) else 0.0)

    corr_df = pd.DataFrame({
        "gene": adata.var_names,
        "correlation": correlations,
        "abs_correlation": np.abs(correlations),
    }).sort_values("abs_correlation", ascending=False)

    corr_df.to_csv(os.path.join(save_path, "pseudotime_correlation.csv"), index=False)
    print(f"[PseudotimeExpr] Top correlated gene: {corr_df.iloc[0]['gene']} (rho={corr_df.iloc[0]['correlation']:.3f})")

    # 상위 유전자 선택
    top_genes = corr_df.head(top_n_genes)["gene"].tolist()

    # 3. Pseudotime bins으로 평균 발현량 계산
    print("[PseudotimeExpr] Binning expression along pseudotime...")
    bin_edges = np.linspace(pt_valid.min(), pt_valid.max(), n_bins + 1)
    bin_labels = np.digitize(pt_valid, bin_edges) - 1
    bin_labels = np.clip(bin_labels, 0, n_bins - 1)

    gene_indices = [list(adata.var_names).index(g) for g in top_genes]
    binned_expr = np.zeros((top_n_genes, n_bins))
    for b in range(n_bins):
        mask = bin_labels == b
        if mask.sum() > 0:
            binned_expr[:, b] = X[mask][:, gene_indices].mean(axis=0)

    # z-score 정규화 (행 기준)
    row_means = binned_expr.mean(axis=1, keepdims=True)
    row_stds = binned_expr.std(axis=1, keepdims=True)
    row_stds[row_stds == 0] = 1
    binned_expr_z = (binned_expr - row_means) / row_stds

    # 4. 시각화
    # 4-1. Pseudotime 히트맵
    print("[PseudotimeExpr] Generating pseudotime heatmap...")
    fig, ax = plt.subplots(figsize=(14, max(6, top_n_genes * 0.25)))
    im = ax.imshow(binned_expr_z, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
    ax.set_yticks(range(len(top_genes)))
    ax.set_yticklabels(top_genes, fontsize=7)
    ax.set_xlabel("Pseudotime bins")
    ax.set_ylabel("Genes")
    ax.set_title(f"Gene Expression along Pseudotime (top {top_n_genes} genes)")
    plt.colorbar(im, ax=ax, label="z-score")
    save_figure(fig, save_path, "pseudotime_heatmap.png")

    # 4-2. 상위 유전자 트렌드 라인 플롯
    print("[PseudotimeExpr] Generating trend plot...")
    n_plot = min(10, top_n_genes)
    fig, axes = plt.subplots(
        n_plot, 1,
        figsize=(10, n_plot * 2),
        sharex=True,
    )
    if n_plot == 1:
        axes = [axes]

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    for i in range(n_plot):
        ax = axes[i]
        gene = top_genes[i]
        rho = corr_df[corr_df["gene"] == gene]["correlation"].values[0]

        ax.plot(bin_centers, binned_expr[i], color="steelblue", linewidth=2)
        ax.fill_between(bin_centers, binned_expr[i], alpha=0.2, color="steelblue")
        ax.set_ylabel(gene, fontsize=9, rotation=0, labelpad=60, ha="right")
        ax.text(0.98, 0.85, f"rho={rho:.3f}", transform=ax.transAxes, fontsize=8, ha="right")
        ax.tick_params(axis="y", labelsize=7)

    axes[-1].set_xlabel("Pseudotime")
    fig.suptitle("Gene Expression Trends along Pseudotime", fontsize=14, y=1.01)
    save_figure(fig, save_path, "pseudotime_trend.png")

    print(f"[PseudotimeExpr] Complete. Results: {save_path}")
    return adata
