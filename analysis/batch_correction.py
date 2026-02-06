"""
배치 보정

BBKNN 기반 배치 효과 제거
"""
import os

import scanpy as sc
import bbknn
import matplotlib.pyplot as plt

from config.defaults import BBKNN_NEIGHBORS_WITHIN_BATCH, NEIGHBORS
from core.neighbors import compute_pca, compute_neighbors, compute_clustering
from plotting.utils import save_figure


def _get_bbknn_neighbors(n_obs: int) -> int:
    """세포 수 기반 neighbors_within_batch 자동 결정"""
    for threshold, value in sorted(BBKNN_NEIGHBORS_WITHIN_BATCH.items(), reverse=True):
        if n_obs > threshold:
            return value
    return BBKNN_NEIGHBORS_WITHIN_BATCH[0]


def batch_correction(
    adata,
    save_path: str,
    batch_key: str = "sample",
    neighbors_within_batch: int = None,
) -> sc.AnnData:
    """
    BBKNN 기반 배치 보정

    Parameters
    ----------
    adata : AnnData
    save_path : str
    batch_key : str
        배치 키 컬럼
    neighbors_within_batch : int, optional
        배치 내 neighbors 수 (None이면 자동)

    Returns
    -------
    AnnData
        배치 보정된 AnnData
    """
    os.makedirs(save_path, exist_ok=True)

    # neighbors_within_batch 자동 설정
    if neighbors_within_batch is None:
        neighbors_within_batch = _get_bbknn_neighbors(adata.n_obs)
        print(f"[BatchCorrection] neighbors_within_batch 자동 설정: {neighbors_within_batch}")

    # 1. PCA
    print("[BatchCorrection] Computing PCA...")
    compute_pca(adata)

    # 2. BBKNN 적용 (복사본에서)
    print("[BatchCorrection] Running BBKNN...")
    adata_bbknn = adata.copy()
    bbknn.bbknn(adata_bbknn, batch_key=batch_key, neighbors_within_batch=neighbors_within_batch)

    # 3. Leiden 클러스터링
    print("[BatchCorrection] Running Leiden clustering...")
    compute_clustering(adata_bbknn, method="leiden")

    # 4. UMAP
    print("[BatchCorrection] Computing UMAP...")
    sc.tl.umap(adata_bbknn)

    # 5. 비보정 데이터 UMAP (비교용)
    compute_neighbors(adata, n_neighbors=10, n_pcs=NEIGHBORS["n_pcs"])
    sc.tl.umap(adata)

    # 6. 시각화: 비보정 vs 보정
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.umap(adata, color="sample", title="Uncorrected UMAP", show=False, ax=axes[0])
    sc.pl.umap(adata_bbknn, color="sample", title="BBKNN Corrected UMAP", show=False, ax=axes[1])
    save_figure(fig, save_path, "umap_batch_correction_comparison.png")

    # 7. 시각화: Leiden 클러스터
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata_bbknn, color="leiden", title="Leiden Clusters", show=False, ax=ax)
    save_figure(fig, save_path, "umap_leiden_clusters.png")

    print(f"[BatchCorrection] Complete. Results: {save_path}")
    return adata_bbknn
