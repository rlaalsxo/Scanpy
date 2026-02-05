import matplotlib
matplotlib.use("Agg")
import os
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np


def _detect_cluster_key(adata):
    """클러스터 키 자동 탐지 (leiden > louvain > clusters > cluster)"""
    for cand in ["leiden", "louvain", "clusters", "cluster"]:
        if cand in adata.obs:
            return cand
    return None


def _ensure_prerequisites(adata):
    """Trajectory 분석 필수 조건 확인 및 자동 생성"""
    # PCA 없으면 생성
    if "X_pca" not in adata.obsm:
        print("[Trajectory] Computing PCA...")
        sc.tl.pca(adata, svd_solver="arpack")

    # neighbors 재계산 (BBKNN 호환성 - 항상 표준 neighbors 사용)
    print("[Trajectory] Computing neighbors for PAGA compatibility...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)

    # 클러스터 키 탐지
    cluster_key = _detect_cluster_key(adata)
    if cluster_key is None:
        print("[Trajectory] Computing leiden clustering...")
        sc.tl.leiden(adata, key_added="leiden")
        cluster_key = "leiden"

    return cluster_key


def _select_root_cell(adata, root_cluster=None, root_gene=None):
    """
    DPT root cell 선택
    - root_cluster: 해당 클러스터에서 DC1이 최소인 세포
    - root_gene: 해당 유전자 발현이 최대인 세포
    - 둘 다 없으면: DC1이 최소인 세포 (자동)
    """
    if root_cluster is not None:
        cluster_key = _detect_cluster_key(adata)
        if cluster_key and root_cluster in adata.obs[cluster_key].values:
            mask = adata.obs[cluster_key] == root_cluster
            dc1 = adata.obsm["X_diffmap"][mask, 0]
            return np.where(mask)[0][np.argmin(dc1)]

    if root_gene is not None and root_gene in adata.var_names:
        gene_idx = adata.var_names.get_loc(root_gene)
        expr = adata.X[:, gene_idx]
        if hasattr(expr, "toarray"):
            expr = expr.toarray().flatten()
        return int(np.argmax(expr))

    return int(np.argmin(adata.obsm["X_diffmap"][:, 0]))


def trajectory_analysis(adata, save_path, species="human",
                        cluster_key=None, root_cluster=None,
                        root_gene=None, paga_threshold=0.05, n_dcs=15):
    """
    Trajectory 분석 (PAGA + Diffusion Pseudotime)

    Parameters
    ----------
    adata : AnnData
    save_path : str
    species : str
    cluster_key : str, optional (None이면 자동 탐지)
    root_cluster : str, optional (DPT 시작 클러스터)
    root_gene : str, optional (DPT 시작 마커 유전자)
    paga_threshold : float
    n_dcs : int (Diffusion Component 수)

    Returns
    -------
    AnnData with trajectory results
    """
    os.makedirs(save_path, exist_ok=True)

    # 1. 필수 조건 확인 및 자동 생성
    detected_key = _ensure_prerequisites(adata)
    if cluster_key is None:
        cluster_key = detected_key
    print(f"[Trajectory] cluster_key: {cluster_key}")

    n_clusters = adata.obs[cluster_key].nunique()

    # 2. PAGA 분석
    print("[Trajectory] Running PAGA...")
    sc.tl.paga(adata, groups=cluster_key)

    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.paga(adata, threshold=paga_threshold, show=False, ax=ax,
               fontsize=10, node_size_scale=1.5)
    fig.tight_layout()
    fig.savefig(os.path.join(save_path, "paga_graph.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # PAGA-initialized UMAP
    sc.tl.umap(adata, init_pos="paga")

    # PAGA + UMAP 오버레이
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sc.pl.paga(adata, threshold=paga_threshold, show=False, ax=axes[0],
               fontsize=10, node_size_scale=1.5, title="PAGA Graph")
    sc.pl.umap(adata, color=cluster_key, show=False, ax=axes[1],
               title="PAGA-initialized UMAP")
    fig.tight_layout()
    fig.savefig(os.path.join(save_path, "paga_umap_overlay.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # 3. Diffusion Map
    print("[Trajectory] Computing diffusion map...")
    sc.tl.diffmap(adata, n_comps=n_dcs)

    # Diffusion Components 시각화
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    for i, ax in enumerate(axes.flatten()):
        if i < min(4, n_dcs - 1):
            sc.pl.embedding(adata, basis="diffmap", color=cluster_key,
                            components=[1, i + 2], ax=ax, show=False,
                            title=f"DC1 vs DC{i + 2}")
    fig.tight_layout()
    fig.savefig(os.path.join(save_path, "diffusion_components.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # 4. Diffusion Pseudotime
    print("[Trajectory] Computing diffusion pseudotime...")
    root_idx = _select_root_cell(adata, root_cluster, root_gene)
    adata.uns["iroot"] = root_idx
    print(f"[Trajectory] Root cell: {root_idx} (cluster: {adata.obs[cluster_key].iloc[root_idx]})")

    sc.tl.dpt(adata, n_dcs=n_dcs)

    # 5. 시각화
    # Pseudotime UMAP
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color="dpt_pseudotime", ax=ax, show=False,
               color_map="viridis", title="Diffusion Pseudotime")
    fig.tight_layout()
    fig.savefig(os.path.join(save_path, "diffusion_pseudotime_umap.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # 클러스터별 pseudotime 분포
    fig_width = max(8, n_clusters * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    sc.pl.violin(adata, keys="dpt_pseudotime", groupby=cluster_key,
                 ax=ax, show=False, rotation=45)
    ax.set_title("Pseudotime Distribution by Cluster")
    fig.tight_layout()
    fig.savefig(os.path.join(save_path, "pseudotime_by_cluster.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Root cell 위치 표시
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color=cluster_key, ax=ax, show=False)
    root_coords = adata.obsm["X_umap"][root_idx]
    ax.scatter(root_coords[0], root_coords[1], c="red", s=200,
               marker="*", edgecolors="black", linewidths=2, zorder=10)
    ax.set_title(f"Root Cell (Cluster: {adata.obs[cluster_key].iloc[root_idx]})")
    fig.tight_layout()
    fig.savefig(os.path.join(save_path, "root_cell_location.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # 6. 요약 테이블
    summary = adata.obs.groupby(cluster_key)["dpt_pseudotime"].agg(["mean", "std", "min", "max"])
    summary = summary.sort_values("mean")
    summary.to_csv(os.path.join(save_path, "trajectory_summary.csv"))

    print(f"[Trajectory] Complete. Results: {save_path}")
    return adata
