"""
Trajectory 분석

PAGA + Diffusion Pseudotime
"""
import os

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

from config.defaults import TRAJECTORY, NEIGHBORS
from core.neighbors import compute_pca, compute_neighbors, compute_clustering
from plotting.utils import save_figure


def _select_root_cell(adata, root_cluster=None, root_gene=None, cluster_key=None):
    """
    DPT root cell 선택

    - root_cluster: 해당 클러스터에서 DC1이 최소인 세포
    - root_gene: 해당 유전자 발현이 최대인 세포
    - 둘 다 없으면: DC1이 최소인 세포 (자동)
    """
    if root_cluster is not None:
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


def trajectory_analysis(
    adata,
    save_path: str,
    species: str = "human",
    root_cluster: str = None,
    root_gene: str = None,
    paga_threshold: float = None,
    n_dcs: int = None,
) -> sc.AnnData:
    """
    Trajectory 분석 (PAGA + Diffusion Pseudotime)

    Parameters
    ----------
    adata : AnnData
    save_path : str
    species : str
    root_cluster : str, optional
    root_gene : str, optional
    paga_threshold : float, optional
    n_dcs : int, optional

    Returns
    -------
    AnnData
    """
    os.makedirs(save_path, exist_ok=True)

    if paga_threshold is None:
        paga_threshold = TRAJECTORY["paga_threshold"]
    if n_dcs is None:
        n_dcs = TRAJECTORY["n_dcs"]

    # 1. 복사본에서 trajectory 분석 수행 (기본 키 사용으로 PAGA 호환성 확보)
    print("[Trajectory] Creating copy for trajectory analysis...")
    adata_traj = adata.copy()

    print("[Trajectory] Computing PCA...")
    sc.tl.pca(adata_traj, n_comps=50, svd_solver="arpack")

    print("[Trajectory] Computing neighbors...")
    sc.pp.neighbors(adata_traj, n_neighbors=NEIGHBORS["n_neighbors"], n_pcs=NEIGHBORS["n_pcs"])

    cluster_key = "leiden_traj"
    print("[Trajectory] Computing clustering...")
    sc.tl.leiden(adata_traj, key_added=cluster_key)

    print(f"[Trajectory] cluster_key: {cluster_key}")
    n_clusters = adata_traj.obs[cluster_key].nunique()

    # 2. PAGA 분석
    print("[Trajectory] Running PAGA...")
    print(f"[DEBUG] n_clusters: {adata_traj.obs[cluster_key].nunique()}")
    print(f"[DEBUG] cluster distribution:\n{adata_traj.obs[cluster_key].value_counts()}")
    print(f"[DEBUG] connectivities shape: {adata_traj.obsp['connectivities'].shape}")
    print(f"[DEBUG] connectivities nnz: {adata_traj.obsp['connectivities'].nnz}")
    sc.tl.paga(adata_traj, groups=cluster_key)

    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.paga(adata_traj, threshold=paga_threshold, show=False, ax=ax, fontsize=10, node_size_scale=1.5)
    save_figure(fig, save_path, "paga_graph.png")

    # PAGA-initialized UMAP
    sc.tl.umap(adata_traj, init_pos="paga")

    # PAGA + UMAP 오버레이
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sc.pl.paga(adata_traj, threshold=paga_threshold, show=False, ax=axes[0], fontsize=10, node_size_scale=1.5, title="PAGA Graph")
    sc.pl.umap(adata_traj, color=cluster_key, show=False, ax=axes[1], title="PAGA-initialized UMAP")
    save_figure(fig, save_path, "paga_umap_overlay.png")

    # 3. Diffusion Map
    print("[Trajectory] Computing diffusion map...")
    sc.tl.diffmap(adata_traj, n_comps=n_dcs)

    # Diffusion Components 시각화
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    for i, ax in enumerate(axes.flatten()):
        if i < min(4, n_dcs - 1):
            sc.pl.embedding(adata_traj, basis="diffmap", color=cluster_key,
                            components=[1, i + 2], ax=ax, show=False,
                            title=f"DC1 vs DC{i + 2}")
    save_figure(fig, save_path, "diffusion_components.png")

    # 4. Diffusion Pseudotime
    print("[Trajectory] Computing diffusion pseudotime...")
    root_idx = _select_root_cell(adata_traj, root_cluster, root_gene, cluster_key)
    adata_traj.uns["iroot"] = root_idx
    print(f"[Trajectory] Root cell: {root_idx} (cluster: {adata_traj.obs[cluster_key].iloc[root_idx]})")

    sc.tl.dpt(adata_traj, n_dcs=n_dcs)

    # 5. 시각화
    # Pseudotime UMAP
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata_traj, color="dpt_pseudotime", ax=ax, show=False, color_map="viridis", title="Diffusion Pseudotime")
    save_figure(fig, save_path, "diffusion_pseudotime_umap.png")

    # 클러스터별 pseudotime 분포
    fig_width = max(8, n_clusters * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    sc.pl.violin(adata_traj, keys="dpt_pseudotime", groupby=cluster_key, ax=ax, show=False, rotation=45)
    ax.set_title("Pseudotime Distribution by Cluster")
    save_figure(fig, save_path, "pseudotime_by_cluster.png")

    # Root cell 위치 표시
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata_traj, color=cluster_key, ax=ax, show=False)
    root_coords = adata_traj.obsm["X_umap"][root_idx]
    ax.scatter(root_coords[0], root_coords[1], c="red", s=200, marker="*", edgecolors="black", linewidths=2, zorder=10)
    ax.set_title(f"Root Cell (Cluster: {adata_traj.obs[cluster_key].iloc[root_idx]})")
    save_figure(fig, save_path, "root_cell_location.png")

    # 6. 요약 테이블
    summary = adata_traj.obs.groupby(cluster_key)["dpt_pseudotime"].agg(["mean", "std", "min", "max"])
    summary = summary.sort_values("mean")
    summary.to_csv(os.path.join(save_path, "trajectory_summary.csv"))

    # 7. 결과를 원본 adata에 복사
    adata.obs[cluster_key] = adata_traj.obs[cluster_key]
    adata.obs["dpt_pseudotime"] = adata_traj.obs["dpt_pseudotime"]
    adata.obsm["X_diffmap"] = adata_traj.obsm["X_diffmap"]
    adata.obsm["X_umap_traj"] = adata_traj.obsm["X_umap"]
    adata.uns["paga"] = adata_traj.uns["paga"]
    adata.uns["iroot"] = adata_traj.uns["iroot"]

    print(f"[Trajectory] Complete. Results: {save_path}")
    return adata
