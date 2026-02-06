"""
공통 유틸리티 모듈

scRNA-seq 파이프라인 전체에서 사용하는 공통 함수들.
전처리 상태 확인 및 생성 로직을 한 곳에서 관리.
"""
import scanpy as sc


def detect_cluster_key(adata):
    """
    클러스터 키 자동 탐지 (leiden > louvain > clusters > cluster)

    Parameters
    ----------
    adata : AnnData

    Returns
    -------
    str or None
        발견된 클러스터 키, 없으면 None
    """
    for cand in ["leiden", "louvain", "clusters", "cluster"]:
        if cand in adata.obs:
            return cand
    return None


def ensure_pca(adata, n_comps=50):
    """
    PCA 확보 (없으면 생성)

    Parameters
    ----------
    adata : AnnData
    n_comps : int
        PCA 컴포넌트 수
    """
    if "X_pca" not in adata.obsm:
        print("[Common] Computing PCA...")
        sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")


def ensure_neighbors(adata, n_neighbors=15, n_pcs=40, key_added=None):
    """
    Neighbors 확보

    Parameters
    ----------
    adata : AnnData
    n_neighbors : int
    n_pcs : int
    key_added : str, optional
        - None: 기본 neighbors (없으면 생성)
        - "neighbors_paga" 등: 별도 키로 저장 (항상 새로 생성)
    """
    ensure_pca(adata)
    if key_added:
        # 별도 키는 항상 새로 생성
        print(f"[Common] Computing neighbors (key={key_added})...")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, key_added=key_added)
    else:
        # 기본 neighbors는 없을 때만 생성
        if "neighbors" not in adata.uns:
            print("[Common] Computing neighbors...")
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)


def ensure_clustering(adata, method="leiden", neighbors_key=None):
    """
    클러스터링 확보 (없으면 생성)

    Parameters
    ----------
    adata : AnnData
    method : str
        "leiden" 또는 "louvain"
    neighbors_key : str, optional
        사용할 neighbors 키

    Returns
    -------
    str
        클러스터 키 이름
    """
    cluster_key = detect_cluster_key(adata)
    if cluster_key:
        return cluster_key

    print(f"[Common] Computing {method} clustering...")
    ensure_neighbors(adata, key_added=neighbors_key)

    if method == "leiden":
        sc.tl.leiden(adata, key_added="leiden", neighbors_key=neighbors_key)
    elif method == "louvain":
        sc.tl.louvain(adata, key_added="louvain", neighbors_key=neighbors_key)
    else:
        raise ValueError(f"Unknown clustering method: {method}")

    return method


def ensure_umap(adata, neighbors_key=None):
    """
    UMAP 확보 (없으면 생성)

    Parameters
    ----------
    adata : AnnData
    neighbors_key : str, optional
        사용할 neighbors 키
    """
    if "X_umap" not in adata.obsm:
        print("[Common] Computing UMAP...")
        ensure_neighbors(adata, key_added=neighbors_key)
        sc.tl.umap(adata, neighbors_key=neighbors_key)
