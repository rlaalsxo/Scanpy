"""
Neighbors 및 Clustering 관련 함수

PCA, neighbors, clustering 로직 통합
"""
import scanpy as sc

from config.defaults import PCA, NEIGHBORS


def compute_pca(adata, n_comps: int = None, svd_solver: str = None):
    """
    PCA 계산

    Parameters
    ----------
    adata : AnnData
    n_comps : int, optional
    svd_solver : str, optional
    """
    if n_comps is None:
        n_comps = PCA["n_comps"]
    if svd_solver is None:
        svd_solver = PCA["svd_solver"]

    sc.tl.pca(adata, n_comps=n_comps, svd_solver=svd_solver)


def compute_neighbors(
    adata,
    n_neighbors: int = None,
    n_pcs: int = None,
    key_added: str = None,
    use_rep: str = None,
):
    """
    Neighbors 계산

    Parameters
    ----------
    adata : AnnData
    n_neighbors : int, optional
    n_pcs : int, optional
    key_added : str, optional
        저장할 키 (None이면 기본 neighbors)
    use_rep : str, optional
        사용할 representation (기본: X_pca)
    """
    if n_neighbors is None:
        n_neighbors = NEIGHBORS["n_neighbors"]
    if n_pcs is None:
        n_pcs = NEIGHBORS["n_pcs"]

    kwargs = {
        "n_neighbors": n_neighbors,
        "n_pcs": n_pcs,
    }
    if key_added:
        kwargs["key_added"] = key_added
    if use_rep:
        kwargs["use_rep"] = use_rep

    sc.pp.neighbors(adata, **kwargs)


def compute_clustering(
    adata,
    method: str = "leiden",
    neighbors_key: str = None,
    key_added: str = None,
    resolution: float = 1.0,
):
    """
    Clustering 계산

    Parameters
    ----------
    adata : AnnData
    method : str
        "leiden" 또는 "louvain"
    neighbors_key : str, optional
        사용할 neighbors 키
    key_added : str, optional
        결과 저장할 키 (None이면 method 이름 사용)
    resolution : float
        클러스터링 해상도
    """
    if key_added is None:
        key_added = method

    kwargs = {
        "key_added": key_added,
        "resolution": resolution,
    }
    if neighbors_key:
        kwargs["neighbors_key"] = neighbors_key

    if method == "leiden":
        sc.tl.leiden(adata, **kwargs)
    elif method == "louvain":
        sc.tl.louvain(adata, **kwargs)
    else:
        raise ValueError(f"Unknown clustering method: {method}")


def detect_cluster_key(adata) -> str:
    """
    클러스터 키 자동 탐지

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


def ensure_pca(adata, n_comps: int = None):
    """PCA 확보 (없으면 생성)"""
    if "X_pca" not in adata.obsm:
        print("[Core] Computing PCA...")
        compute_pca(adata, n_comps=n_comps)


def ensure_neighbors(adata, key_added: str = None, n_neighbors: int = None, n_pcs: int = None):
    """
    Neighbors 확보

    Parameters
    ----------
    adata : AnnData
    key_added : str, optional
        None이면 기본 neighbors (없을 때만 생성)
        지정하면 별도 키로 항상 새로 생성
    """
    ensure_pca(adata)

    if key_added:
        print(f"[Core] Computing neighbors (key={key_added})...")
        compute_neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, key_added=key_added)
    else:
        if "neighbors" not in adata.uns:
            print("[Core] Computing neighbors...")
            compute_neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)


def ensure_umap(adata, neighbors_key: str = None):
    """UMAP 확보 (없으면 생성)"""
    if "X_umap" not in adata.obsm:
        print("[Core] Computing UMAP...")
        ensure_neighbors(adata)
        sc.tl.umap(adata, neighbors_key=neighbors_key)
