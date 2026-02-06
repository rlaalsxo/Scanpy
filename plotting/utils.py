"""
시각화 유틸리티

matplotlib 설정 및 공통 plotting 함수
"""
import matplotlib
matplotlib.use("Agg")

import os
import math
import matplotlib.pyplot as plt

from config.defaults import PLOTTING


def save_figure(fig, save_path: str, filename: str, dpi: int = None, close: bool = True):
    """
    Figure 저장 (tight_layout + savefig + close 통합)

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    save_path : str
        저장 디렉토리
    filename : str
        파일명 (확장자 포함)
    dpi : int, optional
        해상도 (기본: PLOTTING["dpi"])
    close : bool
        저장 후 figure 닫을지 여부
    """
    if dpi is None:
        dpi = PLOTTING["dpi"]

    os.makedirs(save_path, exist_ok=True)
    fig.tight_layout()
    fig.savefig(os.path.join(save_path, filename), dpi=dpi, bbox_inches="tight")

    if close:
        plt.close(fig)


def calc_figsize(
    n_genes: int,
    n_clusters: int,
    base_size: tuple = (3.5, 3.5),
    gene_width: float = 0.7,
    cluster_height: float = 0.30,
    min_size: tuple = (4.5, 4.0),
    max_size: tuple = (18.0, 18.0),
) -> tuple:
    """
    Gene/Cluster 수에 따른 동적 figsize 계산

    Parameters
    ----------
    n_genes : int
    n_clusters : int
    base_size : tuple
    gene_width : float
        gene 하나 늘어날 때 가로 증가량
    cluster_height : float
        cluster 하나 늘어날 때 세로 증가량
    min_size : tuple
    max_size : tuple

    Returns
    -------
    tuple: (width, height)
    """
    # gene 수가 적을 때는 가로를 너무 넓히지 않음
    gene_width_eff = gene_width * 0.6 if n_genes <= 8 else gene_width

    # cluster 수가 적을 때는 세로를 너무 길게 만들지 않음
    cluster_height_eff = cluster_height * 0.7 if n_clusters <= 8 else cluster_height

    # 가로: gene 수에 비례
    width = base_size[0] + gene_width_eff * max(n_genes - 1, 0)
    # 세로: cluster 수에 비례
    height = base_size[1] + cluster_height_eff * max(n_clusters - 1, 0)

    # 최소/최대 사이즈 제한
    width = max(min_size[0], min(max_size[0], width))
    height = max(min_size[1], min(max_size[1], height))

    return (width, height)


def _chunk_list(lst, chunk_size):
    """리스트를 chunk_size 단위로 나눔"""
    for i in range(0, len(lst), chunk_size):
        yield lst[i : i + chunk_size]


def plot_gene_cluster_panels(
    plot_func,
    adata,
    genes: list,
    groupby: str,
    save_path: str,
    base_name: str,
    n_clusters: int,
    genes_per_fig: int = None,
):
    """
    dotplot / matrixplot / stacked_violin / heatmap 공통 처리

    - genes_per_fig 개 단위로 여러 장으로 나눔
    - 각 장마다 gene/cluster 수에 맞게 figsize 계산
    - bbox_inches="tight" 로 저장해서 여백 제거

    Parameters
    ----------
    plot_func : callable
        scanpy plotting function (sc.pl.dotplot, sc.pl.matrixplot, etc.)
    adata : AnnData
    genes : list
    groupby : str
    save_path : str
    base_name : str
    n_clusters : int
    genes_per_fig : int, optional
    """
    if genes_per_fig is None:
        genes_per_fig = PLOTTING["genes_per_fig"]

    n_genes_total = len(genes)
    if n_genes_total == 0:
        return

    n_chunks = math.ceil(n_genes_total / genes_per_fig)

    for idx, genes_chunk in enumerate(_chunk_list(genes, genes_per_fig), start=1):
        n_genes_chunk = len(genes_chunk)
        fig_w, fig_h = calc_figsize(n_genes=n_genes_chunk, n_clusters=n_clusters)

        plot_func(
            adata,
            var_names=genes_chunk,
            groupby=groupby,
            use_raw=False,
            show=False,
            figsize=(fig_w, fig_h),
        )
        fig = plt.gcf()

        if n_chunks == 1:
            fname = f"{base_name}.png"
        else:
            fname = f"{base_name}_{idx}.png"

        save_figure(fig, save_path, fname)
