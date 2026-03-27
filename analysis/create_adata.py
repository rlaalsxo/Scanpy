"""
AnnData 생성 및 전처리

10x / BD Rhapsody 데이터 로딩, QC, 더블렛 제거
"""
import os

import pandas as pd
import numpy as np
import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt

from core.io import load_10x_data, load_bd_data, save_adata
from core.preprocessing import (
    check_standard_symbols, validate_species,
    calculate_qc_metrics, filter_cells_qc, remove_genes, normalize_and_hvg,
)
from core.neighbors import compute_pca, compute_neighbors
from plotting.utils import save_figure


def create_adata(
    parent_dir: str,
    save_path: str,
    species: str,
    platform: str = "10x",
    sample_names: list = None,
    obs_name_style: str = "folder_barcode",
    min_genes: int = None,
    max_pct_mt: float = None,
    min_pct_ribo: float = None,
    perform_analysis: bool = True,
) -> sc.AnnData:
    """
    AnnData 생성 및 전처리

    Parameters
    ----------
    parent_dir : str
        데이터 디렉토리 또는 tar 파일
    save_path : str
    species : str
        "human" 또는 "mouse"
    platform : str
        "10x" 또는 "bd"
    sample_names : list, optional
    obs_name_style : str
    min_genes : int, optional
    max_pct_mt : float, optional
    min_pct_ribo : float, optional
    perform_analysis : bool
    """
    result_path = os.path.join(save_path, "CreateAdata")
    os.makedirs(result_path, exist_ok=True)

    # 1. 데이터 로딩
    if platform == "bd":
        print("[CreateAdata] Loading BD Rhapsody data...")
        adata = load_bd_data(parent_dir, sample_names)
    else:
        print("[CreateAdata] Loading 10x data...")
        adata = load_10x_data(parent_dir, sample_names, obs_name_style)

    # 2. 표준 심볼 판별 및 종 검증
    is_standard = check_standard_symbols(adata, species)
    adata.uns["is_standard_symbol"] = is_standard

    if is_standard:
        print("[CreateAdata] Validating species...")
        validate_species(adata, species)

    # 3. QC metrics 계산
    print("[CreateAdata] Calculating QC metrics...")
    calculate_qc_metrics(adata, species, is_standard=is_standard)

    # 4. QC 기반 필터링
    print("[CreateAdata] Filtering cells...")
    adata = filter_cells_qc(adata, species, min_genes, max_pct_mt, min_pct_ribo, is_standard=is_standard)

    # 5. 불필요한 유전자 제거
    print("[CreateAdata] Removing unwanted genes...")
    adata = remove_genes(adata, species, is_standard=is_standard)

    # 6. 더블렛 제거 전 시각화
    print("[CreateAdata] Pre-doublet visualization...")
    ad_pre = adata.copy()
    normalize_and_hvg(ad_pre, batch_key="sample")
    compute_pca(ad_pre)
    compute_neighbors(ad_pre, n_neighbors=10)
    sc.tl.umap(ad_pre)

    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(ad_pre, color="sample", show=False, ax=ax)
    save_figure(fig, result_path, "umap_before_doublet.png")

    # 7. 더블렛 탐지 및 제거
    print("[CreateAdata] Detecting doublets...")
    dbl_score = pd.Series(index=adata.obs_names, dtype=float)
    dbl_pred = pd.Series(index=adata.obs_names, dtype=bool)

    for s in adata.obs["sample"].unique():
        mask = adata.obs["sample"] == s
        sub = adata[mask]
        scrub = scr.Scrublet(sub.X)
        score, pred = scrub.scrub_doublets()
        dbl_score.loc[sub.obs_names] = score
        dbl_pred.loc[sub.obs_names] = pred

    adata.obs["doublet_score"] = dbl_score
    adata.obs["predicted_doublet"] = dbl_pred
    adata.obs["doublet_info"] = dbl_pred.astype(str)

    # 더블렛 violin plot
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.violin(adata, "n_genes_by_counts", groupby="doublet_info", size=0, rotation=45, show=False, ax=ax)
    save_figure(fig, result_path, "violin_doublet.png")

    # 더블렛 제거
    n_before = adata.n_obs
    adata = adata[adata.obs["doublet_info"] == "False", :]
    print(f"[CreateAdata] Removed {n_before - adata.n_obs} doublets")

    # 8. 후처리 및 시각화
    if perform_analysis:
        print("[CreateAdata] Post-processing...")
        normalize_and_hvg(adata, batch_key="sample")
        compute_pca(adata)
        compute_neighbors(adata, n_neighbors=10)
        sc.tl.umap(adata)

        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.umap(adata, color="sample", show=False, ax=ax)
        save_figure(fig, result_path, "umap_after_doublet.png")

    # 9. 메타데이터 정리
    keep_obs = ["sample"]
    keep_var = ["gene_ids", "feature_types"]
    adata.obs = adata.obs[[c for c in adata.obs.columns if c in keep_obs]]
    adata.var = adata.var[[c for c in adata.var.columns if c in keep_var]]
    adata.uns = {"is_standard_symbol": is_standard}
    adata.obsm = {}
    adata.varm = {}
    adata.obsp = {}

    # 10. 저장
    save_adata(adata, save_path, "scRNA_preprocessing.h5ad")

    print(f"[CreateAdata] Complete. Results: {result_path}")
    return adata
