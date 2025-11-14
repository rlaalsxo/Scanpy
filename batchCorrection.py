import matplotlib
matplotlib.use("Agg")
import os
import matplotlib.pyplot as plt
import scanpy as sc
import bbknn

def BatchCorrection(adata, save_path, neighbors_within_batch=None):
    os.makedirs(save_path, exist_ok=True)

    if neighbors_within_batch is None:
        if adata.n_obs > 100_000: neighbors_within_batch = 20
        elif adata.n_obs > 50_000: neighbors_within_batch = 15
        elif adata.n_obs > 10_000: neighbors_within_batch = 10
        else: neighbors_within_batch = 3
        print(f"[BBKNN] neighbors_within_batch 자동 설정됨: {neighbors_within_batch}")

    sc.tl.pca(adata, svd_solver="arpack")

    adata_bbknn = adata.copy()
    bbknn.bbknn(adata_bbknn, batch_key="sample", neighbors_within_batch=neighbors_within_batch)
    sc.tl.leiden(adata_bbknn, key_added="leiden")
    print("[Leiden] 클러스터링 완료. adata.obs['leiden'] 생성됨.")

    sc.tl.umap(adata_bbknn)

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    sc.pl.umap(adata, color="sample", title="Uncorrected UMAP", show=False, ax=axs[0])
    sc.pl.umap(adata_bbknn, color="sample", title="BBKNN Corrected UMAP", show=False, ax=axs[1])
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "umap_batch_correction_comparison.png"), dpi=300)
    plt.close()

    sc.pl.umap(adata_bbknn, color="leiden", title="Leiden Clusters", show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "umap_leiden_clusters.png"), dpi=300)
    plt.close()

    return adata_bbknn