"""
GSEA Enrichment 분석

DEG 결과 기반 GO/KEGG pathway enrichment
"""
import os

import numpy as np
import pandas as pd
import scanpy as sc
import gseapy
import matplotlib.pyplot as plt

from config.defaults import GSEA, DEG
from config.species import get_species_config
from core.neighbors import detect_cluster_key
from plotting.utils import save_figure


def gsea_enrichment(
    adata,
    save_path: str,
    species: str = "human",
    padj_threshold: float = None,
    top_n_pathways: int = None,
) -> sc.AnnData:
    """
    GSEA Enrichment 분석

    Parameters
    ----------
    adata : AnnData
    save_path : str
    species : str
    padj_threshold : float, optional
    top_n_pathways : int, optional

    Returns
    -------
    AnnData
    """
    os.makedirs(save_path, exist_ok=True)

    if not adata.uns.get("is_standard_symbol", True):
        print("[GSEA] 비표준 유전자 심볼 → GSEA Enrichment 분석 스킵")
        return adata

    if padj_threshold is None:
        padj_threshold = GSEA["padj_threshold"]
    if top_n_pathways is None:
        top_n_pathways = GSEA["top_n_pathways"]

    config = get_species_config(species)
    gene_sets = config.get("gsea_gene_sets", ["GO_Biological_Process_2023", "KEGG_2021_Human"])
    organism = "human" if species == "human" else "mouse"

    # 1. DEG 결과 확인
    deg_key = "wilcoxon"
    if deg_key not in adata.uns:
        raise ValueError("[GSEA] DEG results not found. Run DEG analysis (Step 4) first.")

    groupby = None
    if "cell_type" in adata.obs:
        groupby = detect_cluster_key(adata)
    else:
        groupby = detect_cluster_key(adata)

    if groupby is None:
        raise ValueError("[GSEA] No cluster key found in adata.obs.")

    clusters = adata.obs[groupby].astype("category").cat.categories
    print(f"[GSEA] groupby: {groupby}, clusters: {len(clusters)}")
    print(f"[GSEA] Gene sets: {gene_sets}")

    # 2. 클러스터별 Enrichr 실행
    all_results = []
    for cl in clusters:
        print(f"[GSEA] Processing cluster {cl}...")
        glist = (
            sc.get.rank_genes_groups_df(adata, group=cl, key=deg_key)
            .query("pvals_adj < @padj_threshold")["names"]
            .dropna()
            .astype(str)
            .str.strip()
            .tolist()
        )

        if len(glist) == 0:
            print(f"[GSEA] Cluster {cl}: no significant DEGs, skipping.")
            continue

        try:
            enr = gseapy.enrichr(
                gene_list=glist[:300],
                gene_sets=gene_sets,
                organism=organism,
                background=adata.var_names,
                cutoff=1.0,
                no_plot=True,
                outdir=None,
                verbose=False,
            )
            res = getattr(enr, "results", None)
            if res is not None and len(res) > 0:
                res = res.copy()
                res["cluster"] = cl
                all_results.append(res)
        except Exception as e:
            print(f"[GSEA] Cluster {cl} failed: {e}")

    if not all_results:
        print("[GSEA] No enrichment results found.")
        return adata

    results_df = pd.concat(all_results, ignore_index=True)
    results_df = results_df.sort_values("Adjusted P-value")
    print(f"[GSEA] Result columns: {list(results_df.columns)}")

    # 3. 결과 저장
    results_df.to_csv(os.path.join(save_path, "gsea_results.csv"), index=False)
    print(f"[GSEA] Total enrichment results: {len(results_df)}")

    # 유의한 결과만 필터
    sig_results = results_df[results_df["Adjusted P-value"] < padj_threshold].copy()
    print(f"[GSEA] Significant results (padj < {padj_threshold}): {len(sig_results)}")

    # 4. 시각화
    # 4-1. 클러스터별 상위 pathway barplot
    print("[GSEA] Generating top pathways barplot...")
    fig, axes = plt.subplots(
        len(clusters), 1,
        figsize=(12, max(4, len(clusters) * 3)),
        squeeze=False,
    )
    for idx, cl in enumerate(clusters):
        ax = axes[idx, 0]
        cl_data = sig_results[sig_results["cluster"] == cl].head(top_n_pathways)
        if cl_data.empty:
            ax.text(0.5, 0.5, "No significant pathways", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(f"Cluster {cl}")
            ax.set_yticks([])
            continue
        cl_data = cl_data.sort_values("Adjusted P-value", ascending=True)
        bars = ax.barh(
            range(len(cl_data)),
            -np.log10(cl_data["Adjusted P-value"].values),
        )
        ax.set_yticks(range(len(cl_data)))
        ax.set_yticklabels(cl_data["Term"].values, fontsize=7)
        ax.set_xlabel("-log10(Adjusted P-value)")
        ax.set_title(f"Cluster {cl}")
        ax.invert_yaxis()
    plt.tight_layout()
    save_figure(fig, save_path, "gsea_top_pathways.png")

    # 4-2. 클러스터 × pathway 히트맵
    if not sig_results.empty:
        print("[GSEA] Generating pathway heatmap...")
        top_terms = sig_results.groupby("Term")["Adjusted P-value"].min().nsmallest(top_n_pathways * 2).index
        heatmap_data = sig_results[sig_results["Term"].isin(top_terms)].pivot_table(
            index="Term", columns="cluster", values="Adjusted P-value", aggfunc="min"
        )
        heatmap_data = -np.log10(heatmap_data.fillna(1.0))

        fig, ax = plt.subplots(figsize=(max(8, len(clusters) * 0.8), max(6, len(top_terms) * 0.3)))
        im = ax.imshow(heatmap_data.values, cmap="Reds", aspect="auto")
        ax.set_xticks(range(len(heatmap_data.columns)))
        ax.set_yticks(range(len(heatmap_data.index)))
        ax.set_xticklabels(heatmap_data.columns, rotation=45, ha="right")
        ax.set_yticklabels(heatmap_data.index, fontsize=7)
        ax.set_xlabel("Cluster")
        ax.set_title("Pathway Enrichment (-log10 Adjusted P-value)")
        plt.colorbar(im, ax=ax, label="-log10(padj)")
        save_figure(fig, save_path, "gsea_heatmap.png")

    # 4-3. Dotplot
    if not sig_results.empty:
        print("[GSEA] Generating dotplot...")
        top_terms = sig_results.groupby("Term")["Adjusted P-value"].min().nsmallest(top_n_pathways).index
        dot_data = sig_results[sig_results["Term"].isin(top_terms)].copy()

        # gene count 추출
        if "Overlap" in dot_data.columns:
            dot_data["gene_count"] = dot_data["Overlap"].apply(lambda x: int(x.split("/")[0]))
        elif "Genes" in dot_data.columns:
            dot_data["gene_count"] = dot_data["Genes"].apply(lambda x: len(str(x).split(";")) if pd.notna(x) else 0)
        else:
            dot_data["gene_count"] = 5
        dot_data["neg_log_padj"] = -np.log10(dot_data["Adjusted P-value"].clip(lower=1e-300))

        fig, ax = plt.subplots(figsize=(max(8, len(clusters) * 0.8), max(6, len(top_terms) * 0.4)))
        for _, row in dot_data.iterrows():
            cluster_idx = list(clusters).index(row["cluster"])
            term_idx = list(top_terms).index(row["Term"])
            ax.scatter(
                cluster_idx, term_idx,
                s=row["gene_count"] * 10,
                c=row["neg_log_padj"],
                cmap="Reds",
                vmin=0,
                vmax=dot_data["neg_log_padj"].max(),
                edgecolors="black",
                linewidths=0.5,
            )
        ax.set_xticks(range(len(clusters)))
        ax.set_xticklabels(clusters, rotation=45, ha="right")
        ax.set_yticks(range(len(top_terms)))
        ax.set_yticklabels(top_terms, fontsize=7)
        ax.set_xlabel("Cluster")
        ax.set_title("Pathway Enrichment Dotplot")

        # 컬러바
        sm = plt.cm.ScalarMappable(cmap="Reds", norm=plt.Normalize(0, dot_data["neg_log_padj"].max()))
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label="-log10(padj)")
        save_figure(fig, save_path, "gsea_dotplot.png")

    print(f"[GSEA] Complete. Results: {save_path}")
    return adata
