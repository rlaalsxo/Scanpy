import os
import matplotlib.pyplot as plt
import scanpy as sc
import gseapy
import pandas as pd

def deg_analysis_with_sex_gene_filtering(
    adata,
    save_path,
    groupby="leiden",
    top_n=20,
    padj_threshold=0.05
):
    os.makedirs(save_path, exist_ok=True)

    # 성염색체 유전자 제거
    print("Removing sex chromosome genes...")
    annot = sc.queries.biomart_annotations(
        "hsapiens",
        ["ensembl_gene_id", "external_gene_name", "chromosome_name"],
    ).set_index("external_gene_name")

    chrY_genes = adata.var_names.intersection(annot.index[annot.chromosome_name == "Y"])
    chrX_genes = adata.var_names.intersection(annot.index[annot.chromosome_name == "X"])
    sex_genes = chrY_genes.union(chrX_genes)
    adata = adata[:, [g for g in adata.var_names if g not in sex_genes]]

    print(f"Removed {len(sex_genes)} sex chromosome genes.")

    # DEG 분석
    print("Running DEG analysis...")
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", key_added="wilcoxon")
    deg_df = sc.get.rank_genes_groups_df(adata, group=None, key="wilcoxon")
    deg_df = deg_df[deg_df['pvals_adj'] < padj_threshold]
    top_genes = deg_df.groupby('group') \
                      .apply(lambda x: x.nlargest(top_n, 'logfoldchanges'))['names'] \
                      .unique().tolist()

    print(f"Selected {len(top_genes)} top genes.")

    # Cell Type 예측
    print("Predicting cell types via GSEA...")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    marker_path = os.path.join(script_dir, "input", "human_cell_markers.txt")
    marker_df = pd.read_table(marker_path)
    marker_df["nG"] = marker_df.geneSymbol.str.split(",").str.len()
    marker_df = marker_df[(marker_df["nG"] > 5) & (marker_df["nG"] < 100)]
    marker_df = marker_df[marker_df["cancerType"] == "Normal"]
    marker_df.index = marker_df.cellName
    gene_dict = marker_df.geneSymbol.str.split(",").to_dict()

    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", key_added="gsea_wilcoxon")
    pred = {}

    for cl in adata.obs[groupby].cat.categories:
        glist = (
            sc.get.rank_genes_groups_df(adata, group=cl, key="gsea_wilcoxon")["names"]
            .squeeze().str.strip().tolist()
        )
        enr_res = gseapy.enrichr(
            gene_list=glist[:300],
            organism="Human",
            gene_sets=gene_dict,
            background=adata.var_names,
            cutoff=1,
        )
        if enr_res.results.shape[0] == 0:
            pred[cl] = "Unassigned"
        else:
            enr_res.results.sort_values(by="P-value", inplace=True)
            pred[cl] = enr_res.results["Term"].iloc[0]

    adata.obs["cell_type"] = adata.obs[groupby].map(pred)
    print("✅ Cell type prediction complete.")

    # 시각화 저장
    print("Saving plots...")
    sc.pl.rank_genes_groups(adata, n_genes=top_n, key="wilcoxon", sharey=False, show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "rank_genes_groups.png"), dpi=300)
    plt.close()

    sc.pl.dotplot(adata, var_names=top_genes, groupby=groupby, use_raw=False, show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "dotplot.png"), dpi=300)
    plt.close()

    sc.pl.matrixplot(adata, var_names=top_genes, groupby=groupby, use_raw=False, show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "matrixplot.png"), dpi=300)
    plt.close()

    sc.pl.stacked_violin(adata, var_names=top_genes, groupby=groupby, use_raw=False, show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "stacked_violin.png"), dpi=300)
    plt.close()

    sc.pl.heatmap(adata, var_names=top_genes, groupby=groupby, use_raw=False, show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "heatmap.png"), dpi=300)
    plt.close()

    sc.pl.umap(adata, color="cell_type", show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "predicted_cell_types.png"), dpi=300)
    plt.close()

    print("✅ DEG + Cell Type Prediction 완료")
    return adata
