import matplotlib
matplotlib.use("Agg")
import os
import matplotlib.pyplot as plt
import scanpy as sc
import gseapy
import pandas as pd

def _biomart_org(species):
    if species == "human": return "hsapiens"
    if species == "mouse": return "mmusculus"
    raise ValueError(f"Invalid species: {species}")

def _enrichr_lib(species):
    if species == "human":
        return [
            "CellMarker_Augmented_2021",
            "PanglaoDB_Augmented_2021",
            "Descartes_Cell_Types_and_Tissue_2021",
            "Tabula_Sapiens",
            "Azimuth_Cell_Types_2021"
        ]
    if species == "mouse":
        return [
            "Tabula_Muris",
            "Mouse_Gene_Atlas",
            "PanglaoDB_Augmented_2021"
        ]
    raise ValueError(f"Invalid species: {species}")

def deg_analysis_with_sex_gene_filtering(adata, save_path, species, groupby="leiden", top_n=20, padj_threshold=0.05):
    os.makedirs(save_path, exist_ok=True)

    print("Removing sex chromosome genes...")
    org = _biomart_org(species)
    annot = sc.queries.biomart_annotations(
        org, ["ensembl_gene_id","external_gene_name","chromosome_name"]
    ).set_index("external_gene_name")

    chrY = adata.var_names.intersection(annot.index[annot.chromosome_name=="Y"])
    chrX = adata.var_names.intersection(annot.index[annot.chromosome_name=="X"])
    sex_genes = chrY.union(chrX)
    adata = adata[:, [g for g in adata.var_names if g not in sex_genes]]
    print(f"Removed {len(sex_genes)} sex chromosome genes.")

    print("Running DEG analysis...")
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", key_added="wilcoxon")
    deg_df = sc.get.rank_genes_groups_df(adata, group=None, key="wilcoxon")
    deg_df = deg_df[deg_df["pvals_adj"] < padj_threshold]
    if deg_df.empty:
        print("No significant DEGs found.")
        return adata

    top_genes = (
        deg_df.groupby("group")
        .apply(lambda x: x.nlargest(top_n, "logfoldchanges"))["names"]
        .unique()
        .tolist()
    )
    print(f"Selected {len(top_genes)} top genes.")

    print("Predicting cell types via Enrichr...")
    gene_sets = _enrichr_lib(species)
    organism = "human" if species=="human" else "mouse"

    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", key_added="gsea_wilcoxon")
    pred = {}

    for cl in adata.obs[groupby].cat.categories:
        glist = (
            sc.get.rank_genes_groups_df(adata, group=cl, key="gsea_wilcoxon")["names"]
            .dropna().astype(str).str.strip().tolist()
        )
        if len(glist)==0:
            pred[cl] = "Unassigned"
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
                verbose=False
            )
            res = getattr(enr, "results", None)
            if res is None or res.shape[0]==0:
                pred[cl] = "Unassigned"
                continue
            res = res.sort_values("P-value")
            pred[cl] = res["Term"].iloc[0]
        except:
            pred[cl] = "Unassigned"

    adata.obs["cell_type"] = adata.obs[groupby].map(pred)
    print("Cell type prediction complete.")

    sc.pl.rank_genes_groups(adata, n_genes=top_n, key="wilcoxon", sharey=False, show=False)
    plt.tight_layout(); plt.savefig(os.path.join(save_path,"rank_genes_groups.png"), dpi=300); plt.close()

    sc.pl.dotplot(adata, var_names=top_genes, groupby=groupby, use_raw=False, show=False)
    plt.tight_layout(); plt.savefig(os.path.join(save_path,"dotplot.png"), dpi=300); plt.close()

    sc.pl.matrixplot(adata, var_names=top_genes, groupby=groupby, use_raw=False, show=False)
    plt.tight_layout(); plt.savefig(os.path.join(save_path,"matrixplot.png"), dpi=300); plt.close()

    sc.pl.stacked_violin(adata, var_names=top_genes, groupby=groupby, use_raw=False, show=False)
    plt.tight_layout(); plt.savefig(os.path.join(save_path,"stacked_violin.png"), dpi=300); plt.close()

    sc.pl.heatmap(adata, var_names=top_genes, groupby=groupby, use_raw=False, show=False)
    plt.tight_layout(); plt.savefig(os.path.join(save_path,"heatmap.png"), dpi=300); plt.close()

    sc.pl.umap(adata, color="cell_type", show=False)
    plt.tight_layout(); plt.savefig(os.path.join(save_path,"predicted_cell_types.png"), dpi=300); plt.close()

    print("DEG + Cell Type Prediction 완료")
    return adata
