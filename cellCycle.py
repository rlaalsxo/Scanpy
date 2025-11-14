import matplotlib
matplotlib.use("Agg")
import os
import matplotlib.pyplot as plt
import scanpy as sc

def _load_cycle_genes(species):
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "input", f"{species}_cell_cycle.txt")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Cell cycle file not found: {path}")
    with open(path, "r") as f:
        genes = [x.strip() for x in f.readlines() if x.strip()]
    return genes[:43], genes[43:]   # S, G2M 분리 규칙 유지

def score_cell_cycle(adata, save_path, species="human"):
    os.makedirs(save_path, exist_ok=True)

    s_genes, g2m_genes = _load_cycle_genes(species)
    s_genes = [g for g in s_genes if g in adata.var_names]
    g2m_genes = [g for g in g2m_genes if g in adata.var_names]

    print(f"S genes used: {len(s_genes)}")
    print(f"G2M genes used: {len(g2m_genes)}")

    sc.pp.scale(adata, max_value=10)
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    sc.pl.violin(adata, ["S_score", "G2M_score"], jitter=0.4, groupby="sample", rotation=60, show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "violin_cell_cycle.png"), dpi=300)
    plt.close()

    sc.pl.scatter(adata, x="S_score", y="G2M_score", color="phase", show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "scatter_cell_cycle.png"), dpi=300)
    plt.close()

    return adata
