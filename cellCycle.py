import os
import matplotlib.pyplot as plt
import scanpy as sc

def score_cell_cycle(adata, save_path):

    os.makedirs(save_path, exist_ok=True)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    cell_cycle_gene_path = os.path.join(script_dir, "input", "regev_lab_cell_cycle_genes.txt")

    with open(cell_cycle_gene_path, 'r') as f:
        cell_cycle_genes = [line.strip() for line in f.readlines()]

    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]

    # 존재하는 유전자만 필터링
    s_genes_filtered = [gene for gene in s_genes if gene in adata.var_names]
    g2m_genes_filtered = [gene for gene in g2m_genes if gene in adata.var_names]

    print(f"S genes used: {len(s_genes_filtered)} / {len(s_genes)}")
    print(f"G2M genes used: {len(g2m_genes_filtered)} / {len(g2m_genes)}")

    # 정규화된 데이터로 스케일링
    sc.pp.scale(adata, max_value=10)

    # 세포 주기 점수 계산
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_filtered, g2m_genes=g2m_genes_filtered)

    # 바이올린 플롯 저장
    sc.pl.violin(adata, ["S_score", "G2M_score"], jitter=0.4, groupby="sample", rotation=60, show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "violin_cell_cycle.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # 산점도 플롯 저장
    sc.pl.scatter(adata, x="S_score", y="G2M_score", color="phase", show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "scatter_cell_cycle.png"), dpi=300, bbox_inches="tight")
    plt.close()

    return adata
