"""
Cell-Cell Communication 분석

Squidpy 기반 Ligand-Receptor 분석
"""
import os

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt

from config.defaults import CELL_COMMUNICATION
from core.neighbors import detect_cluster_key
from plotting.utils import save_figure


def _plot_interaction_network(interaction_counts, save_path: str):
    """세포 타입 간 상호작용 네트워크 시각화"""
    import networkx as nx

    G = nx.DiGraph()
    all_types = set(interaction_counts.index) | set(interaction_counts.columns)
    G.add_nodes_from(all_types)

    for source in interaction_counts.index:
        for target in interaction_counts.columns:
            weight = interaction_counts.loc[source, target]
            if weight > 0:
                G.add_edge(source, target, weight=weight)

    if len(G.edges()) == 0:
        print("[CellComm] No edges to plot in network.")
        return

    fig, ax = plt.subplots(figsize=(12, 12))
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    # 노드 크기
    node_sizes = []
    for node in G.nodes():
        total = 0
        if node in interaction_counts.index:
            total += interaction_counts.loc[node, :].sum()
        if node in interaction_counts.columns:
            total += interaction_counts.loc[:, node].sum()
        node_sizes.append(300 + total * 30)

    # 엣지 두께
    edge_weights = [G[u][v]["weight"] for u, v in G.edges()]
    max_weight = max(edge_weights) if edge_weights else 1
    edge_widths = [1 + (w / max_weight) * 5 for w in edge_weights]

    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color="lightblue", alpha=0.8, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=9, ax=ax)
    nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.6, edge_color="gray",
                           connectionstyle="arc3,rad=0.1", arrows=True, arrowsize=15, ax=ax)

    ax.set_title("Cell-Cell Communication Network")
    ax.axis("off")
    save_figure(fig, save_path, "interaction_network.png")


def cellcell_communication(
    adata,
    save_path: str,
    species: str = "human",
    groupby: str = None,
    n_perms: int = None,
    pvalue_threshold: float = None,
    top_n_interactions: int = None,
) -> sc.AnnData:
    """
    Cell-Cell Communication 분석

    Parameters
    ----------
    adata : AnnData
    save_path : str
    species : str
    groupby : str, optional
    n_perms : int, optional
    pvalue_threshold : float, optional
    top_n_interactions : int, optional

    Returns
    -------
    AnnData
    """
    os.makedirs(save_path, exist_ok=True)

    if not adata.uns.get("is_standard_symbol", True):
        print("[CellComm] 비표준 유전자 심볼 → Cell-Cell Communication 분석 스킵")
        return adata

    if n_perms is None:
        n_perms = CELL_COMMUNICATION["n_perms"]
    if pvalue_threshold is None:
        pvalue_threshold = CELL_COMMUNICATION["pvalue_threshold"]
    if top_n_interactions is None:
        top_n_interactions = CELL_COMMUNICATION["top_n_interactions"]

    # 1. groupby 확인
    if groupby is None:
        if "cell_type" in adata.obs:
            groupby = "cell_type"
        else:
            groupby = detect_cluster_key(adata)

    if groupby is None or groupby not in adata.obs:
        raise ValueError(f"No valid groupby key found. Run DEG or clustering first.")

    print(f"[CellComm] groupby: {groupby}")
    n_groups = adata.obs[groupby].nunique()
    print(f"[CellComm] Number of groups: {n_groups}")

    # 2. Ligand-Receptor 분석
    print(f"[CellComm] Running ligrec analysis (n_perms={n_perms})...")
    sq.gr.ligrec(adata, cluster_key=groupby, n_perms=n_perms, threshold=0.01, copy=False, use_raw=False)

    ligrec_key = f"{groupby}_ligrec"
    if ligrec_key not in adata.uns:
        print("[CellComm] ligrec results not found.")
        return adata

    # 3. 유의한 상호작용 추출
    print("[CellComm] Extracting significant interactions...")
    pvals = adata.uns[ligrec_key]["pvalues"]
    means = adata.uns[ligrec_key]["means"]

    interactions_list = []
    for lr_pair in pvals.index:
        for cell_pair in pvals.columns:
            pval = pvals.loc[lr_pair, cell_pair]
            mean_expr = means.loc[lr_pair, cell_pair]
            if pval < pvalue_threshold and not np.isnan(pval):
                source, target = cell_pair
                ligand, receptor = lr_pair if isinstance(lr_pair, tuple) else (lr_pair, lr_pair)
                interactions_list.append({
                    "ligand": ligand,
                    "receptor": receptor,
                    "source": source,
                    "target": target,
                    "pvalue": pval,
                    "mean_expression": mean_expr,
                })

    sig_df = pd.DataFrame(interactions_list)
    if sig_df.empty:
        print("[CellComm] No significant interactions found.")
        # 유의한 상호작용이 없어도 전체 dotplot 생성
        print("[CellComm] Generating dotplot for all interactions...")
        try:
            fig_width = max(10, n_groups * 1.5)
            sq.pl.ligrec(adata, cluster_key=groupby, pvalue_threshold=1.0,
                         remove_empty_interactions=True, show=False, figsize=(fig_width, 12))
            plt.tight_layout()
            plt.savefig(os.path.join(save_path, "ligrec_dotplot_all.png"), dpi=300, bbox_inches="tight")
            plt.close()
        except Exception as e:
            print(f"[CellComm] Dotplot failed: {e}")
        # ligrec 결과 삭제 (h5ad 저장 시 튜플 컬럼명 호환성 문제 방지)
        if ligrec_key in adata.uns:
            del adata.uns[ligrec_key]
        print(f"[CellComm] Complete. Results: {save_path}")
        return adata

    sig_df = sig_df.sort_values("pvalue")
    sig_df.to_csv(os.path.join(save_path, "significant_interactions.csv"), index=False)
    print(f"[CellComm] Found {len(sig_df)} significant interactions")

    # 4. 시각화
    # 4-1. Dotplot
    print("[CellComm] Generating dotplot...")
    try:
        fig_height = max(8, min(top_n_interactions * 0.3, 20))
        fig_width = max(10, n_groups * 1.5)
        sq.pl.ligrec(adata, cluster_key=groupby, pvalue_threshold=pvalue_threshold,
                     remove_empty_interactions=True, show=False, figsize=(fig_width, fig_height))
        plt.tight_layout()
        plt.savefig(os.path.join(save_path, "ligrec_dotplot.png"), dpi=300, bbox_inches="tight")
        plt.close()
    except Exception as e:
        print(f"[CellComm] Dotplot failed: {e}")

    # 4-2. 상호작용 히트맵
    print("[CellComm] Generating interaction heatmap...")
    interaction_counts = sig_df.groupby(["source", "target"]).size().unstack(fill_value=0)

    if not interaction_counts.empty:
        fig, ax = plt.subplots(figsize=(max(8, n_groups * 0.8), max(6, n_groups * 0.6)))
        im = ax.imshow(interaction_counts.values, cmap="Reds", aspect="auto")
        ax.set_xticks(range(len(interaction_counts.columns)))
        ax.set_yticks(range(len(interaction_counts.index)))
        ax.set_xticklabels(interaction_counts.columns, rotation=45, ha="right")
        ax.set_yticklabels(interaction_counts.index)
        ax.set_xlabel("Target")
        ax.set_ylabel("Source")
        ax.set_title("Number of Significant L-R Interactions")
        plt.colorbar(im, ax=ax, label="Count")
        save_figure(fig, save_path, "interaction_heatmap.png")

        # 4-3. 네트워크 그래프
        print("[CellComm] Generating network plot...")
        _plot_interaction_network(interaction_counts, save_path)

    # 4-4. 상위 L-R pairs
    print("[CellComm] Generating top L-R pairs plot...")
    top_pairs = sig_df.groupby(["ligand", "receptor"]).size().sort_values(ascending=False).head(20)

    if not top_pairs.empty:
        fig, ax = plt.subplots(figsize=(10, 8))
        top_pairs.plot(kind="barh", ax=ax)
        ax.set_xlabel("Number of Cell Type Pairs")
        ax.set_ylabel("Ligand-Receptor Pair")
        ax.set_title("Top 20 Ligand-Receptor Pairs")
        ax.invert_yaxis()
        save_figure(fig, save_path, "top_lr_pairs.png")

    # 5. 요약 테이블
    summary = pd.DataFrame({
        "outgoing": sig_df.groupby("source").size(),
        "incoming": sig_df.groupby("target").size(),
    }).fillna(0).astype(int)
    summary["total"] = summary["outgoing"] + summary["incoming"]
    summary = summary.sort_values("total", ascending=False)
    summary.to_csv(os.path.join(save_path, "communication_summary.csv"))

    # 6. ligrec 결과 삭제 (h5ad 저장 시 튜플 컬럼명 호환성 문제 방지)
    if ligrec_key in adata.uns:
        del adata.uns[ligrec_key]

    print(f"[CellComm] Complete. Results: {save_path}")
    return adata
