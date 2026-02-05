import matplotlib
matplotlib.use("Agg")
import os
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np

from trajectory import _detect_cluster_key


def _validate_groupby(adata, groupby):
    """groupby 컬럼 검증 및 대체"""
    if groupby in adata.obs:
        return groupby
    # cell_type 우선, 없으면 클러스터 키
    if "cell_type" in adata.obs:
        return "cell_type"
    detected = _detect_cluster_key(adata)
    if detected:
        return detected
    raise ValueError("No valid groupby column found (cell_type, leiden, louvain, clusters)")


def _plot_interaction_network(interaction_counts, save_path):
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

    # 노드 크기: 총 상호작용 수에 비례
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

    nx.draw_networkx_nodes(G, pos, node_size=node_sizes,
                           node_color="lightblue", alpha=0.8, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=9, ax=ax)
    nx.draw_networkx_edges(G, pos, width=edge_widths,
                           alpha=0.6, edge_color="gray",
                           connectionstyle="arc3,rad=0.1",
                           arrows=True, arrowsize=15, ax=ax)

    ax.set_title("Cell-Cell Communication Network")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(save_path, "interaction_network.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def cellcell_communication(adata, save_path, species="human",
                           groupby="cell_type", n_perms=1000,
                           pvalue_threshold=0.05, top_n_interactions=50):
    """
    Cell-Cell Communication 분석 (squidpy 기반)

    Parameters
    ----------
    adata : AnnData
    save_path : str
    species : str
    groupby : str (cell_type, leiden 등)
    n_perms : int (permutation 수)
    pvalue_threshold : float
    top_n_interactions : int

    Returns
    -------
    AnnData with communication results
    """
    os.makedirs(save_path, exist_ok=True)

    # 1. groupby 검증
    groupby = _validate_groupby(adata, groupby)
    print(f"[CellComm] groupby: {groupby}")

    n_groups = adata.obs[groupby].nunique()
    print(f"[CellComm] Number of groups: {n_groups}")

    # 2. Ligand-Receptor 분석
    print(f"[CellComm] Running ligrec analysis (n_perms={n_perms})...")
    sq.gr.ligrec(
        adata,
        cluster_key=groupby,
        n_perms=n_perms,
        threshold=0.01,
        copy=False,
    )

    ligrec_key = f"{groupby}_ligrec"
    if ligrec_key not in adata.uns:
        print(f"[CellComm] ligrec results not found. Check squidpy installation.")
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
        sq.pl.ligrec(
            adata,
            cluster_key=groupby,
            pvalue_threshold=pvalue_threshold,
            remove_empty_interactions=True,
            show=False,
            figsize=(fig_width, fig_height),
        )
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
        fig.tight_layout()
        fig.savefig(os.path.join(save_path, "interaction_heatmap.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

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
        fig.tight_layout()
        fig.savefig(os.path.join(save_path, "top_lr_pairs.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    # 5. 요약 테이블
    summary = pd.DataFrame({
        "outgoing": sig_df.groupby("source").size(),
        "incoming": sig_df.groupby("target").size(),
    }).fillna(0).astype(int)
    summary["total"] = summary["outgoing"] + summary["incoming"]
    summary = summary.sort_values("total", ascending=False)
    summary.to_csv(os.path.join(save_path, "communication_summary.csv"))

    print(f"[CellComm] Complete. Results: {save_path}")
    return adata
