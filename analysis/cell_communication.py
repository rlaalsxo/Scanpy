"""
Cell-Cell Communication 분석

Squidpy 기반 Ligand-Receptor 분석
"""
import os
import re
from math import ceil

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt

from config.defaults import CELL_COMMUNICATION
from core.neighbors import detect_cluster_key
from plotting.utils import save_figure


def _truncate_name(name: str, max_len: int = 30) -> str:
    """셀 타입명에서 ontology ID 제거 + 길이 제한"""
    name = re.sub(r"\s*(CL|UBERON|EFO):\S+", "", str(name)).strip()
    return name[:max_len] + "..." if len(name) > max_len else name


def _plot_ligrec_dotplot(sig_df, save_path: str, top_n: int = 50, filename: str = "ligrec_dotplot.png"):
    """상위 N개 L-R 상호작용 커스텀 dotplot"""
    top = sig_df.nsmallest(top_n, "pvalue").copy()
    top["lr_pair"] = top["ligand"] + " → " + top["receptor"]
    top["cell_pair"] = top["source"].apply(_truncate_name) + " → " + top["target"].apply(_truncate_name)

    lr_pairs = top["lr_pair"].unique()
    cell_pairs = top["cell_pair"].unique()

    if len(lr_pairs) == 0 or len(cell_pairs) == 0:
        return

    lr_to_idx = {lr: i for i, lr in enumerate(lr_pairs)}
    cp_to_idx = {cp: i for i, cp in enumerate(cell_pairs)}

    fig_w = max(10, len(cell_pairs) * 0.8 + 3)
    fig_h = max(6, len(lr_pairs) * 0.25 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    neg_log_p = -np.log10(top["pvalue"].clip(lower=1e-300))
    mean_vals = top["mean_expression"]
    size_min, size_max = 20, 300
    if mean_vals.max() > mean_vals.min():
        sizes = size_min + (mean_vals - mean_vals.min()) / (mean_vals.max() - mean_vals.min()) * (size_max - size_min)
    else:
        sizes = (size_min + size_max) / 2

    sc_plot = ax.scatter(
        [cp_to_idx[cp] for cp in top["cell_pair"]],
        [lr_to_idx[lr] for lr in top["lr_pair"]],
        s=sizes,
        c=neg_log_p,
        cmap="viridis",
        edgecolors="black",
        linewidths=0.3,
        alpha=0.85,
    )

    ax.set_xticks(range(len(cell_pairs)))
    ax.set_xticklabels(cell_pairs, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(lr_pairs)))
    ax.set_yticklabels(lr_pairs, fontsize=7)
    ax.set_xlabel("Cell Type Pair (Source → Target)")
    ax.set_title(f"Top {len(top)} Ligand-Receptor Interactions")

    cbar = plt.colorbar(sc_plot, ax=ax, label="-log10(p-value)", shrink=0.6)
    cbar.ax.tick_params(labelsize=7)

    # 사이즈 범례
    for s_val, s_label in [(size_min, "Low"), (size_max, "High")]:
        ax.scatter([], [], s=s_val, c="gray", edgecolors="black", linewidths=0.3, label=f"Expr: {s_label}")
    ax.legend(loc="upper left", bbox_to_anchor=(1.15, 0.3), fontsize=7, frameon=False)

    save_figure(fig, save_path, filename)


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
    # 4-1. 커스텀 Dotplot (페이지 분할)
    print("[CellComm] Generating dotplots...")
    sorted_df = sig_df.sort_values("pvalue")
    n_pages = ceil(len(sorted_df) / top_n_interactions)
    for page in range(n_pages):
        chunk = sorted_df.iloc[page * top_n_interactions : (page + 1) * top_n_interactions]
        filename = f"ligrec_dotplot_{page + 1}.png"
        _plot_ligrec_dotplot(chunk, save_path, top_n=len(chunk), filename=filename)
    print(f"[CellComm] Generated {n_pages} dotplot page(s)")

    # 4-2. 상호작용 히트맵
    print("[CellComm] Generating interaction heatmap...")
    interaction_counts = sig_df.groupby(["source", "target"]).size().unstack(fill_value=0)

    if not interaction_counts.empty:
        trunc_cols = [_truncate_name(c) for c in interaction_counts.columns]
        trunc_idx = [_truncate_name(i) for i in interaction_counts.index]
        fig, ax = plt.subplots(figsize=(max(8, n_groups * 0.8), max(6, n_groups * 0.6)))
        im = ax.imshow(interaction_counts.values, cmap="Reds", aspect="auto")
        ax.set_xticks(range(len(trunc_cols)))
        ax.set_yticks(range(len(trunc_idx)))
        ax.set_xticklabels(trunc_cols, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(trunc_idx, fontsize=8)
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
