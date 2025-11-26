import matplotlib
matplotlib.use("Agg")
import os
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np

def _load_cycle_genes(species):
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "input", f"{species}_cell_cycle.txt")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Cell cycle file not found: {path}")
    with open(path, "r") as f:
        genes = [x.strip() for x in f.readlines() if x.strip()]
    # 앞 43개: S phase, 나머지: G2M phase
    return genes[:43], genes[43:]

def _plot_phase_composition(adata, groupby_key, save_path):
    """
    phase(=adata.obs['phase'])를 기준으로
    groupby_key(예: 'sample', 'leiden') 별 cell cycle 구성비를 요약하고,
    counts / fraction 테이블 + stacked bar plot을 저장합니다.
    """
    if "phase" not in adata.obs:
        print("`phase` not found in adata.obs. Skip cell cycle composition.")
        return
    if groupby_key not in adata.obs:
        print(f"`{groupby_key}` not found in adata.obs. Skip composition for this key.")
        return

    os.makedirs(save_path, exist_ok=True)

    # 교차표: row = group (sample/cluster), col = phase
    comp_counts = pd.crosstab(adata.obs[groupby_key], adata.obs["phase"])
    if comp_counts.empty:
        print(f"No data for `{groupby_key}` x `phase`. Skip.")
        return

    comp_frac = comp_counts.div(comp_counts.sum(axis=1), axis=0)

    # 테이블 저장
    comp_counts.to_csv(
        os.path.join(save_path, f"cell_cycle_composition_{groupby_key}_counts.csv")
    )
    comp_frac.to_csv(
        os.path.join(save_path, f"cell_cycle_composition_{groupby_key}_fraction.csv")
    )

    # stacked bar plot
    # phase 순서는 G1, S, G2M 우선, 그 외 phase가 있으면 뒤에 붙임
    phase_order = ["G1", "S", "G2M"]
    cols = [p for p in phase_order if p in comp_frac.columns] + [
        p for p in comp_frac.columns if p not in phase_order
    ]
    comp_frac = comp_frac[cols]

    x = np.arange(comp_frac.shape[0])
    # 그룹 수에 따라 가로 크기를 조금 늘림
    fig_width = max(6.0, comp_frac.shape[0] * 0.6)
    fig, ax = plt.subplots(figsize=(fig_width, 4.0))

    bottom = np.zeros_like(x, dtype=float)
    for ph in comp_frac.columns:
        vals = comp_frac[ph].values
        ax.bar(x, vals, bottom=bottom, label=ph)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(comp_frac.index, rotation=45, ha="right")
    ax.set_ylabel("Fraction of cells")
    ax.set_xlabel(groupby_key)
    ax.set_title(f"Cell cycle phase composition by {groupby_key}")
    ax.legend(title="Phase", bbox_to_anchor=(1.05, 1.0), loc="upper left")

    fig.tight_layout()
    fig.savefig(
        os.path.join(save_path, f"cell_cycle_composition_{groupby_key}_stacked_bar.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def _summarize_cell_cycle_composition(adata, save_path):
    """
    샘플 / 클러스터 기준 cell cycle 구성비 요약 및 플롯 생성.
    - sample: adata.obs['sample']가 있을 때만
    - cluster: 'leiden' → 'louvain' → 'clusters' → 'cluster' 순으로 존재 여부를 확인
    """
    # sample 기준
    if "sample" in adata.obs:
        _plot_phase_composition(adata, "sample", save_path)
    else:
        print("`sample` not found in adata.obs. Skip sample-level composition.")

    # cluster 기준 (자동 탐색)
    cluster_key = None
    for cand in ["leiden", "louvain", "clusters", "cluster"]:
        if cand in adata.obs:
            cluster_key = cand
            break

    if cluster_key is not None:
        _plot_phase_composition(adata, cluster_key, save_path)
    else:
        print(
            "No cluster key found for cell cycle composition "
            "(tried: 'leiden', 'louvain', 'clusters', 'cluster'). "
            "Skip cluster-level composition."
        )


def _compare_umap_cell_cycle_effect(adata, save_path, color="phase"):
    """
    cell cycle 효과 제거 전후 UMAP 비교 플롯 생성.
    - 전(U) : adata에 기존에 존재하는 X_umap (있을 경우)
    - 후(U) : S_score, G2M_score를 regress_out 한 AnnData에서 새로 계산한 UMAP
    - 결과:
        - before & after가 모두 가능하면: 1x2 subplot 그림
        - before가 없으면: after만 단독 그림
    - regressed UMAP은 adata.obsm['X_umap_cell_cycle_regressed'] 로 원본에도 저장
    """
    if "S_score" not in adata.obs or "G2M_score" not in adata.obs:
        print("S_score / G2M_score not found. Skip UMAP comparison.")
        return adata

    os.makedirs(save_path, exist_ok=True)

    has_umap_before = "X_umap" in adata.obsm

    # cell cycle score를 회귀로 제거한 복사본 생성
    adata_cc = adata.copy()
    # S_score, G2M_score를 공변량으로 회귀 제거
    sc.pp.regress_out(adata_cc, ["S_score", "G2M_score"])
    # 회귀 후 PCA → neighbors → UMAP
    sc.tl.pca(adata_cc, svd_solver="arpack")
    sc.pp.neighbors(adata_cc, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata_cc)

    # regressed UMAP을 원본 AnnData에도 저장
    adata.obsm["X_umap_cell_cycle_regressed"] = adata_cc.obsm["X_umap"].copy()

    # 플롯 생성
    if has_umap_before:
        # before / after 한 그림에 비교
        fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))

        sc.pl.umap(
            adata,
            color=color,
            ax=axes[0],
            show=False,
            title="Before regressing cell cycle",
        )
        sc.pl.umap(
            adata_cc,
            color=color,
            ax=axes[1],
            show=False,
            title="After regressing cell cycle",
        )

        fig.tight_layout()
        fig.savefig(
            os.path.join(save_path, "umap_cell_cycle_before_after.png"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig)
    else:
        # before UMAP이 없으면 after만 단독 플롯
        sc.pl.umap(
            adata_cc,
            color=color,
            show=False,
            title="UMAP after regressing cell cycle",
        )
        plt.tight_layout()
        plt.savefig(
            os.path.join(save_path, "umap_after_cell_cycle_regress.png"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()

    return adata


def score_cell_cycle(adata, save_path, species="human"):
    os.makedirs(save_path, exist_ok=True)

    # 1) cell cycle gene set 로딩
    s_genes, g2m_genes = _load_cycle_genes(species)
    s_genes = [g for g in s_genes if g in adata.var_names]
    g2m_genes = [g for g in g2m_genes if g in adata.var_names]

    print(f"S genes used: {len(s_genes)}")
    print(f"G2M genes used: {len(g2m_genes)}")

    # 2) cell cycle score 계산
    sc.pp.scale(adata, max_value=10)
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    # 3) 기본 violin / scatter 플롯
    sc.pl.violin(
        adata,
        ["S_score", "G2M_score"],
        jitter=0.4,
        groupby="sample",
        rotation=60,
        show=False,
    )
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "violin_cell_cycle.png"), dpi=300)
    plt.close()

    sc.pl.scatter(adata, x="S_score", y="G2M_score", color="phase", show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "scatter_cell_cycle.png"), dpi=300)
    plt.close()

    # 4) 샘플 / 클러스터 기준 cell cycle phase 구성비 요약 + 플롯
    _summarize_cell_cycle_composition(adata, save_path)

    # 5) cell cycle 효과 제거 전/후 UMAP 비교
    _compare_umap_cell_cycle_effect(adata, save_path, color="phase")

    return adata