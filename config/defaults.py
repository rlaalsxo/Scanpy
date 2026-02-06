"""
기본 파라미터 설정

모든 하드코딩 값을 한 곳에서 관리
"""

# PCA 설정
PCA = {
    "n_comps": 50,
    "svd_solver": "arpack",
}

# Neighbors 설정
NEIGHBORS = {
    "n_neighbors": 15,
    "n_pcs": 40,
}

# Highly Variable Genes 설정
HVG = {
    "n_top_genes": 2000,
}

# BBKNN 설정 (세포 수 기반 자동 조정)
BBKNN_NEIGHBORS_WITHIN_BATCH = {
    100_000: 20,  # >100k cells
    50_000: 15,   # >50k cells
    10_000: 10,   # >10k cells
    0: 3,         # default
}

# Trajectory 분석 설정
TRAJECTORY = {
    "n_dcs": 15,
    "paga_threshold": 0.05,
}

# DEG 분석 설정
DEG = {
    "top_n": 20,
    "padj_threshold": 0.05,
    "method": "wilcoxon",
}

# Cell Communication 설정
CELL_COMMUNICATION = {
    "n_perms": 1000,
    "pvalue_threshold": 0.05,
    "top_n_interactions": 50,
}

# 시각화 설정
PLOTTING = {
    "dpi": 300,
    "genes_per_fig": 40,
}
