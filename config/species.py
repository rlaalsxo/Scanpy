"""
종(species)별 설정

human/mouse별로 다른 설정값들을 통합 관리
"""
import os

SPECIES_CONFIG = {
    "human": {
        # Gene prefix/patterns
        "mt_prefix": "MT-",
        "ribo_prefix": ("RPS", "RPL"),
        "malat": "MALAT1",
        "hb_pattern": "^HB[^(P|E|S)]",

        # QC thresholds
        "qc": {
            "min_genes": 200,
            "max_pct_mt": 20,
            "min_pct_ribo": 5,
        },

        # BioMart
        "biomart_org": "hsapiens",

        # Enrichr libraries for cell type prediction
        "enrichr_libs": [
            "CellMarker_Augmented_2021",
            "PanglaoDB_Augmented_2021",
            "Descartes_Cell_Types_and_Tissue_2021",
            "Tabula_Sapiens",
            "Azimuth_Cell_Types_2021",
        ],

        # Cell cycle genes file
        "cell_cycle_file": "human_cell_cycle.txt",
    },

    "mouse": {
        # Gene prefix/patterns
        "mt_prefix": "mt-",
        "ribo_prefix": ("Rps", "Rpl"),
        "malat": "Malat1",
        "hb_pattern": "^Hb[^(p|e|s)]",

        # QC thresholds
        "qc": {
            "min_genes": 200,
            "max_pct_mt": 20,
            "min_pct_ribo": 5,
        },

        # BioMart
        "biomart_org": "mmusculus",

        # Enrichr libraries for cell type prediction
        "enrichr_libs": [
            "Tabula_Muris",
            "Mouse_Gene_Atlas",
            "PanglaoDB_Augmented_2021",
        ],

        # Cell cycle genes file
        "cell_cycle_file": "mouse_cell_cycle.txt",
    },
}


def get_species_config(species: str) -> dict:
    """종별 설정 반환"""
    if species not in SPECIES_CONFIG:
        raise ValueError(f"Invalid species: {species}. Choose from {list(SPECIES_CONFIG.keys())}")
    return SPECIES_CONFIG[species]


def get_cell_cycle_genes(species: str) -> tuple:
    """
    Cell cycle genes 로딩 (S phase, G2M phase)

    Returns
    -------
    tuple: (s_genes, g2m_genes)
    """
    config = get_species_config(species)
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(here, "input", config["cell_cycle_file"])

    if not os.path.exists(path):
        raise FileNotFoundError(f"Cell cycle file not found: {path}")

    with open(path, "r") as f:
        genes = [x.strip() for x in f.readlines() if x.strip()]

    # 앞 43개: S phase, 나머지: G2M phase
    return genes[:43], genes[43:]
