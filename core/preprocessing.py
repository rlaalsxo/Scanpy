"""
전처리 관련 함수

QC, 필터링, 정규화, HVG 선택
"""
import numpy as np
import scanpy as sc

from config.defaults import HVG
from config.species import get_species_config


def validate_species(adata, species: str):
    config = get_species_config(species)
    mt_count = int(adata.var_names.str.startswith(config["mt_prefix"]).sum())
    ribo_count = int(adata.var_names.str.startswith(config["ribo_prefix"]).sum())

    if mt_count == 0 and ribo_count == 0:
        other = "human" if species == "mouse" else "mouse"
        other_config = get_species_config(other)
        other_mt = int(adata.var_names.str.startswith(other_config["mt_prefix"]).sum())
        other_ribo = int(adata.var_names.str.startswith(other_config["ribo_prefix"]).sum())

        if other_mt > 0 or other_ribo > 0:
            raise ValueError(
                f"Species is set to '{species}', but gene name patterns "
                f"match '{other}'.\n"
                f"  Current({species}): {mt_count} MT genes, {ribo_count} Ribo genes\n"
                f"  Detected({other}): {other_mt} MT genes, {other_ribo} Ribo genes\n"
                f"Please change species to '{other}'."
            )


def calculate_qc_metrics(adata, species: str):
    """
    QC metrics 계산 (mt, ribo, hb)

    Parameters
    ----------
    adata : AnnData
    species : str
    """
    config = get_species_config(species)

    adata.var["mt"] = adata.var_names.str.startswith(config["mt_prefix"])
    adata.var["ribo"] = adata.var_names.str.startswith(config["ribo_prefix"])
    adata.var["hb"] = adata.var_names.str.contains(config["hb_pattern"])

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        log1p=True,
        percent_top=None,
        inplace=True,
    )


def filter_cells_qc(
    adata,
    species: str,
    min_genes: int = None,
    max_pct_mt: float = None,
    min_pct_ribo: float = None,
):
    """
    QC 기반 세포 필터링

    Parameters
    ----------
    adata : AnnData
    species : str
    min_genes : int, optional
    max_pct_mt : float, optional
    min_pct_ribo : float, optional

    Returns
    -------
    AnnData (필터링된)
    """
    config = get_species_config(species)
    qc = config["qc"]

    if min_genes is None:
        min_genes = qc["min_genes"]
    if max_pct_mt is None:
        max_pct_mt = qc["max_pct_mt"]
    if min_pct_ribo is None:
        min_pct_ribo = qc["min_pct_ribo"]

    print(f"[QC] min_genes={min_genes}, max_pct_mt={max_pct_mt}, min_pct_ribo={min_pct_ribo}")

    n_total = adata.n_obs

    sc.pp.filter_cells(adata, min_genes=min_genes)
    n_after_genes = adata.n_obs

    sc.pp.filter_genes(adata, min_cells=3)

    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt, :]
    n_after_mt = adata.n_obs

    adata = adata[adata.obs["pct_counts_ribo"] > min_pct_ribo, :]
    n_after_ribo = adata.n_obs

    if n_after_ribo == 0:
        raise ValueError(
            f"No cells remaining after QC filtering.\n"
            f"  Total: {n_total}\n"
            f"  After min_genes={min_genes}: {n_after_genes}\n"
            f"  After max_pct_mt={max_pct_mt}: {n_after_mt}\n"
            f"  After min_pct_ribo={min_pct_ribo}: {n_after_ribo}\n"
            f"Please check if species='{species}' is correct."
        )

    return adata


def remove_genes(adata, species: str):
    """
    MALAT1, mt, hb 유전자 제거

    Parameters
    ----------
    adata : AnnData
    species : str

    Returns
    -------
    AnnData
    """
    config = get_species_config(species)

    malat = adata.var_names.str.startswith(config["malat"])
    mito = adata.var_names.str.startswith(config["mt_prefix"])
    hb = adata.var_names.str.contains(config["hb_pattern"])

    remove_mask = np.add(np.add(malat, mito), hb)
    adata = adata[:, ~remove_mask]

    sc.pp.filter_genes(adata, min_cells=3)

    return adata


def normalize_and_hvg(adata, n_top_genes: int = None, batch_key: str = None):
    """
    정규화 및 HVG 선택

    Parameters
    ----------
    adata : AnnData
    n_top_genes : int, optional
    batch_key : str, optional
    """
    if n_top_genes is None:
        n_top_genes = HVG["n_top_genes"]

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    kwargs = {"n_top_genes": n_top_genes}
    if batch_key:
        kwargs["batch_key"] = batch_key

    sc.pp.highly_variable_genes(adata, **kwargs)
