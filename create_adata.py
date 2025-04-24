import warnings
import os
import shutil
import gzip
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import tempfile
import re
import uuid

warnings.filterwarnings("ignore")

def fix_features_file(data_dir):
    features_path = os.path.join(data_dir, "features.tsv.gz")
    if not os.path.exists(features_path):
        return
    with gzip.open(features_path, 'rt') as f:
        df = pd.read_csv(f, sep='\t', header=None)
    if df.shape[1] == 2:
        df[2] = "Gene Expression"
        new_features_path = os.path.join(data_dir, "features_fixed.tsv.gz")
        df.to_csv(new_features_path, sep='\t', header=False, index=False, compression="gzip")
        shutil.move(features_path, features_path + ".bak")
        shutil.move(new_features_path, features_path)

def standardize_filenames(data_dir):
    files = os.listdir(data_dir)
    barcode_std = os.path.join(data_dir, "barcodes.tsv.gz")
    if not os.path.exists(barcode_std):
        barcode_files = [f for f in files if "barcodes.tsv.gz" in f]
        if barcode_files:
            shutil.copy(os.path.join(data_dir, barcode_files[0]), barcode_std)

    features_std = os.path.join(data_dir, "features.tsv.gz")
    genes_std = os.path.join(data_dir, "genes.tsv.gz")
    if not os.path.exists(genes_std):
        genes_files = [f for f in files if "genes.tsv.gz" in f]
        if genes_files:
            shutil.copy(os.path.join(data_dir, genes_files[0]), genes_std)
    if not os.path.exists(features_std):
        features_files = [f for f in files if "features.tsv.gz" in f]
        if features_files:
            shutil.copy(os.path.join(data_dir, features_files[0]), features_std)
        elif os.path.exists(genes_std):
            shutil.copy(genes_std, features_std)

    fix_features_file(data_dir)

    matrix_std = os.path.join(data_dir, "matrix.mtx.gz")
    if not os.path.exists(matrix_std):
        matrix_files = [f for f in files if "matrix.mtx" in f]
        if matrix_files:
            if matrix_files[0].endswith(".gz"):
                shutil.copy(os.path.join(data_dir, matrix_files[0]), matrix_std)
            else:
                shutil.copy(os.path.join(data_dir, matrix_files[0]), os.path.join(data_dir, "matrix.mtx"))

def CreateAdata(
    basic_save_path,
    parent_dir,
    output_filename,
    sample_names=None,
    obs_name_style="folder_barcode",
    min_genes=200,
    max_pct_mt=20,
    min_pct_ribo=5,
    perform_analysis=True
):
    save_path = os.path.join(basic_save_path, "CreateAdata")
    os.makedirs(save_path, exist_ok=True)
    data_list = []

    all_entries = os.listdir(parent_dir)
    full_paths = [os.path.join(parent_dir, f) for f in all_entries]
    is_flat = all(os.path.isfile(p) for p in full_paths)

    def generate_obs_name(style, sample_id, barcode):
        if style == "folder_barcode":
            return f"{sample_id}_{barcode}"
        elif style == "uuid":
            return f"{barcode}_{uuid.uuid4().hex[:6]}"
        elif style == "barcode_only":
            return barcode
        else:
            raise ValueError(f"Invalid obs_name_style: {style}")

    if is_flat:
        sample_dict = {}
        for f in all_entries:
            m = re.match(r"(.+?)_(barcodes|features|genes|matrix).*", f)
            if m:
                sample_name = m.group(1)
                sample_dict.setdefault(sample_name, []).append(f)

        temp_dir = tempfile.mkdtemp()
        for idx, (sample_id, files) in enumerate(sample_dict.items()):
            sample_path = os.path.join(temp_dir, sample_id)
            os.makedirs(sample_path, exist_ok=True)
            for file in files:
                shutil.copy(os.path.join(parent_dir, file), os.path.join(sample_path, file))

            standardize_filenames(sample_path)
            ad = sc.read_10x_mtx(sample_path, var_names='gene_symbols', cache=True)
            ad.obs_names = [generate_obs_name(obs_name_style, sample_id, bc) for bc in ad.obs_names]
            ad.obs["sample"] = sample_names[idx] if sample_names else str(idx)
            data_list.append(ad)
    else:
        folders = [os.path.join(parent_dir, d) for d in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, d))]
        for idx, data_dir in enumerate(folders):
            sample_id = os.path.basename(data_dir)
            standardize_filenames(data_dir)
            ad = sc.read_10x_mtx(data_dir, var_names='gene_symbols', cache=True)
            ad.obs_names = [generate_obs_name(obs_name_style, sample_id, bc) for bc in ad.obs_names]
            ad.obs["sample"] = sample_names[idx] if sample_names else str(idx)
            data_list.append(ad)

    adata = data_list[0].concatenate(*data_list[1:], batch_key="batch") if len(data_list) > 1 else data_list[0]

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P|E|S)]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=True, inplace=True)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=3)

    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt, :]
    adata = adata[adata.obs["pct_counts_ribo"] > min_pct_ribo, :]

    malat1 = adata.var_names.str.startswith("MALAT1")
    mito_genes = adata.var_names.str.startswith("MT-")
    hb_genes = adata.var_names.str.contains("^HB[^(P|E|S)]")
    remove = np.add(np.add(malat1, mito_genes), hb_genes)
    adata = adata[:, np.invert(remove)]

    sc.pp.filter_genes(adata, min_cells=3)

    sc.pp.scrublet(adata, batch_key="sample")
    adata.obs["doublet_info"] = adata.obs["predicted_doublet"].astype(str)

    sc.pl.violin(adata, "n_genes_by_counts", size=0, groupby="doublet_info", rotation=45, show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "violin_doublet.png"), dpi=300)
    plt.close()

    adata = adata[adata.obs["doublet_info"] == "False", :]

    if perform_analysis:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
        sc.tl.pca(adata)
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata)

        sc.pl.umap(adata, color=["sample"], show=False)
        plt.tight_layout()
        plt.savefig(os.path.join(save_path, "umap_after_doublet.png"), dpi=300)
        plt.close()

    keep_obs = ["sample"]
    keep_var = ["gene_ids", "feature_types"]
    adata.obs = adata.obs[[col for col in adata.obs.columns if col in keep_obs]]
    adata.var = adata.var[[col for col in adata.var.columns if col in keep_var]]
    adata.uns = {}; adata.obsm = {}; adata.varm = {}; adata.obsp = {}

    adata.write_h5ad(os.path.join(basic_save_path, output_filename), compression="gzip")
    return adata
