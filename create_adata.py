import matplotlib
matplotlib.use("Agg")
import warnings, os, shutil, gzip, pandas as pd, scanpy as sc, numpy as np
import matplotlib.pyplot as plt, tempfile, re, uuid, tarfile
import scrublet as scr

warnings.filterwarnings("ignore")

def _gene_prefix(species):
    if species == "human":
        return {"mt": "MT-", "ribo": ("RPS", "RPL"), "malat": "MALAT1", "hb": "^HB[^(P|E|S)]"}
    if species == "mouse":
        return {"mt": "mt-", "ribo": ("Rps", "Rpl"), "malat": "Malat1", "hb": "^Hb[^(p|e|s)]"}
    raise ValueError(f"Invalid species: {species}")

def fix_features_file(data_dir):
    path = os.path.join(data_dir, "features.tsv.gz")
    if not os.path.exists(path): return
    with gzip.open(path, 'rt') as f:
        df = pd.read_csv(f, sep='\t', header=None)
    if df.shape[1] == 2:
        df[2] = "Gene Expression"
        tmp = os.path.join(data_dir, "features_fixed.tsv.gz")
        df.to_csv(tmp, sep='\t', header=False, index=False, compression="gzip")
        shutil.move(path, path + ".bak"); shutil.move(tmp, path)

def standardize_filenames(data_dir):
    files = os.listdir(data_dir)
    bc = os.path.join(data_dir, "barcodes.tsv.gz")
    if not os.path.exists(bc):
        cands = [f for f in files if "barcodes.tsv.gz" in f]
        if cands: shutil.copy(os.path.join(data_dir, cands[0]), bc)

    feats = os.path.join(data_dir, "features.tsv.gz")
    genes = os.path.join(data_dir, "genes.tsv.gz")
    if not os.path.exists(genes):
        cands = [f for f in files if "genes.tsv.gz" in f]
        if cands: shutil.copy(os.path.join(data_dir, cands[0]), genes)
    if not os.path.exists(feats):
        cands = [f for f in files if "features.tsv.gz" in f]
        if cands: shutil.copy(os.path.join(data_dir, cands[0]), feats)
        elif os.path.exists(genes): shutil.copy(genes, feats)

    fix_features_file(data_dir)

    mtx = os.path.join(data_dir, "matrix.mtx.gz")
    if not os.path.exists(mtx):
        cands = [f for f in files if "matrix.mtx" in f]
        if cands:
            src = os.path.join(data_dir, cands[0])
            shutil.copy(src, mtx if cands[0].endswith(".gz") else os.path.join(data_dir, "matrix.mtx"))

def CreateAdata(basic_save_path, parent_dir, sample_names=None,
                obs_name_style="folder_barcode", min_genes=200, max_pct_mt=20,
                min_pct_ribo=5, perform_analysis=True, species="human"):

    if parent_dir.endswith(".tar") or parent_dir.endswith(".tar.gz"):
        tmp = tempfile.mkdtemp()
        with tarfile.open(parent_dir, "r:*") as tar: tar.extractall(path=tmp)
        for f in os.listdir(tmp):
            p = os.path.join(tmp, f)
            if tarfile.is_tarfile(p):
                out = os.path.join(tmp, "nested_" + os.path.splitext(f)[0])
                os.makedirs(out, exist_ok=True)
                with tarfile.open(p, "r:*") as t: t.extractall(path=out)
        parent_dir = tmp

    save_path = os.path.join(basic_save_path, "CreateAdata")
    os.makedirs(save_path, exist_ok=True)

    entries = os.listdir(parent_dir)
    full = [os.path.join(parent_dir, x) for x in entries]
    is_flat = all(os.path.isfile(p) for p in full)
    data_list = []

    def obs_name(style, sid, bc):
        if style == "folder_barcode": return f"{sid}_{bc}"
        if style == "uuid": return f"{bc}_{uuid.uuid4().hex[:6]}"
        if style == "barcode_only": return bc
        raise ValueError(f"Invalid obs_name_style: {style}")

    if is_flat:
        samp = {}
        for f in entries:
            m = re.match(r"(.+?)[._-](barcodes|features|genes|matrix).*", f)
            if m: samp.setdefault(m.group(1), []).append(f)
        tmp = tempfile.mkdtemp()
        for idx, (sid, files) in enumerate(samp.items()):
            path = os.path.join(tmp, sid); os.makedirs(path, exist_ok=True)
            for f in files: shutil.copy(os.path.join(parent_dir, f), os.path.join(path, f))
            standardize_filenames(path)
            ad = sc.read_10x_mtx(path, var_names="gene_symbols")
            ad.obs_names = [obs_name(obs_name_style, sid, bc) for bc in ad.obs_names]
            ad.obs["sample"] = sample_names[idx] if sample_names else str(idx)
            data_list.append(ad)
    else:
        folders = [os.path.join(parent_dir, d) for d in entries if os.path.isdir(os.path.join(parent_dir, d))]
        for idx, d in enumerate(folders):
            sid = os.path.basename(d)
            standardize_filenames(d)
            ad = sc.read_10x_mtx(d, var_names="gene_symbols")
            ad.obs_names = [obs_name(obs_name_style, sid, bc) for bc in ad.obs_names]
            ad.obs["sample"] = sample_names[idx] if sample_names else str(idx)
            data_list.append(ad)

    adata = data_list[0].concatenate(*data_list[1:], batch_key="batch") if len(data_list) > 1 else data_list[0]

    pf = _gene_prefix(species)
    adata.var["mt"] = adata.var_names.str.startswith(pf["mt"])
    adata.var["ribo"] = adata.var_names.str.startswith(pf["ribo"])
    adata.var["hb"] = adata.var_names.str.contains(pf["hb"])
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], log1p=True, percent_top=None, inplace=True)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt, :]
    adata = adata[adata.obs["pct_counts_ribo"] > min_pct_ribo, :]

    malat = adata.var_names.str.startswith(pf["malat"])
    mito = adata.var_names.str.startswith(pf["mt"])
    hb = adata.var_names.str.contains(pf["hb"])
    rm = np.add(np.add(malat, mito), hb)
    adata = adata[:, ~rm]

    sc.pp.filter_genes(adata, min_cells=3)

    dbl_score = pd.Series(index=adata.obs_names, dtype=float)
    dbl_pred = pd.Series(index=adata.obs_names, dtype=bool)
    for s in adata.obs["sample"].unique():
        mask = adata.obs["sample"] == s
        sub = adata[mask]
        scrub = scr.Scrublet(sub.X)
        score, pred = scrub.scrub_doublets()
        dbl_score.loc[sub.obs_names] = score
        dbl_pred.loc[sub.obs_names] = pred

    adata.obs["doublet_score"] = dbl_score
    adata.obs["predicted_doublet"] = dbl_pred
    adata.obs["doublet_info"] = dbl_pred.astype(str)

    sc.pl.violin(adata, "n_genes_by_counts", groupby="doublet_info", size=0, rotation=45, show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "violin_doublet.png"), dpi=300); plt.close()

    adata = adata[adata.obs["doublet_info"] == "False", :]

    if perform_analysis:
        sc.pp.normalize_total(adata); sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
        sc.tl.pca(adata); sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata)
        sc.pl.umap(adata, color="sample", show=False)
        plt.tight_layout()
        plt.savefig(os.path.join(save_path, "umap_after_doublet.png"), dpi=300)
        plt.close()

    keep_obs, keep_var = ["sample"], ["gene_ids", "feature_types"]
    adata.obs = adata.obs[[c for c in adata.obs.columns if c in keep_obs]]
    adata.var = adata.var[[c for c in adata.var.columns if c in keep_var]]
    adata.uns = {}; adata.obsm = {}; adata.varm = {}; adata.obsp = {}

    adata.write_h5ad(os.path.join(basic_save_path, "scRNA_preprocessing.h5ad"), compression="gzip")
    return adata