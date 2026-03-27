"""
파일 I/O 유틸리티

10x / BD Rhapsody 데이터 로딩, 파일명 표준화, h5ad 저장
"""
import os
import re
import uuid
import glob
import gzip
import shutil
import tempfile
import tarfile

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix


def _gzip_file(src: str, dst: str):
    with open(src, 'rb') as f_in:
        with gzip.open(dst, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def standardize_filenames(data_dir: str):
    """
    10x 파일명 표준화

    다양한 파일명 형식을 scanpy가 인식하는 형식으로 변환
    비압축 파일(.tsv, .mtx)은 .gz로 압축
    """
    files = os.listdir(data_dir)

    # barcodes.tsv.gz
    bc = os.path.join(data_dir, "barcodes.tsv.gz")
    if not os.path.exists(bc):
        cands = [f for f in files if "barcodes.tsv.gz" in f]
        if cands:
            shutil.copy(os.path.join(data_dir, cands[0]), bc)
        else:
            cands = [f for f in files if "barcodes.tsv" in f and not f.endswith(".gz")]
            if cands:
                _gzip_file(os.path.join(data_dir, cands[0]), bc)

    # features.tsv.gz / genes.tsv.gz
    feats = os.path.join(data_dir, "features.tsv.gz")
    genes = os.path.join(data_dir, "genes.tsv.gz")

    if not os.path.exists(genes):
        cands = [f for f in files if "genes.tsv.gz" in f]
        if cands:
            shutil.copy(os.path.join(data_dir, cands[0]), genes)
        else:
            cands = [f for f in files if "genes.tsv" in f and not f.endswith(".gz")]
            if cands:
                _gzip_file(os.path.join(data_dir, cands[0]), genes)

    if not os.path.exists(feats):
        cands = [f for f in files if "features.tsv.gz" in f]
        if cands:
            shutil.copy(os.path.join(data_dir, cands[0]), feats)
        else:
            cands = [f for f in files if "features.tsv" in f and not f.endswith(".gz")]
            if cands:
                _gzip_file(os.path.join(data_dir, cands[0]), feats)
            elif os.path.exists(genes):
                shutil.copy(genes, feats)

    if not os.path.exists(feats) and os.path.exists(genes):
        shutil.copy(genes, feats)

    # features.tsv.gz 컬럼 수정 (2열 -> 3열)
    _fix_features_file(data_dir)

    # matrix.mtx.gz
    mtx = os.path.join(data_dir, "matrix.mtx.gz")
    if not os.path.exists(mtx):
        cands = [f for f in files if "matrix.mtx" in f]
        if cands:
            src = os.path.join(data_dir, cands[0])
            if cands[0].endswith(".gz"):
                shutil.copy(src, mtx)
            else:
                _gzip_file(src, mtx)


def _fix_features_file(data_dir: str):
    """features.tsv.gz가 2열이면 3열로 수정"""
    path = os.path.join(data_dir, "features.tsv.gz")
    if not os.path.exists(path):
        return

    with gzip.open(path, 'rt') as f:
        df = pd.read_csv(f, sep='\t', header=None)

    if df.shape[1] == 2:
        df[2] = "Gene Expression"
        tmp = os.path.join(data_dir, "features_fixed.tsv.gz")
        df.to_csv(tmp, sep='\t', header=False, index=False, compression="gzip")
        shutil.move(path, path + ".bak")
        shutil.move(tmp, path)


def _extract_tar(path: str) -> str:
    """tar/tar.gz 파일 추출 (중첩 tar 포함). 추출된 디렉토리 반환"""
    tmp = tempfile.mkdtemp()
    with tarfile.open(path, "r:*") as tar:
        tar.extractall(path=tmp)
    for f in os.listdir(tmp):
        p = os.path.join(tmp, f)
        if os.path.isfile(p) and tarfile.is_tarfile(p):
            out = os.path.join(tmp, "nested_" + os.path.splitext(f)[0])
            os.makedirs(out, exist_ok=True)
            with tarfile.open(p, "r:*") as t:
                t.extractall(path=out)
    return tmp


def _get_obs_name(style: str, sample_id: str, barcode: str) -> str:
    """obs_name 생성"""
    if style == "folder_barcode":
        return f"{sample_id}_{barcode}"
    if style == "uuid":
        return f"{barcode}_{uuid.uuid4().hex[:6]}"
    if style == "barcode_only":
        return barcode
    raise ValueError(f"Invalid obs_name_style: {style}")


def load_10x_data(
    parent_dir: str,
    sample_names: list = None,
    obs_name_style: str = "folder_barcode",
) -> sc.AnnData:
    """
    10x Genomics 데이터 로딩

    Parameters
    ----------
    parent_dir : str
        샘플 폴더들이 있는 상위 디렉토리 또는 tar 파일 경로
    sample_names : list, optional
        샘플 이름 리스트 (없으면 인덱스 사용)
    obs_name_style : str
        obs_names 생성 방식 ("folder_barcode", "uuid", "barcode_only")

    Returns
    -------
    AnnData
        합쳐진 AnnData 객체
    """
    # tar 파일 처리
    if parent_dir.endswith(".tar") or parent_dir.endswith(".tar.gz"):
        parent_dir = _extract_tar(parent_dir)

    entries = os.listdir(parent_dir)
    full = [os.path.join(parent_dir, x) for x in entries]
    is_flat = all(os.path.isfile(p) for p in full)

    data_list = []

    if is_flat:
        # 플랫 구조: 파일명에서 샘플 ID 추출
        samp = {}
        for f in entries:
            m = re.match(r"(.+?)[._-](barcodes|features|genes|matrix).*", f)
            if m:
                samp.setdefault(m.group(1), []).append(f)

        if samp:
            tmp = tempfile.mkdtemp()
            for idx, (sid, files) in enumerate(samp.items()):
                path = os.path.join(tmp, sid)
                os.makedirs(path, exist_ok=True)
                for f in files:
                    shutil.copy(os.path.join(parent_dir, f), os.path.join(path, f))

                standardize_filenames(path)
                ad = sc.read_10x_mtx(path, var_names="gene_symbols")
                ad.obs_names = [_get_obs_name(obs_name_style, sid, bc) for bc in ad.obs_names]
                ad.obs["sample"] = sample_names[idx] if sample_names else str(idx)
                data_list.append(ad)
        else:
            standardize_filenames(parent_dir)
            ad = sc.read_10x_mtx(parent_dir, var_names="gene_symbols")
            sid = os.path.basename(parent_dir.rstrip("/"))
            ad.obs_names = [_get_obs_name(obs_name_style, sid, bc) for bc in ad.obs_names]
            ad.obs["sample"] = sample_names[0] if sample_names else "0"
            data_list.append(ad)
    else:
        # 폴더 구조
        folders = [os.path.join(parent_dir, d) for d in entries
                   if os.path.isdir(os.path.join(parent_dir, d))]

        for idx, d in enumerate(folders):
            sid = os.path.basename(d)
            standardize_filenames(d)
            ad = sc.read_10x_mtx(d, var_names="gene_symbols")
            ad.obs_names = [_get_obs_name(obs_name_style, sid, bc) for bc in ad.obs_names]
            ad.obs["sample"] = sample_names[idx] if sample_names else str(idx)
            data_list.append(ad)

    # 합치기
    if len(data_list) > 1:
        adata = data_list[0].concatenate(*data_list[1:], batch_key="batch")
    else:
        adata = data_list[0]

    return adata


def _find_bd_csv(parent_dir: str) -> list:
    """BD Rhapsody 발현 매트릭스 CSV 파일 탐색 (재귀, 압축/비압축 모두)"""
    for ext in ["*.csv", "*.csv.gz"]:
        pat = f"*_MolsPerCell{ext}"
        found = glob.glob(os.path.join(parent_dir, "**", pat), recursive=True)
        if not found:
            found = glob.glob(os.path.join(parent_dir, pat))
        if found:
            return sorted(found)

    raise FileNotFoundError(
        f"BD Rhapsody CSV 파일을 찾을 수 없습니다: {parent_dir}\n"
        "Expected: *_MolsPerCell*.csv 또는 *_MolsPerCell*.csv.gz"
    )


def _is_abseq_column(col: str) -> bool:
    """AbSeq(단백질) 컬럼 여부 판별"""
    return bool(re.search(r"[:|]", col))


def load_bd_data(
    parent_dir: str,
    sample_names: list = None,
) -> sc.AnnData:
    """
    BD Rhapsody 데이터 로딩

    tar/tar.gz, 중첩 폴더 구조 모두 지원 (10x 로더와 동일)

    Parameters
    ----------
    parent_dir : str
        BD Rhapsody CSV 파일이 있는 디렉토리 또는 tar 파일 경로
    sample_names : list, optional
        샘플 이름 리스트

    Returns
    -------
    AnnData
    """
    # tar 파일 처리
    if parent_dir.endswith(".tar") or parent_dir.endswith(".tar.gz"):
        parent_dir = _extract_tar(parent_dir)

    csv_files = _find_bd_csv(parent_dir)
    print(f"[IO] Found {len(csv_files)} BD CSV file(s)")

    data_list = []
    for idx, csv_path in enumerate(csv_files):
        print(f"[IO] Loading: {os.path.basename(csv_path)}")
        df = pd.read_csv(csv_path, index_col=0, compression="infer")

        # AbSeq(단백질) 컬럼 분리
        abseq_cols = [c for c in df.columns if _is_abseq_column(c)]
        mrna_cols = [c for c in df.columns if not _is_abseq_column(c)]

        if abseq_cols:
            print(f"[IO] AbSeq columns detected: {len(abseq_cols)} (separated to obsm['protein'])")

        mrna_df = df[mrna_cols]
        cell_indices = mrna_df.index.astype(str)

        ad = sc.AnnData(
            X=csr_matrix(mrna_df.values.astype(np.float32)),
            obs=pd.DataFrame(index=[f"BD_{ci}" for ci in cell_indices]),
            var=pd.DataFrame(index=mrna_df.columns),
        )
        ad.var_names_make_unique()

        if abseq_cols:
            ad.obsm["protein"] = df[abseq_cols].values.astype(np.float32)

        ad.obs["sample"] = sample_names[idx] if sample_names else str(idx)
        data_list.append(ad)

    # Sample Tag 처리 (재귀 탐색, 압축/비압축 모두)
    tag_files = []
    for ext in ["*.csv", "*.csv.gz"]:
        tag_files = glob.glob(os.path.join(parent_dir, "**", f"*Sample_Tag_Calls{ext}"), recursive=True)
        if not tag_files:
            tag_files = glob.glob(os.path.join(parent_dir, f"*Sample_Tag_Calls{ext}"))
        if tag_files:
            break
    if tag_files and len(data_list) == 1:
        ad = data_list[0]
        tags = pd.read_csv(tag_files[0], index_col=0, compression="infer")

        tag_col = None
        for col in ["Sample_Tag", "Sample_Name"]:
            if col in tags.columns:
                tag_col = col
                break

        if tag_col:
            tags = tags[~tags[tag_col].isin(["Multiplet", "Undetermined"])]
            common_idx = ad.obs.index.intersection(
                pd.Index([f"BD_{i}" for i in tags.index.astype(str)])
            )
            if len(common_idx) > 0:
                ad = ad[common_idx].copy()
                tag_map = {f"BD_{str(i)}": t for i, t in zip(tags.index, tags[tag_col])}
                ad.obs["sample"] = [tag_map[ci] for ci in ad.obs_names]
                n_removed = data_list[0].n_obs - ad.n_obs
                if n_removed > 0:
                    print(f"[IO] Removed {n_removed} Multiplet/Undetermined cells")

                # Sample Tag 기반 obs_names 재생성
                ad.obs_names = [
                    f"{sample}_{ci.split('_', 1)[1]}"
                    for ci, sample in zip(ad.obs_names, ad.obs["sample"])
                ]
                data_list = [ad]

    if len(data_list) > 1:
        adata = data_list[0].concatenate(*data_list[1:], batch_key="batch")
    else:
        adata = data_list[0]

    print(f"[IO] BD data loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def save_adata(adata, save_path: str, filename: str = "adata.h5ad"):
    """
    AnnData 저장

    Parameters
    ----------
    adata : AnnData
    save_path : str
    filename : str
    """
    os.makedirs(save_path, exist_ok=True)
    adata.write_h5ad(os.path.join(save_path, filename), compression="gzip")
    print(f"[IO] Saved: {os.path.join(save_path, filename)}")
