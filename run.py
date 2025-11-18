import os
import argparse
import scanpy as sc
from create_adata import CreateAdata
from batchCorrection import BatchCorrection
from cellCycle import score_cell_cycle
from DEG import deg_analysis_with_sex_gene_filtering

def ensure_umap(adata):
    """UMAP이 없으면 자동으로 PCA/Neighbors/UMAP 생성"""
    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata)
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata)
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)
    return adata

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True, help="입력 디렉토리")
    parser.add_argument("--output_dir", required=True, help="결과 저장 디렉토리")
    parser.add_argument("--output_filename", default="adata_final.h5ad", help="최종 저장 파일 이름")
    # 종 / 조직 정보
    parser.add_argument("--species", required=True, choices=["human", "mouse"],
                        help="종 선택 (human 또는 mouse)")
    # 선택 실행 옵션
    parser.add_argument("--run_3", action="store_true", help="Cell Cycle 단계 실행")
    parser.add_argument("--run_4", action="store_true", help="DEG 분석 실행")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    save_path = args.output_dir
    # Step 1: Create Adata (전처리 - 항상 실행)
    print("📌 Step 1: Create Adata")
    adata = CreateAdata( basic_save_path=save_path, parent_dir=args.input_dir, species=args.species )
    # Step 2: Batch Correction (항상 실행)
    print("📌 Step 2: Batch Correction")
    adata = BatchCorrection( adata, save_path=os.path.join(save_path, "BatchCorrection"))
    # Step 3: Cell Cycle (선택 실행)
    if args.run_3:
        print("📌 Step 3: Cell Cycle Scoring")
        adata = score_cell_cycle(adata, save_path=os.path.join(save_path, "CellCycle"), species=args.species)
    else:
        print("⏭ Step 3 skipped")
    # Step 4: DEG + Cell Type Prediction (선택 실행)
    if args.run_4:
        print("📌 Step 4: DEG + Cell Type Prediction")
        # UMAP 존재 여부 확인 → 없으면 자동 생성
        adata = ensure_umap(adata)
        adata = deg_analysis_with_sex_gene_filtering( adata, save_path=os.path.join(save_path, "DEG"), species=args.species )
    else:
        print("⏭ Step 4 skipped")

    print("✅ All steps completed.")
    # 저장 경로 설정 (DEG 실행 여부에 따라)
    final_dir = os.path.join(save_path, "DEG") if args.run_4 else save_path
    os.makedirs(save_path, exist_ok=True)

    adata.write_h5ad( os.path.join(save_path, args.output_filename), compression="gzip" )

if __name__ == "__main__":
    main()
