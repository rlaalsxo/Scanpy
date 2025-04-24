import os
import argparse
import scanpy as sc
from create_adata import CreateAdata
from batchCorrection import BatchCorrection
from cellCycle import score_cell_cycle
from DEG import deg_analysis_with_sex_gene_filtering

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True, help="입력 디렉토리")
    parser.add_argument("--output_dir", required=True, help="결과 저장 디렉토리")
    parser.add_argument("--output_filename", default="adata_final.h5ad", help="최종 저장 파일 이름")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    save_path = args.output_dir
    print("📌 Step 1: Create Adata")
    adata = CreateAdata(
        basic_save_path=save_path,
        parent_dir=args.input_dir,
        output_filename=args.output_filename,
    )

    print("📌 Step 2: Batch Correction")
    adata = BatchCorrection(adata, save_path=os.path.join(save_path, "BatchCorrection"))

    print("📌 Step 3: Cell Cycle Scoring")
    score_cell_cycle(adata, save_path=os.path.join(save_path, "CellCycle"))

    print("📌 Step 4: DEG + Cell Type Prediction")
    adata = deg_analysis_with_sex_gene_filtering(
        adata,
        save_path=os.path.join(save_path, "DEG")
    )

    print("✅ All steps completed.")
    adata.write_h5ad(os.path.join(save_path, "DEG",args.output_filename), compression="gzip")

if __name__ == "__main__":
    main()
