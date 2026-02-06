import os
import argparse
from create_adata import CreateAdata
from batchCorrection import BatchCorrection
from cellCycle import score_cell_cycle
from DEG import deg_analysis_with_sex_gene_filtering
from trajectory import trajectory_analysis
from cellCommunication import cellcell_communication
from common import ensure_neighbors, ensure_umap, detect_cluster_key


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
    parser.add_argument("--run_5", action="store_true", help="Trajectory 분석 실행")
    parser.add_argument("--run_6", action="store_true", help="Cell-Cell Communication 분석 실행")
    # Trajectory 옵션
    parser.add_argument("--root_cluster", type=str, default=None,
                        help="Trajectory 시작점 클러스터 ID")
    parser.add_argument("--root_gene", type=str, default=None,
                        help="Trajectory 시작점 마커 유전자")
    # Cell Communication 옵션
    parser.add_argument("--n_perms", type=int, default=1000,
                        help="Cell communication permutation 수")
    parser.add_argument("--pvalue_threshold", type=float, default=0.05,
                        help="유의성 임계값")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    save_path = args.output_dir

    # =========================================================================
    # Step 1: Create Adata (전처리 - 항상 실행)
    # =========================================================================
    print("Step 1: Create Adata")
    adata = CreateAdata(basic_save_path=save_path, parent_dir=args.input_dir, species=args.species)

    # =========================================================================
    # Step 2: Batch Correction (항상 실행)
    # - BBKNN neighbors + leiden + UMAP 생성
    # =========================================================================
    print("Step 2: Batch Correction")
    adata = BatchCorrection(adata, save_path=os.path.join(save_path, "BatchCorrection"))

    # =========================================================================
    # 전처리 준비 (필요한 경우에만)
    # =========================================================================
    # Step 5용 PAGA neighbors 준비 (별도 키로 저장, BBKNN neighbors 유지)
    if args.run_5:
        print("Preparing neighbors for Trajectory analysis...")
        ensure_neighbors(adata, key_added="neighbors_paga")

    # =========================================================================
    # Step 3: Cell Cycle (선택 실행)
    # =========================================================================
    if args.run_3:
        print("Step 3: Cell Cycle Scoring")
        adata = score_cell_cycle(adata, save_path=os.path.join(save_path, "CellCycle"), species=args.species)
    else:
        print("Step 3 skipped")

    # =========================================================================
    # Step 4: DEG + Cell Type Prediction (선택 실행)
    # =========================================================================
    if args.run_4:
        print("Step 4: DEG + Cell Type Prediction")
        ensure_umap(adata)  # 이미 있으면 스킵
        adata = deg_analysis_with_sex_gene_filtering(adata, save_path=os.path.join(save_path, "DEG"), species=args.species)
    else:
        print("Step 4 skipped")

    # =========================================================================
    # Step 5: Trajectory Analysis (선택 실행)
    # =========================================================================
    if args.run_5:
        print("Step 5: Trajectory Analysis")
        adata = trajectory_analysis(
            adata,
            save_path=os.path.join(save_path, "Trajectory"),
            species=args.species,
            root_cluster=args.root_cluster,
            root_gene=args.root_gene,
            neighbors_key="neighbors_paga",  # PAGA용 neighbors 사용
        )
    else:
        print("Step 5 skipped")

    # =========================================================================
    # Step 6: Cell-Cell Communication (선택 실행)
    # =========================================================================
    if args.run_6:
        print("Step 6: Cell-Cell Communication")
        # groupby 결정: cell_type (Step 4 실행 시) > leiden
        groupby = "cell_type" if "cell_type" in adata.obs else detect_cluster_key(adata)
        adata = cellcell_communication(
            adata,
            save_path=os.path.join(save_path, "CellCommunication"),
            species=args.species,
            groupby=groupby,
            n_perms=args.n_perms,
            pvalue_threshold=args.pvalue_threshold,
        )
    else:
        print("Step 6 skipped")

    # =========================================================================
    # 최종 저장
    # =========================================================================
    print("All steps completed.")
    os.makedirs(save_path, exist_ok=True)
    adata.write_h5ad(os.path.join(save_path, args.output_filename), compression="gzip")
    print(f"Saved: {os.path.join(save_path, args.output_filename)}")


if __name__ == "__main__":
    main()
