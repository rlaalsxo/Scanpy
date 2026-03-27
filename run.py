"""
scRNA-seq 분석 파이프라인

메인 실행 스크립트
"""
import os
import argparse

from analysis import (
    create_adata,
    batch_correction,
    score_cell_cycle,
    deg_analysis,
    trajectory_analysis,
    cellcell_communication,
    gsea_enrichment,
    pseudotime_expression,
    cell_proportion_analysis,
    marker_feature_plot,
)
from core.neighbors import ensure_umap, detect_cluster_key


def main():
    parser = argparse.ArgumentParser(description="scRNA-seq Analysis Pipeline")

    # 필수 인자
    parser.add_argument("--input_dir", required=True, help="입력 디렉토리")
    parser.add_argument("--output_dir", required=True, help="결과 저장 디렉토리")
    parser.add_argument("--species", required=True, choices=["human", "mouse"], help="종 선택")
    parser.add_argument("--platform", default="10x", choices=["10x", "bd"], help="시퀀싱 플랫폼 (10x 또는 bd)")

    # 선택 인자
    parser.add_argument("--output_filename", default="adata_final.h5ad", help="최종 저장 파일 이름")

    # 선택 실행 옵션
    parser.add_argument("--run_3", action="store_true", help="Cell Cycle 단계 실행")
    parser.add_argument("--run_4", action="store_true", help="DEG 분석 실행")
    parser.add_argument("--run_5", action="store_true", help="Trajectory 분석 실행")
    parser.add_argument("--run_6", action="store_true", help="Cell-Cell Communication 분석 실행")
    parser.add_argument("--run_7", action="store_true", help="GSEA Enrichment 분석 실행")
    parser.add_argument("--run_8", action="store_true", help="Pseudotime Gene Expression 분석 실행")
    parser.add_argument("--run_9", action="store_true", help="Cell Proportion 분석 실행")
    parser.add_argument("--run_10", action="store_true", help="Marker Gene Feature Plot 실행")

    # Trajectory 옵션
    parser.add_argument("--root_cluster", type=str, default=None, help="Trajectory 시작점 클러스터 ID")
    parser.add_argument("--root_gene", type=str, default=None, help="Trajectory 시작점 마커 유전자")

    # Cell Communication 옵션
    parser.add_argument("--n_perms", type=int, default=None, help="Cell communication permutation 수")
    parser.add_argument("--pvalue_threshold", type=float, default=None, help="유의성 임계값")

    # Marker Feature Plot 옵션
    parser.add_argument("--marker_genes", type=str, default=None, help="쉼표 구분 유전자 리스트 (예: CD3D,CD8A,MS4A1)")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    save_path = args.output_dir

    # =========================================================================
    # Step 1: Create Adata (전처리 - 항상 실행)
    # =========================================================================
    print("=" * 60)
    print("Step 1: Create Adata")
    print("=" * 60)
    adata = create_adata(
        parent_dir=args.input_dir,
        save_path=save_path,
        species=args.species,
        platform=args.platform,
    )

    # =========================================================================
    # Step 2: Batch Correction (항상 실행)
    # =========================================================================
    print("=" * 60)
    print("Step 2: Batch Correction")
    print("=" * 60)
    adata = batch_correction(
        adata,
        save_path=os.path.join(save_path, "BatchCorrection"),
    )

    # =========================================================================
    # Step 3: Cell Cycle (선택 실행)
    # =========================================================================
    if args.run_3:
        print("=" * 60)
        print("Step 3: Cell Cycle Scoring")
        print("=" * 60)
        adata = score_cell_cycle(
            adata,
            save_path=os.path.join(save_path, "CellCycle"),
            species=args.species,
        )
    else:
        print("Step 3 skipped")

    # =========================================================================
    # Step 4: DEG + Cell Type Prediction (선택 실행)
    # =========================================================================
    if args.run_4:
        print("=" * 60)
        print("Step 4: DEG + Cell Type Prediction")
        print("=" * 60)
        ensure_umap(adata)
        adata = deg_analysis(
            adata,
            save_path=os.path.join(save_path, "DEG"),
            species=args.species,
        )
    else:
        print("Step 4 skipped")

    # =========================================================================
    # Step 5: Trajectory Analysis (선택 실행)
    # =========================================================================
    if args.run_5:
        print("=" * 60)
        print("Step 5: Trajectory Analysis")
        print("=" * 60)
        adata = trajectory_analysis(
            adata,
            save_path=os.path.join(save_path, "Trajectory"),
            species=args.species,
            root_cluster=args.root_cluster,
            root_gene=args.root_gene,
        )
    else:
        print("Step 5 skipped")

    # =========================================================================
    # Step 6: Cell-Cell Communication (선택 실행)
    # =========================================================================
    if args.run_6:
        print("=" * 60)
        print("Step 6: Cell-Cell Communication")
        print("=" * 60)
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
    # Step 7: GSEA Enrichment (선택 실행, Step 4 필요)
    # =========================================================================
    if args.run_7:
        print("=" * 60)
        print("Step 7: GSEA Enrichment Analysis")
        print("=" * 60)
        if "wilcoxon" not in adata.uns:
            print("[GSEA] Skipped: DEG results not found. Run with --run_4 first.")
        else:
            adata = gsea_enrichment(
                adata,
                save_path=os.path.join(save_path, "GSEA"),
                species=args.species,
            )
    else:
        print("Step 7 skipped")

    # =========================================================================
    # Step 8: Pseudotime Gene Expression (선택 실행, Step 5 필요)
    # =========================================================================
    if args.run_8:
        print("=" * 60)
        print("Step 8: Pseudotime Gene Expression")
        print("=" * 60)
        if "dpt_pseudotime" not in adata.obs:
            print("[PseudotimeExpr] Skipped: Pseudotime not found. Run with --run_5 first.")
        else:
            adata = pseudotime_expression(
                adata,
                save_path=os.path.join(save_path, "PseudotimeExpression"),
            )
    else:
        print("Step 8 skipped")

    # =========================================================================
    # Step 9: Cell Proportion (선택 실행)
    # =========================================================================
    if args.run_9:
        print("=" * 60)
        print("Step 9: Cell Proportion Analysis")
        print("=" * 60)
        adata = cell_proportion_analysis(
            adata,
            save_path=os.path.join(save_path, "CellProportion"),
        )
    else:
        print("Step 9 skipped")

    # =========================================================================
    # Step 10: Marker Gene Feature Plot (선택 실행)
    # =========================================================================
    if args.run_10:
        print("=" * 60)
        print("Step 10: Marker Gene Feature Plot")
        print("=" * 60)
        marker_genes = None
        if args.marker_genes:
            marker_genes = [g.strip() for g in args.marker_genes.split(",")]
        adata = marker_feature_plot(
            adata,
            save_path=os.path.join(save_path, "MarkerFeature"),
            genes=marker_genes,
        )
    else:
        print("Step 10 skipped")

    # =========================================================================
    # 최종 저장
    # =========================================================================
    print("=" * 60)
    print("All steps completed.")
    print("=" * 60)
    adata.write_h5ad(os.path.join(save_path, args.output_filename), compression="gzip")
    print(f"Saved: {os.path.join(save_path, args.output_filename)}")


if __name__ == "__main__":
    main()
