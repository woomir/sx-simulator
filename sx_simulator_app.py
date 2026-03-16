#!/usr/bin/env python3
"""
SX Simulator App — 메인 실행 파일
=================================
Li, Ni, Co, Mn Mixer-Settler 용매추출 시뮬레이션.

사용법: python3 sx_simulator_app.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sx_simulator.extraction_isotherm import print_isotherm_table
from sx_simulator.single_stage import print_stage_result
from sx_simulator.multistage_sx import solve_multistage_countercurrent, print_multistage_result


def main():
    print("=" * 70)
    print("  Li/Ni/Co/Mn Mixer-Settler SX 시뮬레이터")
    print("  Based on Vasilyev Competitive Extraction Model (Vasilyev et al. 2019)")
    print("=" * 70)

    # =================================================================
    # SIMULATION PARAMETERS — 여기를 수정하세요
    # =================================================================

    # 1. Feed (수계) 조건
    C_aq_feed = {
        "Li": 1.5,    # g/L
        "Ni": 5.0,    # g/L
        "Co": 3.0,    # g/L
        "Mn": 2.0,    # g/L
    }
    pH_feed = 3.0
    Q_aq = 100.0          # L/hr

    # 2. 용매 (유기계) 조건
    extractant = "Cyanex 272"   # "Cyanex 272" 또는 "D2EHPA"
    C_ext = 0.5                 # 추출제 농도 (M)
    Q_org = 100.0               # L/hr  (A/O = 1:1)

    # 3. Mixer-Settler 설정
    n_stages = 4

    # 4. pH 제어 모드 (목표 pH)
    #    - 모든 stage 동일 pH:  target_pH = 5.0
    #    - stage별 개별 pH:     target_pH_per_stage = [4.0, 4.5, 5.0, 5.5]
    #    - None이면 고정 NaOH 모드
    target_pH = 5.0             # 모든 stage에서 pH 5.0 유지 목표

    metals = ["Li", "Ni", "Co", "Mn"]

    # =================================================================

    # [1] pH-추출률 Isotherm 테이블
    print("\n\n[1] pH-추출률 Isotherm (열역학 평형 기준)")
    print_isotherm_table(extractant, C_ext, pH_range=(1.0, 10.0), pH_step=0.5)

    # [2] 균일 목표 pH 시뮬레이션
    print(f"\n\n[2] {n_stages}-Stage 역류 시뮬레이션 (목표 pH = {target_pH})")
    print(f"    추출제: {extractant} ({C_ext}M), A/O = {Q_aq/Q_org:.1f}")

    result1 = solve_multistage_countercurrent(
        C_aq_feed=C_aq_feed, pH_feed=pH_feed,
        Q_aq=Q_aq, Q_org=Q_org,
        extractant=extractant, C_ext=C_ext, n_stages=n_stages,
        target_pH=target_pH, metals=metals,
    )
    print_multistage_result(result1, C_aq_feed)

    # Stage별 상세
    for i, sr in enumerate(result1["stages"]):
        print_stage_result(sr, i + 1)

    # [3] Stage별 차등 pH 시뮬레이션 (Mn/Co 선택 추출 예시)
    staged_pHs = [3.5, 4.0, 4.5, 5.0]
    print(f"\n\n[3] {n_stages}-Stage 역류 시뮬레이션 (차등 pH: {staged_pHs})")
    print(f"    → 낮은 pH에서 Mn 선택 추출 → 높은 pH에서 Co 추출")

    result2 = solve_multistage_countercurrent(
        C_aq_feed=C_aq_feed, pH_feed=pH_feed,
        Q_aq=Q_aq, Q_org=Q_org,
        extractant=extractant, C_ext=C_ext, n_stages=n_stages,
        target_pH_per_stage=staged_pHs, metals=metals,
    )
    print_multistage_result(result2, C_aq_feed)

    for i, sr in enumerate(result2["stages"]):
        print_stage_result(sr, i + 1)

    # [4] D2EHPA 비교 시뮬레이션
    print(f"\n\n[4] D2EHPA 비교: {n_stages}-Stage (목표 pH = {target_pH})")
    result3 = solve_multistage_countercurrent(
        C_aq_feed=C_aq_feed, pH_feed=pH_feed,
        Q_aq=Q_aq, Q_org=Q_org,
        extractant="D2EHPA", C_ext=0.64, n_stages=n_stages,
        target_pH=target_pH, metals=metals,
    )
    print_multistage_result(result3, C_aq_feed)

    print("\n\n시뮬레이션 완료!")
    print("=" * 70)


if __name__ == "__main__":
    main()
