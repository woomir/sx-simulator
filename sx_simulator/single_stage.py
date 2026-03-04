"""
Single Stage Module
===================
1개 Mixer-Settler Stage에서의 물질수지 및 pH 변화를 계산하는 모듈.

주요 기능:
1. 수계/유기계 간 금속 분배 (분배 계수 D 기반)
2. pH 변화 추적 (H⁺ 수지 — 금속 추출에 의한 H⁺ 방출 + NaOH 중화)
3. 목표 pH 모드: 원하는 pH를 유지하기 위한 NaOH 자동 계산
4. 추출제 로딩 한계: 유기상 로딩률에 따른 D값 감쇠 (금속 경쟁)
"""

import math
from .extraction_isotherm import (
    distribution_coefficient,
    get_proton_release,
    calc_loading_fraction,
    loading_damping_factor,
)
from .config import MOLAR_MASS, DEFAULT_METALS


def solve_single_stage(
    C_aq_in: dict,
    C_org_in: dict,
    pH_in: float,
    Q_aq: float,
    Q_org: float,
    extractant: str,
    C_ext: float,
    C_NaOH: float = 0.0,
    Q_NaOH: float = 0.0,
    target_pH: float = None,
    metals: list = None,
    temperature: float = None,
) -> dict:
    """
    단일 Mixer-Settler stage의 평형 계산을 수행합니다.

    두 가지 모드를 지원합니다:
    1. **고정 NaOH 모드** (target_pH = None): 주어진 NaOH 유량/농도로 계산
    2. **목표 pH 모드** (target_pH = 값): 지정한 pH를 달성하는 데 필요한 NaOH를 자동 계산

    추출제 로딩 한계가 활성화되어 있으며, 유기상에 이미 로딩된 금속이
    추출제 용량에 가까우면 추출 효율이 자동으로 감소합니다.

    Parameters
    ----------
    C_aq_in : dict   — 수계 입구 금속 농도 {금속: g/L}
    C_org_in : dict  — 유기계 입구 금속 농도 {금속: g/L}
    pH_in : float    — 수계 입구 pH
    Q_aq : float     — 수계 유량 (L/hr)
    Q_org : float    — 유기계 유량 (L/hr)
    extractant : str — 추출제 종류
    C_ext : float    — 추출제 농도 (M)
    C_NaOH : float   — NaOH 농도 (M), 고정 NaOH 모드용
    Q_NaOH : float   — NaOH 유량 (L/hr), 고정 NaOH 모드용
    target_pH : float — 목표 pH (None이면 고정 NaOH 모드)
    metals : list    — 계산할 금속 리스트
    temperature : float — 온도 (°C)

    Returns
    -------
    dict:
        C_aq_out, C_org_out, pH_out, extraction, H_released_mol_hr,
        NaOH_consumed_mol_hr, loading_fraction, iterations
    """
    if metals is None:
        metals = DEFAULT_METALS

    use_target_pH = target_pH is not None

    # 목표 pH 모드: 목표 pH에서의 평형을 직접 계산
    if use_target_pH:
        return _solve_at_target_pH(
            C_aq_in, C_org_in, pH_in, Q_aq, Q_org,
            extractant, C_ext, target_pH, metals, temperature
        )
    else:
        return _solve_with_fixed_NaOH(
            C_aq_in, C_org_in, pH_in, Q_aq, Q_org,
            extractant, C_ext, C_NaOH, Q_NaOH, metals, temperature
        )


def _solve_at_target_pH(C_aq_in, C_org_in, pH_in, Q_aq, Q_org,
                         extractant, C_ext, target_pH, metals, temperature=None):
    """
    목표 pH 모드: 지정된 pH에서 평형을 계산하고, 필요한 NaOH를 역산합니다.
    로딩 한계를 반영하여 반복 수렴합니다.
    """
    C_aq_out = {}
    C_org_out = {}
    extraction = {}
    total_H_released = 0.0

    max_loading_iter = 30
    loading_tol = 1e-4
    damping = 1.0  # 시작 감쇠 = 1 (제한 없음)
    relaxation = 0.5  # under-relaxation 계수

    for load_iter in range(max_loading_iter):
        C_aq_out = {}
        C_org_out = {}
        extraction = {}
        total_H_released = 0.0

        # 목표 pH에서 각 금속의 분배 계수를 계산 (로딩 감쇠 적용)
        for metal in metals:
            D_raw = distribution_coefficient(target_pH, metal, extractant, C_ext,
                                             temperature=temperature)
            D = D_raw * damping  # 로딩에 의한 감쇠
            MW = MOLAR_MASS[metal]

            C_aq_mol_in = C_aq_in.get(metal, 0.0) / MW
            C_org_mol_in = C_org_in.get(metal, 0.0) / MW

            total_metal_flow = C_aq_mol_in * Q_aq + C_org_mol_in * Q_org

            # 물질수지: D = C_org_out / C_aq_out
            denominator = Q_aq + D * Q_org
            if denominator <= 0:
                C_aq_mol_out = 0.0
            else:
                C_aq_mol_out = total_metal_flow / denominator

            C_org_mol_out = D * C_aq_mol_out

            C_aq_out[metal] = C_aq_mol_out * MW
            C_org_out[metal] = C_org_mol_out * MW

            # 추출된 금속량
            metal_extracted = (C_aq_mol_in * Q_aq) - (C_aq_mol_out * Q_aq)
            if C_aq_mol_in > 0 and C_aq_mol_in * Q_aq > 0:
                extraction[metal] = max(0.0, metal_extracted / (C_aq_mol_in * Q_aq) * 100.0)
            else:
                extraction[metal] = 0.0

            # H⁺ 방출량 계산
            n_H = get_proton_release(metal, extractant)
            total_H_released += n_H * max(0.0, metal_extracted)

        # 새로운 로딩률과 감쇠 계산
        new_loading = calc_loading_fraction(C_org_out, extractant, C_ext, metals)
        new_damping = loading_damping_factor(new_loading)

        # under-relaxation으로 진동 방지
        damping_prev = damping
        damping = relaxation * new_damping + (1 - relaxation) * damping

        # 수렴 확인
        if load_iter > 0 and abs(damping - damping_prev) < loading_tol:
            break

    # 최종 로딩률 계산
    final_loading = calc_loading_fraction(C_org_out, extractant, C_ext, metals)

    # NaOH 필요량 역산
    H_in = (10 ** (-pH_in)) * Q_aq
    H_out_target = (10 ** (-target_pH)) * Q_aq
    NaOH_consumed = H_in + total_H_released - H_out_target
    NaOH_consumed = max(0.0, NaOH_consumed)

    return {
        "C_aq_out": C_aq_out,
        "C_org_out": C_org_out,
        "pH_out": round(target_pH, 4),
        "extraction": extraction,
        "H_released_mol_hr": total_H_released,
        "NaOH_consumed_mol_hr": NaOH_consumed,
        "loading_fraction": final_loading,
        "iterations": load_iter + 1,
    }


def _solve_with_fixed_NaOH(C_aq_in, C_org_in, pH_in, Q_aq, Q_org,
                             extractant, C_ext, C_NaOH, Q_NaOH, metals, temperature=None):
    """
    고정 NaOH 모드: 주어진 NaOH로 pH를 반복 계산합니다.
    로딩 한계를 반영합니다.
    """
    C_aq_out = {}
    C_org_out = {}
    extraction = {}
    total_H_released = 0.0
    Q_aq_eff = Q_aq + Q_NaOH

    pH_current = pH_in
    max_iter = 100
    tol = 1e-4

    for iteration in range(max_iter):
        pH_prev = pH_current
        total_H_released = 0.0

        # 현재 반복의 유기상 로딩률 계산
        loading_frac = calc_loading_fraction(C_org_out if C_org_out else C_org_in,
                                              extractant, C_ext, metals)
        damping = loading_damping_factor(loading_frac)

        for metal in metals:
            D_raw = distribution_coefficient(pH_current, metal, extractant, C_ext,
                                             temperature=temperature)
            D = D_raw * damping  # 로딩에 의한 감쇠
            MW = MOLAR_MASS[metal]
            C_aq_mol_in = C_aq_in.get(metal, 0.0) / MW
            C_org_mol_in = C_org_in.get(metal, 0.0) / MW

            total_metal_flow = C_aq_mol_in * Q_aq + C_org_mol_in * Q_org
            denominator = Q_aq_eff + D * Q_org
            C_aq_mol_out = total_metal_flow / denominator if denominator > 0 else 0.0
            C_org_mol_out = D * C_aq_mol_out

            C_aq_out[metal] = C_aq_mol_out * MW
            C_org_out[metal] = C_org_mol_out * MW

            metal_extracted = (C_aq_mol_in * Q_aq) - (C_aq_mol_out * Q_aq_eff)
            if C_aq_mol_in > 0:
                extraction[metal] = max(0.0, metal_extracted / (C_aq_mol_in * Q_aq) * 100.0)
            else:
                extraction[metal] = 0.0

            n_H = get_proton_release(metal, extractant)
            total_H_released += n_H * max(0.0, metal_extracted)

        # pH 계산
        H_in = (10 ** (-pH_in)) * Q_aq
        OH_added = C_NaOH * Q_NaOH
        H_out = H_in + total_H_released - OH_added

        if Q_aq_eff > 0 and H_out > 0:
            pH_current = -math.log10(max(H_out / Q_aq_eff, 1e-14))
        elif H_out <= 0:
            OH_excess = -H_out / Q_aq_eff if Q_aq_eff > 0 else 0
            pH_current = 14.0 + math.log10(OH_excess) if OH_excess > 0 else 7.0
        else:
            pH_current = 7.0

        pH_current = max(0.0, min(14.0, pH_current))
        if abs(pH_current - pH_prev) < tol:
            break

    # 최종 로딩률 계산
    final_loading = calc_loading_fraction(C_org_out, extractant, C_ext, metals)

    return {
        "C_aq_out": C_aq_out,
        "C_org_out": C_org_out,
        "pH_out": round(pH_current, 4),
        "extraction": extraction,
        "H_released_mol_hr": total_H_released,
        "NaOH_consumed_mol_hr": C_NaOH * Q_NaOH,
        "loading_fraction": final_loading,
        "iterations": iteration + 1,
    }


def print_stage_result(result: dict, stage_num: int = 1):
    """단일 stage 결과를 포맷팅하여 출력합니다."""
    print(f"\n{'='*55}")
    print(f"  Stage {stage_num} 결과")
    print(f"{'='*55}")
    print(f"  출구 pH: {result['pH_out']:.2f}")
    print(f"  NaOH 소비량: {result.get('NaOH_consumed_mol_hr', 0):.4f} mol/hr")
    print(f"  H⁺ 방출량: {result['H_released_mol_hr']:.4f} mol/hr")

    print(f"\n  {'금속':>6} | {'수계출구(g/L)':>13} | {'유기계출구(g/L)':>15} | {'추출률(%)':>9}")
    print(f"  {'-'*6}-+-{'-'*13}-+-{'-'*15}-+-{'-'*9}")
    for metal in result["C_aq_out"]:
        c_aq = result["C_aq_out"][metal]
        c_org = result["C_org_out"][metal]
        ext = result["extraction"][metal]
        print(f"  {metal:>6} | {c_aq:>13.4f} | {c_org:>15.4f} | {ext:>9.2f}")
    print(f"{'='*55}")
