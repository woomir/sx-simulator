"""
Single Stage Module
===================
1개 Mixer-Settler Stage에서의 물질수지 및 pH 변화를 계산하는 모듈.

주요 기능:
1. 수계/유기계 간 금속 분배 (분배 계수 D 기반)
2. pH 변화 추적 (H⁺ 수지 — 금속 추출에 의한 H⁺ 방출 + NaOH 중화)
3. 목표 pH 모드: 원하는 pH를 유지하기 위한 NaOH 자동 계산
4. 추출제 로딩 한계: 유기상 로딩률에 따른 D값 감쇠 (금속 경쟁)
5. 추출제 경쟁 모드 (Phase 3): 공유 추출제 풀에서 금속 간 경쟁 반영
6. 수계 종분화 (Phase 3): MOH⁺, MSO₄⁰ 종의 pH buffer 효과 반영
"""

from .extraction_isotherm import (
    distribution_coefficient,
    get_proton_release,
    calc_loading_fraction,
    loading_damping_factor,
    calc_free_NaL,
    compute_competitive_extractions,
    get_aqueous_speciation_state,
)
from .config import (MOLAR_MASS, DEFAULT_METALS,
                     SPECIATION_CONSTANTS)


def _partition_stage_with_damping(
    pH: float,
    C_aq_in: dict,
    C_org_in: dict,
    Q_aq_in: float,
    Q_aq_out: float,
    Q_org: float,
    extractant: str,
    C_ext: float,
    metals: list,
    temperature: float,
    damping: float,
    extractant_params: dict = None,
) -> dict:
    """고정된 damping 값으로 금속 분배와 H+ 방출량을 계산합니다."""
    C_aq_out = {}
    C_org_out = {}
    extraction = {}
    total_H_released = 0.0

    for metal in metals:
        D_raw = distribution_coefficient(
            pH, metal, extractant, C_ext, temperature=temperature,
            extractant_params=extractant_params
        )
        D = D_raw * damping
        MW = MOLAR_MASS[metal]

        C_aq_mol_in = C_aq_in.get(metal, 0.0) / MW
        C_org_mol_in = C_org_in.get(metal, 0.0) / MW
        total_metal_flow = C_aq_mol_in * Q_aq_in + C_org_mol_in * Q_org
        denominator = Q_aq_out + D * Q_org

        C_aq_mol_out = total_metal_flow / denominator if denominator > 0 else 0.0
        C_org_mol_out = D * C_aq_mol_out

        C_aq_out[metal] = C_aq_mol_out * MW
        C_org_out[metal] = C_org_mol_out * MW

        metal_extracted = (C_aq_mol_in * Q_aq_in) - (C_aq_mol_out * Q_aq_out)
        if C_aq_mol_in > 0:
            extraction[metal] = max(0.0, metal_extracted / (C_aq_mol_in * Q_aq_in) * 100.0)
        else:
            extraction[metal] = 0.0

        n_H = get_proton_release(metal, extractant, extractant_params)
        total_H_released += n_H * max(0.0, metal_extracted)

    return {
        "C_aq_out": C_aq_out,
        "C_org_out": C_org_out,
        "extraction": extraction,
        "H_released_mol_hr": total_H_released,
    }


def _solve_competitive_stage_state(
    pH: float,
    C_aq_in: dict,
    C_org_in: dict,
    Q_aq_in: float,
    Q_aq_out: float,
    Q_org: float,
    extractant: str,
    C_ext: float,
    metals: list,
    temperature: float,
    extractant_params: dict = None,
) -> dict:
    """경쟁 추출 모드의 단일 stage 평형 상태를 계산합니다."""
    result = compute_competitive_extractions(
        pH, extractant, C_ext, C_aq_in, C_org_in,
        Q_aq_in, Q_org, metals, temperature, Q_aq_eff=Q_aq_out,
        extractant_params=extractant_params
    )

    total_H_released = 0.0
    for metal in metals:
        MW = MOLAR_MASS[metal]
        C_aq_mol_in = C_aq_in.get(metal, 0.0) / MW
        C_aq_mol_out = result["C_aq_out"].get(metal, 0.0) / MW
        metal_extracted = (C_aq_mol_in * Q_aq_in) - (C_aq_mol_out * Q_aq_out)
        n_H = get_proton_release(metal, extractant, extractant_params)
        total_H_released += n_H * max(0.0, metal_extracted)

    return {
        "C_aq_out": result["C_aq_out"],
        "C_org_out": result["C_org_out"],
        "extraction": result["extraction"],
        "H_released_mol_hr": total_H_released,
        "loading_fraction": calc_loading_fraction(
            result["C_org_out"], extractant, C_ext, metals, extractant_params
        ),
    }


def solve_single_stage(
    C_aq_in: dict,
    C_org_in: dict,
    pH_in: float,
    Q_aq: float,
    Q_org: float,
    extractant: str,
    C_ext: float,
    target_pH: float = None,
    C_NaOH: float = 0.0,
    Q_NaOH: float = 0.0,
    metals: list = None,
    temperature: float = None,
    C_sulfate: float = 0.0,
    use_competition: bool = False,
    use_speciation: bool = False,
    extractant_params: dict = None,
) -> dict:
    """
    단일 Mixer-Settler stage의 평형 계산을 수행합니다.

    두 가지 모드를 지원합니다:
    1. **고정 NaOH 모드** (target_pH = None): 주어진 NaOH 유량/농도로 계산
    2. **목표 pH 모드** (target_pH = 값): 지정한 pH를 달성하는 데 필요한 NaOH를 자동 계산

    고급 옵션 (Phase 3):
    - use_competition: 추출제 경쟁 반영 (공유 추출제 풀)
    - use_speciation: 수계 종분화 반영 (MOH⁺, MSO₄⁰)
    
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
    use_competition : bool — 추출제 경쟁 반영 여부
    use_speciation : bool — 수계 종분화 반영 여부

    Returns
    -------
    dict:
        C_aq_out, C_org_out, pH_out, extraction, H_released_mol_hr,
        NaOH_consumed_mol_hr, loading_fraction, iterations
    """
    if metals is None:
        metals = DEFAULT_METALS

    # 목표 pH 모드: 목표 pH에서의 평형을 직접 계산
    if target_pH is not None:
        return _solve_at_target_pH(
            C_aq_in, C_org_in, pH_in, Q_aq, Q_org,
            extractant, C_ext, target_pH, metals,
            temperature, C_sulfate,
            use_competition, use_speciation, extractant_params
        )
    else:
        return _solve_with_fixed_NaOH(
            C_aq_in, C_org_in, pH_in, Q_aq, Q_org,
            extractant, C_ext, C_NaOH, Q_NaOH, metals,
            temperature, C_sulfate,
            use_competition, use_speciation, extractant_params
        )

def calc_aq_protons(pH: float, Q_aq: float, C_sulfate: float, C_aq: dict = None, use_speciation: bool = False) -> float:
    """
    수계의 총 해리가능 양성자(H+) 몰스피드(mol/hr)를 계산합니다.
    황산(HSO4-) 수용액의 버퍼 효과(pKa=1.99) 포함.
    Phase 3: 종분화(MOH⁺)에 의한 OH⁻ 소비(즉, H⁺ 생성과 동등한 효과) 반영.
    또한 황산염 착물(MSO4^0)이 자유 금속 분율을 낮춰 hydrolysis를 완화하는
    효과를 간접 반영합니다.
    """
    free_H = 10.0 ** (-pH)
    Ka_HSO4 = 10.0 ** (-1.99)
    # [HSO4-] = [H+] * C_SO4_total / (Ka + [H+])
    bound_H = free_H * C_sulfate / (Ka_HSO4 + free_H) if C_sulfate > 0 else 0.0
    
    hydrolysis_H = 0.0
    if use_speciation and C_aq:
        for metal, conc_gL in C_aq.items():
            if metal in SPECIATION_CONSTANTS:
                MW = MOLAR_MASS[metal]
                C_M_total = conc_gL / MW
                speciation_state = get_aqueous_speciation_state(
                    metal, pH, C_sulfate
                )
                conc_MOH = C_M_total * speciation_state["hydroxo_fraction"]
                hydrolysis_H += conc_MOH
                
    return (free_H + bound_H + hydrolysis_H) * Q_aq


def _solve_at_target_pH(C_aq_in, C_org_in, pH_in, Q_aq, Q_org,
                          extractant, C_ext, target_pH, metals,
                          temperature=None, C_sulfate=0.0,
                          use_competition=False, use_speciation=False,
                          extractant_params=None):
    """
    목표 pH 모드: 지정된 pH에서 평형을 계산하고, 필요한 NaOH를 역산합니다.
    - use_competition ON: compute_competitive_extractions 1회 호출
    - use_competition OFF: 기존 sigmoid 감쇠 기반 반복 계산

    주의사항 (구조적 한계):
    - NaOH 역산 구조상 NaOH 유량을 사전에 알 수 없으므로, Q_aq_eff = Q_aq으로
      계산합니다. NaOH 투입에 의한 수계 희석 효과는 반영되지 않습니다.
    - pH > 7 OH- 보정에서도 동일한 이유로 Q_aq를 사용합니다.
    - 정밀한 희석 반영이 필요한 경우 고정 NaOH 모드(_solve_with_fixed_NaOH)를
      사용하십시오.
    """
    total_H_released = 0.0

    if use_competition:
        stage_state = _solve_competitive_stage_state(
            target_pH, C_aq_in, C_org_in, Q_aq, Q_aq, Q_org,
            extractant, C_ext, metals, temperature, extractant_params
        )
        C_aq_out = stage_state["C_aq_out"]
        C_org_out = stage_state["C_org_out"]
        extraction = stage_state["extraction"]
        total_H_released = stage_state["H_released_mol_hr"]
        load_iter = 1
        final_loading = stage_state["loading_fraction"]
    else:
        max_loading_iter = 30
        loading_tol = 1e-4
        damping = 1.0
        relaxation = 0.5

        for load_iter in range(max_loading_iter):
            stage_state = _partition_stage_with_damping(
                target_pH, C_aq_in, C_org_in, Q_aq, Q_aq, Q_org,
                extractant, C_ext, metals, temperature, damping, extractant_params
            )
            C_aq_out = stage_state["C_aq_out"]
            C_org_out = stage_state["C_org_out"]
            extraction = stage_state["extraction"]
            total_H_released = stage_state["H_released_mol_hr"]

            new_loading = calc_loading_fraction(
                C_org_out, extractant, C_ext, metals, extractant_params
            )
            new_damping = loading_damping_factor(new_loading)
            damping_prev = damping
            damping = relaxation * new_damping + (1 - relaxation) * damping
            if load_iter > 0 and abs(damping - damping_prev) < loading_tol:
                break

        final_loading = calc_loading_fraction(
            C_org_out, extractant, C_ext, metals, extractant_params
        )

    # 추출제 사포닌화에 의한 NaOH 간접 소모량 계산
    saponification_OH_consumed = 0.0
    NaL_free = calc_free_NaL(
        target_pH, C_ext, extractant, C_org_out, metals, extractant_params
    )
    saponification_OH_consumed = NaL_free * Q_org

    # NaOH 필요량 역산
    total_H_in = calc_aq_protons(pH_in, Q_aq, C_sulfate, C_aq_in, use_speciation)
    total_H_out_target = calc_aq_protons(target_pH, Q_aq, C_sulfate, C_aq_out, use_speciation)
    
    if target_pH > 7.0:
        total_H_out_target -= Q_aq * (10.0 ** -(14.0 - target_pH))

    NaOH_consumed = total_H_in + total_H_released - total_H_out_target + saponification_OH_consumed
    NaOH_consumed = max(0.0, NaOH_consumed)

    return {
        "C_aq_out": C_aq_out,
        "C_org_out": C_org_out,
        "pH_out": round(target_pH, 4),
        "extraction": extraction,
        "H_released_mol_hr": total_H_released,
        "NaOH_consumed_mol_hr": NaOH_consumed,
        "loading_fraction": final_loading,
        "iterations": load_iter + 1 if not use_competition else 1,
    }


def _solve_with_fixed_NaOH(C_aq_in, C_org_in, pH_in, Q_aq, Q_org,
                             extractant, C_ext, C_NaOH, Q_NaOH, metals,
                             temperature=None, C_sulfate=0.0,
                             use_competition=False, use_speciation=False,
                             extractant_params=None):
    """
    고정 NaOH 모드: 주어진 NaOH로 pH를 반복 계산합니다.
    """
    C_aq_out = {}
    C_org_out = {}
    extraction = {}
    total_H_released = 0.0
    Q_aq_eff = Q_aq + Q_NaOH

    pH_low = -1.0
    pH_high = 14.0
    pH_current = pH_in
    max_iter = 100
    tol = 1e-4

    for iteration in range(max_iter):
        pH_current = (pH_low + pH_high) / 2.0
        total_H_released = 0.0

        if use_competition:
            stage_state = _solve_competitive_stage_state(
                pH_current, C_aq_in, C_org_in, Q_aq, Q_aq_eff, Q_org,
                extractant, C_ext, metals, temperature, extractant_params
            )
            C_aq_out = stage_state["C_aq_out"]
            C_org_out = stage_state["C_org_out"]
            extraction = stage_state["extraction"]
            total_H_released = stage_state["H_released_mol_hr"]
        else:
            # 로딩 감쇠 내부 반복 (target_pH 모드와 동일한 패턴)
            damping = 1.0
            for _load_iter in range(15):
                trial_state = _partition_stage_with_damping(
                    pH_current, C_aq_in, C_org_in, Q_aq, Q_aq_eff, Q_org,
                    extractant, C_ext, metals, temperature, damping, extractant_params
                )
                new_loading = calc_loading_fraction(
                    trial_state["C_org_out"], extractant, C_ext, metals, extractant_params
                )
                new_damping = loading_damping_factor(new_loading)
                if abs(new_damping - damping) < 1e-4:
                    damping = new_damping
                    break
                damping = 0.5 * new_damping + 0.5 * damping

            # 확정된 damping으로 최종 계산
            stage_state = _partition_stage_with_damping(
                pH_current, C_aq_in, C_org_in, Q_aq, Q_aq_eff, Q_org,
                extractant, C_ext, metals, temperature, damping, extractant_params
            )
            C_aq_out = stage_state["C_aq_out"]
            C_org_out = stage_state["C_org_out"]
            extraction = stage_state["extraction"]
            total_H_released = stage_state["H_released_mol_hr"]

        total_H_in = calc_aq_protons(pH_in, Q_aq, C_sulfate, C_aq_in, use_speciation)
        OH_added = C_NaOH * Q_NaOH
        
        NaL_free = calc_free_NaL(
            pH_current, C_ext, extractant, C_org_out, metals, extractant_params
        )
        saponification_OH_consumed = NaL_free * Q_org
        OH_effective = max(0.0, OH_added - saponification_OH_consumed)

        H_out_balance = total_H_in + total_H_released - OH_effective
        
        H_out_actual = calc_aq_protons(pH_current, Q_aq_eff, C_sulfate, C_aq_out, use_speciation)
        if pH_current > 7.0:
            H_out_actual -= Q_aq_eff * (10.0 ** -(14.0 - pH_current))
            
        error = H_out_balance - H_out_actual

        if error > 0:
            pH_high = pH_current
        else:
            pH_low = pH_current

        if (pH_high - pH_low) < tol:
            break

    final_loading = calc_loading_fraction(
        C_org_out, extractant, C_ext, metals, extractant_params
    )

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
