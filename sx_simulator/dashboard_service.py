from typing import Optional

"""
Dashboard Service Module
========================
Streamlit 대시보드에서 사용하는 입력 묶음과 시뮬레이션 조립 로직을 담당합니다.

목표:
- UI 코드와 시뮬레이션 조립 코드를 분리
- 동일한 입력 세트를 메인 계산/비교 계산에서 재사용
- 리팩토링 중에도 기존 모델 거동은 유지
"""

from dataclasses import dataclass, field

from .config import VALIDATED_FIELD_WINDOW
from .multistage_sx import solve_multistage_countercurrent


@dataclass
class SimulationInputs:
    """대시보드가 수집한 현재 시뮬레이션 입력."""

    C_aq_feed: dict
    pH_feed: float
    Q_aq: float
    Q_org: float
    extractant: str
    C_ext: float
    n_stages: int
    temperature: float
    C_sulfate: float
    pH_mode: str
    target_pH: Optional[float] = None
    staged_pHs: Optional[list[float]] = None
    C_NaOH: float = 0.0
    Q_NaOH: float = 0.0
    naoh_mode: str = "aqueous_direct"
    naoh_wt_pct: Optional[float] = None
    saponification_model: str = "physical_v2"
    saponification_fraction: Optional[float] = None
    naoh_strategy: str = "uniform"
    naoh_weights: Optional[list[float]] = None
    metals: tuple[str, ...] = field(default_factory=tuple)


def build_simulation_kwargs(
    inputs: SimulationInputs,
    extractant_params: dict,
) -> dict:
    """엔진 호출용 kwargs를 현재 입력 상태에서 조립합니다."""
    sim_kwargs = dict(
        C_aq_feed=inputs.C_aq_feed,
        pH_feed=inputs.pH_feed,
        Q_aq=inputs.Q_aq,
        Q_org=inputs.Q_org,
        extractant=inputs.extractant,
        C_ext=inputs.C_ext,
        n_stages=inputs.n_stages,
        metals=list(inputs.metals),
        temperature=inputs.temperature,
        C_sulfate=inputs.C_sulfate,

        use_speciation=True,
        extractant_params=extractant_params,
        naoh_mode=inputs.naoh_mode,
        saponification_model=inputs.saponification_model,
    )

    if inputs.saponification_fraction is not None:
        sim_kwargs["saponification_fraction"] = inputs.saponification_fraction

    # 사포니피케이션 모드는 항상 단일 목표 pH (후액)만 지원하므로 staged_pHs를 무시합니다
    if inputs.staged_pHs and inputs.naoh_mode != "saponification":
        sim_kwargs["target_pH_per_stage"] = inputs.staged_pHs
    else:
        sim_kwargs["target_pH"] = inputs.target_pH
        
    if inputs.C_NaOH > 0:
        sim_kwargs["C_NaOH"] = inputs.C_NaOH
    if inputs.naoh_mode == "saponification" and inputs.Q_NaOH > 0:
        sim_kwargs["Q_NaOH"] = inputs.Q_NaOH

    return sim_kwargs


def run_simulation(inputs: SimulationInputs, extractant_params: dict) -> dict:
    """현재 입력으로 메인 시뮬레이션을 실행합니다."""
    return solve_multistage_countercurrent(
        **build_simulation_kwargs(inputs, extractant_params)
    )


def run_compare_simulations(
    inputs: SimulationInputs,
    compare_pH: float,
    extractant_params: dict,
) -> dict:
    """추출제 비교 탭에서 사용하는 공통 비교 계산."""
    results = {}
    for extractant in ("Cyanex 272", "D2EHPA"):
        compare_inputs = SimulationInputs(
            C_aq_feed=inputs.C_aq_feed,
            pH_feed=inputs.pH_feed,
            Q_aq=inputs.Q_aq,
            Q_org=inputs.Q_org,
            extractant=extractant,
            C_ext=inputs.C_ext,
            n_stages=inputs.n_stages,
            temperature=inputs.temperature,
            C_sulfate=inputs.C_sulfate,
            pH_mode="목표 pH (자동 NaOH)",
            target_pH=compare_pH,
            C_NaOH=inputs.C_NaOH,
            naoh_mode=inputs.naoh_mode,
            naoh_wt_pct=inputs.naoh_wt_pct,
            saponification_model=inputs.saponification_model,
            saponification_fraction=inputs.saponification_fraction,
            metals=inputs.metals,
        )
        results[extractant] = run_simulation(compare_inputs, extractant_params)
    return results


def compute_loading_pct(result: dict) -> float:
    """Stage 결과에서 최대 로딩률(%)을 계산합니다."""
    if not result.get("stages"):
        return 0.0
    max_loading = max(sr.get("loading_fraction", 0.0) for sr in result["stages"])
    return max_loading * 100.0


def build_scope_assessment(
    profile_key: str,
    inputs: SimulationInputs,
    result: dict,
) -> dict:
    """
    현재 입력이 Data1~6 현장 검증 범위와 얼마나 겹치는지 평가합니다.

    이는 절대적인 정오 판단이 아니라, 현재 조건의 "검증 커버리지"를
    사용자에게 알려주기 위한 안내 지표입니다.
    """
    highlights = []
    cautions = []
    score = 0

    eval_pH = (
        inputs.target_pH
        if inputs.target_pH is not None
        else result["pH_profile"][-1]
    )
    total_feed_gL = sum(inputs.C_aq_feed.values())
    loading_pct = compute_loading_pct(result)
    total_naoh_mol_hr = result["total_NaOH_mol_hr"]

    pH_min, pH_max = VALIDATED_FIELD_WINDOW["pH"]
    cext_min, cext_max = VALIDATED_FIELD_WINDOW["C_ext_m"]
    temp_ref = VALIDATED_FIELD_WINDOW["temperature_c"]
    temp_margin = VALIDATED_FIELD_WINDOW["temperature_margin_c"]
    validated_stages = VALIDATED_FIELD_WINDOW["n_stages"]

    if profile_key == "field_calibrated":
        highlights.append("현재 프로필은 Isd 108 현장 보정값을 기준으로 합니다.")
        score += 1
    else:
        cautions.append(
            "현재 프로필은 문헌 기본값입니다. 현장 프리셋(Data1~6)과의 직접 정합도는 낮아질 수 있습니다."
        )

    if pH_min <= eval_pH <= pH_max:
        score += 1
    else:
        cautions.append(
            f"현재 운전 pH {eval_pH:.2f}는 현장 검증 범위({pH_min:.1f}~{pH_max:.1f}) 밖입니다."
        )

    if cext_min <= inputs.C_ext <= cext_max:
        score += 1
    else:
        cautions.append(
            f"현재 추출제 농도 {inputs.C_ext:.4f} M는 현장 검증 범위({cext_min:.4f}~{cext_max:.4f} M) 밖입니다."
        )

    if abs(inputs.temperature - temp_ref) <= temp_margin:
        score += 1
    else:
        cautions.append(
            f"현재 온도 {inputs.temperature:.0f}°C는 25°C 기반 검증 범위를 크게 벗어납니다."
        )

    if inputs.n_stages == validated_stages:
        score += 1
    else:
        cautions.append(
            f"현재 {inputs.n_stages}단 조건은 현장 검증 기준({validated_stages}단)과 다릅니다."
        )

    if total_feed_gL > VALIDATED_FIELD_WINDOW["total_metals_warning_gL"]:
        cautions.append(
            f"총 금속 농도 {total_feed_gL:.1f} g/L는 고이온강도 구간입니다. 활동도/비이상성 미반영 영향이 커질 수 있습니다."
        )

    if loading_pct > 85.0:
        cautions.append(
            f"최대 로딩률 {loading_pct:.1f}%로 포화에 가깝습니다. 경쟁 추출 휴리스틱 민감도가 커지는 조건입니다."
        )

    if inputs.pH_mode == "목표 pH (자동 NaOH)":
        if inputs.naoh_mode == "saponification":
            highlights.append(
                "NaOH를 수계 직접 중화가 아니라 fresh organic sap inventory로 해석합니다."
            )
            if inputs.naoh_wt_pct is not None:
                highlights.append(
                    f"입력 NaOH 농도는 약 {inputs.naoh_wt_pct:.1f} wt% 기준으로 환산됩니다."
                )
        elif inputs.C_NaOH > 0:
            highlights.append(
                f"목표 pH 희석 추정이 활성화되어 있습니다. 가정 NaOH 농도 {inputs.C_NaOH:.1f} M 기준으로 수계 유량 증가를 함께 계산합니다."
            )
            cautions.append(
                "희석 추정은 입력 NaOH 농도 가정에 민감합니다. 실제 현장 농도와 다르면 후액 농도 예측도 달라질 수 있습니다."
            )
        elif eval_pH > 7.0:
            cautions.append(
                "목표 pH 모드는 NaOH 유량을 역산하므로 pH>7 영역에서 실제 희석 효과를 완전 반영하지 못합니다."
            )
    elif inputs.naoh_mode == "saponification":
        if inputs.saponification_model == "legacy_equivalent_target":
            highlights.append(
                "고정 NaOH + 사포니피케이션은 raw-feed 기준 fresh organic sap condition을 equivalent target pH로 환산해 계산합니다."
            )
            cautions.append(
                "이 경로는 현장 replay 안정성을 위한 field-calibrated legacy fallback입니다."
            )
        else:
            highlights.append(
                "고정 NaOH + 사포니피케이션은 fresh organic sap inventory를 stage 간 counter-current로 전달하는 physical v2 경로를 사용합니다."
            )

    if total_naoh_mol_hr > max(25.0, 0.5 * sum(inputs.C_aq_feed.values())):
        highlights.append(
            "알칼리 수요가 큰 조건입니다. 실제 설비 적용 전에는 고정 NaOH 모드 또는 현장 시험으로 재확인하는 편이 안전합니다."
        )

    if score >= 5 and not cautions:
        level = "high"
        title = "현장 검증 범위 안"
    elif score >= 3:
        level = "medium"
        title = "부분 검증 범위"
    else:
        level = "low"
        title = "검증 범위 밖"

    return {
        "level": level,
        "title": title,
        "highlights": highlights,
        "cautions": cautions,
    }
