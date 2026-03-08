"""
Extraction Isotherm Module
==========================
pH에 따른 금속별 추출률(%) 및 분배 계수(D)를 계산하는 핵심 모듈.

시그모이드 함수로 pH-추출률 관계를 근사:
  E(pH, T) = E_max / (1 + exp(-k(T) · (pH - pH50(T))))

온도 보정:
  pH50(T) = pH50(T_ref) + beta · (T - T_ref) + alpha · log(C_ext/C_ref)
  k(T)        = k_ref * exp(gamma * (T - T_ref))

경쟁 추출 (Phase 3 — Mechanistic Light):
  D_adjusted(M) = D_sigmoid(M) × ([HA]_free / C_ext)^n_ext(M)
  Vasilyev et al. (2019)의 공유 추출제 풀 개념을 시그모이드와 결합
"""

import math
from .config import (EXTRACTANT_PARAMS, MOLAR_MASS, DEFAULT_METALS,
                     T_REF, DEFAULT_TEMPERATURE, MAX_LOADING_FRACTION,
                     EXTRACTION_PRIORITY)


def _resolve_extractant_params(extractant_params: dict = None) -> dict:
    """명시적 파라미터 세트가 없으면 기본 전역 설정을 사용합니다."""
    return extractant_params if extractant_params is not None else EXTRACTANT_PARAMS


def _get_metal_params(
    metal: str,
    extractant: str,
    extractant_params: dict = None,
) -> dict:
    """추출제/금속 조합의 파라미터를 반환하고 입력을 검증합니다."""
    params_store = _resolve_extractant_params(extractant_params)
    if extractant not in params_store:
        raise ValueError(f"지원하지 않는 추출제: {extractant}")
    if metal not in params_store[extractant]:
        raise ValueError(f"지원하지 않는 금속: {metal} (추출제: {extractant})")
    return params_store[extractant][metal]


def get_effective_pH50(metal: str, extractant: str, C_ext: float,
                       temperature: float = None,
                       extractant_params: dict = None) -> float:
    """
    추출제 농도와 온도에 따라 보정된 pH50 값을 반환합니다.

    Parameters
    ----------
    metal : str
        금속 종류 ("Li", "Ni", "Co", "Mn")
    extractant : str
        추출제 종류 ("Cyanex 272", "D2EHPA")
    C_ext : float
        추출제 농도 (M, mol/L)
    temperature : float, optional
        온도 (°C). None이면 T_REF(25°C) 사용

    Returns
    -------
    float
        보정된 pH50 값
    """
    params = _get_metal_params(metal, extractant, extractant_params)
    pH50_ref = params["pH50"]
    alpha = params["alpha"]
    C_ref = params["C_ref"]
    beta = params.get("beta", 0.0)

    if C_ext <= 0:
        raise ValueError(f"추출제 농도는 양수여야 합니다: {C_ext}")

    T = temperature if temperature is not None else T_REF

    # 추출제 농도 보정 + 온도 보정
    pH50_eff = pH50_ref - alpha * math.log10(C_ext / C_ref) + beta * (T - T_REF)
    return pH50_eff


def get_effective_k(metal: str, extractant: str,
                    temperature: float = None,
                    extractant_params: dict = None) -> float:
    """
    온도에 따라 보정된 시그모이드 기울기(k) 값을 반환합니다.

    k(T) = k_ref * exp(gamma * (T - T_ref))

    Parameters
    ----------
    metal : str
    extractant : str
    temperature : float, optional
        온도 (°C). None이면 T_REF(25°C) 사용

    Returns
    -------
    float
        보정된 k 값
    """
    params = _get_metal_params(metal, extractant, extractant_params)
    k_ref = params["k"]
    gamma = params.get("gamma", 0.0)

    T = temperature if temperature is not None else T_REF
    k_eff = k_ref * math.exp(gamma * (T - T_REF))
    return k_eff


def extraction_efficiency(pH: float, metal: str, extractant: str, C_ext: float,
                          temperature: float = None,
                          extractant_params: dict = None) -> float:
    """
    주어진 pH와 온도에서 특정 금속의 추출률(%)을 계산합니다.

    Parameters
    ----------
    pH : float
        수계상의 현재 pH
    metal : str
        금속 종류 ("Li", "Ni", "Co", "Mn")
    extractant : str
        추출제 종류 ("Cyanex 272", "D2EHPA")
    C_ext : float
        추출제 농도 (M, mol/L)
    temperature : float, optional
        온도 (°C). None이면 T_REF(25°C) 사용

    Returns
    -------
    float
        추출률 (0 ~ E_max, %)
    """
    params = _get_metal_params(metal, extractant, extractant_params)
    E_max = params["E_max"]

    pH50_eff = get_effective_pH50(metal, extractant, C_ext, temperature, extractant_params)
    k_eff = get_effective_k(metal, extractant, temperature, extractant_params)

    # 시그모이드 함수
    exponent = -k_eff * (pH - pH50_eff)
    # overflow 방지
    if exponent > 500:
        return 0.0
    elif exponent < -500:
        return E_max

    E = E_max / (1.0 + math.exp(exponent))
    return E


def distribution_coefficient(pH: float, metal: str, extractant: str,
                              C_ext: float,
                              temperature: float = None,
                              extractant_params: dict = None) -> float:
    """
    분배 계수 D를 계산합니다.

    시그모이드 모델 내부의 pseudo-D로, 물질수지 식
    C_aq_out = total_flow / (Q_aq + D * Q_org) 과 쌍으로 사용됩니다.

    D = E / (100 - E)

    Parameters
    ----------
    pH : float
    metal : str
    extractant : str
    C_ext : float
    temperature : float, optional

    Returns
    -------
    float
        분배 계수 D
    """
    E = extraction_efficiency(pH, metal, extractant, C_ext, temperature, extractant_params)

    if E >= 99.99:
        return 1e6
    elif E <= 0.01:
        return 1e-6

    E_frac = E / 100.0
    D = E_frac / (1.0 - E_frac)
    return D


def calc_loading_fraction(C_org: dict, extractant: str, C_ext: float,
                          metals: list = None,
                          extractant_params: dict = None) -> float:
    """
    유기상의 현재 로딩률을 계산합니다.

    로딩률 = Σ(n_ext(M) × C_M,org / MW_M) / C_ext

    Parameters
    ----------
    C_org : dict
        유기상 금속 농도 {metal: g/L}
    extractant : str
    C_ext : float
        추출제 농도 (M, mol/L)
    metals : list, optional

    Returns
    -------
    float
        로딩률 (0.0 ~ 1.0+)
    """
    if metals is None:
        metals = DEFAULT_METALS
    if C_ext <= 0:
        return 1.0

    total_ext_consumed = 0.0
    for metal in metals:
        C_org_m = C_org.get(metal, 0.0)
        if C_org_m <= 0:
            continue
        MW = MOLAR_MASS[metal]
        n_ext = _get_metal_params(metal, extractant, extractant_params).get("n_ext", 2)
        # 금속 몰 농도 (mol/L) × 추출제 화학양론
        total_ext_consumed += n_ext * (C_org_m / MW)

    loading = total_ext_consumed / C_ext
    return loading


def loading_damping_factor(loading_fraction: float) -> float:
    """
    로딩률에 따른 추출 효율 감쇠 계수를 반환합니다.

    감쇠 계수는 0 ~ 1 범위이며:
    - loading << MAX_LOADING: ~1.0 (감쇠 없음)
    - loading → MAX_LOADING: → 0.0 (추출 불가)

    시그모이드 감쇠 사용 (급격한 차단 방지):
    f = 1 / (1 + exp(12 * (L/L_max - 0.85)))

    Parameters
    ----------
    loading_fraction : float
        현재 로딩률 (0.0 ~ 1.0+)

    Returns
    -------
    float
        감쇠 계수 (0.0 ~ 1.0)
    """
    L_max = MAX_LOADING_FRACTION
    if loading_fraction <= 0:
        return 1.0
    if loading_fraction >= L_max:
        return 0.0

    # 시그모이드 감쇠: 로딩률이 85%를 넘으면 급격히 감소
    ratio = loading_fraction / L_max
    exponent = 12.0 * (ratio - 0.85)
    if exponent > 500:
        return 0.0
    elif exponent < -500:
        return 1.0
    factor = 1.0 / (1.0 + math.exp(exponent))
    return factor


# =========================================================================
# ESI 모델 함수 (Lu et al., 2024)
# =========================================================================

def calc_free_NaL(pH: float, C_ext: float, extractant: str,
                  C_org: dict = None, metals: list = None,
                  extractant_params: dict = None) -> float:
    """
    유기상의 자유 사포닌화 추출제 농도 [NaL]_free를 계산합니다.
    사포닌화 계산 시 ESI 모델 대신 화학양론비(n_ext)를 사용하여 근사합니다.

    [NaL]/[HL] = K'a × 10^pH
    C_ext_total = [HL] + [NaL] + consumed
    → [NaL]_free = (C_ext - consumed) × K'a×10^pH / (1 + K'a×10^pH)
    """
    if extractant != "Cyanex 272":
        return 0.0

    K_a = 7.674e-8

    # 금속에 의해 소비된 추출제 몰 농도 계산
    consumed = 0.0
    if C_org and metals:
        for metal in metals:
            n_ext = _get_metal_params(metal, extractant, extractant_params).get("n_ext", 2)
            C_M_org_gL = C_org.get(metal, 0.0)
            MW = MOLAR_MASS.get(metal, 1.0)
            C_M_org_mol = C_M_org_gL / MW
            consumed += n_ext * C_M_org_mol

    C_avail = max(0.0, C_ext - consumed)

    # [NaL]/[HL] = K'a × 10^pH
    ratio = K_a * (10 ** pH)
    NaL_free = C_avail * ratio / (1.0 + ratio)

    return max(0.0, NaL_free)


# =========================================================================
# 추출제 경쟁 모델 (Phase 3 — Vasilyev-inspired)
# =========================================================================

def compute_competitive_extractions(
    pH: float, extractant: str, C_ext: float,
    C_aq_in: dict, C_org_in: dict,
    Q_aq: float, Q_org: float,
    metals: list = None,
    temperature: float = None,
    Q_aq_eff: float = None,
    extractant_params: dict = None,
) -> dict:
    """
    추출제 경쟁을 반영한 분배 계수 및 추출 결과를 계산합니다.

    Vasilyev et al. (2019)의 "공유 추출제 풀" 개념을 시그모이드 모델과 결합:
      D_adjusted(M) = D_sigmoid(M) × ([HA]_free / C_ext)^n_ext(M)

    pH₅₀ 순서대로 금속을 순차 처리하며, 앞선 금속이 소비한 추출제를
    차감하여 뒤 금속의 D를 보정합니다.

    Parameters
    ----------
    pH : float          — 현재 pH
    extractant : str    — 추출제 종류
    C_ext : float       — 추출제 농도 (M)
    C_aq_in : dict      — 수계 입구 금속 농도 {metal: g/L}
    C_org_in : dict     — 유기계 입구 금속 농도 {metal: g/L}
    Q_aq : float        — 수계 입구 유량 (L/hr)
    Q_org : float       — 유기계 유량 (L/hr)
    metals : list       — 금속 리스트
    temperature : float — 온도 (°C)
    Q_aq_eff : float    — 수계 출구 유량 (L/hr) (NaOH 등 추가 시, None이면 Q_aq 사용)

    Returns
    -------
    dict:
        {
          "C_aq_out": {metal: g/L},
          "C_org_out": {metal: g/L},
          "extraction": {metal: %},
          "D_adjusted": {metal: float},
          "HA_free_fraction": float,  # 잔여 추출제 비율
        }
    """
    if metals is None:
        metals = DEFAULT_METALS
    if Q_aq_eff is None:
        Q_aq_eff = Q_aq
    params_store = _resolve_extractant_params(extractant_params)

    # 추출 우선순위: pH₅₀가 낮은 금속 먼저
    priority = EXTRACTION_PRIORITY.get(extractant, metals)
    ordered_metals = [m for m in priority if m in metals]
    # priority에 없는 금속을 뒤에 추가
    for m in metals:
        if m not in ordered_metals:
            ordered_metals.append(m)

    # 유기상에 이미 로딩된 금속 총 몰수 계산
    total_M_org_in_mol = 0.0
    for m in metals:
        C_org_m = C_org_in.get(m, 0.0)
        if C_org_m > 0:
            total_M_org_in_mol += (C_org_m / MOLAR_MASS[m])
            
    # [v2.0 고로딩 다핵 착물 방어 알고리즘]
    # 전체 금속 로딩률 기반으로 동적 n_eff_global 결정 (최대 60% 완화)
    # 유기상 금속 구성 기반 가중 평균 n_ext로 최대 적재량 계산
    if total_M_org_in_mol > 0:
        weighted_n_ext = 0.0
        for m in metals:
            C_org_m = C_org_in.get(m, 0.0)
            if C_org_m > 0:
                mol_m = C_org_m / MOLAR_MASS[m]
                n_ext_m = params_store[extractant][m].get("n_ext", 2)
                weighted_n_ext += n_ext_m * (mol_m / total_M_org_in_mol)
        weighted_n_ext = max(weighted_n_ext, 1.0)
    else:
        weighted_n_ext = 2.0
    max_theoretical_M_mol = C_ext / weighted_n_ext
    if max_theoretical_M_mol > 0:
        loading_ratio = total_M_org_in_mol / max_theoretical_M_mol
    else:
        loading_ratio = 0.0
        
    loading_ratio_clamped = min(1.0, loading_ratio)
    
    # n_eff_multiplier: 1.0 (낮은 로딩) -> 0.35 (극한 로딩)
    n_eff_multiplier = max(0.35, 1.0 - (loading_ratio_clamped ** 2))

    ext_consumed = 0.0
    for m in metals:
        C_org_m = C_org_in.get(m, 0.0)
        if C_org_m > 0:
            MW = MOLAR_MASS[m]
            n_ext = params_store[extractant][m].get("n_ext", 2)
            n_eff = n_ext * n_eff_multiplier
            ext_consumed += n_eff * (C_org_m / MW)

    HA_free = max(0.0, C_ext - ext_consumed)

    C_aq_out = {}
    C_org_out = {}
    extraction = {}
    D_adjusted = {}

    for metal in ordered_metals:
        MW = MOLAR_MASS[metal]
        n_ext = params_store[extractant][metal].get("n_ext", 2)
        n_eff = n_ext * n_eff_multiplier

        # 시그모이드 기반 D
        D_sigmoid = distribution_coefficient(
            pH, metal, extractant, C_ext, temperature=temperature,
            extractant_params=params_store
        )

        # 잔여 추출제 비율로 보정
        if C_ext > 0 and HA_free > 0:
            ratio = HA_free / C_ext
            ratio_eff = max(1e-4, ratio)
            competition_factor = ratio_eff ** n_eff
            competition_factor = max(0.0, min(1.0, competition_factor))
        else:
            competition_factor = 0.0

        D = D_sigmoid * competition_factor
        D_adjusted[metal] = D

        # 물질수지 계산
        C_aq_mol_in = C_aq_in.get(metal, 0.0) / MW
        C_org_mol_in = C_org_in.get(metal, 0.0) / MW
        total_metal_flow = C_aq_mol_in * Q_aq + C_org_mol_in * Q_org

        denominator = Q_aq_eff + D * Q_org
        if denominator <= 0:
            C_aq_mol_out = 0.0
        else:
            C_aq_mol_out = total_metal_flow / denominator
        C_org_mol_out = D * C_aq_mol_out

        C_aq_out[metal] = C_aq_mol_out * MW
        C_org_out[metal] = C_org_mol_out * MW

        # 추출률
        if C_aq_mol_in > 0:
            metal_extracted_mol = C_aq_mol_in - C_aq_mol_out * (Q_aq_eff / Q_aq)
            extraction[metal] = max(0.0, metal_extracted_mol / C_aq_mol_in * 100.0)
        else:
            extraction[metal] = 0.0

        # 이 금속이 유기상에서 순증가한 추출제 소비량 차감
        delta_org_mol = C_org_mol_out - C_org_in.get(metal, 0.0) / MW
        if delta_org_mol > 0:
            HA_free = max(0.0, HA_free - n_eff * delta_org_mol)

    return {
        "C_aq_out": C_aq_out,
        "C_org_out": C_org_out,
        "extraction": extraction,
        "D_adjusted": D_adjusted,
        "HA_free_fraction": HA_free / C_ext if C_ext > 0 else 0.0,
    }


def get_proton_release(metal: str, extractant: str,
                       extractant_params: dict = None) -> int:
    """
    금속 1몰 추출 시 방출되는 H⁺ 몰수를 반환합니다.

    Parameters
    ----------
    metal : str
    extractant : str

    Returns
    -------
    int
        H⁺ 방출 몰수 (n_H)
    """
    return _get_metal_params(metal, extractant, extractant_params)["n_H"]


def compute_all_extractions(pH: float, extractant: str, C_ext: float,
                             metals: list = None,
                             temperature: float = None,
                             extractant_params: dict = None) -> dict:
    """
    주어진 pH와 온도에서 모든 금속의 추출률을 한 번에 계산합니다.

    Parameters
    ----------
    pH : float
    extractant : str
    C_ext : float
    metals : list, optional
    temperature : float, optional

    Returns
    -------
    dict
        {금속: 추출률(%)} 딕셔너리
    """
    if metals is None:
        metals = DEFAULT_METALS

    results = {}
    for metal in metals:
        results[metal] = extraction_efficiency(
            pH, metal, extractant, C_ext, temperature, extractant_params
        )
    return results


def print_isotherm_table(extractant: str, C_ext: float,
                          pH_range: tuple = (1.0, 10.0), pH_step: float = 0.5,
                          metals: list = None, temperature: float = None,
                          extractant_params: dict = None):
    """
    pH 범위에 대한 추출률 테이블을 출력합니다.

    Parameters
    ----------
    extractant : str
    C_ext : float
    pH_range : tuple
    pH_step : float
    metals : list, optional
    temperature : float, optional
    """
    if metals is None:
        metals = DEFAULT_METALS

    T = temperature if temperature is not None else T_REF

    header = f"{'pH':>6}" + "".join(f"{m:>10}" for m in metals)
    print(f"\n추출률 Isotherm 테이블 ({extractant}, {C_ext}M, {T:.0f}°C)")
    print("=" * len(header))
    print(header)
    print("-" * len(header))

    pH = pH_range[0]
    while pH <= pH_range[1] + 0.001:
        row = f"{pH:6.1f}"
        for metal in metals:
            E = extraction_efficiency(
                pH, metal, extractant, C_ext, temperature, extractant_params
            )
            row += f"{E:10.2f}"
        print(row)
        pH += pH_step

    print("=" * len(header))
