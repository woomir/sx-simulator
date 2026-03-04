"""
Extraction Isotherm Module
==========================
pH에 따른 금속별 추출률(%) 및 분배 계수(D)를 계산하는 핵심 모듈.

시그모이드 함수로 pH-추출률 관계를 근사:
  E(pH) = E_max / (1 + exp(-k(T) * (pH - pH50_eff(T))))

여기서 pH50_eff는 추출제 농도와 온도에 의해 보정:
  pH50_eff(T) = pH50_ref - alpha * log10(C_ext / C_ref) + beta * (T - T_ref)
  k(T)        = k_ref * exp(gamma * (T - T_ref))
"""

import math
from .config import EXTRACTANT_PARAMS, MOLAR_MASS, DEFAULT_METALS, T_REF, DEFAULT_TEMPERATURE


def get_effective_pH50(metal: str, extractant: str, C_ext: float,
                       temperature: float = None) -> float:
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
    params = EXTRACTANT_PARAMS[extractant][metal]
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
                    temperature: float = None) -> float:
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
    params = EXTRACTANT_PARAMS[extractant][metal]
    k_ref = params["k"]
    gamma = params.get("gamma", 0.0)

    T = temperature if temperature is not None else T_REF
    k_eff = k_ref * math.exp(gamma * (T - T_REF))
    return k_eff


def extraction_efficiency(pH: float, metal: str, extractant: str, C_ext: float,
                          temperature: float = None) -> float:
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
    params = EXTRACTANT_PARAMS[extractant][metal]
    E_max = params["E_max"]

    pH50_eff = get_effective_pH50(metal, extractant, C_ext, temperature)
    k_eff = get_effective_k(metal, extractant, temperature)

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
                              C_ext: float, ao_ratio: float = 1.0,
                              temperature: float = None) -> float:
    """
    분배 계수 D를 계산합니다.

    D = C_org / C_aq = E / (100 - E)  (순수 열역학적 D)

    Parameters
    ----------
    pH : float
    metal : str
    extractant : str
    C_ext : float
    ao_ratio : float
    temperature : float, optional

    Returns
    -------
    float
        분배 계수 D
    """
    E = extraction_efficiency(pH, metal, extractant, C_ext, temperature)

    if E >= 99.99:
        return 1e6
    elif E <= 0.01:
        return 1e-6

    E_frac = E / 100.0
    D = E_frac / (1.0 - E_frac)
    return D


def get_proton_release(metal: str, extractant: str) -> int:
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
    return EXTRACTANT_PARAMS[extractant][metal]["n_H"]


def compute_all_extractions(pH: float, extractant: str, C_ext: float,
                             metals: list = None,
                             temperature: float = None) -> dict:
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
        results[metal] = extraction_efficiency(pH, metal, extractant, C_ext, temperature)
    return results


def print_isotherm_table(extractant: str, C_ext: float,
                          pH_range: tuple = (1.0, 10.0), pH_step: float = 0.5,
                          metals: list = None, temperature: float = None):
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
            E = extraction_efficiency(pH, metal, extractant, C_ext, temperature)
            row += f"{E:10.2f}"
        print(row)
        pH += pH_step

    print("=" * len(header))
