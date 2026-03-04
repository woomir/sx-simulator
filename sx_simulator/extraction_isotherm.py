"""
Extraction Isotherm Module
==========================
pH에 따른 금속별 추출률(%) 및 분배 계수(D)를 계산하는 핵심 모듈.

시그모이드 함수로 pH-추출률 관계를 근사:
  E(pH) = E_max / (1 + exp(-k * (pH - pH50_eff)))

여기서 pH50_eff는 추출제 농도에 의해 보정:
  pH50_eff = pH50_ref - alpha * log10(C_ext / C_ref)
"""

import math
from .config import EXTRACTANT_PARAMS, MOLAR_MASS, DEFAULT_METALS


def get_effective_pH50(metal: str, extractant: str, C_ext: float) -> float:
    """
    추출제 농도에 따라 보정된 pH50 값을 반환합니다.

    Parameters
    ----------
    metal : str
        금속 종류 ("Li", "Ni", "Co", "Mn")
    extractant : str
        추출제 종류 ("Cyanex 272", "D2EHPA")
    C_ext : float
        추출제 농도 (M, mol/L)

    Returns
    -------
    float
        보정된 pH50 값
    """
    params = EXTRACTANT_PARAMS[extractant][metal]
    pH50_ref = params["pH50"]
    alpha = params["alpha"]
    C_ref = params["C_ref"]

    if C_ext <= 0:
        raise ValueError(f"추출제 농도는 양수여야 합니다: {C_ext}")

    # 추출제 농도 증가 → pH50 감소 (더 낮은 pH에서도 추출 가능)
    pH50_eff = pH50_ref - alpha * math.log10(C_ext / C_ref)
    return pH50_eff


def extraction_efficiency(pH: float, metal: str, extractant: str, C_ext: float) -> float:
    """
    주어진 pH에서 특정 금속의 추출률(%)을 계산합니다.

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

    Returns
    -------
    float
        추출률 (0 ~ E_max, %)
    """
    params = EXTRACTANT_PARAMS[extractant][metal]
    k = params["k"]
    E_max = params["E_max"]

    pH50_eff = get_effective_pH50(metal, extractant, C_ext)

    # 시그모이드 함수
    exponent = -k * (pH - pH50_eff)
    # overflow 방지
    if exponent > 500:
        return 0.0
    elif exponent < -500:
        return E_max

    E = E_max / (1.0 + math.exp(exponent))
    return E


def distribution_coefficient(pH: float, metal: str, extractant: str,
                              C_ext: float, ao_ratio: float = 1.0) -> float:
    """
    분배 계수 D를 계산합니다.

    D = E / (100 - E) * (1 / ao_ratio) 가 아니라,
    정의에 따라 D = C_org / C_aq 이므로:
      E = D * ao_ratio / (1 + D * ao_ratio) * 100
      → D = E / (100 - E) / ao_ratio (ao_ratio = Q_org / Q_aq 일 때)

    여기서는 순수한 열역학적 D (AO ratio 무관)를 반환합니다:
      D = E / (100 - E)  (단, E를 분율로 변환)

    Parameters
    ----------
    pH : float
    metal : str
    extractant : str
    C_ext : float
    ao_ratio : float
        A/O 부피비 (여기서는 사용하지 않음, 참고용)

    Returns
    -------
    float
        분배 계수 D
    """
    E = extraction_efficiency(pH, metal, extractant, C_ext)

    if E >= 99.99:
        return 1e6  # 사실상 무한대
    elif E <= 0.01:
        return 1e-6  # 사실상 0

    # E(%) → 분율 변환
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
                             metals: list = None) -> dict:
    """
    주어진 pH에서 모든 금속의 추출률을 한 번에 계산합니다.

    Parameters
    ----------
    pH : float
    extractant : str
    C_ext : float
    metals : list, optional
        계산할 금속 리스트 (기본: ["Li", "Ni", "Co", "Mn"])

    Returns
    -------
    dict
        {금속: 추출률(%)} 딕셔너리
    """
    if metals is None:
        metals = DEFAULT_METALS

    results = {}
    for metal in metals:
        results[metal] = extraction_efficiency(pH, metal, extractant, C_ext)
    return results


def print_isotherm_table(extractant: str, C_ext: float,
                          pH_range: tuple = (1.0, 10.0), pH_step: float = 0.5,
                          metals: list = None):
    """
    pH 범위에 대한 추출률 테이블을 출력합니다.

    Parameters
    ----------
    extractant : str
    C_ext : float
    pH_range : tuple
        (pH_min, pH_max)
    pH_step : float
    metals : list, optional
    """
    if metals is None:
        metals = DEFAULT_METALS

    # 헤더
    header = f"{'pH':>6}" + "".join(f"{m:>10}" for m in metals)
    print(f"\n추출률 Isotherm 테이블 ({extractant}, {C_ext}M)")
    print("=" * len(header))
    print(header)
    print("-" * len(header))

    pH = pH_range[0]
    while pH <= pH_range[1] + 0.001:
        row = f"{pH:6.1f}"
        for metal in metals:
            E = extraction_efficiency(pH, metal, extractant, C_ext)
            row += f"{E:10.2f}"
        print(row)
        pH += pH_step

    print("=" * len(header))
