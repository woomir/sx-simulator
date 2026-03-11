"""
SX Simulator Configuration
==========================
금속별/추출제별 Isotherm 파라미터 데이터베이스.

추출 반응 메커니즘:
  M²⁺(aq) + n(HA)₂(org) → M(A₂)_n(org) + 2nH⁺(aq)

시그모이드 함수로 pH-추출률 관계를 근사:
  E(pH) = E_max / (1 + exp(-k * (pH - pH50)))

파라미터:
  - pH50     : 추출률이 50%에 도달하는 pH
  - k        : 시그모이드 기울기 (전이 급격도)
  - E_max    : 최대 추출률 (%)
  - n_H      : 금속 1몰 추출 시 방출되는 H⁺ 몰수
  - n_ext    : 금속 1몰 추출 시 소비되는 추출제 dimer (HA)₂ 몰수
  - alpha    : 추출제 농도 감도 계수 (pH50 보정용)
  - C_ref    : 기준 추출제 농도 (M)

온도 의존성 파라미터 (반-경험적 모델):
  - beta     : pH50 온도 계수 (°C⁻¹) — pH50(T) = pH50(T_ref) + beta*(T - T_ref)
  - gamma    : k 온도 계수 (°C⁻¹)    — k(T) = k(T_ref) · exp(gamma*(T - T_ref))
  양이온 추출은 발열 반응이 많아 beta < 0 (온도↑ → pH50↓ → 추출 용이)

추출제 로딩 한계:
  각 금속 추출 시 추출제 (HA)₂ dimer를 소비합니다:
    M²⁺ + n_ext(HA)₂ → MA_n(org) + 2nH⁺
  유기상에 로딩된 금속 총량이 추출제 용량에 근접하면 추출 효율이 감소합니다.

참고 문헌:
  - ALTA-2024-LBR-Paper-OLI (Miller et al.)
  - Wang et al. (2002, 2004, 2006)
  - Vasilyev et al. (2019) — stoichiometry 참고
"""

import copy
import math

# =============================================================================
# 금속 원자량 (g/mol) — 농도 변환용
# =============================================================================
MOLAR_MASS = {
    "Li": 6.941,
    "Ni": 58.693,
    "Co": 58.933,
    "Mn": 54.938,
    "Ca": 40.078,
    "Mg": 24.305,
    "Zn": 65.380,
    "Na": 22.990,
}

# =============================================================================
# 추출제별 금속 Isotherm 파라미터
# =============================================================================
# 각 추출제에 대해 금속별 파라미터를 딕셔너리로 정리합니다.
# pH50, k 값은 ALTA 2024 Figure 1 및 관련 문헌에서 추정한 초기값입니다.
# 실제 데이터로 보정(fitting) 작업이 필요합니다.
#
# 아래 `LITERATURE_BASE_PARAMS`는 문헌/기본 추정값을 보관합니다.
# 현장 적용용 조정치는 별도 override로 관리하고, 실제 기본 런타임 값은
# `FIELD_CALIBRATED_PARAMS`와 `EXTRACTANT_PARAMS`를 통해 제공합니다.

LITERATURE_BASE_PARAMS = {
    # -----------------------------------------------------------------
    # Cyanex 272 (Bis(2,4,4-trimethylpentyl)phosphinic acid)
    # -----------------------------------------------------------------
    # ALTA 2024 Figure 1(a),(b) 기반 추정:
    # - Co 추출은 Ni보다 낮은 pH에서 시작 (Co가 dimer, Ni가 trimer 복합체)
    # - Mn은 Co보다도 낮은 pH에서 추출 시작
    # - Li는 매우 높은 pH에서만 추출 (pH > 7)
    "Cyanex 272": {
        "Mn": {
            "pH50": 3.5,       # Mn은 가장 먼저 추출됨
            "k": 3.0,          # 비교적 급격한 전이
            "E_max": 99.5,     # 최대 추출률 (%)
            "n_H": 2,          # Mn²⁺ → 2H⁺ 방출
            "n_ext": 2,        # Mn²⁺ + 2(HA)₂ → MnA₂(HA)₂ + 2H⁺
            "alpha": 0.8,      # 추출제 농도 감도
            "C_ref": 0.5,      # 기준 농도 0.5M
            "beta": -0.04,     # 온도↑ → pH50↓ (발열 반응)
            "gamma": 0.005,    # 온도↑ → 전이 약간 급격해짐
        },
        "Co": {
            "pH50": 4.0,       # Mn 다음으로 추출
            "k": 3.5,          # 꽤 급격한 전이
            "E_max": 99.5,
            "n_H": 2,          # Co²⁺ → 2H⁺ 방출 (dimer 복합체)
            "n_ext": 2,        # Co²⁺ + 2(HA)₂ → CoA₂(HA)₂ + 2H⁺
            "alpha": 0.8,
            "C_ref": 0.5,
            "beta": -0.05,     # Co 추출은 온도에 민감
            "gamma": 0.008,
        },
        "Ni": {
            "pH50": 5.8,       # 문헌 기본값
            "k": 2.5,          # 상대적으로 완만한 전이
            "E_max": 99.0,
            "n_H": 2,          # Ni²⁺ → 2H⁺ 방출
            "n_ext": 2,        # Ni²⁺ + 2(HA)₂ → NiA₂(HA)₂ + 2H⁺
            "alpha": 0.7,
            "C_ref": 0.5,
            "beta": -0.03,     # Ni 추출은 온도에 상대적으로 덜 민감
            "gamma": 0.006,
        },
        "Li": {
            "pH50": 8.0,       # 매우 높은 pH에서만 추출
            "k": 2.0,          # 완만한 전이
            "E_max": 95.0,
            "n_H": 1,          # Li⁺ → 1H⁺ 방출
            "n_ext": 1,        # Li⁺ + (HA)₂ → LiA·HA + H⁺
            "alpha": 0.5,
            "C_ref": 0.5,
            "beta": -0.02,     # Li 추출은 온도 영향 적음
            "gamma": 0.003,
        },
        # Ca, Mg, Zn — Mohapatra et al. (2007), Pereira et al.
        "Zn": {
            "pH50": 1.85,      # Zn은 Mn보다도 낮은 pH에서 추출
            "k": 3.5,
            "E_max": 99.5,
            "n_H": 2,          # Zn²⁺ → 2H⁺ 방출
            "n_ext": 2,
            "alpha": 0.9,
            "C_ref": 0.5,
            "beta": -0.05,
            "gamma": 0.007,
        },
        "Ca": {
            "pH50": 5.5,       # Ni보다 약간 낮은 pH에서 추출
            "k": 2.0,
            "E_max": 90.0,
            "n_H": 2,          # Ca²⁺ → 2H⁺ 방출
            "n_ext": 2,
            "alpha": 0.6,
            "C_ref": 0.5,
            "beta": -0.02,
            "gamma": 0.003,
        },
        "Mg": {
            "pH50": 5.7,       # Ca와 유사, 약간 높은 pH
            "k": 2.0,
            "E_max": 85.0,
            "n_H": 2,          # Mg²⁺ → 2H⁺ 방출
            "n_ext": 2,
            "alpha": 0.6,
            "C_ref": 0.5,
            "beta": -0.02,
            "gamma": 0.003,
        },
    },

    # -----------------------------------------------------------------
    # D2EHPA (Di(2-ethylhexyl) phosphoric acid)
    # -----------------------------------------------------------------
    # ALTA 2024 Figure 1(c),(d),(e) 기반 추정:
    # - D2EHPA는 Cyanex 272보다 전반적으로 낮은 pH에서 추출 시작
    # - Co: pH < 5에서 급격히 증가, pH ≥ 6.5에서 완전 추출
    # - Li: pH < 5.5에서 추출 안 됨, pH > 5.5 이상에서 추출 시작
    # - TBP 첨가 시 pH isotherm이 높은 pH쪽으로 이동
    "D2EHPA": {
        "Mn": {
            "pH50": 2.5,       # D2EHPA에서 Mn은 매우 낮은 pH에서 추출
            "k": 3.5,
            "E_max": 99.5,
            "n_H": 2,
            "n_ext": 2,
            "alpha": 0.9,
            "C_ref": 0.64,     # 논문 기준 0.64M
            "beta": -0.03,
            "gamma": 0.004,
        },
        "Co": {
            "pH50": 3.5,       # pH < 5에서 급격 증가
            "k": 3.0,
            "E_max": 99.5,
            "n_H": 2,
            "n_ext": 2,
            "alpha": 0.9,
            "C_ref": 0.64,
            "beta": -0.04,
            "gamma": 0.006,
        },
        "Ni": {
            "pH50": 4.5,
            "k": 2.5,
            "E_max": 99.0,
            "n_H": 2,
            "n_ext": 2,
            "alpha": 0.8,
            "C_ref": 0.64,
            "beta": -0.025,
            "gamma": 0.005,
        },
        "Li": {
            "pH50": 6.5,       # pH > 5.5에서 추출 시작
            "k": 2.0,
            "E_max": 90.0,     # Li 최대 추출률이 상대적으로 낮음
            "n_H": 1,
            "n_ext": 1,
            "alpha": 0.6,
            "C_ref": 0.64,
            "beta": -0.02,
            "gamma": 0.003,
        },
        # Ca, Mg, Zn — Mohapatra et al. (2007), Pereira et al.
        "Zn": {
            "pH50": 1.5,       # D2EHPA에서 Zn은 가장 낮은 pH에서 추출
            "k": 4.0,
            "E_max": 99.5,
            "n_H": 2,
            "n_ext": 2,
            "alpha": 1.0,
            "C_ref": 0.64,
            "beta": -0.04,
            "gamma": 0.006,
        },
        "Ca": {
            "pH50": 3.5,       # pH 3.0~3.5에서 72% 추출
            "k": 2.0,
            "E_max": 90.0,
            "n_H": 2,
            "n_ext": 2,
            "alpha": 0.7,
            "C_ref": 0.64,
            "beta": -0.02,
            "gamma": 0.004,
        },
        "Mg": {
            "pH50": 4.5,       # 문헌 기본값
            "k": 2.0,
            "E_max": 85.0,
            "n_H": 2,
            "n_ext": 2,
            "alpha": 0.7,
            "C_ref": 0.64,
            "beta": -0.02,
            "gamma": 0.004,
        },
    },
}

# =============================================================================
# 모델 보정 프로필
# =============================================================================
# 문헌 기본값과 현장 보정값을 코드상에서 분리해 둡니다.
# - literature_default: 문헌/기본 추정값
# - field_calibrated: 현장 Isd 108 데이터에 맞춘 오버레이 적용값

SITE_PARAMETER_OVERRIDES = {
    ("Cyanex 272", "Ni"): {
        "pH50": 6.55,  # 현장 Data1~3 기준 저 pH Ni 과추출 완화
        "k": 2.0,      # Data1/2와 Data3 사이 전이 폭 균형화
    },
    ("Cyanex 272", "Li"): {
        "pH50": 8.10,  # Data1/2 과소추출과 Data3 과대추출 사이 절충
        "k": 1.8,
    },
    ("D2EHPA", "Li"): {
        "pH50": 6.3,  # field 기준 Li 잔여오차 완화
    },
    ("D2EHPA", "Co"): {
        "E_max": 99.0,  # trace Co 과추출 완화
    },
    ("D2EHPA", "Mg"): {
        "pH50": 3.2,  # 현장 데이터 기반 하향 보정
    },
}


def _build_parameter_profile(base_params: dict, overrides: dict | None = None) -> dict:
    """기본 파라미터에 선택적 override를 적용한 깊은 복사본을 반환합니다."""
    params = copy.deepcopy(base_params)
    if not overrides:
        return params

    for (extractant, metal), values in overrides.items():
        params[extractant][metal].update(values)
    return params


FIELD_CALIBRATED_PARAMS = _build_parameter_profile(
    LITERATURE_BASE_PARAMS,
    SITE_PARAMETER_OVERRIDES,
)

# 하위 호환용 기본 별칭: 엔진이 명시적 파라미터 세트를 받지 않으면 이 값을 사용합니다.
EXTRACTANT_PARAMS = copy.deepcopy(FIELD_CALIBRATED_PARAMS)


def get_parameter_profile(profile_name: str = "field_calibrated") -> dict:
    """
    시뮬레이터 파라미터 프로필의 깊은 복사본을 반환합니다.

    Parameters
    ----------
    profile_name : str
        - "field_calibrated": 현장 보정값 (기본값)
        - "literature_default": 문헌 기본값으로 일부 pH50 복원
    """
    if profile_name == "field_calibrated":
        return copy.deepcopy(FIELD_CALIBRATED_PARAMS)

    if profile_name == "literature_default":
        return copy.deepcopy(LITERATURE_BASE_PARAMS)

    raise ValueError(f"지원하지 않는 파라미터 프로필: {profile_name}")

# =============================================================================
# 추출 선택성 순서 (참고)
# =============================================================================
# Cyanex 272:  Zn > Mn > Co >> Ni > Ca ≈ Mg >> Li  (낮은 pH에서 높은 pH 순)
# D2EHPA:      Zn > Mn > Ca > Co > Mg > Ni >> Li

# =============================================================================
# 기본 시뮬레이션 설정
# =============================================================================
DEFAULT_METALS = ["Li", "Ni", "Co", "Mn", "Ca", "Mg", "Zn"]
DEFAULT_EXTRACTANT = "Cyanex 272"
DEFAULT_MODEL_TYPE = "sigmoid"
DEFAULT_TEMPERATURE = 25.0  # °C (기준 온도, T_ref)
T_REF = 25.0                # 온도 보정 기준점 (°C)
DEFAULT_STAGES = 4
DEFAULT_AO_RATIO = 1.0      # A/O = 1:1

# 수렴 관련 설정 (다단 역류 시뮬레이션용)
CONVERGENCE_TOLERANCE = 1e-6
MAX_ITERATIONS = 500

# 추출제 로딩 한계 설정
MAX_LOADING_FRACTION = 0.95   # 최대 로딩률 (추출제 용량의 95%까지 사용 가능)

# =============================================================================
# 추출제 경쟁 모델 설정 (Phase 3 — Mechanistic Light)
# =============================================================================
# Vasilyev et al. (2019)의 메커니즘을 단순화한 경쟁 추출 모델
# 핵심: 공유 추출제 풀에서 금속 간 경쟁 반영
# D_adjusted(M) = D_sigmoid(M) × ([HA]_free / C_ext)^n_ext(M)

DEFAULT_USE_COMPETITION = False   # 기본값: 경쟁 OFF (기존 시그모이드 모드)

# 추출 우선순위: pH₅₀가 낮은 금속이 먼저 추출제를 소비
EXTRACTION_PRIORITY = {
    "Cyanex 272": ["Zn", "Mn", "Co", "Ni", "Ca", "Mg", "Li"],
    "D2EHPA":     ["Zn", "Mn", "Ca", "Co", "Mg", "Ni", "Li"],
}

# =============================================================================
# 사포니피케이션 근사 설정
# =============================================================================
# 현장에서는 NaOH가 수계 직접 주입이 아니라 pre-saponified organic 준비에 쓰이므로,
# 유기상 자유 추출제의 Na-form 평형과 pH50 이동을 별도 근사합니다.
SAPONIFICATION_MAX_FRACTION = 0.95

# 자유 추출제(미로딩 site) 중 Na-form으로 남는 평형 분율을 계산할 때 사용하는
# apparent organic pKa. Cyanex 272 값은 기존 calc_free_NaL 근사와 정합되게 맞춥니다.
SAPONIFICATION_ORGANIC_PKA = {
    "Cyanex 272": 7.1149,
    "D2EHPA": 5.1000,
}

# 사포니피케이션으로 준비된 유기상이 자유산을 직접 중화하는 정도는
# 금속 교환보다 제한적이므로, direct neutralization은 부분 접근성으로 감쇠합니다.
SAPONIFICATION_DIRECT_NEUTRALIZATION_FACTOR = {
    "Cyanex 272": 0.25,
    "D2EHPA": 0.15,
}

# 현장 raw-feed fixed-saponification 데이터를 기준으로, sap condition을
# 등가 cascade target pH로 환산하는 간이 회귀식입니다.
SAPONIFICATION_EQUIVALENT_TARGET_PH_COEFFS = {
    "Cyanex 272": {
        "intercept": 1.90000000,
        "sap_acid_ratio": 0.0,
        "oa_ratio": 0.84745763,
        "max_pH": 7.2,
    },
    "D2EHPA": {
        "intercept": -7.60855021,
        "sap_acid_ratio": 0.14088917,
        "oa_ratio": 2.57995412,
        "max_pH": 6.5,
    },
}

# sap fraction이 높아질수록 같은 pH에서 추출이 쉬워지는 효과를 pH50 shift로 근사합니다.
# shift = coefficient * saponification_fraction
SAPONIFICATION_P_H50_SHIFT = {
    "Cyanex 272": {
        "default": 0.35,
        "Li": 0.40,
        "Ni": 0.65,
        "Co": 0.40,
        "Mn": 0.30,
        "Ca": 0.30,
        "Mg": 0.30,
        "Zn": 0.20,
    },
    "D2EHPA": {
        "default": 0.30,
        "Li": 0.70,
        "Ni": 0.40,
        "Co": 0.35,
        "Mn": 0.25,
        "Ca": 0.30,
        "Mg": 0.25,
        "Zn": 0.20,
    },
}

# =============================================================================
# 수계 종분화 평형 상수 (Phase 3 — Aqueous Speciation)
# =============================================================================
# M²⁺ + OH⁻ ⇌ MOH⁺   (K_MOH = [MOH⁺] / ([M²⁺]·[OH⁻]))
# M²⁺ + SO₄²⁻ ⇌ MSO₄⁰  (K_MSO4 = [MSO₄⁰] / ([M²⁺]·[SO₄²⁻]))
#
# 참고: Baes & Mesmer (1976), Smith & Martell (1976)

DEFAULT_USE_SPECIATION = False    # 기본값: 종분화 OFF

SPECIATION_CONSTANTS = {
    "Co": {"K_MOH": 10**(4.3),  "K_MSO4": 10**(2.36)},   # β₁/Kw: 10^(-9.7)/10^(-14) = 10^4.3
    "Ni": {"K_MOH": 10**(4.1),  "K_MSO4": 10**(2.29)},   # β₁/Kw: 10^(-9.9)/10^(-14) = 10^4.1
    "Mn": {"K_MOH": 10**(3.4),  "K_MSO4": 10**(2.26)},   # β₁/Kw: 10^(-10.6)/10^(-14) = 10^3.4
    "Li": {"K_MOH": 10**(0.4),  "K_MSO4": 10**(0.64)},   # β₁/Kw: 10^(-13.6)/10^(-14) = 10^0.4
    # Baes & Mesmer (1976), Smith & Martell (1976)
    "Zn": {"K_MOH": 10**(5.06), "K_MSO4": 10**(1.48)},   # β₁/Kw: 10^(-8.94)/10^(-14) = 10^5.06
    "Mg": {"K_MOH": 10**(2.30), "K_MSO4": 10**(2.13)},   # β₁/Kw: 10^(-11.70)/10^(-14) = 10^2.30
    "Ca": {"K_MOH": 10**(1.43), "K_MSO4": 10**(0.36)},   # β₁/Kw: 10^(-12.57)/10^(-14) = 10^1.43
}

# 자유 금속 기반 sulfate 보정은 현재 pH50 현장 보정을 보조하는 상대 보정으로만 사용합니다.
# 기준 황산염 농도(typical field sulfate)에서의 자유 금속 분율 대비 현재 분율의
# 비율을 사용하고, 먼저 기준 상대계수를 계산한 뒤 extractant/metal별로 선택 적용합니다.
FREE_METAL_CORRECTION_REFERENCE_SULFATE_M = 1.0
FREE_METAL_CORRECTION_MIN = 0.75
FREE_METAL_CORRECTION_MAX = 1.15

# 현재는 D2EHPA-Co 조합에만 강화된 sulfate 상대 보정을 적용합니다.
# raw_factor = relative_free_metal_factor
# correction = clamp(raw_factor ** power, min, max)
SULFATE_D_CORRECTION_RULES = {
    ("D2EHPA", "Co"): {"power": 4.0, "min": 0.2, "max": 1.0, "start_m": 0.0},
    ("D2EHPA", "Ni"): {"power": 2.5, "min": 0.35, "max": 1.0, "start_m": 1.55},
    ("D2EHPA", "Li"): {"power": 1.75, "min": 0.55, "max": 1.0, "start_m": 1.55},
}

# =============================================================================
# 현장 검증 범위 (Data1~6 기준)
# =============================================================================
VALIDATED_FIELD_WINDOW = {
    "pH": (3.9, 7.0),
    "C_ext_m": (0.6053, 0.6308),
    "temperature_c": 25.0,
    "temperature_margin_c": 5.0,
    "n_stages": 5,
    "total_metals_warning_gL": 80.0,
}

# Force Streamlit Cloud to reload this file
