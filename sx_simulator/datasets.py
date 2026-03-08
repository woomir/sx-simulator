"""
Shared SX field datasets
========================
대시보드 프리셋과 검증 스크립트가 동일한 현장 데이터 정의를 공유하도록
묶어둔 모듈입니다.
"""

import copy

from .config import MOLAR_MASS


FIELD_DATASETS = {
    "Data1 (CoSX-D9-data)": {
        "ext": "Cyanex 272",
        "C_ext": 0.6308,
        "pH": 6.10,
        "n_stages": 5,
        "T": 25.0,
        "feed_flow": 25.0,
        "naoh_flow": 11.0,
        "org_flow": 118.0,
        "input": {"Li": 5.283, "Ni": 33.894, "Co": 3.096, "Mn": 0.067, "Ca": 0.002, "Mg": 0.026, "Zn": 0.0},
        "output": {"Li": 3.322, "Ni": 7.166, "Co": 0.001, "Mn": 0.0, "Ca": 0.0, "Mg": 0.0, "Zn": 0.0},
    },
    "Data2 (CoSX-D5-data)": {
        "ext": "Cyanex 272",
        "C_ext": 0.6308,
        "pH": 5.89,
        "n_stages": 5,
        "T": 25.0,
        "feed_flow": 25.0,
        "naoh_flow": 9.4,
        "org_flow": 118.0,
        "input": {"Li": 5.260, "Ni": 33.876, "Co": 3.053, "Mn": 0.071, "Ca": 0.001, "Mg": 0.030, "Zn": 0.0},
        "output": {"Li": 3.478, "Ni": 8.951, "Co": 0.0, "Mn": 0.0, "Ca": 0.0, "Mg": 0.0, "Zn": 0.0},
    },
    "Data3 (CoSX-D2-data)": {
        "ext": "Cyanex 272",
        "C_ext": 0.6308,
        "pH": 7.00,
        "n_stages": 5,
        "T": 25.0,
        "feed_flow": 20.0,
        "naoh_flow": 9.4,
        "org_flow": 118.0,
        "input": {"Li": 5.852, "Ni": 28.476, "Co": 4.106, "Mn": 0.091, "Ca": 0.002, "Mg": 0.035, "Zn": 0.0},
        "output": {"Li": 3.151, "Ni": 0.044, "Co": 0.0, "Mn": 0.0, "Ca": 0.0, "Mg": 0.0, "Zn": 0.0},
    },
    "Data4 (IMSX-D5-Data)": {
        "ext": "D2EHPA",
        "C_ext": 0.6053,
        "pH": 4.00,
        "n_stages": 5,
        "T": 25.0,
        "feed_flow": 30.0,
        "naoh_flow": 4.2,
        "org_flow": 124.0,
        "input": {"Li": 9.767, "Ni": 28.889, "Co": 16.178, "Mn": 15.204, "Ca": 0.336, "Mg": 0.149, "Zn": 0.006},
        "output": {"Li": 6.335, "Ni": 18.186, "Co": 0.230, "Mn": 0.0, "Ca": 0.0, "Mg": 0.001, "Zn": 0.0},
    },
    "Data5 (IMSX-D9-Data)": {
        "ext": "D2EHPA",
        "C_ext": 0.6053,
        "pH": 4.20,
        "n_stages": 5,
        "T": 25.0,
        "feed_flow": 25.0,
        "naoh_flow": 3.2,
        "org_flow": 100.0,
        "input": {"Li": 0.013, "Ni": 71.000, "Co": 8.700, "Mn": 8.900, "Ca": 0.250, "Mg": 3.000, "Zn": 0.600},
        "output": {"Li": 5.400, "Ni": 16.000, "Co": 0.160, "Mn": 0.0, "Ca": 0.0, "Mg": 0.0, "Zn": 0.0},
    },
    "Data6 (IMSX-D?)": {
        "ext": "D2EHPA",
        "C_ext": 0.6053,
        "pH": 3.90,
        "n_stages": 5,
        "T": 25.0,
        "feed_flow": 30.0,
        "naoh_flow": 3.5,
        "org_flow": 100.0,
        "input": {"Li": 8.390, "Ni": 27.917, "Co": 15.400, "Mn": 12.883, "Ca": 0.413, "Mg": 0.129, "Zn": 0.012},
        "output": {"Li": 6.591, "Ni": 21.034, "Co": 0.848, "Mn": 0.0, "Ca": 0.003, "Mg": 0.0, "Zn": 0.0},
    },
}

KNOWN_PRESET_NOTES = {
    "Data5 (IMSX-D9-Data)": (
        "Data5의 Li 후액 값은 입력 Li 농도보다 높아 물질수지상 일관되지 않습니다. "
        "검증 리포트에서는 Li 지표를 제외해서 해석하는 편이 안전합니다."
    ),
}

KNOWN_INVALID_POINTS = {
    ("Data5", "Li"): (
        "원본 검증표의 Li 후액 값(5.400 g/L)이 입력 Li 농도(0.013 g/L)보다 커서 "
        "물질수지상 일관되지 않으므로 해당 지표는 검증 집계에서 제외합니다."
    ),
}


def _strip_expected_outputs(case: dict) -> dict:
    """대시보드 프리셋용으로 기대 출력값을 제외한 복사본을 반환합니다."""
    preset = copy.deepcopy(case)
    preset.pop("output", None)
    return preset


def _build_verification_name(base_name: str, n_stages: int) -> str:
    """기존 검증 스크립트와 동일한 표시 이름을 생성합니다."""
    return base_name[:-1] + f", {n_stages} Stages)"


DASHBOARD_PRESETS = {
    name: _strip_expected_outputs(case)
    for name, case in FIELD_DATASETS.items()
}

VERIFICATION_EXPERIMENTS = [
    {"name": _build_verification_name(name, case["n_stages"]), **copy.deepcopy(case)}
    for name, case in FIELD_DATASETS.items()
]


def calc_sulfate_from_feed(feed: dict) -> float:
    """대시보드와 검증 스크립트가 공유하는 총 황산염 농도 계산식."""
    return (
        feed["Li"] / (2 * MOLAR_MASS["Li"])
        + feed["Ni"] / MOLAR_MASS["Ni"]
        + feed["Co"] / MOLAR_MASS["Co"]
        + feed["Mn"] / MOLAR_MASS["Mn"]
        + feed["Ca"] / MOLAR_MASS["Ca"]
        + feed["Mg"] / MOLAR_MASS["Mg"]
        + feed["Zn"] / MOLAR_MASS["Zn"]
    )
