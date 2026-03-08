"""
Shared SX field datasets
========================
대시보드 프리셋과 검증 스크립트가 동일한 현장 데이터 정의를 공유하도록
묶어둔 모듈입니다.
"""

import copy

from .config import MOLAR_MASS

DEFAULT_ASSUMED_NAOH_CONCENTRATION_M = 5.0

VALIDATION_BASES = (
    "legacy_premixed_target_pH",
    "raw_feed_target_pH",
)

LEGACY_REGRESSION_MAE_THRESHOLDS = {
    "Li": 0.90,
    "Ni": 0.80,
    "Co": 0.08,
    "Mn": 0.01,
    "Ca": 0.01,
    "Mg": 0.01,
    "Zn": 0.01,
}

LEGACY_REGRESSION_DATASET_MAE_THRESHOLDS = {
    "Data1": 0.35,
    "Data2": 0.20,
    "Data3": 0.15,
    "Data4": 0.40,
    "Data5": 0.08,
    "Data6": 0.30,
}

LEGACY_REGRESSION_ABS_ERROR_THRESHOLDS = {
    ("Data1", "Ni"): 2.20,
    ("Data2", "Ni"): 1.20,
    ("Data3", "Li"): 0.90,
    ("Data4", "Li"): 2.30,
    ("Data4", "Co"): 0.15,
    ("Data5", "Ni"): 0.20,
    ("Data5", "Co"): 0.20,
    ("Data6", "Ni"): 1.00,
    ("Data6", "Co"): 0.10,
}


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


def dilute_feed_by_flow(feed: dict, inlet_flow: float, outlet_flow: float) -> dict:
    """보존 성분 농도를 유량 증가에 맞춰 희석 환산합니다."""
    if inlet_flow <= 0 or outlet_flow <= 0:
        return {metal: 0.0 for metal in feed}
    factor = inlet_flow / outlet_flow
    return {metal: conc * factor for metal, conc in feed.items()}


def prepare_verification_case(
    case: dict,
    basis: str = "legacy_premixed_target_pH",
    assumed_naoh_concentration_m: float | None = None,
) -> dict:
    """검증 데이터 한 건을 지정한 basis 기준의 엔진 입력으로 변환합니다."""
    if basis not in VALIDATION_BASES:
        raise ValueError(f"Unknown verification basis: {basis}")

    assumed_naoh_m = (
        assumed_naoh_concentration_m
        if assumed_naoh_concentration_m is not None
        else case.get(
            "assumed_naoh_concentration_m",
            DEFAULT_ASSUMED_NAOH_CONCENTRATION_M,
        )
    )

    prepared = {
        "name": case["name"],
        "basis": basis,
        "dataset_tag": case["name"].split()[0],
        "input_original": copy.deepcopy(case["input"]),
        "expected_output": copy.deepcopy(case["output"]),
        "listed_naoh_flow_l_hr": case["naoh_flow"],
        "assumed_naoh_concentration_m": assumed_naoh_m,
        "sim_kwargs": {},
    }

    if basis == "legacy_premixed_target_pH":
        q_aq = case["feed_flow"] + case["naoh_flow"]
        c_aq_feed = dilute_feed_by_flow(case["input"], case["feed_flow"], q_aq)
        sulfate = calc_sulfate_from_feed(c_aq_feed)
        prepared.update(
            {
                "input_model_basis": c_aq_feed,
                "dilution_factor": case["feed_flow"] / q_aq,
                "sulfate_m": sulfate,
                "sim_kwargs": {
                    "C_aq_feed": c_aq_feed,
                    "pH_feed": case["pH"],
                    "Q_aq": q_aq,
                    "Q_org": case["org_flow"],
                    "extractant": case["ext"],
                    "C_ext": case["C_ext"],
                    "n_stages": case["n_stages"],
                    "temperature": case["T"],
                    "C_sulfate": sulfate,
                    "target_pH": case["pH"],
                },
            }
        )
        return prepared

    c_aq_feed = copy.deepcopy(case["input"])
    sulfate = calc_sulfate_from_feed(c_aq_feed)
    prepared.update(
        {
            "input_model_basis": c_aq_feed,
            "dilution_factor": 1.0,
            "sulfate_m": sulfate,
            "sim_kwargs": {
                "C_aq_feed": c_aq_feed,
                "pH_feed": case["pH"],
                "Q_aq": case["feed_flow"],
                "Q_org": case["org_flow"],
                "extractant": case["ext"],
                "C_ext": case["C_ext"],
                "n_stages": case["n_stages"],
                "temperature": case["T"],
                "C_sulfate": sulfate,
                "target_pH": case["pH"],
                "C_NaOH": assumed_naoh_m,
            },
        }
    )
    return prepared
