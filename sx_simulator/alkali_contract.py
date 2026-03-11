"""
Alkali interface contracts
==========================
알칼리 적용 위치/제어 방식/목표 종류를 명시적으로 정리하는 경량 계약 계층입니다.
"""

from dataclasses import dataclass
from typing import Literal


AlkaliApplicationMode = Literal[
    "aqueous_direct",
    "fresh_organic_saponification",
]

AlkaliControlMode = Literal[
    "fixed_input",
    "solve_to_target",
]

AlkaliTargetKind = Literal[
    "none",
    "uniform_stage_pH",
    "stage_pH_profile",
    "raffinate_pH",
]

SaponificationModel = Literal[
    "physical_v2",
    "legacy_equivalent_target",
]


@dataclass(frozen=True)
class AlkaliContract:
    application: AlkaliApplicationMode
    control: AlkaliControlMode
    target_kind: AlkaliTargetKind
    saponification_model: SaponificationModel = "physical_v2"


def build_alkali_contract(
    *,
    naoh_mode: str,
    target_pH: float | None,
    target_pH_per_stage: list[float] | None,
    saponification_model: str = "physical_v2",
) -> AlkaliContract:
    """legacy kwargs를 명시적 알칼리 계약으로 정규화합니다."""
    if naoh_mode not in {"aqueous_direct", "saponification"}:
        raise ValueError(f"Unsupported naoh_mode: {naoh_mode}")
    if saponification_model not in {"physical_v2", "legacy_equivalent_target"}:
        raise ValueError(
            f"Unsupported saponification_model: {saponification_model}"
        )

    if naoh_mode == "aqueous_direct":
        if target_pH_per_stage is not None:
            return AlkaliContract(
                application="aqueous_direct",
                control="solve_to_target",
                target_kind="stage_pH_profile",
            )
        if target_pH is not None:
            return AlkaliContract(
                application="aqueous_direct",
                control="solve_to_target",
                target_kind="uniform_stage_pH",
            )
        return AlkaliContract(
            application="aqueous_direct",
            control="fixed_input",
            target_kind="none",
        )

    if target_pH_per_stage is not None:
        raise ValueError(
            "fresh organic saponification does not support stage_pH_profile targets."
        )

    if target_pH is not None:
        if saponification_model == "legacy_equivalent_target":
            return AlkaliContract(
                application="fresh_organic_saponification",
                control="solve_to_target",
                target_kind="uniform_stage_pH",
                saponification_model="legacy_equivalent_target",
            )
        return AlkaliContract(
            application="fresh_organic_saponification",
            control="solve_to_target",
            target_kind="raffinate_pH",
            saponification_model="physical_v2",
        )

    return AlkaliContract(
        application="fresh_organic_saponification",
        control="fixed_input",
        target_kind="none",
        saponification_model=saponification_model,
    )
