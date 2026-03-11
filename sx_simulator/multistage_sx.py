"""
Multistage SX Module
====================
역류(Counter-current) 다단 Mixer-Settler 시뮬레이션 엔진.

    Stage 1    Stage 2    Stage 3    ...    Stage N
    Aq →       Aq →       Aq →              Aq → (후액 Raffinate)
        ← Org      ← Org      ← Org  ...       ← Org (신선 용매)

목표 pH 모드: 각 stage에서 원하는 pH를 유지하도록 NaOH 자동 계산
"""

import copy
from .alkali_contract import build_alkali_contract
from .single_stage import (
    estimate_equivalent_saponification_target_pH,
    estimate_saponified_extractant_mol_flow,
    solve_single_stage,
)
from .config import (
    DEFAULT_METALS,
    CONVERGENCE_TOLERANCE,
    MAX_ITERATIONS,
    CONVERGENCE_TRACE_ORGANIC_ABS_TOL_G_L,
    CONVERGENCE_SAP_ABS_TOL_MOL_HR,
)


def _clone_org_profile(org_profile: list, metals: list) -> list:
    """유기상 stage profile을 금속 키 기준으로 안전하게 복사합니다."""
    return [
        {metal: stage.get(metal, 0.0) for metal in metals}
        for stage in org_profile
    ]


def _clone_scalar_profile(profile: list[float]) -> list[float]:
    """float profile의 얕은 복사본을 반환합니다."""
    return [float(value) for value in profile]


def _build_naoh_distribution(
    test_q_naoh: float,
    n_stages: int,
    naoh_weights: list = None,
    naoh_strategy: str = "uniform",
) -> list:
    """다단 고정 NaOH 분배량을 계산합니다."""
    if naoh_weights is not None and len(naoh_weights) == n_stages:
        weights = naoh_weights
    elif naoh_strategy == "front_loaded":
        weights = [0.5 ** i for i in range(n_stages)]
    else:
        weights = [1.0] * n_stages

    weight_sum = sum(weights)
    if weight_sum <= 0:
        return [0.0] * n_stages
    return [test_q_naoh * (weight / weight_sum) for weight in weights]


def _get_relaxation_alpha(
    iteration: int,
    last_max_diff: float = None,
    target_mode: bool = False,
    relaxation_scale: float = 1.0,
) -> float:
    """
    외부 fixed-point 반복의 under-relaxation 계수를 반환합니다.

    target pH 모드는 강한 비선형성을 보이므로 기본적으로 더 보수적인 damping을 사용합니다.
    """
    if target_mode:
        if iteration < 40:
            alpha = 0.10
        elif iteration < 120:
            alpha = 0.18
        elif iteration < 250:
            alpha = 0.28
        else:
            alpha = 0.40
    else:
        if iteration < 30:
            alpha = 0.18
        elif iteration < 100:
            alpha = 0.30
        elif iteration < 250:
            alpha = 0.45
        else:
            alpha = 0.60

    if last_max_diff is not None:
        if last_max_diff > 10.0:
            alpha = min(alpha, 0.05)
        elif last_max_diff > 3.0:
            alpha = min(alpha, 0.08)
        elif last_max_diff > 1.0:
            alpha = min(alpha, 0.12)
        elif last_max_diff > 0.3:
            alpha = min(alpha, 0.18)
        elif last_max_diff < 0.03:
            alpha = min(0.75, alpha * 1.15)

    alpha *= relaxation_scale
    return max(0.03, min(0.85, alpha))


def _compute_max_relative_diff(
    current_org_out: list,
    previous_org_out: list,
    metals: list,
    current_sap_out: list[float] | None = None,
    previous_sap_out: list[float] | None = None,
) -> float:
    """유기상 프로파일 변화의 최대 상대차를 계산합니다."""
    max_diff = 0.0
    for stage_idx in range(len(current_org_out)):
        for metal in metals:
            current_value = current_org_out[stage_idx][metal]
            previous_value = previous_org_out[stage_idx].get(metal, 0.0)
            abs_diff = abs(current_value - previous_value)
            if abs_diff <= CONVERGENCE_TRACE_ORGANIC_ABS_TOL_G_L:
                continue
            ref_val = max(abs(previous_value), abs(current_value), 1e-10)
            diff = abs_diff / ref_val
            if diff > max_diff:
                max_diff = diff
        if current_sap_out is not None and previous_sap_out is not None:
            current_value = current_sap_out[stage_idx]
            previous_value = previous_sap_out[stage_idx]
            abs_diff = abs(current_value - previous_value)
            if abs_diff <= CONVERGENCE_SAP_ABS_TOL_MOL_HR:
                continue
            ref_val = max(abs(previous_value), abs(current_value), 1e-10)
            diff = abs_diff / ref_val
            if diff > max_diff:
                max_diff = diff
    return max_diff


def _interpolate_stage_targets(stage_targets: list, pH_feed: float, level: float) -> list:
    """
    target pH fallback continuation용 stage target 보간값을 생성합니다.
    level=1.0이면 원래 target과 동일합니다.
    """
    interpolated = []
    for target in stage_targets:
        if target is None:
            interpolated.append(None)
        else:
            interpolated.append(pH_feed + level * (target - pH_feed))
    return interpolated


def _resolve_stage_saponified_input_mol_hr(
    stage_idx: int,
    n_stages: int,
    sap_profile_seed: list[float] | None,
    fresh_sap_input_mol_hr: float,
) -> float:
    """현재 stage 유기상에 유입되는 sap inventory를 seed profile에서 읽어옵니다."""
    if stage_idx >= n_stages - 1:
        return max(0.0, fresh_sap_input_mol_hr)
    if sap_profile_seed is None or len(sap_profile_seed) <= stage_idx + 1:
        return 0.0
    return max(0.0, float(sap_profile_seed[stage_idx + 1]))


def solve_multistage_countercurrent(
    C_aq_feed: dict,
    pH_feed: float,
    Q_aq: float,
    Q_org: float,
    extractant: str,
    C_ext: float,
    n_stages: int,
    target_pH: float = None,
    target_pH_per_stage: list = None,
    C_NaOH: float = 0.0,
    Q_NaOH: float = 0.0,
    naoh_mode: str = "aqueous_direct",
    saponification_model: str = "physical_v2",
    saponification_fraction: float | None = None,
    naoh_strategy: str = "uniform",
    naoh_weights: list = None,
    C_org_fresh: dict = None,
    metals: list = None,
    tolerance: float = None,
    max_iter: int = None,
    temperature: float = None,
    C_sulfate: float = 0.0,
    use_competition: bool = False,
    use_speciation: bool = False,
    extractant_params: dict = None,
) -> dict:
    """
    역류 다단 Mixer-Settler 시뮬레이션.

    pH 제어 모드:
    1. target_pH (float): 모든 stage에 동일한 목표 pH 적용
    2. target_pH_per_stage (list): stage별 개별 목표 pH 지정
    3. 둘 다 None이면 고정 NaOH 모드
    """
    if metals is None:
        metals = DEFAULT_METALS
    if tolerance is None:
        tolerance = CONVERGENCE_TOLERANCE
    if max_iter is None:
        max_iter = MAX_ITERATIONS
    if C_org_fresh is None:
        C_org_fresh = {m: 0.0 for m in metals}

    alkali_contract = build_alkali_contract(
        naoh_mode=naoh_mode,
        target_pH=target_pH,
        target_pH_per_stage=target_pH_per_stage,
        saponification_model=saponification_model,
    )
    physical_saponification_mode = (
        alkali_contract.application == "fresh_organic_saponification"
        and alkali_contract.saponification_model == "physical_v2"
    )
    legacy_equivalent_saponification_mode = (
        alkali_contract.application == "fresh_organic_saponification"
        and alkali_contract.saponification_model == "legacy_equivalent_target"
    )

    # 목표 pH 결정
    if alkali_contract.target_kind == "stage_pH_profile":
        stage_target_pHs = target_pH_per_stage
    elif alkali_contract.target_kind == "uniform_stage_pH":
        # 단일 목표 pH → 모든 stage에 동일 pH 적용
        stage_target_pHs = [target_pH] * n_stages
    else:
        stage_target_pHs = [None] * n_stages

    equivalent_saponification_target_mode = False
    equivalent_saponification_target_pH = None
    if (
        alkali_contract.application == "fresh_organic_saponification"
        and alkali_contract.control == "fixed_input"
        and legacy_equivalent_saponification_mode
    ):
        equivalent_saponification_target_pH = estimate_equivalent_saponification_target_pH(
            C_aq_in=C_aq_feed,
            pH_in=pH_feed,
            Q_aq=Q_aq,
            Q_org=Q_org,
            extractant=extractant,
            C_ext=C_ext,
            C_NaOH=C_NaOH,
            Q_NaOH=Q_NaOH,
            C_sulfate=C_sulfate,
            use_speciation=use_speciation,
            saponification_fraction=saponification_fraction,
        )
        stage_target_pHs = [equivalent_saponification_target_pH] * n_stages
        equivalent_saponification_target_mode = True

    def _run_exact_stage_pass(
        org_profile_seed: list,
        Q_NaOH_dist: list,
        stage_targets: list,
        sap_profile_seed: list | None = None,
        fresh_sap_input_mol_hr: float = 0.0,
    ) -> list:
        """주어진 유기상 seed를 사용해 relaxation 없는 stage 결과를 계산합니다."""
        stage_results = []
        C_aq_current = copy.deepcopy(C_aq_feed)
        pH_current = pH_feed
        Q_aq_current = Q_aq
        C_sulfate_current = C_sulfate

        for stage_idx in range(n_stages):
            C_org_in = (
                org_profile_seed[stage_idx + 1]
                if stage_idx < n_stages - 1
                else copy.deepcopy(C_org_fresh)
            )
            stage_target = stage_targets[stage_idx]
            if stage_target is not None:
                target_naoh_mode = (
                    "aqueous_direct"
                    if equivalent_saponification_target_mode
                    else naoh_mode
                )
                result = solve_single_stage(
                    C_aq_in=C_aq_current, C_org_in=C_org_in,
                    pH_in=pH_current, Q_aq=Q_aq_current, Q_org=Q_org,
                    extractant=extractant, C_ext=C_ext,
                    target_pH=stage_target, metals=metals,
                    C_NaOH=0.0 if equivalent_saponification_target_mode else C_NaOH,
                    Q_NaOH=0.0,
                    naoh_mode=target_naoh_mode,
                    saponification_fraction=(
                        None
                        if equivalent_saponification_target_mode
                        else saponification_fraction
                    ),
                    temperature=temperature, C_sulfate=C_sulfate_current,
                    use_competition=use_competition, use_speciation=use_speciation,
                    extractant_params=extractant_params,
                )
            else:
                q_naoh = Q_NaOH_dist[stage_idx]
                if physical_saponification_mode:
                    stage_sap_in_mol_hr = _resolve_stage_saponified_input_mol_hr(
                        stage_idx,
                        n_stages,
                        sap_profile_seed,
                        fresh_sap_input_mol_hr,
                    )
                    result = solve_single_stage(
                        C_aq_in=C_aq_current, C_org_in=C_org_in,
                        pH_in=pH_current, Q_aq=Q_aq_current, Q_org=Q_org,
                        extractant=extractant, C_ext=C_ext,
                        C_NaOH=0.0, Q_NaOH=0.0, naoh_mode=naoh_mode,
                        saponified_extractant_mol_hr=stage_sap_in_mol_hr,
                        saponification_fraction=saponification_fraction,
                        metals=metals,
                        temperature=temperature, C_sulfate=C_sulfate_current,
                        use_competition=use_competition, use_speciation=use_speciation,
                        extractant_params=extractant_params,
                    )
                else:
                    result = solve_single_stage(
                        C_aq_in=C_aq_current, C_org_in=C_org_in,
                        pH_in=pH_current, Q_aq=Q_aq_current, Q_org=Q_org,
                        extractant=extractant, C_ext=C_ext,
                        C_NaOH=C_NaOH, Q_NaOH=q_naoh, naoh_mode=naoh_mode,
                        saponification_fraction=saponification_fraction,
                        metals=metals,
                        temperature=temperature, C_sulfate=C_sulfate_current,
                        use_competition=use_competition, use_speciation=use_speciation,
                        extractant_params=extractant_params,
                    )
                if alkali_contract.application == "aqueous_direct":
                    Q_aq_current += q_naoh

            stage_results.append(result)
            C_aq_current = copy.deepcopy(result["C_aq_out"])
            pH_current = result["pH_out"]
            C_sulfate_current = result.get("C_sulfate_out_M", C_sulfate_current)
            if stage_target is not None:
                Q_aq_current = result.get("Q_aq_out_L_hr", Q_aq_current)

        return stage_results

    def _build_result_payload(
        stage_results: list,
        converged: bool,
        iteration_count: int,
        best_residual: float,
        solver_strategy: str,
        has_target_mode: bool,
        total_q_naoh_input: float,
        fresh_sap_input_mol_hr: float,
    ) -> dict:
        raffinate = stage_results[-1]["C_aq_out"]
        pH_profile = [result["pH_out"] for result in stage_results]
        overall_extraction = {}
        for metal in metals:
            C_feed = C_aq_feed.get(metal, 0.0)
            C_raff = raffinate.get(metal, 0.0)
            overall_extraction[metal] = (1.0 - C_raff / C_feed) * 100.0 if C_feed > 0 else 0.0

        naoh_profile = [result.get("NaOH_consumed_mol_hr", 0.0) for result in stage_results]
        naoh_flow_profile = [
            result.get("Q_NaOH_estimated_L_hr", 0.0)
            for result in stage_results
        ]
        if naoh_mode == "saponification" and not has_target_mode and C_NaOH > 0.0:
            naoh_flow_profile = [value / C_NaOH for value in naoh_profile]

        actual_naoh_consumed_mol_hr = sum(naoh_profile)
        total_naoh_mol_hr = (
            fresh_sap_input_mol_hr
            if naoh_mode == "saponification" and (not has_target_mode or equivalent_saponification_target_mode)
            else actual_naoh_consumed_mol_hr
        )

        return {
            "stages": stage_results,
            "raffinate": raffinate,
            "loaded_organic": stage_results[0]["C_org_out"],
            "pH_profile": pH_profile,
            "NaOH_profile": naoh_profile,
            "NaOH_flow_profile": naoh_flow_profile,
            "overall_extraction": overall_extraction,
            "converged": converged,
            "iterations": iteration_count,
            "best_residual": best_residual,
            "solver_strategy": solver_strategy,
            "total_NaOH_mol_hr": total_naoh_mol_hr,
            "total_NaOH_consumed_mol_hr": actual_naoh_consumed_mol_hr,
            "naoh_mode": naoh_mode,
            "alkali_application": alkali_contract.application,
            "alkali_control_mode": alkali_contract.control,
            "alkali_target_kind": alkali_contract.target_kind,
            "saponification_model": alkali_contract.saponification_model,
            "legacy_fallback_used": equivalent_saponification_target_mode,
            "fresh_saponification_capacity_mol_hr": fresh_sap_input_mol_hr,
            "equivalent_saponification_target_pH": equivalent_saponification_target_pH,
        }

    def _run_with_q_naoh(
        test_q_naoh: float,
        stage_targets: list = None,
        initial_org_out: list = None,
        initial_sap_out: list | None = None,
        relaxation_scale: float = 1.0,
        solver_strategy: str = "direct",
    ) -> dict:
        if stage_targets is None:
            stage_targets = stage_target_pHs

        has_target_mode = any(target is not None for target in stage_targets)
        fresh_sap_input_mol_hr = (
            estimate_saponified_extractant_mol_flow(
                C_ext, Q_org, C_NaOH, test_q_naoh, saponification_fraction
            )
            if alkali_contract.application == "fresh_organic_saponification"
            and (equivalent_saponification_target_mode or not has_target_mode)
            else 0.0
        )
        Q_NaOH_dist = (
            [test_q_naoh] * n_stages
            if naoh_mode == "saponification"
            else _build_naoh_distribution(
                test_q_naoh, n_stages, naoh_weights, naoh_strategy
            )
        )
        org_out = (
            _clone_org_profile(initial_org_out, metals)
            if initial_org_out is not None
            else [{metal: 0.0 for metal in metals} for _ in range(n_stages)]
        )
        sap_out = (
            _clone_scalar_profile(initial_sap_out)
            if initial_sap_out is not None
            else [0.0] * n_stages
        )
        converged = False
        last_max_diff = None
        best_residual = float("inf")
        best_org_out = _clone_org_profile(org_out, metals)
        best_sap_out = _clone_scalar_profile(sap_out)
        best_iteration = 0

        for iteration in range(max_iter):
            org_out_prev = copy.deepcopy(org_out)
            sap_out_prev = _clone_scalar_profile(sap_out)
            stage_results = []
            C_aq_current = copy.deepcopy(C_aq_feed)
            pH_current = pH_feed
            Q_aq_current = Q_aq
            C_sulfate_current = C_sulfate
            target_stage_converged = True
            alpha = _get_relaxation_alpha(
                iteration,
                last_max_diff=last_max_diff,
                target_mode=has_target_mode,
                relaxation_scale=relaxation_scale,
            )

            for s in range(n_stages):
                C_org_in = org_out[s + 1] if s < n_stages - 1 else copy.deepcopy(C_org_fresh)
                stage_num = s + 1
                t_pH = stage_targets[s]

                if t_pH is not None:
                    target_naoh_mode = (
                        "aqueous_direct"
                        if equivalent_saponification_target_mode
                        else naoh_mode
                    )
                    result = solve_single_stage(
                        C_aq_in=C_aq_current, C_org_in=C_org_in,
                        pH_in=pH_current, Q_aq=Q_aq_current, Q_org=Q_org,
                        extractant=extractant, C_ext=C_ext,
                        target_pH=t_pH, metals=metals,
                        C_NaOH=0.0 if equivalent_saponification_target_mode else C_NaOH,
                        Q_NaOH=0.0,
                        naoh_mode=target_naoh_mode,
                        saponification_fraction=(
                            None
                            if equivalent_saponification_target_mode
                            else saponification_fraction
                        ),
                        temperature=temperature,
                        C_sulfate=C_sulfate_current,
                        use_competition=use_competition,
                        use_speciation=use_speciation,
                        extractant_params=extractant_params,
                    )
                    target_stage_converged = (
                        target_stage_converged
                        and result.get("target_pH_dilution_converged", True)
                    )
                else:
                    q_naoh = Q_NaOH_dist[stage_num - 1]
                    if physical_saponification_mode:
                        stage_sap_in_mol_hr = _resolve_stage_saponified_input_mol_hr(
                            s,
                            n_stages,
                            sap_out_prev,
                            fresh_sap_input_mol_hr,
                        )
                        result = solve_single_stage(
                            C_aq_in=C_aq_current, C_org_in=C_org_in,
                            pH_in=pH_current, Q_aq=Q_aq_current, Q_org=Q_org,
                            extractant=extractant, C_ext=C_ext,
                            C_NaOH=0.0, Q_NaOH=0.0, naoh_mode=naoh_mode,
                            saponified_extractant_mol_hr=stage_sap_in_mol_hr,
                            saponification_fraction=saponification_fraction,
                            metals=metals,
                            temperature=temperature,
                            C_sulfate=C_sulfate_current,
                            use_competition=use_competition,
                            use_speciation=use_speciation,
                            extractant_params=extractant_params,
                        )
                    else:
                        result = solve_single_stage(
                            C_aq_in=C_aq_current, C_org_in=C_org_in,
                            pH_in=pH_current, Q_aq=Q_aq_current, Q_org=Q_org,
                            extractant=extractant, C_ext=C_ext,
                            C_NaOH=C_NaOH, Q_NaOH=q_naoh, naoh_mode=naoh_mode,
                            saponification_fraction=saponification_fraction,
                            metals=metals,
                            temperature=temperature,
                            C_sulfate=C_sulfate_current,
                            use_competition=use_competition,
                            use_speciation=use_speciation,
                            extractant_params=extractant_params,
                        )
                    if alkali_contract.application == "aqueous_direct":
                        Q_aq_current += q_naoh  # NaOH 유량 누적

                stage_results.append(result)

                org_calc = result["C_org_out"]
                org_relaxed = {}
                for m in metals:
                    relaxed_value = alpha * org_calc[m] + (1.0 - alpha) * org_out_prev[s].get(m, 0.0)
                    org_relaxed[m] = max(0.0, relaxed_value)
                org_out[s] = org_relaxed
                if alkali_contract.application == "fresh_organic_saponification" and t_pH is None:
                    sap_calc = result.get("saponified_extractant_out_mol_hr", 0.0)
                    sap_relaxed = alpha * sap_calc + (1.0 - alpha) * sap_out_prev[s]
                    sap_out[s] = max(0.0, sap_relaxed)

                C_aq_current = copy.deepcopy(result["C_aq_out"])
                pH_current = result["pH_out"]
                C_sulfate_current = result.get("C_sulfate_out_M", C_sulfate_current)
                if t_pH is not None:
                    Q_aq_current = result.get("Q_aq_out_L_hr", Q_aq_current)

            max_diff = _compute_max_relative_diff(
                org_out,
                org_out_prev,
                metals,
                current_sap_out=sap_out if alkali_contract.application == "fresh_organic_saponification" and not has_target_mode else None,
                previous_sap_out=sap_out_prev if alkali_contract.application == "fresh_organic_saponification" and not has_target_mode else None,
            )
            last_max_diff = max_diff

            if max_diff < best_residual:
                best_residual = max_diff
                best_org_out = _clone_org_profile(org_out, metals)
                best_sap_out = _clone_scalar_profile(sap_out)
                best_iteration = iteration + 1

            if max_diff < tolerance and target_stage_converged:
                converged = True
                break

        final_org_seed = org_out if converged else best_org_out
        final_sap_seed = sap_out if converged else best_sap_out
        stage_results = _run_exact_stage_pass(
            final_org_seed,
            Q_NaOH_dist,
            stage_targets,
            sap_profile_seed=final_sap_seed,
            fresh_sap_input_mol_hr=fresh_sap_input_mol_hr,
        )
        max_refine_passes = max(6, 2 * n_stages)
        for _ in range(max_refine_passes):
            refined_org_seed = [result["C_org_out"] for result in stage_results]
            refined_sap_seed = [
                result.get("saponified_extractant_out_mol_hr", 0.0)
                for result in stage_results
            ]
            refined_results = _run_exact_stage_pass(
                refined_org_seed,
                Q_NaOH_dist,
                stage_targets,
                sap_profile_seed=refined_sap_seed,
                fresh_sap_input_mol_hr=fresh_sap_input_mol_hr,
            )
            refined_diff = _compute_max_relative_diff(
                [result["C_org_out"] for result in refined_results],
                [result["C_org_out"] for result in stage_results],
                metals,
                current_sap_out=(
                    [result.get("saponified_extractant_out_mol_hr", 0.0) for result in refined_results]
                    if alkali_contract.application == "fresh_organic_saponification" and not has_target_mode
                    else None
                ),
                previous_sap_out=(
                    [result.get("saponified_extractant_out_mol_hr", 0.0) for result in stage_results]
                    if alkali_contract.application == "fresh_organic_saponification" and not has_target_mode
                    else None
                ),
            )
            stage_results = refined_results
            if refined_diff < max(tolerance, 1e-6):
                break
        return _build_result_payload(
            stage_results=stage_results,
            converged=converged,
            iteration_count=(iteration + 1) if converged else best_iteration,
            best_residual=best_residual,
            solver_strategy=solver_strategy,
            has_target_mode=has_target_mode,
            total_q_naoh_input=test_q_naoh,
            fresh_sap_input_mol_hr=fresh_sap_input_mol_hr,
        )

    if (
        physical_saponification_mode
        and alkali_contract.control == "solve_to_target"
        and alkali_contract.target_kind == "raffinate_pH"
    ):
        if C_NaOH <= 0.0:
            raise ValueError(
                "fresh organic saponification solve_to_target requires C_NaOH > 0."
            )

        q_capacity_max = (0.95 * C_ext * Q_org) / C_NaOH if C_NaOH > 0.0 else 0.0
        q_low = 0.0
        q_high = max(Q_NaOH, q_capacity_max, 1e-9)
        low_result = _run_with_q_naoh(
            q_low,
            solver_strategy="fresh_sap_raffinate_target_low",
        )
        low_residual = low_result["pH_profile"][-1] - target_pH
        if low_residual >= 0.0:
            low_result["objective_converged"] = abs(low_residual) < 1e-3
            low_result["objective_residual_pH"] = low_residual
            low_result["solver_strategy"] = "fresh_sap_raffinate_target"
            return low_result

        high_result = _run_with_q_naoh(
            q_high,
            solver_strategy="fresh_sap_raffinate_target_high",
            initial_org_out=[stage["C_org_out"] for stage in low_result["stages"]],
            initial_sap_out=[
                stage.get("saponified_extractant_out_mol_hr", 0.0)
                for stage in low_result["stages"]
            ],
            relaxation_scale=0.85,
        )
        high_residual = high_result["pH_profile"][-1] - target_pH
        best_result = (
            low_result if abs(low_residual) <= abs(high_residual) else high_result
        )
        best_residual = min(abs(low_residual), abs(high_residual))

        if high_residual <= 0.0:
            high_result["objective_converged"] = False
            high_result["objective_residual_pH"] = high_residual
            high_result["solver_strategy"] = "fresh_sap_raffinate_target_unreachable"
            return high_result if abs(high_residual) <= abs(low_residual) else low_result

        seed_org = [stage["C_org_out"] for stage in high_result["stages"]]
        seed_sap = [
            stage.get("saponified_extractant_out_mol_hr", 0.0)
            for stage in high_result["stages"]
        ]
        objective_converged = False
        for iteration in range(24):
            q_mid = 0.5 * (q_low + q_high)
            mid_result = _run_with_q_naoh(
                q_mid,
                solver_strategy="fresh_sap_raffinate_target",
                initial_org_out=seed_org,
                initial_sap_out=seed_sap,
                relaxation_scale=0.80,
            )
            mid_residual = mid_result["pH_profile"][-1] - target_pH
            if abs(mid_residual) < best_residual:
                best_result = mid_result
                best_residual = abs(mid_residual)
            seed_org = [stage["C_org_out"] for stage in mid_result["stages"]]
            seed_sap = [
                stage.get("saponified_extractant_out_mol_hr", 0.0)
                for stage in mid_result["stages"]
            ]
            if abs(mid_residual) < 1e-3 or (q_high - q_low) < 1e-4:
                best_result = mid_result
                objective_converged = abs(mid_residual) < 1e-3
                break
            if mid_residual < 0.0:
                q_low = q_mid
            else:
                q_high = q_mid

        final_residual = best_result["pH_profile"][-1] - target_pH
        best_result["objective_converged"] = objective_converged or abs(final_residual) < 1e-3
        best_result["objective_residual_pH"] = final_residual
        best_result["solver_strategy"] = "fresh_sap_raffinate_target"
        return best_result

    initial_org_seed = None
    if (
        physical_saponification_mode
        and alkali_contract.control == "fixed_input"
        and Q_NaOH > 0.0
    ):
        legacy_seed_result = solve_multistage_countercurrent(
            C_aq_feed=C_aq_feed,
            pH_feed=pH_feed,
            Q_aq=Q_aq,
            Q_org=Q_org,
            extractant=extractant,
            C_ext=C_ext,
            n_stages=n_stages,
            target_pH=None,
            target_pH_per_stage=None,
            C_NaOH=C_NaOH,
            Q_NaOH=Q_NaOH,
            naoh_mode=naoh_mode,
            saponification_model="legacy_equivalent_target",
            saponification_fraction=saponification_fraction,
            naoh_strategy=naoh_strategy,
            naoh_weights=naoh_weights,
            C_org_fresh=C_org_fresh,
            metals=metals,
            tolerance=tolerance,
            max_iter=max_iter,
            temperature=temperature,
            C_sulfate=C_sulfate,
            use_competition=use_competition,
            use_speciation=use_speciation,
            extractant_params=extractant_params,
        )
        initial_org_seed = [
            stage["C_org_out"] for stage in legacy_seed_result["stages"]
        ]

    primary_result = _run_with_q_naoh(
        Q_NaOH,
        initial_org_out=initial_org_seed,
        solver_strategy="direct",
    )
    if primary_result["converged"]:
        primary_result["objective_converged"] = primary_result["converged"]
        primary_result["objective_residual_pH"] = None
        return primary_result

    has_target_mode = any(target is not None for target in stage_target_pHs)

    if has_target_mode:
        continuation_levels = [0.35, 0.65, 0.85, 1.0]
        continuation_seed = None
        best_result = primary_result

        for level in continuation_levels:
            continuation_targets = _interpolate_stage_targets(
                stage_target_pHs, pH_feed, level
            )
            continuation_result = _run_with_q_naoh(
                Q_NaOH,
                stage_targets=continuation_targets,
                initial_org_out=continuation_seed,
                initial_sap_out=None,
                relaxation_scale=0.85 if level < 1.0 else 0.75,
                solver_strategy=f"target_continuation_{level:.2f}",
            )
            continuation_seed = [
                stage_result["C_org_out"] for stage_result in continuation_result["stages"]
            ]

            if continuation_result["converged"]:
                continuation_result["objective_converged"] = continuation_result["converged"]
                continuation_result["objective_residual_pH"] = None
                return continuation_result
            if continuation_result["best_residual"] < best_result["best_residual"]:
                best_result = continuation_result

        best_result["objective_converged"] = best_result["converged"]
        best_result["objective_residual_pH"] = None
        return best_result

    fallback_result = _run_with_q_naoh(
        Q_NaOH,
        initial_org_out=[stage["C_org_out"] for stage in primary_result["stages"]],
        initial_sap_out=[
            stage.get("saponified_extractant_out_mol_hr", 0.0)
            for stage in primary_result["stages"]
        ] if alkali_contract.application == "fresh_organic_saponification" and not has_target_mode else None,
        relaxation_scale=0.70,
        solver_strategy="fixed_naoh_fallback",
    )
    if fallback_result["converged"] or fallback_result["best_residual"] < primary_result["best_residual"]:
        fallback_result["objective_converged"] = fallback_result["converged"]
        fallback_result["objective_residual_pH"] = None
        return fallback_result
    primary_result["objective_converged"] = primary_result["converged"]
    primary_result["objective_residual_pH"] = None
    return primary_result


def print_multistage_result(result: dict, C_aq_feed: dict = None):
    """다단 시뮬레이션 결과를 포맷팅하여 출력합니다."""
    n_stages = len(result["stages"])
    metals = list(result["raffinate"].keys())

    print(f"\n{'='*70}")
    print(f"  다단 역류 Mixer-Settler SX 시뮬레이션 결과")
    print(f"{'='*70}")
    print(f"  Stage 수: {n_stages} | 수렴: {'Yes' if result['converged'] else 'No'}"
          f" | 반복: {result['iterations']}")

    # pH 프로파일
    print(f"\n  Stage별 pH 및 NaOH 소비량:")
    print(f"  {'Stage':>7} | {'pH':>6} | {'NaOH(mol/hr)':>13} | pH 바")
    print(f"  {'-'*7}-+-{'-'*6}-+-{'-'*13}-+-{'-'*20}")
    for i in range(n_stages):
        pH = result["pH_profile"][i]
        naoh = result["NaOH_profile"][i]
        bar = "█" * int(pH * 3)
        print(f"  Stage {i+1} | {pH:>5.2f} | {naoh:>13.2f} | {bar}")
    print(f"  {'':>7}   {'합계':>6}   {result['total_NaOH_mol_hr']:>13.2f}")

    # 전체 추출률
    print(f"\n  {'금속':>6} | {'피드(g/L)':>10} | {'후액(g/L)':>10} | {'추출률(%)':>10}")
    print(f"  {'-'*6}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}")
    for m in metals:
        c_feed = C_aq_feed.get(m, 0.0) if C_aq_feed else 0.0
        c_raff = result["raffinate"][m]
        ext = result["overall_extraction"][m]
        print(f"  {m:>6} | {c_feed:>10.4f} | {c_raff:>10.4f} | {ext:>10.2f}")

    print(f"\n  후액 최종 pH: {result['pH_profile'][-1]:.2f}")
    print(f"  총 NaOH 소비: {result['total_NaOH_mol_hr']:.2f} mol/hr")
    print(f"{'='*70}")
