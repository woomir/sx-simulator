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
from .single_stage import solve_single_stage
from .config import DEFAULT_METALS, CONVERGENCE_TOLERANCE, MAX_ITERATIONS


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
    naoh_strategy: str = "uniform",
    naoh_weights: list = None,
    C_org_fresh: dict = None,
    metals: list = None,
    tolerance: float = None,
    max_iter: int = None,
    temperature: float = None,
    model_type: str = "sigmoid",
    C_sulfate: float = 0.0,
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

    # 목표 pH 결정
    use_bisection = False
    if target_pH_per_stage is not None:
        stage_target_pHs = target_pH_per_stage
    elif target_pH is not None:
        # 단일 목표 pH만 주어졌을 때: 전체 Q_NaOH 최적화
        stage_target_pHs = [None] * n_stages
        use_bisection = True
    else:
        stage_target_pHs = [None] * n_stages

    def _run_with_q_naoh(test_q_naoh: float) -> dict:
        # NaOH 분배 (고정 NaOH 모드용)
        if naoh_weights is not None and len(naoh_weights) == n_stages:
            weights = naoh_weights
        elif naoh_strategy == "front_loaded":
            weights = [0.5 ** i for i in range(n_stages)]
        else: # "uniform"
            weights = [1.0] * n_stages
            
        w_sum = sum(weights)
        if w_sum > 0:
            Q_NaOH_dist = [test_q_naoh * (w / w_sum) for w in weights]
        else:
            Q_NaOH_dist = [0.0] * n_stages

        org_out = [{m: 0.0 for m in metals} for _ in range(n_stages)]
        converged = False

        for iteration in range(max_iter):
            org_out_prev = copy.deepcopy(org_out)
            stage_results = []
            C_aq_current = copy.deepcopy(C_aq_feed)
            pH_current = pH_feed

            for s in range(n_stages):
                C_org_in = org_out[s + 1] if s < n_stages - 1 else copy.deepcopy(C_org_fresh)
                stage_num = s + 1
                t_pH = stage_target_pHs[s]

                if t_pH is not None:
                    result = solve_single_stage(
                        C_aq_in=C_aq_current, C_org_in=C_org_in,
                        pH_in=pH_current, Q_aq=Q_aq, Q_org=Q_org,
                        extractant=extractant, C_ext=C_ext,
                        target_pH=t_pH, metals=metals,
                        temperature=temperature,
                        model_type=model_type,
                        C_sulfate=C_sulfate,
                    )
                else:
                    q_naoh = Q_NaOH_dist[stage_num - 1]
                    result = solve_single_stage(
                        C_aq_in=C_aq_current, C_org_in=C_org_in,
                        pH_in=pH_current, Q_aq=Q_aq, Q_org=Q_org,
                        extractant=extractant, C_ext=C_ext,
                        C_NaOH=C_NaOH, Q_NaOH=q_naoh, metals=metals,
                        temperature=temperature,
                        model_type=model_type,
                        C_sulfate=C_sulfate,
                    )

                stage_results.append(result)
                org_out[s] = copy.deepcopy(result["C_org_out"])
                C_aq_current = copy.deepcopy(result["C_aq_out"])
                pH_current = result["pH_out"]

            max_diff = 0.0
            for s in range(n_stages):
                for m in metals:
                    max_diff = max(max_diff, abs(org_out[s][m] - org_out_prev[s][m]))
            if max_diff < tolerance:
                converged = True
                break

        raffinate = stage_results[-1]["C_aq_out"]
        pH_profile = [r["pH_out"] for r in stage_results]
        overall_extraction = {}
        for m in metals:
            C_feed = C_aq_feed.get(m, 0.0)
            C_raff = raffinate.get(m, 0.0)
            overall_extraction[m] = (1.0 - C_raff / C_feed) * 100.0 if C_feed > 0 else 0.0

        return {
            "stages": stage_results,
            "raffinate": raffinate,
            "loaded_organic": stage_results[0]["C_org_out"],
            "pH_profile": pH_profile,
            "NaOH_profile": [r.get("NaOH_consumed_mol_hr", 0) for r in stage_results],
            "overall_extraction": overall_extraction,
            "converged": converged,
            "iterations": iteration + 1,
            "total_NaOH_mol_hr": sum(r.get("NaOH_consumed_mol_hr", 0) for r in stage_results),
        }

    if not use_bisection:
        return _run_with_q_naoh(Q_NaOH)

    # Bisection search to find Q_NaOH that results in target_pH
    low_q = 0.0
    high_q = float(Q_aq * 2.0) # Start with a reasonable upper bound
    
    res_low = _run_with_q_naoh(low_q)
    if res_low["pH_profile"][-1] >= target_pH:
        return res_low
        
    res_high = _run_with_q_naoh(high_q)
    for _ in range(5):
        if res_high["pH_profile"][-1] >= target_pH:
            break
        high_q *= 2.0
        res_high = _run_with_q_naoh(high_q)
        
    best_res = res_high
    for _ in range(30):
        mid_q = (low_q + high_q) / 2.0
        best_res = _run_with_q_naoh(mid_q)
        mid_ph = best_res["pH_profile"][-1]
        
        if abs(mid_ph - target_pH) < 0.005:
            break
            
        if mid_ph < target_pH:
            low_q = mid_q
        else:
            high_q = mid_q
            
    return best_res


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
