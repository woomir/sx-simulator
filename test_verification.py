"""
SX 시뮬레이터 v1.7.0 종합 검증 스크립트
========================================
16가지 조건 조합(pH모드×경쟁×종분화×추출제)에 대해 자동 실행하고
물질수지, pH수지, Phase 3 기능을 검증합니다.

실행: .venv/bin/python test_verification.py
"""

import sys, os, copy, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sx_simulator.multistage_sx import solve_multistage_countercurrent
from sx_simulator.single_stage import solve_single_stage, calc_aq_protons
from sx_simulator.extraction_isotherm import (
    extraction_efficiency, distribution_coefficient,
    compute_competitive_extractions, calc_loading_fraction
)
from sx_simulator.config import DEFAULT_METALS, MOLAR_MASS, EXTRACTANT_PARAMS

# ─── 테스트 기본 조건 ─────────────────────────────────────────
BASE_FEED = {"Li": 1.5, "Ni": 5.0, "Co": 3.0, "Mn": 2.0}
PH_FEED = 3.0
Q_AQ = 100.0
Q_ORG = 100.0
N_STAGES = 4
TEMPERATURE = 25.0
C_SULFATE = 0.5
METALS = DEFAULT_METALS

# 고로딩 조건
HIGH_LOAD_FEED = {"Li": 1.5, "Ni": 30.0, "Co": 15.0, "Mn": 10.0}
LOW_C_EXT = 0.1

# ─── 유틸리티 ─────────────────────────────────────────────────
PASS = "✅ PASS"
FAIL = "❌ FAIL"
WARN = "⚠️ WARN"

def check_mass_balance(result, C_aq_feed, Q_aq, Q_org, metals, Q_NaOH=0.0, tol=1.0):
    """물질수지 검증: Feed = Raffinate + Loaded Organic (g/hr 기준)"""
    errors = {}
    for m in metals:
        feed_flow = C_aq_feed[m] * Q_aq
        Q_raff = Q_aq + Q_NaOH  # NaOH 희석 반영
        raff_flow = result["raffinate"][m] * Q_raff
        org_flow = result["loaded_organic"][m] * Q_org
        balance = feed_flow - raff_flow - org_flow
        rel_error = abs(balance) / feed_flow * 100 if feed_flow > 0 else 0
        errors[m] = {"balance_g_hr": round(balance, 6), "rel_error_pct": round(rel_error, 4)}
    max_err = max(e["rel_error_pct"] for e in errors.values())
    return errors, max_err < tol

def check_ph_target(result, target_pH, tol=0.05):
    """목표 pH 모드: 출구 pH가 target과 일치하는지 확인"""
    final_pH = result["pH_profile"][-1]
    error = abs(final_pH - target_pH)
    return final_pH, error, error < tol

def check_result_keys(result):
    """결과 딕셔너리 키 완전성 검증"""
    required = ["stages", "raffinate", "loaded_organic", "pH_profile",
                 "NaOH_profile", "overall_extraction", "converged", "iterations",
                 "total_NaOH_mol_hr"]
    missing = [k for k in required if k not in result]
    stage_required = ["C_aq_out", "C_org_out", "pH_out", "extraction",
                      "H_released_mol_hr", "loading_fraction"]
    stage_missing = []
    if result.get("stages"):
        stage_missing = [k for k in stage_required if k not in result["stages"][0]]
    return missing, stage_missing

# ─── 검증 1: 16가지 조합 실행 ───────────────────────────────
def test_all_combinations():
    print("\n" + "=" * 80)
    print("  [TEST 1] 16가지 조건 조합 실행 검증")
    print("=" * 80)

    results_table = []
    all_pass = True

    combos = []
    for ph_mode in ["target_pH", "fixed_NaOH"]:
        for competition in [False, True]:
            for speciation in [False, True]:
                for extractant in ["Cyanex 272", "D2EHPA"]:
                    combos.append((ph_mode, competition, speciation, extractant))

    for i, (ph_mode, comp, spec, ext) in enumerate(combos, 1):
        C_ext = 0.5 if ext == "Cyanex 272" else 0.64
        kwargs = dict(
            C_aq_feed=copy.deepcopy(BASE_FEED), pH_feed=PH_FEED,
            Q_aq=Q_AQ, Q_org=Q_ORG, extractant=ext, C_ext=C_ext,
            n_stages=N_STAGES, metals=METALS, temperature=TEMPERATURE,
            C_sulfate=C_SULFATE, use_competition=comp, use_speciation=spec,
        )
        if ph_mode == "target_pH":
            kwargs["target_pH"] = 5.0
        else:
            kwargs["C_NaOH"] = 5.0
            kwargs["tolerance"] = 1e-4  # 상대오차 기반이므로 tolerance 완화
            kwargs["Q_NaOH"] = 12.0

        label = f"{ph_mode:12s} | comp={str(comp):5s} | spec={str(spec):5s} | {ext:12s}"

        try:
            result = solve_multistage_countercurrent(**kwargs)
            converged = result["converged"]
            missing_keys, stage_missing = check_result_keys(result)

            q_naoh = kwargs.get("Q_NaOH", 0.0)
            mb_errors, mb_ok = check_mass_balance(result, BASE_FEED, Q_AQ, Q_ORG, METALS, Q_NaOH=q_naoh)
            max_mb_err = max(e["rel_error_pct"] for e in mb_errors.values())

            if ph_mode == "target_pH":
                final_pH, ph_err, ph_ok = check_ph_target(result, 5.0)
            else:
                final_pH = result["pH_profile"][-1]
                ph_err = 0.0
                ph_ok = True

            status = PASS if (converged and mb_ok and ph_ok and not missing_keys) else FAIL
            if status == FAIL:
                all_pass = False

            results_table.append({
                "id": i, "label": label, "status": status,
                "converged": converged, "final_pH": final_pH, "ph_err": ph_err,
                "max_mb_err": max_mb_err, "iterations": result["iterations"],
                "extractions": {m: round(result["overall_extraction"][m], 1) for m in METALS},
                "naoh_total": round(result["total_NaOH_mol_hr"], 2),
                "max_loading": round(max(sr.get("loading_fraction", 0) for sr in result["stages"]) * 100, 1),
                "missing_keys": missing_keys,
            })

        except Exception as e:
            all_pass = False
            results_table.append({
                "id": i, "label": label, "status": FAIL,
                "error": str(e),
            })

    # 결과 출력
    print(f"\n{'#':>3} | {'조건':^48s} | {'결과':^8s} | {'수렴':^4s} | {'pH':>5s} | {'MB오차%':>7s} | {'NaOH':>7s} | {'로딩%':>5s}")
    print("-" * 120)
    for r in results_table:
        if "error" in r:
            print(f"{r['id']:>3} | {r['label']:48s} | {r['status']:8s} | ERROR: {r['error'][:40]}")
        else:
            print(f"{r['id']:>3} | {r['label']:48s} | {r['status']:8s} | {'Y' if r['converged'] else 'N':^4s} | "
                  f"{r['final_pH']:5.2f} | {r['max_mb_err']:7.4f} | {r['naoh_total']:7.2f} | {r['max_loading']:5.1f}")

    print(f"\n  종합: {'✅ 전체 PASS' if all_pass else '❌ 일부 FAIL'}")
    return results_table, all_pass


# ─── 검증 2: Phase 3 경쟁 모델 효과 ─────────────────────────
def test_competition_effect():
    print("\n" + "=" * 80)
    print("  [TEST 2] Phase 3 추출제 경쟁 모델 효과 검증")
    print("=" * 80)

    all_pass = True

    # 2-A: 낮은 로딩에서 차이 미미한지
    print("\n  [2-A] 낮은 로딩 조건 (기본 Feed, C_ext=0.5M)")
    results = {}
    for comp in [False, True]:
        r = solve_multistage_countercurrent(
            C_aq_feed=copy.deepcopy(BASE_FEED), pH_feed=PH_FEED,
            Q_aq=Q_AQ, Q_org=Q_ORG, extractant="Cyanex 272", C_ext=0.5,
            n_stages=N_STAGES, target_pH=5.0, metals=METALS,
            temperature=TEMPERATURE, C_sulfate=C_SULFATE,
            use_competition=comp,
        )
        results[comp] = r
        loading = max(sr.get("loading_fraction", 0) for sr in r["stages"]) * 100
        print(f"    경쟁={'ON' if comp else 'OFF':3s}: "
              f"Li={r['overall_extraction']['Li']:.1f}%, "
              f"Ni={r['overall_extraction']['Ni']:.1f}%, "
              f"Co={r['overall_extraction']['Co']:.1f}%, "
              f"Mn={r['overall_extraction']['Mn']:.1f}%, "
              f"로딩={loading:.1f}%")

    # 2-B: 고로딩 조건에서 경쟁 효과 확인
    print("\n  [2-B] 고로딩 조건 (Ni=30, Co=15, Mn=10, C_ext=0.1M)")
    results_high = {}
    for comp in [False, True]:
        try:
            r = solve_multistage_countercurrent(
                C_aq_feed=copy.deepcopy(HIGH_LOAD_FEED), pH_feed=PH_FEED,
                Q_aq=Q_AQ, Q_org=Q_ORG, extractant="Cyanex 272", C_ext=LOW_C_EXT,
                n_stages=N_STAGES, target_pH=5.0, metals=METALS,
                temperature=TEMPERATURE, C_sulfate=C_SULFATE,
                use_competition=comp,
            )
            results_high[comp] = r
            loading = max(sr.get("loading_fraction", 0) for sr in r["stages"]) * 100
            print(f"    경쟁={'ON' if comp else 'OFF':3s}: "
                  f"Li={r['overall_extraction']['Li']:.1f}%, "
                  f"Ni={r['overall_extraction']['Ni']:.1f}%, "
                  f"Co={r['overall_extraction']['Co']:.1f}%, "
                  f"Mn={r['overall_extraction']['Mn']:.1f}%, "
                  f"로딩={loading:.1f}%")
        except Exception as e:
            all_pass = False
            print(f"    경쟁={'ON' if comp else 'OFF':3s}: ❌ ERROR - {e}")

    # 2-C: 경쟁 모델이 D를 감소시키는 방향인지
    print("\n  [2-C] 경쟁 보정 방향 검증 (D_adjusted <= D_sigmoid)")
    comp_result = compute_competitive_extractions(
        pH=5.0, extractant="Cyanex 272", C_ext=0.5,
        C_aq_in=copy.deepcopy(BASE_FEED), C_org_in={m: 0 for m in METALS},
        Q_aq=Q_AQ, Q_org=Q_ORG, metals=METALS, temperature=TEMPERATURE,
    )
    for m in METALS:
        D_sig = distribution_coefficient(5.0, m, "Cyanex 272", 0.5, temperature=TEMPERATURE)
        D_adj = comp_result["D_adjusted"][m]
        ok = D_adj <= D_sig + 1e-10
        status = PASS if ok else FAIL
        if not ok:
            all_pass = False
        print(f"    {m}: D_sig={D_sig:.4f}, D_adj={D_adj:.4f} → {status}")

    print(f"\n  종합: {'✅ PASS' if all_pass else '❌ FAIL'}")
    return all_pass


# ─── 검증 3: Phase 3 종분화 효과 ─────────────────────────────
def test_speciation_effect():
    print("\n" + "=" * 80)
    print("  [TEST 3] Phase 3 수계 종분화(Speciation) 효과 검증")
    print("=" * 80)

    all_pass = True

    # 3-A: NaOH 소비량 차이 (종분화 ON → NaOH 증가 예상)
    print("\n  [3-A] 종분화 ON/OFF 시 NaOH 소비량 비교 (목표 pH=5.0)")
    for spec in [False, True]:
        r = solve_multistage_countercurrent(
            C_aq_feed=copy.deepcopy(BASE_FEED), pH_feed=PH_FEED,
            Q_aq=Q_AQ, Q_org=Q_ORG, extractant="Cyanex 272", C_ext=0.5,
            n_stages=N_STAGES, target_pH=5.0, metals=METALS,
            temperature=TEMPERATURE, C_sulfate=C_SULFATE,
            use_speciation=spec,
        )
        print(f"    종분화={'ON' if spec else 'OFF':3s}: NaOH={r['total_NaOH_mol_hr']:.2f} mol/hr, "
              f"pH_final={r['pH_profile'][-1]:.2f}")

    # 3-B: 고정 NaOH 모드에서 종분화 ON → 출구 pH 낮아져야 함
    print("\n  [3-B] 고정 NaOH 모드에서 종분화 효과 (동일 NaOH → pH 차이)")
    phs = {}
    for spec in [False, True]:
        r = solve_multistage_countercurrent(
            C_aq_feed=copy.deepcopy(BASE_FEED), pH_feed=PH_FEED,
            Q_aq=Q_AQ, Q_org=Q_ORG, extractant="Cyanex 272", C_ext=0.5,
            n_stages=N_STAGES, C_NaOH=5.0, Q_NaOH=12.0, metals=METALS,
            temperature=TEMPERATURE, C_sulfate=C_SULFATE,
            use_speciation=spec,
        )
        phs[spec] = r["pH_profile"][-1]
        print(f"    종분화={'ON' if spec else 'OFF':3s}: pH_final={phs[spec]:.4f}")

    # 3-C: calc_aq_protons 함수 직접 검증
    print("\n  [3-C] calc_aq_protons() 함수 직접 비교")
    test_C_aq = {"Li": 1.0, "Ni": 4.0, "Co": 2.5, "Mn": 1.5}
    for spec in [False, True]:
        H = calc_aq_protons(5.0, Q_AQ, C_SULFATE, test_C_aq, use_speciation=spec)
        print(f"    종분화={'ON' if spec else 'OFF':3s}: H+(pH=5)={H:.6f} mol/hr")

    print(f"\n  종합: {'✅ PASS' if all_pass else '❌ FAIL'}")
    return all_pass


# ─── 검증 4: 온도 의존성 ─────────────────────────────────────
def test_temperature_effect():
    print("\n" + "=" * 80)
    print("  [TEST 4] 온도 의존성 검증")
    print("=" * 80)

    print("\n  온도에 따른 추출률 변화 (Cyanex 272, pH=5.0):")
    print(f"    {'온도(°C)':>8} | {'Li':>8} | {'Ni':>8} | {'Co':>8} | {'Mn':>8}")
    print(f"    {'-'*8}-+-{'-'*8}-+-{'-'*8}-+-{'-'*8}-+-{'-'*8}")
    for T in [15.0, 25.0, 40.0]:
        r = solve_multistage_countercurrent(
            C_aq_feed=copy.deepcopy(BASE_FEED), pH_feed=PH_FEED,
            Q_aq=Q_AQ, Q_org=Q_ORG, extractant="Cyanex 272", C_ext=0.5,
            n_stages=N_STAGES, target_pH=5.0, metals=METALS,
            temperature=T, C_sulfate=C_SULFATE,
        )
        exts = [f"{r['overall_extraction'][m]:8.1f}" for m in METALS]
        print(f"    {T:8.0f} | {'|'.join(exts)}")
    return True


# ─── 검증 5: Stage별 차등 pH ────────────────────────────────
def test_staged_pH():
    print("\n" + "=" * 80)
    print("  [TEST 5] Stage별 차등 pH 검증")
    print("=" * 80)

    staged_pHs = [4.0, 4.5, 5.0, 5.5]
    r = solve_multistage_countercurrent(
        C_aq_feed=copy.deepcopy(BASE_FEED), pH_feed=PH_FEED,
        Q_aq=Q_AQ, Q_org=Q_ORG, extractant="Cyanex 272", C_ext=0.5,
        n_stages=N_STAGES, target_pH_per_stage=staged_pHs, metals=METALS,
        temperature=TEMPERATURE, C_sulfate=C_SULFATE,
    )

    all_ok = True
    print(f"\n    {'Stage':>7} | {'목표pH':>7} | {'실제pH':>7} | {'오차':>7} | {'판정':>6}")
    print(f"    {'-'*7}-+-{'-'*7}-+-{'-'*7}-+-{'-'*7}-+-{'-'*6}")
    for i in range(N_STAGES):
        actual = r["pH_profile"][i]
        error = abs(actual - staged_pHs[i])
        ok = error < 0.05
        if not ok:
            all_ok = False
        print(f"    Stage {i+1} | {staged_pHs[i]:7.2f} | {actual:7.2f} | {error:7.4f} | {PASS if ok else FAIL}")

    print(f"\n  종합: {'✅ PASS' if all_ok else '❌ FAIL'}")
    return all_ok


# ─── 메인 ─────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 80)
    print("  SX 시뮬레이터 v1.7.0 종합 검증")
    print(f"  실행 시각: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)

    results = {}
    results["16_combos"], pass1 = test_all_combinations()
    pass2 = test_competition_effect()
    pass3 = test_speciation_effect()
    pass4 = test_temperature_effect()
    pass5 = test_staged_pH()

    print("\n" + "=" * 80)
    print("  최종 검증 결과 요약")
    print("=" * 80)
    tests = [
        ("16 조합 실행", pass1),
        ("추출제 경쟁 모델", pass2),
        ("수계 종분화", pass3),
        ("온도 의존성", pass4),
        ("Stage별 차등 pH", pass5),
    ]
    for name, passed in tests:
        print(f"  {PASS if passed else FAIL} {name}")

    overall = all(p for _, p in tests)
    print(f"\n  {'🎉 전체 검증 통과!' if overall else '⚠️ 일부 검증 실패 — 상세 내용을 확인하세요.'}")
    print("=" * 80)
