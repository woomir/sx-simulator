"""
compare_phase3.py — Vasilyev 경쟁 추출 모델 단일 모드 확인 스크립트.

기존에는 Phase 3 ON/OFF를 비교했으나, Vasilyev 모델이 기본 엔진으로
통합되었으므로 단일 모드로 실행하여 결과를 확인합니다.
"""
import sys
import os

sys.path.append(os.path.abspath('.'))

from sx_simulator.multistage_sx import solve_multistage_countercurrent
from sx_simulator.config import DEFAULT_METALS
from test_verification import experiments


def run_sim():
    results = {}
    for exp in experiments:
        Q_aq = exp['feed_flow'] + exp['naoh_flow']
        Q_org = exp['org_flow']
        dilution_factor = exp['feed_flow'] / Q_aq

        diluted_inputs = {}
        for m in DEFAULT_METALS:
            diluted_inputs[m] = exp['input'].get(m, 0.0) * dilution_factor

        sim_kwargs = dict(
            C_aq_feed=diluted_inputs,
            pH_feed=exp['pH'],
            Q_aq=Q_aq,
            Q_org=Q_org,
            extractant=exp['ext'],
            C_ext=exp['C_ext'],
            n_stages=exp['n_stages'],
            metals=DEFAULT_METALS,
            temperature=exp['T'],
            C_sulfate=0.3,
            use_speciation=True,
            target_pH=exp['pH']
        )

        try:
            res = solve_multistage_countercurrent(**sim_kwargs)
            results[exp['name']] = res['raffinate']
        except Exception as e:
            print(f"Failed for {exp['name']}: {e}")
            results[exp['name']] = None
    return results


def main():
    results = run_sim()

    for exp in experiments:
        name = exp['name']
        print(f"--- {name} ---")
        if results[name] is None:
            print("Simulation failed, skipping.")
            continue

        print(f"| Metal | Exp Output | Sim (Vasilyev) |")
        print(f"|---|---|---|")
        for m in DEFAULT_METALS:
            if m == 'Li' and name.startswith('Data5'):
                continue
            exp_val = exp['output'].get(m, 0.0)
            sim_val = results[name][m]
            if m in ['Ni', 'Co', 'Mg', 'Mn', 'Li']:
                print(f"| {m} | {exp_val:.3f} | {sim_val:.3f} |")
        print()


if __name__ == "__main__":
    main()
