import sys
import os
import pandas as pd

sys.path.append(os.path.abspath('.'))

from sx_simulator.multistage_sx import solve_multistage_countercurrent
from sx_simulator.config import DEFAULT_METALS
from test_verification import experiments

def run_sim(use_phase3):
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
            use_competition=use_phase3,
            use_speciation=use_phase3,
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
    res_off = run_sim(False)
    res_on = run_sim(True)
    
    for exp in experiments:
        name = exp['name']
        print(f"--- {name} ---")
        if res_off[name] is None or res_on[name] is None:
            print("Simulation failed, skipping.")
            continue
            
        print(f"| Metal | Exp Output | Sim (Phase3 OFF) | Sim (Phase3 ON) | Diff (ON - OFF) |")
        print(f"|---|---|---|---|---|")
        for m in DEFAULT_METALS:
            if m == 'Li' and name.startswith('Data5'):
                continue
            exp_val = exp['output'].get(m, 0.0)
            off_val = res_off[name][m]
            on_val = res_on[name][m]
            diff = on_val - off_val
            
            # Print only if there's a meaningful change or it's a major metal
            if abs(diff) > 0.01 or m in ['Ni', 'Co', 'Mg', 'Mn', 'Li']:
                print(f"| {m} | {exp_val:.3f} | {off_val:.3f} | {on_val:.3f} | {diff:+.3f} |")
        print()

if __name__ == "__main__":
    main()
