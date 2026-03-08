import sys
import os

sys.path.append(os.path.abspath('.'))

from sx_simulator.multistage_sx import solve_multistage_countercurrent
from sx_simulator.config import DEFAULT_METALS
from sx_simulator.datasets import (
    VERIFICATION_EXPERIMENTS as experiments,
    KNOWN_INVALID_POINTS,
    calc_sulfate_from_feed,
)

def main():
    overall_errors = {metal: [] for metal in DEFAULT_METALS}

    for exp in experiments:
        print(f"--- Running {exp['name']} ---")
        
        Q_aq = exp['feed_flow'] + exp['naoh_flow']
        Q_org = exp['org_flow']
        dilution_factor = exp['feed_flow'] / Q_aq
        
        diluted_inputs = {}
        for m in DEFAULT_METALS:
            diluted_inputs[m] = exp['input'].get(m, 0.0) * dilution_factor
        sulfate = calc_sulfate_from_feed(diluted_inputs)
            
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
            C_sulfate=sulfate,
            use_competition=True,
            use_speciation=True,
            target_pH=exp['pH']
        )
        
        try:
            result = solve_multistage_countercurrent(**sim_kwargs)
            res_df = result['raffinate']
            abs_errors = []
            worst_row = None
            
            print(f"Dilution Factor: {dilution_factor:.4f}, O/A Ratio: {Q_org/Q_aq:.4f}, Sulfate: {sulfate:.4f} M")
            print(f"| Metal | In_Original | In_Diluted | Sim_Raffinate | Exp_Output | Sim_Error |")
            print(f"|---|---|---|---|---|---|")
            dataset_tag = exp["name"].split()[0]
            for m in DEFAULT_METALS:
                in_orig = exp['input'].get(m, 0.0)
                in_diluted = diluted_inputs[m]
                sim_out = res_df[m]
                exp_out = exp['output'].get(m, 0.0)
                error = sim_out - exp_out

                exclusion_reason = KNOWN_INVALID_POINTS.get((dataset_tag, m))
                if exclusion_reason:
                    print(f"| {m} | {in_orig:.3f} | {in_diluted:.3f} | {sim_out:.3f} | {exp_out:.3f} | 제외 |")
                    print(f"  Note: {exclusion_reason}")
                    continue

                abs_errors.append(abs(error))
                overall_errors[m].append(abs(error))
                if worst_row is None or abs(error) > abs(worst_row[1]):
                    worst_row = (m, error)
                print(f"| {m} | {in_orig:.3f} | {in_diluted:.3f} | {sim_out:.3f} | {exp_out:.3f} | {error:+.3f} |")
            mae = sum(abs_errors) / len(abs_errors) if abs_errors else 0.0
            if worst_row is not None:
                print(f"Summary: converged={result['converged']} iterations={result['iterations']} dataset_MAE={mae:.3f} g/L worst={worst_row[0]} {worst_row[1]:+.3f} g/L")
        except Exception as e:
            print(f"Simulation failed: {e}")
        print("\n")

    print("--- Overall Metal MAE Summary ---")
    for metal in DEFAULT_METALS:
        if overall_errors[metal]:
            metal_mae = sum(overall_errors[metal]) / len(overall_errors[metal])
            print(f"{metal}: MAE={metal_mae:.3f} g/L over {len(overall_errors[metal])} datasets")

if __name__ == "__main__":
    main()
