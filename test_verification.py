import sys
import os
import pandas as pd

sys.path.append(os.path.abspath('.'))

from sx_simulator.multistage_sx import solve_multistage_countercurrent
from sx_simulator.config import DEFAULT_METALS

experiments = [
    {
        "name": "Data1 (CoSX-D9-data, 5 Stages)",
        "ext": "Cyanex 272", "C_ext": 0.6308, "pH": 6.10, "n_stages": 5, "T": 25.0,
        "feed_flow": 25.0, "naoh_flow": 11.0, "org_flow": 118.0,
        "input": {"Li": 5.283, "Ni": 33.894, "Co": 3.096, "Mn": 0.067, "Ca": 0.002, "Mg": 0.026, "Zn": 0.0},
        "output": {"Li": 3.322, "Ni": 7.166, "Co": 0.001, "Mn": 0.0, "Ca": 0.0, "Mg": 0.0, "Zn": 0.0}
    },
    {
        "name": "Data2 (CoSX-D5-data, 5 Stages)",
        "ext": "Cyanex 272", "C_ext": 0.6308, "pH": 5.89, "n_stages": 5, "T": 25.0,
        "feed_flow": 25.0, "naoh_flow": 9.4, "org_flow": 118.0,
        "input": {"Li": 5.260, "Ni": 33.876, "Co": 3.053, "Mn": 0.071, "Ca": 0.001, "Mg": 0.030, "Zn": 0.0},
        "output": {"Li": 3.478, "Ni": 8.951, "Co": 0.0, "Mn": 0.0, "Ca": 0.0, "Mg": 0.0, "Zn": 0.0}
    },
    {
        "name": "Data3 (CoSX-D2-data, 5 Stages)",
        "ext": "Cyanex 272", "C_ext": 0.6308, "pH": 7.00, "n_stages": 5, "T": 25.0,
        "feed_flow": 20.0, "naoh_flow": 9.4, "org_flow": 118.0,
        "input": {"Li": 5.852, "Ni": 28.476, "Co": 4.106, "Mn": 0.091, "Ca": 0.002, "Mg": 0.035, "Zn": 0.0},
        "output": {"Li": 3.151, "Ni": 0.044, "Co": 0.0, "Mn": 0.0, "Ca": 0.0, "Mg": 0.0, "Zn": 0.0}
    },
    {
        "name": "Data4 (IMSX-D5-Data, 5 Stages)",
        "ext": "D2EHPA", "C_ext": 0.6053, "pH": 4.00, "n_stages": 5, "T": 25.0,
        "feed_flow": 30.0, "naoh_flow": 4.2, "org_flow": 124.0,
        "input": {"Li": 9.767, "Ni": 28.889, "Co": 16.178, "Mn": 15.204, "Ca": 0.336, "Mg": 0.149, "Zn": 0.006},
        "output": {"Li": 6.335, "Ni": 18.186, "Co": 0.230, "Mn": 0.0, "Ca": 0.0, "Mg": 0.001, "Zn": 0.0}
    },
    {
        "name": "Data5 (IMSX-D9-Data, 5 Stages)",
        "ext": "D2EHPA", "C_ext": 0.6053, "pH": 4.20, "n_stages": 5, "T": 25.0,
        "feed_flow": 25.0, "naoh_flow": 3.2, "org_flow": 100.0,
        "input": {"Li": 0.013, "Ni": 71.000, "Co": 8.700, "Mn": 8.900, "Ca": 0.250, "Mg": 3.000, "Zn": 0.600},
        "output": {"Li": 5.400, "Ni": 16.000, "Co": 0.160, "Mn": 0.0, "Ca": 0.0, "Mg": 0.0, "Zn": 0.0}
    },
    {
        "name": "Data6 (IMSX-D?, 5 Stages)",
        "ext": "D2EHPA", "C_ext": 0.6053, "pH": 3.90, "n_stages": 5, "T": 25.0,
        "feed_flow": 30.0, "naoh_flow": 3.5, "org_flow": 100.0,
        "input": {"Li": 8.390, "Ni": 27.917, "Co": 15.400, "Mn": 12.883, "Ca": 0.413, "Mg": 0.129, "Zn": 0.012},
        "output": {"Li": 6.591, "Ni": 21.034, "Co": 0.848, "Mn": 0.0, "Ca": 0.003, "Mg": 0.0, "Zn": 0.0}
    }
]

def main():
    for exp in experiments:
        print(f"--- Running {exp['name']} ---")
        
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
            C_sulfate=0.3, # approximate, just needed to prevent None
            use_competition=True,
            use_speciation=True,
            target_pH=exp['pH']
        )
        
        try:
            result = solve_multistage_countercurrent(**sim_kwargs)
            res_df = result['raffinate']
            
            print(f"Dilution Factor: {dilution_factor:.4f}, O/A Ratio: {Q_org/Q_aq:.4f}")
            print(f"| Metal | In_Original | In_Diluted | Sim_Raffinate | Exp_Output | Sim_Error |")
            print(f"|---|---|---|---|---|---|")
            for m in DEFAULT_METALS:
                in_orig = exp['input'].get(m, 0.0)
                in_diluted = diluted_inputs[m]
                sim_out = res_df[m]
                exp_out = exp['output'].get(m, 0.0)
                error = sim_out - exp_out
                
                if m == 'Li' and exp['name'].startswith('Data5'):
                    continue
                    
                print(f"| {m} | {in_orig:.3f} | {in_diluted:.3f} | {sim_out:.3f} | {exp_out:.3f} | {error:+.3f} |")
        except Exception as e:
            print(f"Simulation failed: {e}")
        print("\n")

if __name__ == "__main__":
    main()
