import sys
import os

sys.path.append(os.path.abspath('.'))

from sx_simulator.multistage_sx import solve_multistage_countercurrent
from sx_simulator.config import DEFAULT_METALS, MOLAR_MASS

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

KNOWN_INVALID_POINTS = {
    ("Data5", "Li"): (
        "원본 검증표의 Li 후액 값(5.400 g/L)이 입력 Li 농도(0.013 g/L)보다 커서 "
        "물질수지상 일관되지 않으므로 해당 지표는 검증 집계에서 제외합니다."
    ),
}


def calc_sulfate_from_feed(feed: dict) -> float:
    """대시보드와 동일한 방식으로 총 황산염 농도를 계산합니다."""
    return (
        feed["Li"] / (2 * MOLAR_MASS["Li"])
        + feed["Ni"] / MOLAR_MASS["Ni"]
        + feed["Co"] / MOLAR_MASS["Co"]
        + feed["Mn"] / MOLAR_MASS["Mn"]
        + feed["Ca"] / MOLAR_MASS["Ca"]
        + feed["Mg"] / MOLAR_MASS["Mg"]
        + feed["Zn"] / MOLAR_MASS["Zn"]
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
