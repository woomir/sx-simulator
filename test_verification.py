import sys
import os

sys.path.append(os.path.abspath('.'))

from sx_simulator.multistage_sx import solve_multistage_countercurrent
from sx_simulator.config import DEFAULT_METALS
from sx_simulator.datasets import (
    VERIFICATION_EXPERIMENTS as experiments,
    VALIDATION_BASES,
    KNOWN_INVALID_POINTS,
    prepare_verification_case,
)

def run_validation_basis(basis: str):
    print(f"=== Verification Basis: {basis} ===")
    overall_errors = {metal: [] for metal in DEFAULT_METALS}

    for exp in experiments:
        prepared = prepare_verification_case(exp, basis=basis)
        sim_kwargs = {
            **prepared["sim_kwargs"],
            "metals": DEFAULT_METALS,
            "use_competition": True,
            "use_speciation": True,
        }

        print(f"--- Running {prepared['name']} ---")
        try:
            result = solve_multistage_countercurrent(**sim_kwargs)
            res_df = result['raffinate']
            abs_errors = []
            worst_row = None

            q_aq = sim_kwargs["Q_aq"]
            q_org = sim_kwargs["Q_org"]
            print(
                f"Dilution Factor: {prepared['dilution_factor']:.4f}, "
                f"O/A Ratio: {q_org/q_aq:.4f}, Sulfate: {prepared['sulfate_m']:.4f} M"
            )
            if basis == "raw_feed_target_pH":
                estimated_naoh_flow = sum(result.get("NaOH_flow_profile", []))
                print(
                    f"Assumed NaOH: {prepared['assumed_naoh_concentration_m']:.2f} M, "
                    f"Listed NaOH Flow: {prepared['listed_naoh_flow_l_hr']:.3f} L/hr, "
                    f"Estimated NaOH Flow: {estimated_naoh_flow:.3f} L/hr"
                )
            print(f"| Metal | In_Original | In_ModelBasis | Sim_Raffinate | Exp_Output | Sim_Error |")
            print(f"|---|---|---|---|---|---|")
            dataset_tag = prepared["dataset_tag"]
            for m in DEFAULT_METALS:
                in_orig = prepared['input_original'].get(m, 0.0)
                in_model_basis = prepared['input_model_basis'].get(m, 0.0)
                sim_out = res_df[m]
                exp_out = prepared['expected_output'].get(m, 0.0)
                error = sim_out - exp_out

                exclusion_reason = KNOWN_INVALID_POINTS.get((dataset_tag, m))
                if exclusion_reason:
                    print(f"| {m} | {in_orig:.3f} | {in_model_basis:.3f} | {sim_out:.3f} | {exp_out:.3f} | 제외 |")
                    print(f"  Note: {exclusion_reason}")
                    continue

                abs_errors.append(abs(error))
                overall_errors[m].append(abs(error))
                if worst_row is None or abs(error) > abs(worst_row[1]):
                    worst_row = (m, error)
                print(f"| {m} | {in_orig:.3f} | {in_model_basis:.3f} | {sim_out:.3f} | {exp_out:.3f} | {error:+.3f} |")
            mae = sum(abs_errors) / len(abs_errors) if abs_errors else 0.0
            if worst_row is not None:
                print(f"Summary: converged={result['converged']} iterations={result['iterations']} dataset_MAE={mae:.3f} g/L worst={worst_row[0]} {worst_row[1]:+.3f} g/L")
        except Exception as e:
            print(f"Simulation failed: {e}")
        print("\n")

    print(f"--- Overall Metal MAE Summary ({basis}) ---")
    for metal in DEFAULT_METALS:
        if overall_errors[metal]:
            metal_mae = sum(overall_errors[metal]) / len(overall_errors[metal])
            print(f"{metal}: MAE={metal_mae:.3f} g/L over {len(overall_errors[metal])} datasets")
    print("")

    return overall_errors


def main():
    for basis in VALIDATION_BASES:
        run_validation_basis(basis)

if __name__ == "__main__":
    main()
