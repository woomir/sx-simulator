import sys
import os

sys.path.append(os.path.abspath('.'))

from sx_simulator.multistage_sx import solve_multistage_countercurrent
from sx_simulator.config import DEFAULT_METALS
from sx_simulator.datasets import (
    VERIFICATION_EXPERIMENTS as experiments,
    VALIDATION_BASES,
    KNOWN_INVALID_POINTS,
    LEGACY_REGRESSION_ABS_ERROR_THRESHOLDS,
    LEGACY_REGRESSION_DATASET_MAE_THRESHOLDS,
    LEGACY_REGRESSION_MAE_THRESHOLDS,
    prepare_verification_case,
)


def append_unique(items: list[str], message: str) -> None:
    if message not in items:
        items.append(message)


def run_validation_basis(basis: str):
    print(f"=== Verification Basis: {basis} ===")
    overall_errors = {metal: [] for metal in DEFAULT_METALS}
    dataset_metrics = []
    failures = []

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
            metal_errors = {}
            estimated_naoh_flow = None

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
                metal_errors[m] = error
                if worst_row is None or abs(error) > abs(worst_row[1]):
                    worst_row = (m, error)
                print(f"| {m} | {in_orig:.3f} | {in_model_basis:.3f} | {sim_out:.3f} | {exp_out:.3f} | {error:+.3f} |")
            mae = sum(abs_errors) / len(abs_errors) if abs_errors else 0.0
            converged = result.get("converged", False)
            iterations = result.get("iterations")
            if worst_row is not None:
                print(f"Summary: converged={converged} iterations={iterations} dataset_MAE={mae:.3f} g/L worst={worst_row[0]} {worst_row[1]:+.3f} g/L")

            if not converged:
                append_unique(failures, f"{dataset_tag}: simulation did not converge")

            dataset_metrics.append(
                {
                    "name": prepared["name"],
                    "dataset_tag": dataset_tag,
                    "converged": converged,
                    "iterations": iterations,
                    "mae": mae,
                    "worst": worst_row,
                    "errors": metal_errors,
                    "estimated_naoh_flow_l_hr": estimated_naoh_flow,
                    "listed_naoh_flow_l_hr": prepared.get("listed_naoh_flow_l_hr"),
                    "exception": None,
                }
            )
        except Exception as e:
            print(f"Simulation failed: {e}")
            append_unique(failures, f"{prepared['dataset_tag']}: simulation raised {e}")
            dataset_metrics.append(
                {
                    "name": prepared["name"],
                    "dataset_tag": prepared["dataset_tag"],
                    "converged": False,
                    "iterations": None,
                    "mae": None,
                    "worst": None,
                    "errors": {},
                    "estimated_naoh_flow_l_hr": None,
                    "listed_naoh_flow_l_hr": prepared.get("listed_naoh_flow_l_hr"),
                    "exception": str(e),
                }
            )
        print("\n")

    print(f"--- Overall Metal MAE Summary ({basis}) ---")
    overall_mae = {}
    for metal in DEFAULT_METALS:
        if overall_errors[metal]:
            metal_mae = sum(overall_errors[metal]) / len(overall_errors[metal])
            overall_mae[metal] = metal_mae
            print(f"{metal}: MAE={metal_mae:.3f} g/L over {len(overall_errors[metal])} datasets")
    print("")

    return {
        "basis": basis,
        "overall_errors": overall_errors,
        "overall_mae": overall_mae,
        "dataset_metrics": dataset_metrics,
        "failures": failures,
    }


def evaluate_legacy_regression(metrics: dict) -> list[str]:
    failures = list(metrics["failures"])

    for dataset in metrics["dataset_metrics"]:
        dataset_tag = dataset["dataset_tag"]
        exception = dataset.get("exception")
        if exception:
            append_unique(failures, f"{dataset_tag}: simulation raised {exception}")
            continue

        if not dataset.get("converged", False):
            append_unique(failures, f"{dataset_tag}: simulation did not converge")

        dataset_mae = dataset.get("mae")
        dataset_threshold = LEGACY_REGRESSION_DATASET_MAE_THRESHOLDS.get(dataset_tag)
        if dataset_threshold is not None and dataset_mae is not None and dataset_mae > dataset_threshold:
            append_unique(
                failures,
                f"{dataset_tag}: dataset MAE {dataset_mae:.3f} g/L exceeds threshold {dataset_threshold:.3f} g/L",
            )

        for (expected_dataset, metal), abs_threshold in LEGACY_REGRESSION_ABS_ERROR_THRESHOLDS.items():
            if expected_dataset != dataset_tag:
                continue
            error = dataset["errors"].get(metal)
            if error is None:
                append_unique(
                    failures,
                    f"{dataset_tag} {metal}: regression point missing from evaluated errors",
                )
                continue
            abs_error = abs(error)
            if abs_error > abs_threshold:
                append_unique(
                    failures,
                    f"{dataset_tag} {metal}: abs error {abs_error:.3f} g/L exceeds threshold {abs_threshold:.3f} g/L",
                )

    for metal, mae in metrics["overall_mae"].items():
        mae_threshold = LEGACY_REGRESSION_MAE_THRESHOLDS.get(metal)
        if mae_threshold is not None and mae > mae_threshold:
            append_unique(
                failures,
                f"{metal}: overall MAE {mae:.3f} g/L exceeds threshold {mae_threshold:.3f} g/L",
            )

    return failures


def main():
    results = {}
    for basis in VALIDATION_BASES:
        results[basis] = run_validation_basis(basis)

    failures = evaluate_legacy_regression(results["legacy_premixed_target_pH"])
    if failures:
        print("=== LEGACY REGRESSION STATUS: FAILED ===")
        for failure in failures:
            print(f"- {failure}")
        print("Raw-feed basis is diagnostic only and does not affect pass/fail yet.")
        raise SystemExit(1)

    print("=== LEGACY REGRESSION STATUS: PASSED ===")
    print("Raw-feed basis is diagnostic only and does not affect pass/fail yet.")

if __name__ == "__main__":
    main()
