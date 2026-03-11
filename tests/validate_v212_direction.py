import os
import sys
import math
import random
import statistics
import time
from collections import defaultdict

sys.path.append(os.path.abspath("."))

from sx_simulator.config import DEFAULT_METALS, get_parameter_profile
from sx_simulator.datasets import (
    FIELD_DATASETS,
    KNOWN_INVALID_POINTS,
    VERIFICATION_EXPERIMENTS,
    calc_sulfate_from_feed,
    prepare_verification_case,
)
from sx_simulator.multistage_sx import solve_multistage_countercurrent


FIELD_BASES = (
    "legacy_premixed_target_pH",
    "raw_feed_target_pH",
    "raw_feed_fixed_saponification",
    "raw_feed_physical_saponification_v2",
)

RANDOM_MODES = (
    "aqueous_direct_fixed",
    "aqueous_direct_target",
    "fresh_saponification_fixed",
    "fresh_saponification_target_physical",
    "fresh_saponification_target_legacy",
)

TRACE_EXPECTED_FLOOR_G_L = 0.1


def percentile(values: list[float], ratio: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    index = max(0, min(len(ordered) - 1, math.ceil(len(ordered) * ratio) - 1))
    return ordered[index]


def mean_or_zero(values: list[float]) -> float:
    return statistics.fmean(values) if values else 0.0


def mean_or_none(values: list[float]) -> float | None:
    return statistics.fmean(values) if values else None


def relative_error_pct(simulated: float, experimental: float) -> float | None:
    if abs(experimental) < TRACE_EXPECTED_FLOOR_G_L:
        return None
    return abs(simulated - experimental) / abs(experimental) * 100.0


def build_extractant_ranges(extractant: str) -> dict:
    matching = [case for case in FIELD_DATASETS.values() if case["ext"] == extractant]
    if not matching:
        raise ValueError(f"No field dataset found for extractant: {extractant}")

    def numeric_range(values: list[float], expand_low: float = 0.9, expand_high: float = 1.1) -> tuple[float, float]:
        low = min(values)
        high = max(values)
        if math.isclose(low, high):
            return max(0.0, low * expand_low), high * expand_high if high > 0.0 else 0.0
        return max(0.0, low * expand_low), high * expand_high

    feed_ranges = {}
    for metal in DEFAULT_METALS:
        values = [case["input"].get(metal, 0.0) for case in matching]
        if max(values) <= 0.0:
            feed_ranges[metal] = (0.0, 0.0)
        else:
            feed_ranges[metal] = numeric_range(values, 0.8, 1.2)

    pH_low, pH_high = numeric_range([case["pH"] for case in matching], 0.95, 1.10)
    feed_pH_low, feed_pH_high = numeric_range([case.get("pH_feed", case["pH"]) for case in matching], 0.9, 1.05)
    q_aq_low, q_aq_high = numeric_range([case["feed_flow"] for case in matching], 0.8, 1.2)
    q_org_low, q_org_high = numeric_range([case["org_flow"] for case in matching], 0.8, 1.2)
    c_ext_low, c_ext_high = numeric_range([case["C_ext"] for case in matching], 0.9, 1.1)
    naoh_flow_low, naoh_flow_high = numeric_range([case["naoh_flow"] for case in matching], 0.6, 1.4)
    naoh_m_low, naoh_m_high = numeric_range(
        [prepare_verification_case({"name": "tmp", **case}, "raw_feed_fixed_saponification")["assumed_naoh_concentration_m"] for case in matching],
        0.8,
        1.2,
    )

    return {
        "feed": feed_ranges,
        "target_pH": (max(2.0, pH_low - 0.2), min(8.0, pH_high + 0.6)),
        "feed_pH": (max(0.5, feed_pH_low - 0.2), min(7.0, feed_pH_high + 0.2)),
        "Q_aq": (max(5.0, q_aq_low), max(10.0, q_aq_high)),
        "Q_org": (max(5.0, q_org_low), max(10.0, q_org_high)),
        "C_ext": (max(0.05, c_ext_low), max(0.1, c_ext_high)),
        "Q_NaOH": (max(0.2, naoh_flow_low), max(0.5, naoh_flow_high)),
        "C_NaOH": (max(0.5, naoh_m_low), max(1.0, naoh_m_high)),
    }


def sample_case(extractant: str, rng: random.Random, idx: int) -> dict:
    ranges = build_extractant_ranges(extractant)
    c_aq_feed = {}
    for metal, (low, high) in ranges["feed"].items():
        c_aq_feed[metal] = 0.0 if high <= 0.0 else rng.uniform(low, high)

    q_aq = rng.uniform(*ranges["Q_aq"])
    q_org = rng.uniform(*ranges["Q_org"])
    case = {
        "case_id": f"{extractant}-{idx:02d}",
        "extractant": extractant,
        "C_aq_feed": c_aq_feed,
        "pH_feed": rng.uniform(*ranges["feed_pH"]),
        "Q_aq": q_aq,
        "Q_org": q_org,
        "C_ext": rng.uniform(*ranges["C_ext"]),
        "n_stages": rng.randint(2, 6),
        "target_pH": rng.uniform(*ranges["target_pH"]),
        "C_NaOH": rng.uniform(*ranges["C_NaOH"]),
        "Q_NaOH": rng.uniform(*ranges["Q_NaOH"]),
        "temperature": 25.0,
    }
    case["C_sulfate"] = calc_sulfate_from_feed(case["C_aq_feed"])
    return case


def run_case(case: dict, mode: str, extractant_params: dict) -> dict:
    kwargs = {
        "C_aq_feed": case["C_aq_feed"],
        "pH_feed": case["pH_feed"],
        "Q_aq": case["Q_aq"],
        "Q_org": case["Q_org"],
        "extractant": case["extractant"],
        "C_ext": case["C_ext"],
        "n_stages": case["n_stages"],
        "temperature": case["temperature"],
        "C_sulfate": case["C_sulfate"],
        "metals": DEFAULT_METALS,
        "use_competition": True,
        "use_speciation": True,
        "extractant_params": extractant_params,
    }

    if mode == "aqueous_direct_fixed":
        kwargs.update(
            {
                "naoh_mode": "aqueous_direct",
                "C_NaOH": case["C_NaOH"],
                "Q_NaOH": case["Q_NaOH"],
            }
        )
    elif mode == "aqueous_direct_target":
        kwargs.update(
            {
                "naoh_mode": "aqueous_direct",
                "target_pH": case["target_pH"],
                "C_NaOH": case["C_NaOH"],
            }
        )
    elif mode == "fresh_saponification_fixed":
        kwargs.update(
            {
                "naoh_mode": "saponification",
                "C_NaOH": case["C_NaOH"],
                "Q_NaOH": case["Q_NaOH"],
                "saponification_model": "physical_v2",
            }
        )
    elif mode == "fresh_saponification_target_physical":
        kwargs.update(
            {
                "naoh_mode": "saponification",
                "target_pH": case["target_pH"],
                "C_NaOH": case["C_NaOH"],
                "Q_NaOH": case["Q_NaOH"],
                "saponification_model": "physical_v2",
            }
        )
    elif mode == "fresh_saponification_target_legacy":
        kwargs.update(
            {
                "naoh_mode": "saponification",
                "target_pH": case["target_pH"],
                "C_NaOH": case["C_NaOH"],
                "Q_NaOH": case["Q_NaOH"],
                "saponification_model": "legacy_equivalent_target",
            }
        )
    else:
        raise ValueError(f"Unknown mode: {mode}")

    started_at = time.perf_counter()
    result = solve_multistage_countercurrent(**kwargs)
    runtime_s = time.perf_counter() - started_at

    q_aq_out = result["stages"][-1].get("Q_aq_out_L_hr", case["Q_aq"])
    max_balance_error_g = 0.0
    for metal in DEFAULT_METALS:
        mass_in = case["C_aq_feed"].get(metal, 0.0) * case["Q_aq"]
        mass_out = result["raffinate"].get(metal, 0.0) * q_aq_out + result["loaded_organic"].get(metal, 0.0) * case["Q_org"]
        max_balance_error_g = max(max_balance_error_g, abs(mass_in - mass_out))

    has_negative = False
    has_non_finite = False
    for stage in result["stages"]:
        for bucket in ("C_aq_out", "C_org_out"):
            for value in stage.get(bucket, {}).values():
                if not math.isfinite(value):
                    has_non_finite = True
                if value < -1e-9:
                    has_negative = True

    return {
        "result": result,
        "runtime_s": runtime_s,
        "max_balance_error_g": max_balance_error_g,
        "has_negative": has_negative,
        "has_non_finite": has_non_finite,
    }


def evaluate_field_by_extractant(extractant: str, extractant_params: dict) -> dict:
    matching = [case for case in VERIFICATION_EXPERIMENTS if case["ext"] == extractant]
    summary = {}

    for basis in FIELD_BASES:
        dataset_metrics = []
        point_abs_errors = defaultdict(list)
        point_relative_errors = defaultdict(list)
        aggregate_abs_errors = []
        aggregate_relative_errors = []
        for case in matching:
            prepared = prepare_verification_case(case, basis=basis)
            sim_kwargs = {
                **prepared["sim_kwargs"],
                "metals": DEFAULT_METALS,
                "use_competition": True,
                "use_speciation": True,
                "extractant_params": extractant_params,
            }
            result = solve_multistage_countercurrent(**sim_kwargs)
            dataset_tag = prepared["dataset_tag"]
            abs_errors = []
            relative_errors = []
            for metal in DEFAULT_METALS:
                if (dataset_tag, metal) in KNOWN_INVALID_POINTS:
                    continue
                sim_out = result["raffinate"].get(metal, 0.0)
                exp_out = prepared["expected_output"].get(metal, 0.0)
                error = sim_out - exp_out
                abs_error = abs(error)
                abs_errors.append(abs_error)
                aggregate_abs_errors.append(abs_error)
                point_abs_errors[metal].append(abs_error)

                rel_error = relative_error_pct(sim_out, exp_out)
                if rel_error is not None:
                    relative_errors.append(rel_error)
                    aggregate_relative_errors.append(rel_error)
                    point_relative_errors[metal].append(rel_error)

            dataset_metrics.append(
                {
                    "dataset": dataset_tag,
                    "mae_g_l": mean_or_zero(abs_errors),
                    "mean_relative_error_pct": mean_or_zero(relative_errors),
                    "converged": bool(result.get("converged", False)),
                    "iterations": int(result.get("iterations", 0)),
                }
            )

        summary[basis] = {
            "dataset_metrics": dataset_metrics,
            "aggregate_mae_g_l": mean_or_zero(aggregate_abs_errors),
            "aggregate_relative_error_pct": mean_or_zero(aggregate_relative_errors),
            "li_mae_g_l": mean_or_zero(point_abs_errors["Li"]),
            "ni_mae_g_l": mean_or_zero(point_abs_errors["Ni"]),
            "co_mae_g_l": mean_or_zero(point_abs_errors["Co"]),
            "li_relative_error_pct": mean_or_none(point_relative_errors["Li"]),
            "ni_relative_error_pct": mean_or_none(point_relative_errors["Ni"]),
            "co_relative_error_pct": mean_or_none(point_relative_errors["Co"]),
        }

    return summary


def evaluate_random_modes_by_extractant(
    extractant: str,
    extractant_params: dict,
    seed: int,
    cases_per_extractant: int = 12,
) -> dict:
    rng = random.Random(seed)
    cases = [sample_case(extractant, rng, idx + 1) for idx in range(cases_per_extractant)]
    mode_summary = {}

    for mode in RANDOM_MODES:
        runs = []
        for case in cases:
            try:
                outcome = run_case(case, mode, extractant_params)
                result = outcome["result"]
                runs.append(
                    {
                        "converged": bool(result.get("converged", False)),
                        "iterations": int(result.get("iterations", 0)),
                        "runtime_s": outcome["runtime_s"],
                        "max_balance_error_g": outcome["max_balance_error_g"],
                        "has_negative": outcome["has_negative"],
                        "has_non_finite": outcome["has_non_finite"],
                    }
                )
            except Exception:
                runs.append(
                    {
                        "converged": False,
                        "iterations": 0,
                        "runtime_s": 0.0,
                        "max_balance_error_g": float("inf"),
                        "has_negative": True,
                        "has_non_finite": True,
                    }
                )

        converged_count = sum(1 for run in runs if run["converged"])
        mode_summary[mode] = {
            "total_cases": len(runs),
            "converged_cases": converged_count,
            "convergence_rate": converged_count / len(runs) if runs else 0.0,
            "negative_cases": sum(1 for run in runs if run["has_negative"]),
            "non_finite_cases": sum(1 for run in runs if run["has_non_finite"]),
            "balance_fail_cases": sum(1 for run in runs if run["max_balance_error_g"] > 1e-3),
            "p95_iterations": percentile([run["iterations"] for run in runs], 0.95),
            "p95_runtime_s": percentile([run["runtime_s"] for run in runs], 0.95),
            "max_balance_error_g": max(run["max_balance_error_g"] for run in runs),
        }

    return mode_summary


def evaluate_target_saponification_q_sensitivity(
    extractant: str,
    extractant_params: dict,
    seed: int,
    cases_per_extractant: int = 12,
) -> dict:
    rng = random.Random(seed)
    cases = [sample_case(extractant, rng, idx + 1) for idx in range(cases_per_extractant)]
    invariant_cases = 0
    max_raffinate_diff = 0.0
    max_pH_profile_diff = 0.0

    for case in cases:
        low_q_case = dict(case)
        high_q_case = dict(case)
        low_q_case["Q_NaOH"] = 0.5
        high_q_case["Q_NaOH"] = 50.0

        low_result = run_case(low_q_case, "fresh_saponification_target_legacy", extractant_params)["result"]
        high_result = run_case(high_q_case, "fresh_saponification_target_legacy", extractant_params)["result"]

        case_raffinate_diff = max(
            abs(low_result["raffinate"].get(metal, 0.0) - high_result["raffinate"].get(metal, 0.0))
            for metal in DEFAULT_METALS
        )
        case_pH_diff = max(
            abs(a - b)
            for a, b in zip(low_result["pH_profile"], high_result["pH_profile"])
        )

        max_raffinate_diff = max(max_raffinate_diff, case_raffinate_diff)
        max_pH_profile_diff = max(max_pH_profile_diff, case_pH_diff)
        if case_raffinate_diff < 1e-12 and case_pH_diff < 1e-12:
            invariant_cases += 1

    return {
        "total_cases": len(cases),
        "invariant_cases": invariant_cases,
        "invariant_ratio": invariant_cases / len(cases) if cases else 0.0,
        "max_raffinate_diff_gL": max_raffinate_diff,
        "max_pH_profile_diff": max_pH_profile_diff,
    }


def print_field_summary(extractant: str, summary: dict) -> None:
    print(f"\n=== FIELD SUMMARY: {extractant} ===")
    for basis in FIELD_BASES:
        basis_summary = summary[basis]
        li_rel_text = (
            f"{basis_summary['li_relative_error_pct']:.1f}%"
            if basis_summary["li_relative_error_pct"] is not None
            else "n/a(trace)"
        )
        ni_rel_text = (
            f"{basis_summary['ni_relative_error_pct']:.1f}%"
            if basis_summary["ni_relative_error_pct"] is not None
            else "n/a(trace)"
        )
        co_rel_text = (
            f"{basis_summary['co_relative_error_pct']:.1f}%"
            if basis_summary["co_relative_error_pct"] is not None
            else "n/a(trace)"
        )
        print(
            f"- {basis}: aggregate_rel_error={basis_summary['aggregate_relative_error_pct']:.1f}%, "
            f"aggregate_MAE={basis_summary['aggregate_mae_g_l']:.3f} g/L, "
            f"Li_rel={li_rel_text}, Ni_rel={ni_rel_text}, Co_rel={co_rel_text}"
        )
        for dataset in basis_summary["dataset_metrics"]:
            print(
                f"  - {dataset['dataset']}: mean_rel_error={dataset['mean_relative_error_pct']:.1f}%, "
                f"MAE={dataset['mae_g_l']:.3f} g/L, converged={dataset['converged']}, "
                f"iterations={dataset['iterations']}"
            )


def print_random_summary(extractant: str, summary: dict) -> None:
    print(f"\n=== RANDOM SUMMARY: {extractant} ===")
    for mode in RANDOM_MODES:
        mode_summary = summary[mode]
        print(
            f"- {mode}: convergence={mode_summary['converged_cases']}/{mode_summary['total_cases']} "
            f"({mode_summary['convergence_rate']:.0%}), balance_fail={mode_summary['balance_fail_cases']}, "
            f"negative={mode_summary['negative_cases']}, non_finite={mode_summary['non_finite_cases']}, "
            f"p95_iter={mode_summary['p95_iterations']:.0f}, p95_runtime={mode_summary['p95_runtime_s']:.3f}s, "
            f"max_balance_error={mode_summary['max_balance_error_g']:.6f} g"
        )


def print_sensitivity_summary(extractant: str, summary: dict) -> None:
    print(f"\n=== TARGET-PH + SAPONIFICATION Q_NaOH SENSITIVITY: {extractant} ===")
    print(
        f"- invariant_cases={summary['invariant_cases']}/{summary['total_cases']} "
        f"({summary['invariant_ratio']:.0%}), "
        f"max_raffinate_diff={summary['max_raffinate_diff_gL']:.12f} g/L, "
        f"max_pH_profile_diff={summary['max_pH_profile_diff']:.12f}"
    )


def main() -> None:
    extractant_params = get_parameter_profile("field_calibrated")
    extractants = ("Cyanex 272", "D2EHPA")

    for idx, extractant in enumerate(extractants, start=1):
        field_summary = evaluate_field_by_extractant(extractant, extractant_params)
        random_summary = evaluate_random_modes_by_extractant(
            extractant,
            extractant_params,
            seed=20260311 + idx * 100,
        )
        sensitivity_summary = evaluate_target_saponification_q_sensitivity(
            extractant,
            extractant_params,
            seed=20260311 + idx * 1000,
        )

        print_field_summary(extractant, field_summary)
        print_random_summary(extractant, random_summary)
        print_sensitivity_summary(extractant, sensitivity_summary)


if __name__ == "__main__":
    main()
