#!/usr/bin/env python3
"""
Jantunen et al. (2022) 논문 파라미터 vs 현장 보정 파라미터 비교
================================================================
Jantunen et al. (2022, Metals 12, 1445)의 D2EHPA/Cyanex 272
pH-추출률 데이터를 기반으로 피팅한 pH50 파라미터와 현재 현장 보정 파라미터를
Data1-6에 적용하여 비교합니다.

논문 조건:
  - D2EHPA: 0.8M + 5vol% TBP → C_ref=0.64M 환산 (TBP pH 이동 효과 포함)
  - Cyanex 272: 0.8M → C_ref=0.5M 환산

Usage:
    python scripts/compare_jantunen_params.py
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sx_simulator.config import (
    LITERATURE_BASE_PARAMS,
    _build_parameter_profile,
    get_parameter_profile,
)
from sx_simulator.datasets import (
    VERIFICATION_EXPERIMENTS,
    KNOWN_INVALID_POINTS,
    prepare_verification_case,
)
from sx_simulator.multistage_sx import solve_multistage_countercurrent

# ─────────────────────────────────────────────────────────────────────
# Jantunen et al. (2022) 논문 기반 pH50 override
# ─────────────────────────────────────────────────────────────────────
# Fig.2b (D2EHPA 0.8M + 5vol% TBP): Mn, Co, Ni, Li pH50 읽기
# Fig.3b (Cyanex 272 0.8M): Mn, Co, Ni pH50 읽기
# Fig.4 (Cyanex 272 selectivity): Li pH50 추정
#
# 논문 원래 농도에서 C_ref 환산 후 pH50 값 (alpha 보정 적용 완료)
# D2EHPA:   C_paper=0.8M → C_ref=0.64M, shift = alpha * ln(0.64/0.8)
# Cyanex 272: C_paper=0.8M → C_ref=0.5M, shift = alpha * ln(0.5/0.8)

JANTUNEN_2022_OVERRIDES = {
    # D2EHPA — 논문 Fig.2b 기반 (TBP 포함 조건 주의)
    ("D2EHPA", "Mn"): {"pH50": 2.07},
    ("D2EHPA", "Co"): {"pH50": 4.04},
    ("D2EHPA", "Ni"): {"pH50": 4.78},
    ("D2EHPA", "Li"): {"pH50": 7.10},
    # Cyanex 272 — 논문 Fig.3b, Fig.4 기반
    ("Cyanex 272", "Mn"): {"pH50": 3.98},
    ("Cyanex 272", "Co"): {"pH50": 4.60},
    ("Cyanex 272", "Ni"): {"pH50": 6.50},
    ("Cyanex 272", "Li"): {"pH50": 7.45},
}

VERIFICATION_BASIS = "legacy_premixed_target_pH"
METALS_OF_INTEREST = ["Li", "Ni", "Co", "Mn"]
REPORT_ALL_METALS = ["Li", "Ni", "Co", "Mn", "Ca", "Mg", "Zn"]


def _dataset_tag(name: str) -> str:
    """'Data1 (CoSX-D9-data, 5 Stages)' → 'Data1'."""
    return name.split()[0]


def _is_invalid_point(dataset_tag: str, metal: str) -> bool:
    return (dataset_tag, metal) in KNOWN_INVALID_POINTS


# ─────────────────────────────────────────────────────────────────────
# 파라미터 비교 테이블 출력
# ─────────────────────────────────────────────────────────────────────
def print_parameter_comparison(field_params: dict, paper_params: dict):
    print("=" * 78)
    print("  파라미터 비교: 현장 보정 vs Jantunen et al. (2022)")
    print("=" * 78)

    for extractant in ["D2EHPA", "Cyanex 272"]:
        print(f"\n  [{extractant}]")
        print(f"  {'Metal':>5} | {'Field pH50':>10} | {'Paper pH50':>10} | {'Delta':>7} | {'Field k':>7} | {'Paper k':>7}")
        print(f"  {'-'*5}-+-{'-'*10}-+-{'-'*10}-+-{'-'*7}-+-{'-'*7}-+-{'-'*7}")
        for metal in METALS_OF_INTEREST:
            fp = field_params[extractant][metal]
            pp = paper_params[extractant][metal]
            delta = pp["pH50"] - fp["pH50"]
            marker = " **" if abs(delta) >= 0.5 else ""
            print(
                f"  {metal:>5} | {fp['pH50']:>10.2f} | {pp['pH50']:>10.2f} | "
                f"{delta:>+7.2f}{marker} | {fp['k']:>7.1f} | {pp['k']:>7.1f}"
            )
    print()


# ─────────────────────────────────────────────────────────────────────
# 시뮬레이션 실행
# ─────────────────────────────────────────────────────────────────────
def run_all_experiments(params: dict, label: str) -> list[dict]:
    """6개 데이터셋에 대해 시뮬레이션을 실행하고 결과를 반환합니다."""
    results = []
    for exp in VERIFICATION_EXPERIMENTS:
        tag = _dataset_tag(exp["name"])
        prep = prepare_verification_case(exp, basis=VERIFICATION_BASIS)
        sim_kwargs = prep["sim_kwargs"]

        try:
            sim = solve_multistage_countercurrent(
                **sim_kwargs, extractant_params=params
            )
            converged = sim.get("converged", False)
            raffinate = sim["raffinate"]
        except Exception as e:
            print(f"  [!] {label} / {tag}: 시뮬레이션 실패 — {e}")
            converged = False
            raffinate = {m: float("nan") for m in REPORT_ALL_METALS}

        results.append({
            "name": exp["name"],
            "tag": tag,
            "extractant": exp["ext"],
            "converged": converged,
            "expected": prep["expected_output"],
            "predicted": raffinate,
        })
    return results


# ─────────────────────────────────────────────────────────────────────
# 상세 결과 테이블
# ─────────────────────────────────────────────────────────────────────
def print_detailed_results(field_results: list, paper_results: list):
    print("=" * 100)
    print("  상세 비교: Dataset x Metal별 후액 농도 (g/L)")
    print("=" * 100)
    print(
        f"  {'Dataset':>6} {'Ext':>6} {'Conv':>4} | "
        f"{'Metal':>5} | {'Actual':>8} | {'Field':>8} {'Err':>7} | "
        f"{'Paper':>8} {'Err':>7} | {'Winner':>6}"
    )
    print(f"  {'-'*6} {'-'*6} {'-'*4}-+-{'-'*5}-+-{'-'*8}-+-{'-'*8}-{'-'*7}-+-{'-'*8}-{'-'*7}-+-{'-'*6}")

    for fr, pr in zip(field_results, paper_results):
        tag = fr["tag"]
        ext = fr["extractant"][:3]
        f_conv = "Y" if fr["converged"] else "N"
        p_conv = "Y" if pr["converged"] else "N"
        conv_str = f"{f_conv}/{p_conv}"

        for metal in REPORT_ALL_METALS:
            actual = fr["expected"].get(metal, 0.0)
            field_pred = fr["predicted"].get(metal, float("nan"))
            paper_pred = pr["predicted"].get(metal, float("nan"))

            if _is_invalid_point(tag, metal):
                print(
                    f"  {tag:>6} {ext:>6} {conv_str:>4} | "
                    f"{metal:>5} | {actual:>8.3f} | {'---':>8} {'(skip)':>7} | "
                    f"{'---':>8} {'(skip)':>7} | {'N/A':>6}"
                )
                continue

            f_err = field_pred - actual
            p_err = paper_pred - actual

            if abs(f_err) < abs(p_err):
                winner = "Field"
            elif abs(p_err) < abs(f_err):
                winner = "Paper"
            else:
                winner = "Tie"

            print(
                f"  {tag:>6} {ext:>6} {conv_str:>4} | "
                f"{metal:>5} | {actual:>8.3f} | {field_pred:>8.3f} {f_err:>+7.3f} | "
                f"{paper_pred:>8.3f} {p_err:>+7.3f} | {winner:>6}"
            )
        print(f"  {'-'*6} {'-'*6} {'-'*4}-+-{'-'*5}-+-{'-'*8}-+-{'-'*8}-{'-'*7}-+-{'-'*8}-{'-'*7}-+-{'-'*6}")


# ─────────────────────────────────────────────────────────────────────
# MAE 집계
# ─────────────────────────────────────────────────────────────────────
def compute_mae_summary(
    field_results: list, paper_results: list
) -> dict:
    """금속별, 데이터셋별, 추출제별, 전체 MAE를 계산합니다."""
    # 수집 구조: {key: [(field_err, paper_err), ...]}
    by_metal: dict[str, list] = {m: [] for m in REPORT_ALL_METALS}
    by_dataset: dict[str, list] = {}
    by_extractant: dict[str, list] = {}
    all_errors: list = []

    for fr, pr in zip(field_results, paper_results):
        tag = fr["tag"]
        ext = fr["extractant"]
        if tag not in by_dataset:
            by_dataset[tag] = []
        if ext not in by_extractant:
            by_extractant[ext] = []

        for metal in REPORT_ALL_METALS:
            if _is_invalid_point(tag, metal):
                continue
            actual = fr["expected"].get(metal, 0.0)
            f_pred = fr["predicted"].get(metal, float("nan"))
            p_pred = pr["predicted"].get(metal, float("nan"))

            f_ae = abs(f_pred - actual)
            p_ae = abs(p_pred - actual)
            pair = (f_ae, p_ae)

            by_metal[metal].append(pair)
            by_dataset[tag].append(pair)
            by_extractant[ext].append(pair)
            all_errors.append(pair)

    def _mae(pairs: list) -> tuple[float, float]:
        if not pairs:
            return (float("nan"), float("nan"))
        f_sum = sum(p[0] for p in pairs)
        p_sum = sum(p[1] for p in pairs)
        n = len(pairs)
        return (f_sum / n, p_sum / n)

    return {
        "by_metal": {m: _mae(by_metal[m]) for m in REPORT_ALL_METALS},
        "by_dataset": {d: _mae(by_dataset[d]) for d in sorted(by_dataset)},
        "by_extractant": {e: _mae(by_extractant[e]) for e in sorted(by_extractant)},
        "overall": _mae(all_errors),
        "n_total": len(all_errors),
    }


def print_mae_summary(mae: dict):
    print("\n" + "=" * 68)
    print("  MAE 요약: 현장 보정 vs Jantunen et al. (2022)")
    print("=" * 68)

    # 금속별
    print(f"\n  [금속별 MAE (g/L)]")
    print(f"  {'Metal':>5} | {'Field MAE':>10} | {'Paper MAE':>10} | {'Winner':>6} | N")
    print(f"  {'-'*5}-+-{'-'*10}-+-{'-'*10}-+-{'-'*6}-+-{'-'*3}")
    for metal in REPORT_ALL_METALS:
        f_mae, p_mae = mae["by_metal"][metal]
        n = len([1 for fr, pr in zip(field_results, paper_results)
                 for m in [metal] if not _is_invalid_point(fr["tag"], m)])
        winner = "Field" if f_mae < p_mae else ("Paper" if p_mae < f_mae else "Tie")
        print(f"  {metal:>5} | {f_mae:>10.4f} | {p_mae:>10.4f} | {winner:>6} | {n}")

    # 데이터셋별
    print(f"\n  [데이터셋별 MAE (g/L)]")
    print(f"  {'Dataset':>7} | {'Field MAE':>10} | {'Paper MAE':>10} | {'Winner':>6} | N")
    print(f"  {'-'*7}-+-{'-'*10}-+-{'-'*10}-+-{'-'*6}-+-{'-'*3}")
    for dataset in sorted(mae["by_dataset"]):
        f_mae, p_mae = mae["by_dataset"][dataset]
        n = len(mae["by_dataset"][dataset]) if isinstance(mae["by_dataset"][dataset], tuple) else 0
        # re-count from errors
        winner = "Field" if f_mae < p_mae else ("Paper" if p_mae < f_mae else "Tie")
        print(f"  {dataset:>7} | {f_mae:>10.4f} | {p_mae:>10.4f} | {winner:>6}")

    # 추출제 그룹별
    print(f"\n  [추출제별 MAE (g/L)]")
    print(f"  {'Extractant':>12} | {'Field MAE':>10} | {'Paper MAE':>10} | {'Winner':>6}")
    print(f"  {'-'*12}-+-{'-'*10}-+-{'-'*10}-+-{'-'*6}")
    for ext in sorted(mae["by_extractant"]):
        f_mae, p_mae = mae["by_extractant"][ext]
        winner = "Field" if f_mae < p_mae else ("Paper" if p_mae < f_mae else "Tie")
        print(f"  {ext:>12} | {f_mae:>10.4f} | {p_mae:>10.4f} | {winner:>6}")

    # 전체
    f_all, p_all = mae["overall"]
    winner = "Field" if f_all < p_all else ("Paper" if p_all < f_all else "Tie")
    print(f"\n  전체 MAE: Field={f_all:.4f}, Paper={p_all:.4f}  →  {winner} ({mae['n_total']}개 비교점)")
    print("=" * 68)


# ─────────────────────────────────────────────────────────────────────
# 수렴 상태 요약
# ─────────────────────────────────────────────────────────────────────
def print_convergence_summary(field_results: list, paper_results: list):
    print("\n  [수렴 상태]")
    print(f"  {'Dataset':>7} | {'Field':>6} | {'Paper':>6}")
    print(f"  {'-'*7}-+-{'-'*6}-+-{'-'*6}")
    for fr, pr in zip(field_results, paper_results):
        f_status = "OK" if fr["converged"] else "FAIL"
        p_status = "OK" if pr["converged"] else "FAIL"
        print(f"  {fr['tag']:>7} | {f_status:>6} | {p_status:>6}")
    f_total = sum(1 for r in field_results if r["converged"])
    p_total = sum(1 for r in paper_results if r["converged"])
    print(f"  {'합계':>7} | {f_total}/{len(field_results)} | {p_total}/{len(paper_results)}")


# ─────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("\n" + "=" * 78)
    print("  Jantunen et al. (2022) 논문 파라미터 vs 현장 보정 파라미터 비교")
    print("  Basis: legacy_premixed_target_pH")
    print("=" * 78)

    # 파라미터 프로필 빌드
    field_params = get_parameter_profile("field_calibrated")
    paper_params = _build_parameter_profile(LITERATURE_BASE_PARAMS, JANTUNEN_2022_OVERRIDES)

    print_parameter_comparison(field_params, paper_params)

    # 시뮬레이션 실행
    print("\n  현장 보정 모델 실행 중...")
    field_results = run_all_experiments(field_params, "Field")
    print("  논문 기반 모델 실행 중...")
    paper_results = run_all_experiments(paper_params, "Paper")

    # 결과 출력
    print_convergence_summary(field_results, paper_results)
    print_detailed_results(field_results, paper_results)

    mae = compute_mae_summary(field_results, paper_results)
    print_mae_summary(mae)
