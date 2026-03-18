#!/usr/bin/env python3
"""
브랜치 간 SX 시뮬레이터 정확도 비교
=====================================
main (4PL Sigmoid) vs colleague/ahh-modifications (5PL Richards)

Git worktree + sys.path 스왑으로 양쪽 시뮬레이터를 순차 로드하여
동일 현장 데이터에 대한 예측 정확도를 비교합니다.
"""

import sys
import os
import re
import csv
import copy
import math
import subprocess
import importlib
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = Path(__file__).resolve().parent.parent
WORKTREE_DIR = PROJECT_ROOT / ".worktrees" / "colleague"
DATA_DIR = PROJECT_ROOT / "data"
REPORTS_DIR = PROJECT_ROOT / "reports"

DEFAULT_METALS = ["Li", "Ni", "Co", "Mn", "Ca", "Mg", "Zn"]


# ═════════════════════════════════════════════════════════════════════
# Part A: Worktree 관리 + 모듈 로더
# ═════════════════════════════════════════════════════════════════════

def setup_worktree() -> Path:
    """colleague 브랜치를 .worktrees/colleague 에 체크아웃."""
    if WORKTREE_DIR.exists():
        print(f"  [INFO] Worktree 이미 존재: {WORKTREE_DIR}")
        return WORKTREE_DIR

    WORKTREE_DIR.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "git", "worktree", "add",
        str(WORKTREE_DIR),
        "colleague/ahh-modifications",
    ]
    print(f"  [INFO] Worktree 생성: {' '.join(cmd)}")
    subprocess.run(cmd, cwd=str(PROJECT_ROOT), check=True, capture_output=True)
    return WORKTREE_DIR


def _purge_sx_modules():
    """sys.modules에서 sx_simulator.* 모듈을 전부 제거."""
    to_remove = [key for key in sys.modules if key.startswith("sx_simulator")]
    for key in to_remove:
        del sys.modules[key]


def load_simulator(branch_root: Path, branch_name: str) -> dict:
    """지정된 브랜치 루트에서 시뮬레이터를 로드.

    Returns:
        {
            "solve": solve_multistage_countercurrent,
            "estimate_naoh_m": estimate_naoh_molarity_from_wt_pct,
            "calc_sulfate": calc_sulfate_from_feed,
            "has_use_competition": bool,
            "branch_name": str,
        }
    """
    _purge_sx_modules()

    branch_str = str(branch_root)
    sys.path.insert(0, branch_str)

    try:
        from sx_simulator.multistage_sx import solve_multistage_countercurrent
        from sx_simulator.datasets import (
            estimate_naoh_molarity_from_wt_pct,
            calc_sulfate_from_feed,
        )

        # use_competition 파라미터 존재 여부 확인
        import inspect
        sig = inspect.signature(solve_multistage_countercurrent)
        has_uc = "use_competition" in sig.parameters

        print(f"  [OK] {branch_name} 시뮬레이터 로드 완료 (use_competition={has_uc})")

        return {
            "solve": solve_multistage_countercurrent,
            "estimate_naoh_m": estimate_naoh_molarity_from_wt_pct,
            "calc_sulfate": calc_sulfate_from_feed,
            "has_use_competition": has_uc,
            "branch_name": branch_name,
        }
    finally:
        sys.path.remove(branch_str)


def run_sim(sim_bundle: dict, **kwargs):
    """use_competition 호환성 래퍼."""
    if sim_bundle["has_use_competition"]:
        kwargs.setdefault("use_competition", True)
    else:
        kwargs.pop("use_competition", None)
    return sim_bundle["solve"](**kwargs)


# ═════════════════════════════════════════════════════════════════════
# Part B: CSV 파싱 (field_validation_main.py 기반)
# ═════════════════════════════════════════════════════════════════════

def _parse_numeric(s: str) -> float:
    """콤마 구분 숫자, "<1" → 0.0 처리."""
    s = s.strip()
    if not s:
        return None
    if s.startswith("<"):
        return 0.0
    s = s.replace(",", "")
    try:
        return float(s)
    except ValueError:
        return None


def parse_nicy_csv(path: str) -> list[dict]:
    """Ni_Cy.csv 파싱 → 22개 날짜별 dict 리스트."""
    with open(path, "r", encoding="utf-8-sig") as f:
        raw_lines = f.readlines()

    reader = csv.reader(raw_lines)
    rows = list(reader)

    date_row = rows[1]
    dates = [d.strip() for d in date_row[2:] if d.strip()]
    n_dates = len(dates)

    def get_values(row_idx: int) -> list:
        row = rows[row_idx] if row_idx < len(rows) else []
        vals = []
        for i in range(2, 2 + n_dates):
            if i < len(row):
                vals.append(_parse_numeric(row[i]))
            else:
                vals.append(None)
        return vals

    org_ml_min = get_values(6)
    aq_ml_min = get_values(7)
    naoh_ml_min = get_values(8)
    naoh_wt_pct = get_values(9)

    feed_metal_names = ["Ca", "Co", "Cr", "Li", "Mg", "Mn", "Ni", "Zn"]
    feed_metals = {}
    for j, metal in enumerate(feed_metal_names):
        feed_metals[metal] = get_values(11 + j)

    feed_ph = get_values(19)

    raff1_metals = {}
    for j, metal in enumerate(feed_metal_names):
        raff1_metals[metal] = get_values(20 + j)

    raff1_ph = get_values(28)

    raff3_metals = {}
    for j, metal in enumerate(feed_metal_names):
        raff3_metals[metal] = get_values(29 + j)

    raff3_ph = get_values(37)

    results = []
    for i in range(n_dates):
        entry = {
            "date": dates[i],
            "C_ext": 0.473107913,
            "n_stages": 3,
            "extractant": "Cyanex 272",
            "org_ml_min": org_ml_min[i],
            "aq_ml_min": aq_ml_min[i],
            "naoh_ml_min": naoh_ml_min[i],
            "naoh_wt_pct": naoh_wt_pct[i] if naoh_wt_pct[i] is not None else 10.0,
            "feed_metals_mg_l": {m: feed_metals[m][i] for m in feed_metal_names if feed_metals[m][i] is not None},
            "feed_ph": feed_ph[i],
            "raff1_metals_mg_l": {m: raff1_metals[m][i] for m in feed_metal_names if raff1_metals[m][i] is not None},
            "raff1_ph": raff1_ph[i],
            "raff3_metals_mg_l": {m: raff3_metals[m][i] for m in feed_metal_names if raff3_metals[m][i] is not None},
            "raff3_ph": raff3_ph[i],
        }
        entry["Q_org_lhr"] = (entry["org_ml_min"] or 0) * 0.06
        entry["Q_aq_lhr"] = (entry["aq_ml_min"] or 0) * 0.06
        entry["Q_naoh_lhr"] = (entry["naoh_ml_min"] or 0) * 0.06 if entry["naoh_ml_min"] is not None else None
        entry["feed_g_l"] = {m: v / 1000.0 for m, v in entry["feed_metals_mg_l"].items()}
        entry["raff1_Ni_g_l"] = entry["raff1_metals_mg_l"].get("Ni", 0) / 1000.0
        entry["raff3_Ni_g_l"] = entry["raff3_metals_mg_l"].get("Ni", 0) / 1000.0
        if entry["Q_aq_lhr"] > 0:
            entry["OA_ratio"] = entry["Q_org_lhr"] / entry["Q_aq_lhr"]
        else:
            entry["OA_ratio"] = 0
        results.append(entry)

    return results


def parse_bmcy_csv(path: str) -> list[dict]:
    """BM-CY.csv 파싱 → 11일분 데이터."""
    with open(path, "r", encoding="utf-8-sig") as f:
        raw_lines = f.readlines()

    reader = csv.reader(raw_lines)
    rows = list(reader)

    metal_cols = {"Ca": 2, "Co": 3, "Cr": 4, "Li": 5, "Mg": 6, "Mn": 7, "Ni": 8, "Zn": 9}
    ph_col = 10
    cond_offset = 13

    conditions = {}
    for row_idx in range(1, min(13, len(rows))):
        row = rows[row_idx]
        if len(row) <= cond_offset:
            continue
        day_label = row[cond_offset].strip() if len(row) > cond_offset else ""
        if not day_label or "일차" not in day_label:
            continue
        day_num = int(re.search(r"(\d+)", day_label).group(1))
        org_ml_min = _parse_numeric(row[cond_offset + 4]) if len(row) > cond_offset + 4 else None
        naoh_ml_min = _parse_numeric(row[cond_offset + 6]) if len(row) > cond_offset + 6 else None
        feed_ml_min = _parse_numeric(row[cond_offset + 8]) if len(row) > cond_offset + 8 else None
        conditions[day_num] = {
            "org_ml_min": org_ml_min,
            "naoh_ml_min": naoh_ml_min,
            "feed_ml_min": feed_ml_min,
        }

    results = []
    for day in range(1, 12):
        day_data = {"day": day, "stages": {}}
        base_row = 1 + (day - 1) * 6

        for s in range(6):
            row_idx = base_row + s
            if row_idx >= len(rows):
                break
            row = rows[row_idx]
            sample_name = row[1].strip() if len(row) > 1 else ""

            metals = {}
            for m, col in metal_cols.items():
                if col < len(row):
                    val = _parse_numeric(row[col])
                    metals[m] = val if val is not None else 0.0
                else:
                    metals[m] = 0.0

            ph_val = _parse_numeric(row[ph_col]) if ph_col < len(row) else None

            if "BM Feed" in sample_name:
                day_data["feed"] = {"metals_mg_l": metals, "pH": ph_val}
            else:
                stage_match = re.search(r"(\d+)단", sample_name)
                if stage_match:
                    stage_num = int(stage_match.group(1))
                    day_data["stages"][stage_num] = {"metals_mg_l": metals, "pH": ph_val}

        cond = conditions.get(day, {})
        day_data["org_ml_min"] = cond.get("org_ml_min", 118.0)
        day_data["feed_ml_min"] = cond.get("feed_ml_min", 25.0)
        day_data["naoh_ml_min"] = cond.get("naoh_ml_min", 0)

        day_data["Q_org_lhr"] = day_data["org_ml_min"] * 0.06
        day_data["Q_aq_lhr"] = day_data["feed_ml_min"] * 0.06
        day_data["Q_naoh_lhr"] = day_data["naoh_ml_min"] * 0.06 if day_data["naoh_ml_min"] else 0

        day_data["OA_ratio"] = day_data["Q_org_lhr"] / day_data["Q_aq_lhr"] if day_data["Q_aq_lhr"] > 0 else 0

        if "feed" in day_data:
            day_data["feed_g_l"] = {m: v / 1000.0 for m, v in day_data["feed"]["metals_mg_l"].items()}
            day_data["feed_ph"] = day_data["feed"]["pH"]
        else:
            day_data["feed_g_l"] = {}
            day_data["feed_ph"] = 3.0

        if 1 in day_data["stages"]:
            day_data["raff1_Ni_g_l"] = day_data["stages"][1]["metals_mg_l"].get("Ni", 0) / 1000.0
            day_data["raff1_Ni_mg_l"] = day_data["stages"][1]["metals_mg_l"].get("Ni", 0)

        day_data["C_ext"] = 0.63081055

        results.append(day_data)

    return results


def _build_feed_dict(feed_g_l: dict) -> dict:
    """시뮬레이터 입력용 금속 딕셔너리 (Cr 제외, 7종만)."""
    return {m: feed_g_l.get(m, 0.0) for m in DEFAULT_METALS}


# ═════════════════════════════════════════════════════════════════════
# Part C: 시뮬레이션 실행
# ═════════════════════════════════════════════════════════════════════

def run_report1(nicy_data: list, sim_bundle: dict) -> list:
    """Report 1: NaOH 사포니피케이션 모드 (Ni-Cy 3단, NaOH 유량 있는 14건)."""
    results = []
    for entry in nicy_data:
        if entry["date"] == "25.03.11":
            continue
        if entry["Q_naoh_lhr"] is None:
            continue

        feed = _build_feed_dict(entry["feed_g_l"])
        sulfate = sim_bundle["calc_sulfate"](feed)
        naoh_m = sim_bundle["estimate_naoh_m"](entry["naoh_wt_pct"])

        try:
            sim = run_sim(
                sim_bundle,
                C_aq_feed=feed,
                pH_feed=entry["feed_ph"],
                Q_aq=entry["Q_aq_lhr"],
                Q_org=entry["Q_org_lhr"],
                extractant="Cyanex 272",
                C_ext=entry["C_ext"],
                n_stages=3,
                C_NaOH=naoh_m,
                Q_NaOH=entry["Q_naoh_lhr"],
                naoh_mode="saponification",
                saponification_model="physical_v2",
                use_competition=True,
                C_sulfate=sulfate,
            )
            sim_raff1_ni = sim["raffinate"]["Ni"]
            sim_raff3_ni = sim["stages"][0]["C_aq_out"]["Ni"]
            sim_raff1_ph = sim["pH_profile"][-1]
            sim_raff3_ph = sim["pH_profile"][0]
            converged = sim["converged"]
        except Exception as e:
            print(f"  [WARN] R1 {sim_bundle['branch_name']} {entry['date']}: {e}")
            sim_raff1_ni = sim_raff3_ni = sim_raff1_ph = sim_raff3_ph = None
            converged = False

        results.append({
            "date": entry["date"],
            "feed_Ni": entry["feed_g_l"].get("Ni", 0),
            "feed_ph": entry["feed_ph"],
            "OA_ratio": round(entry["OA_ratio"], 1),
            "Q_NaOH_lhr": entry["Q_naoh_lhr"],
            "naoh_wt_pct": entry["naoh_wt_pct"],
            "sim_raff1_Ni": sim_raff1_ni,
            "act_raff1_Ni": entry["raff1_Ni_g_l"],
            "sim_raff1_pH": sim_raff1_ph,
            "act_raff1_pH": entry["raff1_ph"],
            "sim_raff3_Ni": sim_raff3_ni,
            "act_raff3_Ni": entry["raff3_Ni_g_l"],
            "sim_raff3_pH": sim_raff3_ph,
            "act_raff3_pH": entry["raff3_ph"],
            "converged": converged,
        })

    return results


def run_report2_nicy(nicy_data: list, sim_bundle: dict) -> list:
    """Report 2 (Ni-Cy part): pH 평형 모드 — 20건."""
    results = []
    for entry in nicy_data:
        if entry["date"] == "25.03.11":
            continue
        if entry["raff3_ph"] is None:
            continue

        feed = _build_feed_dict(entry["feed_g_l"])
        sulfate = sim_bundle["calc_sulfate"](feed)

        raff1_ph = entry["raff1_ph"]
        raff3_ph = entry["raff3_ph"]
        mid_ph = (raff1_ph + raff3_ph) / 2.0
        target_phs = [raff3_ph, mid_ph, raff1_ph]

        try:
            sim = run_sim(
                sim_bundle,
                C_aq_feed=feed,
                pH_feed=entry["feed_ph"],
                Q_aq=entry["Q_aq_lhr"],
                Q_org=entry["Q_org_lhr"],
                extractant="Cyanex 272",
                C_ext=entry["C_ext"],
                n_stages=3,
                target_pH_per_stage=target_phs,
                naoh_mode="aqueous_direct",
                C_NaOH=5.0,
                use_competition=True,
                C_sulfate=sulfate,
            )
            sim_raff_ni = sim["raffinate"]["Ni"]
            sim_ex = (1.0 - sim_raff_ni / feed["Ni"]) * 100.0 if feed["Ni"] > 0 else 0
            converged = sim["converged"]
        except Exception as e:
            print(f"  [WARN] R2-NiCy {sim_bundle['branch_name']} {entry['date']}: {e}")
            sim_raff_ni = sim_ex = None
            converged = False

        act_raff_ni = entry["raff1_Ni_g_l"]
        act_ex = (1.0 - act_raff_ni / feed["Ni"]) * 100.0 if feed["Ni"] > 0 else 0

        results.append({
            "label": entry["date"],
            "feed_Ni": feed["Ni"],
            "sim_raff": sim_raff_ni,
            "act_raff": act_raff_ni,
            "sim_ex": round(sim_ex, 1) if sim_ex is not None else None,
            "act_ex": round(act_ex, 1),
            "converged": converged,
        })

    return results


def run_report2_bmcy(bmcy_data: list, sim_bundle: dict) -> list:
    """Report 2 (BM-CY part): pH 평형, Ni-only — 11건."""
    results = []
    for day_data in bmcy_data:
        day = day_data["day"]

        feed = {"Ni": day_data["feed_g_l"].get("Ni", 0)}
        for m in DEFAULT_METALS:
            if m != "Ni":
                feed[m] = 0.0

        sulfate = sim_bundle["calc_sulfate"](feed)

        stage_phs = []
        for s in range(5, 0, -1):
            if s in day_data["stages"]:
                stage_phs.append(day_data["stages"][s].get("pH", 5.5))
            else:
                stage_phs.append(5.5)

        try:
            sim = run_sim(
                sim_bundle,
                C_aq_feed=feed,
                pH_feed=day_data["feed_ph"],
                Q_aq=day_data["Q_aq_lhr"],
                Q_org=day_data["Q_org_lhr"],
                extractant="Cyanex 272",
                C_ext=day_data["C_ext"],
                n_stages=5,
                target_pH_per_stage=stage_phs,
                naoh_mode="aqueous_direct",
                C_NaOH=5.0,
                use_competition=True,
                C_sulfate=sulfate,
            )
            sim_raff_ni = sim["raffinate"]["Ni"]
            sim_ex = (1.0 - sim_raff_ni / feed["Ni"]) * 100.0 if feed["Ni"] > 0 else 0
            converged = sim["converged"]
        except Exception as e:
            print(f"  [WARN] R2-BMCY {sim_bundle['branch_name']} day {day}: {e}")
            sim_raff_ni = sim_ex = None
            converged = False

        act_raff_ni = day_data.get("raff1_Ni_g_l", 0)
        act_ex = (1.0 - act_raff_ni / feed["Ni"]) * 100.0 if feed["Ni"] > 0 else 0

        results.append({
            "label": f"{day}일차",
            "day": day,
            "feed_Ni": feed["Ni"],
            "sim_raff": sim_raff_ni,
            "act_raff": act_raff_ni,
            "sim_ex": round(sim_ex, 1) if sim_ex is not None else None,
            "act_ex": round(act_ex, 1),
            "converged": converged,
        })

    return results


def run_report3(bmcy_data: list, sim_bundle: dict) -> list:
    """Report 3: BM-CY 다금속 (5단, 11건)."""
    results = []
    for day_data in bmcy_data:
        day = day_data["day"]
        feed = _build_feed_dict(day_data["feed_g_l"])
        sulfate = sim_bundle["calc_sulfate"](feed)

        stage_phs = []
        for s in range(5, 0, -1):
            if s in day_data["stages"]:
                stage_phs.append(day_data["stages"][s].get("pH", 5.5))
            else:
                stage_phs.append(5.5)

        try:
            sim = run_sim(
                sim_bundle,
                C_aq_feed=feed,
                pH_feed=day_data["feed_ph"],
                Q_aq=day_data["Q_aq_lhr"],
                Q_org=day_data["Q_org_lhr"],
                extractant="Cyanex 272",
                C_ext=day_data["C_ext"],
                n_stages=5,
                target_pH_per_stage=stage_phs,
                naoh_mode="aqueous_direct",
                C_NaOH=5.0,
                use_competition=True,
                C_sulfate=sulfate,
            )
            converged = sim["converged"]
            all_ex = {}
            all_raff = {}
            for m in DEFAULT_METALS:
                all_raff[m] = sim["raffinate"].get(m, 0)
                if feed.get(m, 0) > 0:
                    all_ex[m] = (1.0 - all_raff[m] / feed[m]) * 100.0
                else:
                    all_ex[m] = 0.0
        except Exception as e:
            print(f"  [WARN] R3 {sim_bundle['branch_name']} day {day}: {e}")
            converged = False
            all_ex = {m: None for m in DEFAULT_METALS}
            all_raff = {m: None for m in DEFAULT_METALS}

        # 실측 추출률 (plant 1단 Raff' = sim raffinate)
        act_ex = {}
        if 1 in day_data["stages"]:
            raff1 = day_data["stages"][1]["metals_mg_l"]
            for m in DEFAULT_METALS:
                feed_val = feed.get(m, 0)
                raff_val = raff1.get(m, 0) / 1000.0  # mg/L → g/L
                if feed_val > 0:
                    act_ex[m] = (1.0 - raff_val / feed_val) * 100.0
                else:
                    act_ex[m] = 0.0
        else:
            act_ex = {m: None for m in DEFAULT_METALS}

        results.append({
            "label": f"{day}일차",
            "day": day,
            "feed": feed,
            "sim_ex": all_ex,
            "act_ex": act_ex,
            "sim_raff": all_raff,
            "converged": converged,
        })

    return results


# ═════════════════════════════════════════════════════════════════════
# Part D: 통계 함수
# ═════════════════════════════════════════════════════════════════════

def calc_mae(sim_vals: list, act_vals: list) -> float:
    pairs = [(s, a) for s, a in zip(sim_vals, act_vals) if s is not None and a is not None]
    if not pairs:
        return float("nan")
    return sum(abs(s - a) for s, a in pairs) / len(pairs)


def calc_mean_err(sim_vals: list, act_vals: list) -> float:
    pairs = [(s, a) for s, a in zip(sim_vals, act_vals) if s is not None and a is not None]
    if not pairs:
        return float("nan")
    return sum(s - a for s, a in pairs) / len(pairs)


def calc_mape(sim_vals: list, act_vals: list) -> float:
    pairs = [(s, a) for s, a in zip(sim_vals, act_vals) if s is not None and a is not None and abs(a) > 0.001]
    if not pairs:
        return float("nan")
    return sum(abs(s - a) / abs(a) for s, a in pairs) / len(pairs) * 100.0


def _fmt(v, prec=2) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "N/A"
    return f"{v:.{prec}f}"


def _winner(main_val, col_val, lower_better=True) -> str:
    """두 값 비교하여 승자 표시."""
    if main_val is None or col_val is None:
        return "-"
    if isinstance(main_val, float) and math.isnan(main_val):
        return "Colleague"
    if isinstance(col_val, float) and math.isnan(col_val):
        return "Main"
    if lower_better:
        if main_val < col_val:
            return "**Main**"
        elif col_val < main_val:
            return "**Colleague**"
    else:
        if main_val > col_val:
            return "**Main**"
        elif col_val > main_val:
            return "**Colleague**"
    return "Draw"


# ═════════════════════════════════════════════════════════════════════
# Part E: Markdown 보고서 생성
# ═════════════════════════════════════════════════════════════════════

def generate_report(
    main_commit: str,
    col_commit: str,
    r1_main: list, r1_col: list,
    r2_nicy_main: list, r2_nicy_col: list,
    r2_bmcy_main: list, r2_bmcy_col: list,
    r3_main: list, r3_col: list,
) -> str:
    """Markdown 보고서 생성."""

    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    lines = []
    W = lines.append

    W(f"# SX Simulator 브랜치 정확도 비교\n")
    W(f"> 생성: {now}\n")

    # ── 개요 ──
    W("## 개요\n")
    W("| 항목 | Main | Colleague |")
    W("|------|------|-----------|")
    W(f"| Branch | `main` | `colleague/ahh-modifications` |")
    W(f"| Commit | `{main_commit[:7]}` | `{col_commit[:7]}` |")
    W("| 모델 | 4PL Sigmoid | 5PL Richards |")
    W("| Competition | 토글 (ON 사용) | 항상 ON |")
    W("| Loading damping | 있음 | 제거 |")
    W("")

    # ── Report 1 통계 ──
    def _r1_stats(results):
        s1 = [r["sim_raff1_Ni"] for r in results]
        a1 = [r["act_raff1_Ni"] for r in results]
        s3 = [r["sim_raff3_Ni"] for r in results]
        a3 = [r["act_raff3_Ni"] for r in results]
        sp1 = [r["sim_raff1_pH"] for r in results]
        ap1 = [r["act_raff1_pH"] for r in results]
        return {
            "n": len(results),
            "raff1_ni_mae": calc_mae(s1, a1),
            "raff3_ni_mape": calc_mape(s3, a3),
            "raff1_ph_mae": calc_mae(sp1, ap1),
            "raff1_ni_mean_err": calc_mean_err(s1, a1),
        }

    st1_m = _r1_stats(r1_main)
    st1_c = _r1_stats(r1_col)

    # ── Report 2 통계 ──
    def _r2_stats(nicy, bmcy):
        all_sim_raff = [r["sim_raff"] for r in nicy] + [r["sim_raff"] for r in bmcy]
        all_act_raff = [r["act_raff"] for r in nicy] + [r["act_raff"] for r in bmcy]
        all_sim_ex = [r["sim_ex"] for r in nicy] + [r["sim_ex"] for r in bmcy]
        all_act_ex = [r["act_ex"] for r in nicy] + [r["act_ex"] for r in bmcy]
        return {
            "total_n": len(nicy) + len(bmcy),
            "raff_ni_mae": calc_mae(all_sim_raff, all_act_raff),
            "ex_mae": calc_mae(all_sim_ex, all_act_ex),
            "ex_mean_err": calc_mean_err(all_sim_ex, all_act_ex),
            "nicy_raff_mae": calc_mae(
                [r["sim_raff"] for r in nicy], [r["act_raff"] for r in nicy]
            ),
            "nicy_ex_mae": calc_mae(
                [r["sim_ex"] for r in nicy], [r["act_ex"] for r in nicy]
            ),
            "bmcy_raff_mae": calc_mae(
                [r["sim_raff"] for r in bmcy], [r["act_raff"] for r in bmcy]
            ),
            "bmcy_ex_mae": calc_mae(
                [r["sim_ex"] for r in bmcy], [r["act_ex"] for r in bmcy]
            ),
        }

    st2_m = _r2_stats(r2_nicy_main, r2_bmcy_main)
    st2_c = _r2_stats(r2_nicy_col, r2_bmcy_col)

    # ── Report 3 통계 ──
    def _r3_metal_mae(results, metal):
        sim = [r["sim_ex"].get(metal) for r in results]
        act = [r["act_ex"].get(metal) for r in results]
        return calc_mae(sim, act)

    r3_metals_mae_m = {m: _r3_metal_mae(r3_main, m) for m in DEFAULT_METALS}
    r3_metals_mae_c = {m: _r3_metal_mae(r3_col, m) for m in DEFAULT_METALS}

    # ── 요약 ──
    W("## 요약\n")
    W("| 지표 | Main | Colleague | 승자 |")
    W("|------|------|-----------|------|")
    W(f"| R1: Raff1 Ni MAE (g/L) | {_fmt(st1_m['raff1_ni_mae'],2)} | {_fmt(st1_c['raff1_ni_mae'],2)} | {_winner(st1_m['raff1_ni_mae'], st1_c['raff1_ni_mae'])} |")
    W(f"| R1: Raff3 Ni MAPE (%) | {_fmt(st1_m['raff3_ni_mape'],1)} | {_fmt(st1_c['raff3_ni_mape'],1)} | {_winner(st1_m['raff3_ni_mape'], st1_c['raff3_ni_mape'])} |")
    W(f"| R1: Raff1 pH MAE | {_fmt(st1_m['raff1_ph_mae'],2)} | {_fmt(st1_c['raff1_ph_mae'],2)} | {_winner(st1_m['raff1_ph_mae'], st1_c['raff1_ph_mae'])} |")
    W(f"| R2: Raff Ni MAE (g/L) | {_fmt(st2_m['raff_ni_mae'],2)} | {_fmt(st2_c['raff_ni_mae'],2)} | {_winner(st2_m['raff_ni_mae'], st2_c['raff_ni_mae'])} |")
    W(f"| R2: 추출률 MAE (%) | {_fmt(st2_m['ex_mae'],1)} | {_fmt(st2_c['ex_mae'],1)} | {_winner(st2_m['ex_mae'], st2_c['ex_mae'])} |")

    # R3 핵심 금속
    for m in ["Ni", "Li", "Co"]:
        W(f"| R3: {m} 추출률 MAE (%) | {_fmt(r3_metals_mae_m[m],1)} | {_fmt(r3_metals_mae_c[m],1)} | {_winner(r3_metals_mae_m[m], r3_metals_mae_c[m])} |")
    W("")

    # ── Report 1 상세 ──
    W("---")
    W("## Report 1: NaOH 사포니피케이션 (Ni-Cy 3단, 14건)\n")
    W("### 통계 비교\n")
    W("| 지표 | Main | Colleague |")
    W("|------|------|-----------|")
    W(f"| 데이터 수 | {st1_m['n']} | {st1_c['n']} |")
    W(f"| Raff1 Ni MAE (g/L) | {_fmt(st1_m['raff1_ni_mae'],3)} | {_fmt(st1_c['raff1_ni_mae'],3)} |")
    W(f"| Raff1 Ni Mean Error (g/L) | {_fmt(st1_m['raff1_ni_mean_err'],3)} | {_fmt(st1_c['raff1_ni_mean_err'],3)} |")
    W(f"| Raff3 Ni MAPE (%) | {_fmt(st1_m['raff3_ni_mape'],1)} | {_fmt(st1_c['raff3_ni_mape'],1)} |")
    W(f"| Raff1 pH MAE | {_fmt(st1_m['raff1_ph_mae'],3)} | {_fmt(st1_c['raff1_ph_mae'],3)} |")
    W("")

    W("### 상세 데이터\n")
    W("| Date | Feed Ni | Act Raff1 Ni | Main Raff1 Ni | Col Raff1 Ni | Act Raff3 Ni | Main Raff3 Ni | Col Raff3 Ni | Act pH | Main pH | Col pH |")
    W("|------|---------|-------------|---------------|--------------|-------------|---------------|--------------|--------|---------|--------|")
    for rm, rc in zip(r1_main, r1_col):
        W(f"| {rm['date']} "
          f"| {_fmt(rm['feed_Ni'],1)} "
          f"| {_fmt(rm['act_raff1_Ni'],2)} "
          f"| {_fmt(rm['sim_raff1_Ni'],2)} "
          f"| {_fmt(rc['sim_raff1_Ni'],2)} "
          f"| {_fmt(rm['act_raff3_Ni'],2)} "
          f"| {_fmt(rm['sim_raff3_Ni'],2)} "
          f"| {_fmt(rc['sim_raff3_Ni'],2)} "
          f"| {_fmt(rm['act_raff1_pH'],2)} "
          f"| {_fmt(rm['sim_raff1_pH'],2)} "
          f"| {_fmt(rc['sim_raff1_pH'],2)} |")
    W("")

    # ── Report 2 상세 ──
    W("---")
    W("## Report 2: pH 평형 모드 (Ni-Cy 20건 + BM-CY 11건)\n")
    W("### 통계 비교\n")
    W("| 지표 | Main | Colleague |")
    W("|------|------|-----------|")
    W(f"| 전체 데이터 수 | {st2_m['total_n']} | {st2_c['total_n']} |")
    W(f"| 전체 Raff Ni MAE (g/L) | {_fmt(st2_m['raff_ni_mae'],2)} | {_fmt(st2_c['raff_ni_mae'],2)} |")
    W(f"| 전체 추출률 MAE (%) | {_fmt(st2_m['ex_mae'],1)} | {_fmt(st2_c['ex_mae'],1)} |")
    W(f"| 전체 추출률 Mean Error (%) | {_fmt(st2_m['ex_mean_err'],1)} | {_fmt(st2_c['ex_mean_err'],1)} |")
    W(f"| Ni-Cy Raff MAE (g/L) | {_fmt(st2_m['nicy_raff_mae'],2)} | {_fmt(st2_c['nicy_raff_mae'],2)} |")
    W(f"| Ni-Cy 추출률 MAE (%) | {_fmt(st2_m['nicy_ex_mae'],1)} | {_fmt(st2_c['nicy_ex_mae'],1)} |")
    W(f"| BM-CY Raff MAE (g/L) | {_fmt(st2_m['bmcy_raff_mae'],2)} | {_fmt(st2_c['bmcy_raff_mae'],2)} |")
    W(f"| BM-CY 추출률 MAE (%) | {_fmt(st2_m['bmcy_ex_mae'],1)} | {_fmt(st2_c['bmcy_ex_mae'],1)} |")
    W("")

    W("### Ni-Cy 상세\n")
    W("| Label | Feed Ni | Act Raff | Main Raff | Col Raff | Act Ex% | Main Ex% | Col Ex% |")
    W("|-------|---------|---------|-----------|----------|---------|----------|---------|")
    for rm, rc in zip(r2_nicy_main, r2_nicy_col):
        W(f"| {rm['label']} "
          f"| {_fmt(rm['feed_Ni'],3)} "
          f"| {_fmt(rm['act_raff'],3)} "
          f"| {_fmt(rm['sim_raff'],3)} "
          f"| {_fmt(rc['sim_raff'],3)} "
          f"| {_fmt(rm['act_ex'],1)} "
          f"| {_fmt(rm['sim_ex'],1)} "
          f"| {_fmt(rc['sim_ex'],1)} |")
    W("")

    W("### BM-CY 상세\n")
    W("| Label | Feed Ni | Act Raff | Main Raff | Col Raff | Act Ex% | Main Ex% | Col Ex% |")
    W("|-------|---------|---------|-----------|----------|---------|----------|---------|")
    for rm, rc in zip(r2_bmcy_main, r2_bmcy_col):
        W(f"| {rm['label']} "
          f"| {_fmt(rm['feed_Ni'],3)} "
          f"| {_fmt(rm['act_raff'],3)} "
          f"| {_fmt(rm['sim_raff'],3)} "
          f"| {_fmt(rc['sim_raff'],3)} "
          f"| {_fmt(rm['act_ex'],1)} "
          f"| {_fmt(rm['sim_ex'],1)} "
          f"| {_fmt(rc['sim_ex'],1)} |")
    W("")

    # ── Report 3 상세 ──
    W("---")
    W("## Report 3: BM-CY 다금속 (5단, 11건)\n")
    W("### 금속별 추출률 MAE (%)\n")
    W("| Metal | Main | Colleague | 승자 |")
    W("|-------|------|-----------|------|")
    for m in DEFAULT_METALS:
        W(f"| {m} | {_fmt(r3_metals_mae_m[m],1)} | {_fmt(r3_metals_mae_c[m],1)} | {_winner(r3_metals_mae_m[m], r3_metals_mae_c[m])} |")
    W("")

    W("### Ni 추출률 상세\n")
    W("| Label | Act Ni Ex% | Main Ni Ex% | Col Ni Ex% |")
    W("|-------|-----------|-------------|------------|")
    for rm, rc in zip(r3_main, r3_col):
        act_ni = rm["act_ex"].get("Ni")
        sim_ni_m = rm["sim_ex"].get("Ni")
        sim_ni_c = rc["sim_ex"].get("Ni")
        W(f"| {rm['label']} | {_fmt(act_ni,1)} | {_fmt(sim_ni_m,1)} | {_fmt(sim_ni_c,1)} |")
    W("")

    W("### Li 추출률 상세\n")
    W("| Label | Act Li Ex% | Main Li Ex% | Col Li Ex% |")
    W("|-------|-----------|-------------|------------|")
    for rm, rc in zip(r3_main, r3_col):
        act_li = rm["act_ex"].get("Li")
        sim_li_m = rm["sim_ex"].get("Li")
        sim_li_c = rc["sim_ex"].get("Li")
        W(f"| {rm['label']} | {_fmt(act_li,1)} | {_fmt(sim_li_m,1)} | {_fmt(sim_li_c,1)} |")
    W("")

    W("---")
    W(f"*Report generated by `scripts/branch_accuracy_comparison.py` on {now}*\n")

    return "\n".join(lines)


# ═════════════════════════════════════════════════════════════════════
# Part F: main()
# ═════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("  SX Simulator 브랜치 정확도 비교")
    print("=" * 60)

    # 1. Worktree 생성
    print("\n[1/7] Worktree 설정...")
    worktree_path = setup_worktree()

    # 2. CSV 파싱 (1회)
    print("\n[2/7] CSV 데이터 파싱...")
    nicy_data = parse_nicy_csv(str(DATA_DIR / "Ni_Cy.csv"))
    bmcy_data = parse_bmcy_csv(str(DATA_DIR / "BM-CY.csv"))
    print(f"  Ni-Cy: {len(nicy_data)}건, BM-CY: {len(bmcy_data)}건")

    # 3. Main 시뮬레이션
    print("\n[3/7] Main 브랜치 시뮬레이터 로드...")
    main_sim = load_simulator(PROJECT_ROOT, "main")

    # 커밋 해시 추출
    main_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=str(PROJECT_ROOT), capture_output=True, text=True
    ).stdout.strip()
    col_commit = subprocess.run(
        ["git", "rev-parse", "colleague/ahh-modifications"],
        cwd=str(PROJECT_ROOT), capture_output=True, text=True
    ).stdout.strip()

    print("\n[4/7] Main 시뮬레이션 실행...")
    r1_main = run_report1(nicy_data, main_sim)
    r2_nicy_main = run_report2_nicy(nicy_data, main_sim)
    r2_bmcy_main = run_report2_bmcy(bmcy_data, main_sim)
    r3_main = run_report3(bmcy_data, main_sim)
    print(f"  R1: {len(r1_main)}건, R2-NiCy: {len(r2_nicy_main)}건, R2-BMCY: {len(r2_bmcy_main)}건, R3: {len(r3_main)}건")

    # 4. sys.modules 퍼지
    print("\n[5/7] 모듈 퍼지 + Colleague 시뮬레이터 로드...")
    col_sim = load_simulator(worktree_path, "colleague")

    # 5. Colleague 시뮬레이션
    print("\n[6/7] Colleague 시뮬레이션 실행...")
    r1_col = run_report1(nicy_data, col_sim)
    r2_nicy_col = run_report2_nicy(nicy_data, col_sim)
    r2_bmcy_col = run_report2_bmcy(bmcy_data, col_sim)
    r3_col = run_report3(bmcy_data, col_sim)
    print(f"  R1: {len(r1_col)}건, R2-NiCy: {len(r2_nicy_col)}건, R2-BMCY: {len(r2_bmcy_col)}건, R3: {len(r3_col)}건")

    # 6. Markdown 보고서 생성
    print("\n[7/7] 보고서 생성...")
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    report_md = generate_report(
        main_commit, col_commit,
        r1_main, r1_col,
        r2_nicy_main, r2_nicy_col,
        r2_bmcy_main, r2_bmcy_col,
        r3_main, r3_col,
    )

    report_path = REPORTS_DIR / "branch_accuracy_comparison.md"
    report_path.write_text(report_md, encoding="utf-8")
    print(f"\n  보고서 저장: {report_path}")

    # 7. .gitignore에 .worktrees/ 추가
    gitignore_path = PROJECT_ROOT / ".gitignore"
    gitignore_text = gitignore_path.read_text(encoding="utf-8") if gitignore_path.exists() else ""
    if ".worktrees/" not in gitignore_text:
        with open(gitignore_path, "a", encoding="utf-8") as f:
            f.write("\n# Branch comparison worktrees\n.worktrees/\n")
        print("  .gitignore에 .worktrees/ 추가됨")

    print("\n" + "=" * 60)
    print("  완료!")
    print("=" * 60)


if __name__ == "__main__":
    main()
