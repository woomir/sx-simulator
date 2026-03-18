#!/usr/bin/env python3
"""
Jantunen et al. (2022) Figure 3(b) 검증 그래프
==============================================
CSV에서 추출한 pH-E% 데이터를 Plotly 인터랙티브 차트로 재현하고
논문 텍스트의 검증 기준값과 교차 검증합니다.

Usage:
    python scripts/plot_jantunen_fig3b.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_FILE = PROJECT_ROOT / "data" / "jantunen2022_fig3b_cyanex272_OA2.5.csv"
OUTPUT_FILE = PROJECT_ROOT / "reports" / "jantunen_fig3b_verification.html"

# ─────────────────────────────────────────────────────────────────────
# 검증 기준값 (논문 텍스트에서 발췌)
# ─────────────────────────────────────────────────────────────────────
# pH50 값은 Figure 3(a) O/A=1 기준 (graphical interpolation, p.6)
REFERENCE_PH50 = {"Mn": 3.82, "Co": 4.44}

# pH 5.27에서의 기대 추출률 (p.7, O/A=2.5 기준)
VALIDATION_AT_PH527 = {
    "Mn": 99.8,
    "Co": 99.0,
    "Ni": 13.9,
    "Li": 3.6,
}

# ─────────────────────────────────────────────────────────────────────
# 마커 스타일 (원본 논문 재현)
# ─────────────────────────────────────────────────────────────────────
METAL_STYLE = {
    "Li": {"color": "orange", "symbol": "square", "name": "Li"},
    "Na": {"color": "grey", "symbol": "triangle-down", "name": "Na"},
    "Mn": {"color": "deeppink", "symbol": "triangle-up", "name": "Mn"},
    "Co": {"color": "magenta", "symbol": "circle", "name": "Co"},
    "Ni": {"color": "green", "symbol": "diamond", "name": "Ni"},
}

METAL_ORDER = ["Li", "Na", "Mn", "Co", "Ni"]


def load_data(path: Path) -> pd.DataFrame:
    """CSV 파일을 로드합니다 (# 주석 행 무시)."""
    df = pd.read_csv(path, comment="#")
    df.columns = df.columns.str.strip()
    return df


def interpolate_ph50(df: pd.DataFrame, metal: str) -> float | None:
    """선형 보간으로 E=50%에 해당하는 pH를 추정합니다."""
    sub = df[df["metal"] == metal].sort_values("pH")
    ph = sub["pH"].values
    e = sub["E_pct"].values

    if len(ph) < 2 or e.max() < 50:
        return None

    for i in range(len(e) - 1):
        if e[i] <= 50 <= e[i + 1]:
            # 선형 보간
            frac = (50.0 - e[i]) / (e[i + 1] - e[i])
            return float(ph[i] + frac * (ph[i + 1] - ph[i]))
    return None


def build_figure(df: pd.DataFrame) -> go.Figure:
    """Plotly Figure를 생성합니다."""
    fig = go.Figure()

    # 각 금속별 데이터 트레이스
    for metal in METAL_ORDER:
        style = METAL_STYLE[metal]
        sub = df[df["metal"] == metal].sort_values("pH")

        fig.add_trace(go.Scatter(
            x=sub["pH"],
            y=sub["E_pct"],
            mode="markers+lines",
            name=style["name"],
            marker=dict(
                symbol=style["symbol"],
                size=10,
                color=style["color"],
                line=dict(width=1, color="black"),
            ),
            line=dict(color=style["color"], width=1.5, dash="dot"),
            hovertemplate=(
                f"<b>{style['name']}</b><br>"
                "pH = %{x:.2f}<br>"
                "E = %{y:.1f}%<br>"
                "<extra></extra>"
            ),
        ))

    # pH50 기준선 (O/A=1 참고값)
    for metal, ph50 in REFERENCE_PH50.items():
        fig.add_shape(
            type="line",
            x0=ph50, x1=ph50, y0=0, y1=50,
            line=dict(color=METAL_STYLE[metal]["color"], width=1.5, dash="dash"),
        )
        fig.add_shape(
            type="line",
            x0=0, x1=ph50, y0=50, y1=50,
            line=dict(color="lightgrey", width=1, dash="dot"),
        )
        fig.add_annotation(
            x=ph50, y=-4,
            text=f"pH50={ph50}<br>(O/A=1)",
            showarrow=False,
            font=dict(size=9, color=METAL_STYLE[metal]["color"]),
        )

    # pH 5.27 검증선
    fig.add_shape(
        type="line",
        x0=5.27, x1=5.27, y0=0, y1=105,
        line=dict(color="blue", width=1.5, dash="dashdot"),
    )
    fig.add_annotation(
        x=5.27, y=105,
        text="pH 5.27 (검증점)",
        showarrow=False,
        font=dict(size=10, color="blue"),
        yanchor="bottom",
    )

    # 검증값 주석
    annotations_text = []
    for metal, expected in VALIDATION_AT_PH527.items():
        annotations_text.append(f"{metal}: {expected}%")
    fig.add_annotation(
        x=5.27, y=75,
        text="<b>pH 5.27 기대값:</b><br>" + "<br>".join(annotations_text),
        showarrow=True,
        arrowhead=2,
        ax=80, ay=-30,
        bordercolor="blue",
        borderwidth=1,
        borderpad=4,
        bgcolor="rgba(255,255,255,0.9)",
        font=dict(size=10),
    )

    # 레이아웃
    fig.update_layout(
        title=dict(
            text=(
                "Jantunen et al. (2022) Figure 3(b) — "
                "0.8 M Cyanex 272, O/A = 2.5<br>"
                "<sup>T = 25±1°C, t_eq = 15 min | "
                "Feed: Li 2.52, Ni 2.05, Co 16.5, Mn 2.09 g/L</sup>"
            ),
            font=dict(size=16),
        ),
        xaxis=dict(
            title="pH",
            range=[0, 9],
            dtick=1,
            gridcolor="lightgrey",
        ),
        yaxis=dict(
            title="E [%]",
            range=[-5, 110],
            dtick=10,
            gridcolor="lightgrey",
        ),
        legend=dict(
            x=0.02, y=0.98,
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="black",
            borderwidth=1,
        ),
        plot_bgcolor="white",
        width=900,
        height=600,
        hovermode="closest",
    )

    return fig


def run_validation(df: pd.DataFrame):
    """검증 결과를 콘솔에 출력합니다."""
    print("\n" + "=" * 65)
    print("  Jantunen et al. (2022) Fig.3(b) 데이터 검증")
    print("  0.8 M Cyanex 272, O/A = 2.5")
    print("=" * 65)

    # 1) pH50 교차 검증
    print("\n  [1] pH50 교차 검증 (선형 보간)")
    print(f"  {'Metal':>5} | {'보간 pH50':>10} | {'참고 pH50':>10} | {'비고'}")
    print(f"  {'-'*5}-+-{'-'*10}-+-{'-'*10}-+-{'-'*20}")

    for metal in ["Mn", "Co", "Ni", "Li"]:
        ph50_interp = interpolate_ph50(df, metal)
        ph50_ref = REFERENCE_PH50.get(metal)

        if ph50_interp is not None:
            interp_str = f"{ph50_interp:.2f}"
        else:
            interp_str = "N/A"

        if ph50_ref is not None:
            ref_str = f"{ph50_ref:.2f}"
            note = "(O/A=1 참고)"
        else:
            ref_str = "—"
            note = ""

        if ph50_interp is not None and ph50_ref is not None:
            delta = ph50_interp - ph50_ref
            note = f"Δ={delta:+.2f} (O/A=2.5 이동)"

        print(f"  {metal:>5} | {interp_str:>10} | {ref_str:>10} | {note}")

    # 2) pH 5.27 검증
    print(f"\n  [2] pH 5.27 추출률 검증")
    print(f"  {'Metal':>5} | {'데이터 E%':>10} | {'기대 E%':>10} | {'차이':>8} | {'판정'}")
    print(f"  {'-'*5}-+-{'-'*10}-+-{'-'*10}-+-{'-'*8}-+-{'-'*6}")

    all_pass = True
    for metal, expected in VALIDATION_AT_PH527.items():
        sub = df[(df["metal"] == metal) & (np.isclose(df["pH"], 5.27, atol=0.05))]
        if len(sub) > 0:
            actual = sub["E_pct"].iloc[0]
            diff = actual - expected
            ok = abs(diff) < 3.0  # ±3% 허용
            status = "OK" if ok else "WARN"
            if not ok:
                all_pass = False
            print(f"  {metal:>5} | {actual:>10.1f} | {expected:>10.1f} | {diff:>+8.1f} | {status:>6}")
        else:
            print(f"  {metal:>5} | {'없음':>10} | {expected:>10.1f} | {'—':>8} | {'MISS':>6}")
            all_pass = False

    # 3) 데이터 요약
    print(f"\n  [3] 데이터 요약")
    for metal in METAL_ORDER:
        sub = df[df["metal"] == metal]
        print(f"  {metal:>5}: {len(sub)}개 포인트, "
              f"pH [{sub['pH'].min():.1f}–{sub['pH'].max():.1f}], "
              f"E [{sub['E_pct'].min():.0f}–{sub['E_pct'].max():.0f}%]")

    print(f"\n  총 데이터 포인트: {len(df)}개")

    if all_pass:
        print("\n  ✓ pH 5.27 검증 통과 (모든 값 ±3% 이내)")
    else:
        print("\n  ✗ 일부 검증 불일치 — 데이터 재확인 필요")

    print("=" * 65)


def main():
    if not DATA_FILE.exists():
        print(f"[오류] 데이터 파일을 찾을 수 없습니다: {DATA_FILE}")
        sys.exit(1)

    df = load_data(DATA_FILE)
    print(f"  데이터 로드 완료: {len(df)}개 포인트 ({DATA_FILE.name})")

    # 검증
    run_validation(df)

    # Plotly 그래프 생성
    fig = build_figure(df)

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(OUTPUT_FILE), include_plotlyjs=True)
    print(f"\n  → 인터랙티브 그래프 저장: {OUTPUT_FILE}")
    print(f"  → 브라우저에서 열어 원본 Figure 3(b)와 시각적으로 비교하세요.")


if __name__ == "__main__":
    main()
