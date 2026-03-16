from typing import Optional

"""
Dashboard Tab Renderers
=======================
Streamlit 대시보드의 주요 동적 탭 렌더링 로직을 분리합니다.

목표:
- sx_dashboard.py의 길이를 줄이고 책임을 분리
- 시뮬레이션 계산과 UI 렌더링을 느슨하게 연결
- 탭 단위 수정 시 충돌 범위 축소
"""

import math
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from .dashboard_service import run_compare_simulations
from .extraction_isotherm import (
    calc_loading_fraction,
    distribution_coefficient,
    extraction_efficiency,
    get_effective_k,
    get_effective_pH50,
    pka_dissociation_factor,
)
from .config import EXTRACTANT_PKA
from .fitting import fit_sigmoid, sigmoid_model


def _metal_color_list(metals: list, metal_colors: dict) -> list:
    return [metal_colors.get(metal, "#636e72") for metal in metals]


def render_results_tab(
    result: dict,
    metals: list,
    loading_pct: float,
    scope_assessment: dict,
    profile_label: str,
    selected_preset_note: Optional[str],
    C_aq_feed: dict,
    pH_feed: float,
    n_stages: int,
    metal_colors: dict,
    scope_renderer,
    inputs,
) -> None:
    if not result["converged"]:
        st.error(
            "현재 조건에서 다단 계산이 수렴하지 않았습니다. 결과 수치를 의사결정에 직접 사용하기 전에 입력 조건과 전략을 다시 확인하세요."
        )

    st.subheader("🎯 추출 결과 요약")
    color_list = _metal_color_list(metals, metal_colors)

    met_cols_row1 = st.columns(4)
    for i, metal in enumerate(metals[:4]):
        ext_pct = result["overall_extraction"][metal]
        raff = result["raffinate"][metal]
        met_cols_row1[i].metric(
            label=f"{metal} 추출률",
            value=f"{ext_pct:.1f}%",
            delta=f"후액: {raff:.3f} g/L",
        )
    if len(metals) > 4:
        met_cols_row2 = st.columns(4)
        for i, metal in enumerate(metals[4:]):
            ext_pct = result["overall_extraction"][metal]
            raff = result["raffinate"][metal]
            met_cols_row2[i].metric(
                label=f"{metal} 추출률",
                value=f"{ext_pct:.1f}%",
                delta=f"후액: {raff:.3f} g/L",
            )

    col_info1, col_info2, col_info3, col_info4, col_info5 = st.columns(5)
    with col_info1:
        st.metric("후액 최종 pH", f"{result['pH_profile'][-1]:.2f}")
    with col_info2:
        st.metric("총 NaOH 소비", f"{result['total_NaOH_mol_hr']:.1f} mol/hr")
    
    sap_capacity_mol_hr = inputs.C_ext * inputs.Q_org
    sap_pct = (result['total_NaOH_mol_hr'] / sap_capacity_mol_hr * 100.0) if sap_capacity_mol_hr > 0 else 0.0
    with col_info3:
        st.metric("Saponification %", f"{sap_pct:.1f}%", help="입력된 유기계 유량(Q_org)과 추출제 농도(C_ext) 대비, 실제 계산/투입된 NaOH의 몰 비율입니다.")

    with col_info4:
        if loading_pct > 85:
            st.metric(
                "최대 로딩률",
                f"{loading_pct:.1f}%",
                delta="⚠️ 포화 근접",
                delta_color="inverse",
            )
        else:
            st.metric("최대 로딩률", f"{loading_pct:.1f}%")
    with col_info5:
        if result["converged"]:
            st.metric("수렴 상태", "수렴", delta=f"{result['iterations']} iter")
        else:
            st.metric(
                "수렴 상태",
                "미수렴",
                delta=f"{result['iterations']} iter",
                delta_color="inverse",
            )

    st.markdown("---")
    st.subheader("🧭 현재 시뮬레이션 해석 가이드")
    st.caption(f"활성 프로필: **{profile_label}**")
    scope_renderer(st, scope_assessment)
    if selected_preset_note:
        st.warning(f"프리셋 주의: {selected_preset_note}")

    st.markdown("---")
    st.subheader("📊 금속별 추출률")
    ext_data = pd.DataFrame(
        {
            "금속": metals,
            "추출률 (%)": [result["overall_extraction"][m] for m in metals],
            "피드 (g/L)": [C_aq_feed[m] for m in metals],
            "후액 (g/L)": [result["raffinate"][m] for m in metals],
        }
    )
    fig_ext = px.bar(
        ext_data,
        x="금속",
        y="추출률 (%)",
        color="금속",
        color_discrete_sequence=color_list,
        text="추출률 (%)",
    )
    fig_ext.update_traces(texttemplate="%{text:.1f}%", textposition="outside")
    fig_ext.update_layout(
        height=400,
        showlegend=False,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        yaxis=dict(range=[0, 110]),
    )
    st.plotly_chart(fig_ext, use_container_width=True)

    st.markdown("---")
    col_ph, col_naoh = st.columns(2)
    stages_x = [0] + list(range(1, n_stages + 1))

    with col_ph:
        st.subheader("📉 Stage별 pH 프로파일")
        pH_y = [float(pH_feed)] + result["pH_profile"]
        fig_ph = go.Figure()
        fig_ph.add_trace(
            go.Scatter(
                x=stages_x,
                y=pH_y,
                mode="lines+markers",
                name="pH",
                line=dict(color="#e94560", width=3),
                marker=dict(size=12, symbol="circle"),
            )
        )
        fig_ph.update_layout(
            xaxis_title="Stage",
            yaxis_title="pH",
            height=350,
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            xaxis=dict(
                tickmode="array",
                tickvals=stages_x,
                ticktext=["Feed"] + [f"St.{i}" for i in range(1, n_stages + 1)],
            ),
        )
        st.plotly_chart(fig_ph, use_container_width=True)

    with col_naoh:
        st.subheader("🧪 Stage별 NaOH 소비량")
        fig_naoh = go.Figure()
        naoh_y = [0.0] + result["NaOH_profile"]
        fig_naoh.add_trace(
            go.Bar(
                x=stages_x,
                y=naoh_y,
                marker_color="#0f3460",
                text=[f"{v:.1f}" if v > 0 else "" for v in naoh_y],
                textposition="outside",
            )
        )
        fig_naoh.update_layout(
            xaxis_title="Stage",
            yaxis_title="NaOH (mol/hr)",
            height=350,
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            xaxis=dict(
                tickmode="array",
                tickvals=stages_x,
                ticktext=["Feed"] + [f"St.{i}" for i in range(1, n_stages + 1)],
            ),
        )
        st.plotly_chart(fig_naoh, use_container_width=True)

    st.subheader("🔬 Stage별 수계 금속 농도 변화")
    conc_data = {"Stage": ["Feed"] + [f"Stage {i+1}" for i in range(n_stages)]}
    for metal in metals:
        conc_data[metal] = [C_aq_feed[metal]]
        for stage_result in result["stages"]:
            conc_data[metal].append(stage_result["C_aq_out"][metal])

    df_conc = pd.DataFrame(conc_data)
    fig_conc = go.Figure()
    for i, metal in enumerate(metals):
        fig_conc.add_trace(
            go.Scatter(
                x=df_conc["Stage"],
                y=df_conc[metal],
                mode="lines+markers",
                name=metal,
                line=dict(color=color_list[i], width=2),
                marker=dict(size=8),
            )
        )
    fig_conc.update_layout(
        xaxis_title="Stage",
        yaxis_title="농도 (g/L)",
        height=400,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    st.plotly_chart(fig_conc, use_container_width=True)


def render_isotherm_tab(
    metals: list,
    extractant: str,
    C_ext: float,
    temperature: float,
    target_pH: Optional[float],
    active_extractant_params: dict,
    metal_colors: dict,
    ph_plot_range: list,
) -> None:
    st.subheader(f"📈 pH-추출률 Isotherm ({extractant}, {C_ext}M, {temperature:.0f}°C)")

    color_list = _metal_color_list(metals, metal_colors)
    fig_iso = go.Figure()
    for i, metal in enumerate(metals):
        E_values = [
            extraction_efficiency(
                pH,
                metal,
                extractant,
                C_ext,
                temperature,
                extractant_params=active_extractant_params,
            )
            for pH in ph_plot_range
        ]
        fig_iso.add_trace(
            go.Scatter(
                x=ph_plot_range,
                y=E_values,
                mode="lines",
                name=metal,
                line=dict(color=color_list[i], width=3),
            )
        )
    fig_iso.update_layout(
        xaxis_title="pH",
        yaxis_title="추출률 (%)",
        height=500,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        yaxis=dict(range=[0, 105]),
    )
    if target_pH:
        fig_iso.add_vline(
            x=target_pH,
            line_dash="dash",
            line_color="red",
            annotation_text=f"목표 pH = {target_pH}",
        )
    st.plotly_chart(fig_iso, use_container_width=True)

    st.subheader("📋 금속별 pH₅₀ (보정값)")
    ph50_data = []
    for metal in metals:
        ph50_eff = get_effective_pH50(
            metal,
            extractant,
            C_ext,
            temperature,
            extractant_params=active_extractant_params,
        )
        k_eff = get_effective_k(
            metal,
            extractant,
            temperature,
            extractant_params=active_extractant_params,
        )
        params = active_extractant_params[extractant][metal]
        ph50_data.append(
            {
                "금속": metal,
                "pH₅₀ (기준)": params["pH50"],
                "pH₅₀ (보정)": round(ph50_eff, 2),
                "k (기준)": params["k"],
                "k (보정)": round(k_eff, 2),
                "E_max (%)": params["E_max"],
                "n_H": params["n_H"],
            }
        )
    st.dataframe(pd.DataFrame(ph50_data), use_container_width=True, hide_index=True)


def render_mccabe_thiele_tab(
    metals: list,
    result: dict,
    extractant: str,
    C_ext: float,
    temperature: float,
    C_aq_feed: dict,
    n_stages: int,
    active_extractant_params: dict,
) -> None:
    st.subheader("📉 McCabe-Thiele 다이어그램")
    st.markdown(
        "선택한 금속에 대한 조작선(Operating Line)과 평형 곡선(Equilibrium Curve)을 표시합니다."
    )

    mt_metal = st.selectbox(
        "다이어그램을 그릴 금속 선택",
        metals,
        index=metals.index("Co") if "Co" in metals else 0,
    )

    ref_pH = result["pH_profile"][-1]
    st.caption(f"💡 평형 곡선은 최종 후액 기준 pH ({ref_pH:.2f}) 로 계산되었습니다.")

    x_max = C_aq_feed[mt_metal] * 1.2
    x_vals = np.linspace(0, x_max, 100)

    eq_y_vals = []
    for x in x_vals:
        D = distribution_coefficient(
            ref_pH,
            mt_metal,
            extractant,
            C_ext,
            temperature,
            extractant_params=active_extractant_params,
        )
        sim_org = {mt_metal: D * x}
        loading = calc_loading_fraction(
            sim_org,
            extractant,
            C_ext,
            [mt_metal],
            extractant_params=active_extractant_params,
        )
        damping = max(0.0, 1.0 - loading) if loading > 0.85 else 1.0
        D = D * damping
        eq_y_vals.append(D * x)

    fig_mt = go.Figure()
    fig_mt.add_trace(
        go.Scatter(
            x=x_vals,
            y=eq_y_vals,
            mode="lines",
            name="Equilibrium Curve (평형선)",
            line=dict(color="blue", width=2),
        )
    )

    x_in_feed = C_aq_feed[mt_metal]
    y_out_org = result["loaded_organic"][mt_metal]
    x_out_raff = result["raffinate"][mt_metal]
    y_in_org = 0.0

    fig_mt.add_trace(
        go.Scatter(
            x=[x_out_raff, x_in_feed],
            y=[y_in_org, y_out_org],
            mode="lines",
            name="Operating Line (조작선)",
            line=dict(color="red", width=2, dash="dash"),
        )
    )

    stage_x = [x_out_raff]
    stage_y = [y_in_org]
    current_x = x_out_raff

    for stage_index in range(n_stages - 1, -1, -1):
        actual_y = result["stages"][stage_index]["C_org_out"][mt_metal]
        stage_x.append(current_x)
        stage_y.append(actual_y)

        if stage_index == 0:
            current_x = C_aq_feed[mt_metal]
        else:
            current_x = result["stages"][stage_index - 1]["C_aq_out"][mt_metal]

        stage_x.append(current_x)
        stage_y.append(actual_y)

    fig_mt.add_trace(
        go.Scatter(
            x=stage_x,
            y=stage_y,
            mode="lines+markers",
            name="Stages (실제 단수)",
            line=dict(color="green", width=2, shape="hv"),
            marker=dict(size=8, symbol="circle"),
        )
    )

    fig_mt.update_layout(
        title=f"McCabe-Thiele for {mt_metal}",
        xaxis_title=f"Aqueous {mt_metal} (g/L)",
        yaxis_title=f"Organic {mt_metal} (g/L)",
        height=600,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        legend=dict(x=0.01, y=0.99),
    )
    st.plotly_chart(fig_mt, use_container_width=True)


def render_compare_tab(
    metals: list,
    simulation_inputs,
    active_extractant_params: dict,
    temperature: float,
    C_ext: float,
    n_stages: int,
    C_sulfate: float,
    ph_plot_range: list,
) -> None:
    st.subheader("🔬 Cyanex 272 vs D2EHPA 비교")
    compare_pH = st.slider("비교 목표 pH", 2.0, 8.0, 5.0, 0.1, key="compare_ph")
    st.info(
        f"비교 조건: 공통 C_ext {C_ext:.4f} M, {n_stages}단, {temperature:.0f}°C, SO₄²⁻ {C_sulfate:.3f} M"
    )

    results_compare = run_compare_simulations(
        simulation_inputs,
        compare_pH,
        active_extractant_params,
    )

    compare_data = []
    for extractant in ["Cyanex 272", "D2EHPA"]:
        for metal in metals:
            compare_data.append(
                {
                    "추출제": extractant,
                    "금속": metal,
                    "추출률 (%)": results_compare[extractant]["overall_extraction"][metal],
                }
            )

    df_compare = pd.DataFrame(compare_data)
    fig_compare = px.bar(
        df_compare,
        x="금속",
        y="추출률 (%)",
        color="추출제",
        barmode="group",
        color_discrete_sequence=["#e94560", "#0f3460"],
        text="추출률 (%)",
    )
    fig_compare.update_traces(texttemplate="%{text:.1f}%", textposition="outside")
    fig_compare.update_layout(
        height=450,
        yaxis=dict(range=[0, 115]),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
    )
    st.plotly_chart(fig_compare, use_container_width=True)

    col_c1, col_c2 = st.columns(2)
    with col_c1:
        st.metric(
            "Cyanex 272 NaOH",
            f"{results_compare['Cyanex 272']['total_NaOH_mol_hr']:.1f} mol/hr",
        )
    with col_c2:
        st.metric(
            "D2EHPA NaOH",
            f"{results_compare['D2EHPA']['total_NaOH_mol_hr']:.1f} mol/hr",
        )

    st.markdown("---")
    st.subheader("📈 Isotherm 비교")
    selected_metal = st.selectbox("금속 선택", metals, index=2)

    fig_iso_cmp = go.Figure()
    for extractant, color in [("Cyanex 272", "#e94560"), ("D2EHPA", "#0f3460")]:
        ext_c = 0.5 if extractant == "Cyanex 272" else 0.64
        E_vals = [
            extraction_efficiency(
                pH,
                selected_metal,
                extractant,
                ext_c,
                temperature,
                extractant_params=active_extractant_params,
            )
            for pH in ph_plot_range
        ]
        fig_iso_cmp.add_trace(
            go.Scatter(
                x=ph_plot_range,
                y=E_vals,
                mode="lines",
                name=extractant,
                line=dict(color=color, width=3),
            )
        )
    fig_iso_cmp.add_vline(
        x=compare_pH,
        line_dash="dash",
        line_color="gray",
        annotation_text=f"비교 pH = {compare_pH}",
    )
    fig_iso_cmp.update_layout(
        title=f"{selected_metal} 추출 Isotherm 비교",
        xaxis_title="pH",
        yaxis_title="추출률 (%)",
        height=400,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
    )
    st.plotly_chart(fig_iso_cmp, use_container_width=True)


def render_detail_tab(
    result: dict,
    metals: list,
    C_aq_feed: dict,
    Q_aq: float,
    Q_org: float,
) -> None:
    st.subheader("📋 시뮬레이션 상세 결과")
    st.info(
        f"수렴: {'✅ Yes' if result['converged'] else '❌ No'} | "
        f"반복: {result['iterations']} | "
        f"후액 pH: {result['pH_profile'][-1]:.2f} | "
        f"NaOH 총: {result['total_NaOH_mol_hr']:.1f} mol/hr"
    )

    for i, stage_result in enumerate(result["stages"]):
        with st.expander(
            f"Stage {i+1} (pH = {stage_result['pH_out']:.2f}, NaOH = {stage_result.get('NaOH_consumed_mol_hr', 0):.2f} mol/hr)"
        ):
            rows = []
            for metal in metals:
                rows.append(
                    {
                        "금속": metal,
                        "수계 입구 (g/L)": C_aq_feed[metal]
                        if i == 0
                        else result["stages"][i - 1]["C_aq_out"][metal],
                        "수계 출구 (g/L)": stage_result["C_aq_out"][metal],
                        "유기계 출구 (g/L)": stage_result["C_org_out"][metal],
                        "Stage 추출률 (%)": stage_result["extraction"][metal],
                    }
                )
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown("---")
    st.subheader("⚖️ 물질수지 검증")
    balance_data = []
    for metal in metals:
        feed_flow = C_aq_feed[metal] * Q_aq
        raff_flow = result["raffinate"][metal] * Q_aq
        org_flow = result["loaded_organic"][metal] * Q_org
        balance = feed_flow - raff_flow - org_flow
        balance_data.append(
            {
                "금속": metal,
                "피드 (g/hr)": round(feed_flow, 3),
                "후액 (g/hr)": round(raff_flow, 3),
                "유기계 (g/hr)": round(org_flow, 3),
                "잔차 (g/hr)": round(balance, 3),
            }
    )
    st.dataframe(pd.DataFrame(balance_data), use_container_width=True, hide_index=True)


def render_formula_tab(
    extractant: str,
    C_ext: float,
    temperature: float,
    T_REF: float,
    metals: list,
    active_extractant_params: dict,
    target_pH: Optional[float],
    staged_pHs: Optional[list],
    pH_mode: str,
    Q_aq: float,
    Q_org: float,
    C_sulfate: float,
    n_stages: int,
    metal_colors: dict,
) -> None:
    st.subheader(
        f"📐 시스템 수식 및 메커니즘 ({extractant}, {C_ext}M, {temperature:.0f}°C)"
    )
    st.markdown(
        "본 시뮬레이터는 **시그모이드 근사(OLI/MSE 기반 파라미터) + Vasilyev 경쟁 보정**의 "
        "하이브리드 모델을 사용합니다. 파라미터를 변경하면 수식도 실시간으로 업데이트됩니다."
    )
    st.info(
        "🔬 **모델 구조:**\n\n"
        "- **1층 — 기본 D 계산:** OLI/MSE 열역학 데이터에서 피팅한 시그모이드 함수 (pH50, k, E_max)\n"
        "- **2층 — 경쟁 보정:** Vasilyev et al. (2019)의 공유 추출제 풀 모델로 금속 간 경쟁 반영"
    )
    color_list = _metal_color_list(metals, metal_colors)

    st.markdown("---")
    st.markdown("### 1️⃣ 기본 추출률 (시그모이드 근사 — OLI/MSE 파라미터)")
    st.markdown(
        "각 금속의 pH에 따른 **기본 추출률**을 시그모이드(Sigmoid) 함수로 근사합니다. "
        "이 파라미터(pH50, k, E_max)는 OLI Systems의 MSE 열역학 프레임워크에서 "
        "도출한 추출 곡선을 피팅하여 얻은 값입니다:"
    )
    st.latex(
        r"E_M(\text{pH}, T) = \frac{E_{\max}}{1 + \exp\!\left(-k(T) \cdot (\text{pH} - \text{pH}_{50,\text{eff}}(T))\right)}"
    )

    st.markdown(
        "여기서 $\\text{pH}_{50,\\text{eff}}$는 추출제 농도와 **온도**에 따라 보정됩니다:"
    )
    st.latex(
        r"\text{pH}_{50,\text{eff}}(T) = \text{pH}_{50,\text{ref}} - \alpha \cdot \log_{10}\!\left(\frac{C_{\text{ext}}}{C_{\text{ref}}}\right) + \beta \cdot (T - T_{\text{ref}})"
    )

    st.markdown("시그모이드 기울기 $k$도 온도에 따라 보정됩니다:")
    st.latex(
        r"k(T) = k_{\text{ref}} \cdot \exp\!\left(\gamma \cdot (T - T_{\text{ref}})\right)"
    )

    if abs(temperature - T_REF) > 0.5:
        st.info(
            f"🌡️ 온도 보정 활성: T = {temperature:.0f}°C (T_ref = {T_REF:.0f}°C, ΔT = {temperature - T_REF:+.0f}°C)"
        )

    st.markdown("#### 금속별 적용 수식")
    for metal in metals:
        params = active_extractant_params[extractant][metal]
        ph50_eff = get_effective_pH50(
            metal,
            extractant,
            C_ext,
            temperature,
            extractant_params=active_extractant_params,
        )
        k_eff = get_effective_k(
            metal,
            extractant,
            temperature,
            extractant_params=active_extractant_params,
        )
        col_eq, col_plot = st.columns([3, 2])

        with col_eq:
            st.markdown(f"**{metal}** ({extractant}, {temperature:.0f}°C):")
            if abs(C_ext - params["C_ref"]) > 0.001 or abs(temperature - T_REF) > 0.5:
                parts = [f"{params['pH50']:.1f}"]
                if abs(C_ext - params["C_ref"]) > 0.001:
                    log_ratio = math.log10(C_ext / params["C_ref"])
                    parts.append(
                        rf"- {params['alpha']:.1f} \times \log_{{10}}\!\left(\frac{{{C_ext}}}{{{params['C_ref']}}}\right)"
                    )
                if abs(temperature - T_REF) > 0.5:
                    beta = params.get("beta", 0.0)
                    parts.append(
                        rf"+ ({beta:.3f}) \times ({temperature - T_REF:+.0f})"
                    )
                st.latex(
                    rf"\text{{pH}}_{{50,\text{{eff}}}}^{{\text{{{metal}}}}} = "
                    + " ".join(parts)
                    + rf" = \mathbf{{{ph50_eff:.2f}}}"
                )
            else:
                st.latex(
                    rf"\text{{pH}}_{{50,\text{{eff}}}}^{{\text{{{metal}}}}} = "
                    rf"{params['pH50']:.1f} \quad (C_{{\text{{ext}}}} = C_{{\text{{ref}}}}, T = T_{{\text{{ref}}}})"
                )

            if abs(temperature - T_REF) > 0.5:
                gamma = params.get("gamma", 0.0)
                st.latex(
                    rf"k^{{\text{{{metal}}}}}(T) = {params['k']:.1f} \times "
                    rf"\exp({gamma:.3f} \times {temperature - T_REF:+.0f}) = \mathbf{{{k_eff:.2f}}}"
                )

            st.latex(
                rf"E_{{\text{{{metal}}}}}(\text{{pH}}) = "
                rf"\frac{{{params['E_max']:.1f}}}"
                rf"{{1 + \exp\!\left(-{k_eff:.2f} \cdot "
                rf"(\text{{pH}} - {ph50_eff:.2f})\right)}}"
            )

            st.markdown("**2층 — Vasilyev 경쟁 보정: 실제 D값 계산**")
            st.latex(
                r"D_{M}^{\text{adj}} = D_{M}^{\text{sig}} \times \left(\frac{[\overline{\text{HA}}]_{\text{free}}}{C_{\text{ext}}}\right)^{n_{\text{eff}}}"
            )
            st.caption(
                "Vasilyev et al. (2019): 위 시그모이드로 구한 기본 D에 잔여 추출제 비율을 곱하여 금속 간 경쟁을 반영합니다."
            )
            pHs = [x * 0.1 for x in range(10, 101)]
            Es = [
                extraction_efficiency(
                    pH,
                    metal,
                    extractant,
                    C_ext,
                    temperature,
                    extractant_params=active_extractant_params,
                )
                for pH in pHs
            ]
            fig_mini = go.Figure()
            fig_mini.add_trace(
                go.Scatter(
                    x=pHs,
                    y=Es,
                    mode="lines",
                    line=dict(color=color_list[metals.index(metal)], width=2),
                )
            )
            if target_pH:
                E_at_target = extraction_efficiency(
                    target_pH,
                    metal,
                    extractant,
                    C_ext,
                    temperature,
                    extractant_params=active_extractant_params,
                )
                fig_mini.add_trace(
                    go.Scatter(
                        x=[target_pH],
                        y=[E_at_target],
                        mode="markers",
                        marker=dict(size=12, color="red", symbol="x"),
                        name=f"pH={target_pH}",
                    )
                )
                fig_mini.add_annotation(
                    x=target_pH,
                    y=E_at_target,
                    text=f"{E_at_target:.1f}%",
                    showarrow=True,
                    arrowhead=2,
                    yshift=15,
                )
            fig_mini.update_layout(
                height=180,
                margin=dict(l=30, r=10, t=10, b=30),
                xaxis_title="pH",
                yaxis_title="E(%)",
                showlegend=False,
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                yaxis=dict(range=[0, 105]),
            )
            st.plotly_chart(fig_mini, use_container_width=True)

        st.markdown("")

    st.markdown("---")
    st.markdown("### 2️⃣ 분배 계수 (Distribution Coefficient)")
    st.latex(r"D_M = \frac{C_{M,\text{org}}}{C_{M,\text{aq}}} = \frac{E_M / 100}{1 - E_M / 100}")

    if target_pH:
        st.markdown(f"**목표 pH = {target_pH}에서의 분배 계수:**")
        d_data = []
        for metal in metals:
            E_val = extraction_efficiency(
                target_pH,
                metal,
                extractant,
                C_ext,
                temperature,
                extractant_params=active_extractant_params,
            )
            D_val = distribution_coefficient(
                target_pH,
                metal,
                extractant,
                C_ext,
                temperature=temperature,
                extractant_params=active_extractant_params,
            )
            d_data.append(
                {
                    "금속": metal,
                    f"E(pH={target_pH})": f"{E_val:.2f}%",
                    f"D(pH={target_pH})": f"{D_val:.4f}",
                }
            )
        st.dataframe(pd.DataFrame(d_data), use_container_width=True, hide_index=True)

    st.markdown("---")
    st.markdown("### 2️⃣-B 추출제 해리 평형 (pKa)")
    pKa_val = EXTRACTANT_PKA.get(extractant)
    if pKa_val is not None:
        st.markdown(
            f"**{extractant}**는 약산 추출제로, 수계와 접촉 시 해리됩니다:"
        )
        st.latex(r"\overline{\text{HA}}_{\text{(org)}} \rightleftharpoons \text{H}^+_{\text{(aq)}} + \text{A}^-_{\text{(aq)}}")
        st.latex(rf"\text{{pKa}} = {pKa_val:.2f}")
        st.markdown(
            "비해리 분율(HA 형태 유지 비율):"
        )
        st.latex(
            r"f_{\text{pKa}} = \frac{1}{1 + 10^{(\text{pH} - \text{pKa})}}"
        )
        if target_pH:
            f_val = pka_dissociation_factor(target_pH, extractant)
            st.info(
                f"📌 목표 pH = {target_pH}에서 **{extractant}**의 비해리(HA) 분율: "
                f"**f_pKa = {f_val:.4f}** ({f_val*100:.1f}%)"
            )
        st.caption(
            f"현재 모델에서는 시그모이드 파라미터(pH50, k)가 OLI/MSE 열역학으로부터 피팅되었기 때문에 "
            f"pKa 효과가 이미 내재되어 있습니다. 이 섹션은 참고용 정보입니다."
        )
    else:
        st.info(f"{extractant}의 pKa 데이터가 등록되어 있지 않습니다.")

    st.markdown("---")
    st.markdown("### 3️⃣ 단일 Stage 물질수지 (Mass Balance)")
    st.markdown("각 Mixer-Settler stage에서 수계(Aqueous)와 유기계(Organic)의 평형:")
    st.latex(r"D_M = \frac{C_{M,\text{org}}^{\text{out}}}{C_{M,\text{aq}}^{\text{out}}}")
    st.latex(
        r"C_{M,\text{aq}}^{\text{in}} \cdot Q_{\text{aq}} + C_{M,\text{org}}^{\text{in}} \cdot Q_{\text{org}}"
        r" = C_{M,\text{aq}}^{\text{out}} \cdot Q_{\text{aq}} + C_{M,\text{org}}^{\text{out}} \cdot Q_{\text{org}}"
    )
    st.markdown("이를 정리하면:")
    st.latex(
        r"C_{M,\text{aq}}^{\text{out}} = "
        r"\frac{C_{M,\text{aq}}^{\text{in}} \cdot Q_{\text{aq}} + C_{M,\text{org}}^{\text{in}} \cdot Q_{\text{org}}}"
        r"{Q_{\text{aq}} + D_M \cdot Q_{\text{org}}}"
    )
    st.latex(r"C_{M,\text{org}}^{\text{out}} = D_M \cdot C_{M,\text{aq}}^{\text{out}}")

    st.markdown(
        f"현재 설정: $Q_{{\\text{{aq}}}}$ = {Q_aq} L/hr, $Q_{{\\text{{org}}}}$ = {Q_org} L/hr, "
        f"A/O = {Q_aq/Q_org:.2f}"
    )

    st.markdown("---")
    st.markdown("### 4️⃣ pH 수지 (Proton Balance)")
    st.markdown("양이온 추출형 추출제의 금속-양성자 교환 반응:")
    st.latex(
        r"M^{n+}_{(\text{aq})} + n(\text{HA})_{2,(\text{org})} \rightarrow M(\text{A}_2)_{n,(\text{org})} + 2n\text{H}^+_{(\text{aq})}"
    )

    st.markdown("**금속별 H⁺ 방출:**")
    for metal in metals:
        n_H = active_extractant_params[extractant][metal]["n_H"]
        charge = "+" if n_H == 1 else f"{n_H}+"
        st.latex(rf"\text{{{metal}}}^{{{charge}}} \rightarrow {n_H} \text{{H}}^+ \text{{ 방출}}")

    st.markdown("**수계 H⁺ 수지:**")
    st.latex(
        r"[\text{H}^+]_{\text{out}} \cdot Q_{\text{aq}} = "
        r"[\text{H}^+]_{\text{in}} \cdot Q_{\text{aq}} + "
        r"\sum_M n_{\text{H},M} \cdot \Delta C_M \cdot Q_{\text{aq}} - "
        r"C_{\text{NaOH}} \cdot Q_{\text{NaOH}}"
    )
    st.latex(r"\text{pH}_{\text{out}} = -\log_{10}\!\left([\text{H}^+]_{\text{out}}\right)")

    st.markdown("**수계 금속 종분화(Speciation) 효과**")
    st.latex(r"[\text{H}^+]_{\text{hydrolysis}} \approx \sum_M K_{MOH} \cdot [\text{M}^{n+}] \cdot [\text{OH}^-]")
    st.caption(
        "고 pH 역외 접근 시 금속 수산화물 착물(MOH⁺ 등)이 생성되며 양성자를 방출(OH⁻ 소비)해 직접적인 완충제로 작용함을 수지에 반영했습니다."
    )

    if pH_mode == "목표 pH (자동 NaOH)":
        st.info(
            f"**목표 pH 모드**: pH = {target_pH if not staged_pHs else staged_pHs} 를 유지하기 위해 필요한 NaOH를 자동 역산합니다.\n\n"
            f"NaOH 필요량 = $[H^+]_{{in}} \\cdot Q_{{aq}} + \\Sigma(n_H \\cdot \\Delta C_M \\cdot Q_{{aq}}) - [H^+]_{{target}} \\cdot Q_{{aq}}$"
        )

    st.markdown("---")
    st.markdown("### 5️⃣ 시스템 입력 보정 로직 (황산염 농도 자동 연산)")
    st.markdown(
        "피드 용액에 포함된 양이온 금속들이 형태상 금속 황산염($M\\text{SO}_4$, 또는 $M_2(\\text{SO}_4)_3$)으로 투입된다고 가정하고 화학양론비에 따라 **초기 총 음이온 농도**를 연산합니다."
    )
    st.latex(
        r"C_{\text{sulfate}} (\text{M}) = \sum_{M} \left( \frac{C_{M,\text{feed}}}{\text{MW}_M} \times \frac{\text{valency}_M}{2} \right)"
    )
    st.info(
        f"계산된 피드 용액 내 황산염 음이온($\\text{{SO}}_4^{{2-}}$) 농도: **{C_sulfate:.4f} M**"
    )

    st.markdown("---")
    st.markdown("### 6️⃣ 다단 연속 추출 순환 해법 (McCabe-Thiele Iterator)")
    st.markdown(
        "**교류 방식 (Counter-current) 반복 알고리즘:** 다단(Multi-stage) Mixer-Settler 시스템의 작동 원리에 따라 수계(Aqueous)는 $\\text{Stage } 1 \\rightarrow N$ 방향으로, 유기계(Organic)는 $\\text{Stage } N \\rightarrow 1$ 방향으로 엇갈리며 흐릅니다."
    )
    st.markdown(
        f"시뮬레이터는 이 동적 흐름 모델에서 {n_stages}개의 각 Stage가 평형에 도달할 때까지 **교대 수렴(Alternating Convergence)** 반복 루프를 구동합니다. 이전 Stage와 다음 Stage에서 넘어오는 물질들을 교차 연산하여 물질수지가 모두 일치할 때까지 수십 차례 왕복 연산합니다."
    )
    st.latex(r"C_{M,\text{org}}^{\text{stage } i} = f_{\text{eq}}\!\left( C_{M,\text{aq}}^{\text{stage } i}, \text{pH}_i \right)")
    st.latex(
        r"Q_{\text{aq}} \cdot C_{M,\text{aq}}^{\text{stage } i} + Q_{\text{org}} \cdot C_{M,\text{org}}^{\text{stage } i} = Q_{\text{aq}} \cdot C_{M,\text{aq}}^{\text{stage } i-1} + Q_{\text{org}} \cdot C_{M,\text{org}}^{\text{stage } i+1}"
    )

    st.markdown("---")
    st.markdown("### 7️⃣ 추출 농도 물리적 한계 (Vasilyev 경쟁 모델: 다핵 착물 보정)")
    st.markdown(
        "금속 로딩량이 치솟아 남은 자유 추출제($\\overline{\\text{HA}}$)가 모자랄 때, $D_M$ 지표가 급격하게 저하되도록 억제하는 Vasilyev 경쟁 모델의 핵심 수식입니다. (v2.0 다핵 착물 $n_{\\text{eff}}$ 완화 적용)"
    )
    st.latex(r"[\overline{\text{HA}}]_{\text{free}} = C_{\text{ext}} - \sum_{M} n_{\text{eff},M} \cdot C_{M,\text{org}}")
    st.latex(
        r"D_{M}^{\text{adj}} = D_{M}^{\text{sig}} \times \left( \max\left(10^{-4}, \frac{[\overline{\text{HA}}]_{\text{free}}}{C_{\text{ext}}}\right) \right)^{n_{\text{eff}}}"
    )
    st.caption(
        "Vasilyev et al. (2019): 고농도 피드 유입 시 모든 금속이 100% 추출되는 비현실적 결과를 방지하고 가용 유기제 내에서 실제적인 상호 조율점(Crowding Out)을 찾습니다."
    )

    st.markdown("---")
    st.markdown("### 8️⃣ 시뮬레이션 계산 흐름 (Step-by-Step)")
    st.markdown(
        "시뮬레이터가 결과를 도출하는 **전체 계산 순서**를 단계별로 설명합니다."
    )
    st.markdown(
        """
**Step 1. 입력 준비**
- 사용자가 설정한 피드 농도, pH, 유량, 추출제 종류/농도, 온도, 단수(Stage 수) 등을 읽어옵니다.
- 피드 용액의 금속 황산염 조성에서 SO₄²⁻ 농도를 자동 계산합니다.

**Step 2. 각 Stage에서의 평형 계산 (핵심)**

각 Stage에서 다음 순서로 금속별 분배를 계산합니다:

> **① 시그모이드로 기본 D 계산 (OLI/MSE 파라미터)**
> - 각 금속의 pH50, k, E_max를 온도와 추출제 농도로 보정
> - 보정된 시그모이드로 기본 분배계수 $D_{sig}$ 계산
>
> **② 잔여 추출제 풀 계산 (Vasilyev 모델)**
> - 유기상에 이미 로딩된 금속이 소비한 추출제를 차감
> - $[\overline{HA}]_{free} = C_{ext} - \sum n_{eff,M} \cdot C_{M,org}$
>
> **③ pH₅₀ 우선순위로 금속 순차 처리**
> - pH₅₀가 낮은 금속(추출 잘 되는 금속)부터 먼저 처리
> - 예: Zn → Mn → Ca → Mg → Co → Ni → Li 순서
>
> **④ 경쟁 보정된 D값으로 물질수지 계산**
> - $D_{adj} = D_{sig} \times (HA_{free}/C_{ext})^{n_{eff}}$
> - 이 D값으로 수계/유기계 출구 농도를 결정
> - 해당 금속이 소비한 추출제를 풀에서 차감 → 다음 금속에 영향

**Step 3. pH 수지 계산**
- 금속 추출 시 방출된 H⁺ 총량 계산
- 수계 종분화(MOH⁺) 효과, 황산 버퍼(HSO₄⁻) 효과 반영
- NaOH 중화량 차감 → 출구 pH 결정

**Step 4. 다단 역류 수렴 반복**
- Stage 1→N 방향으로 수계, Stage N→1 방향으로 유기계가 역류
- 각 Stage의 유기상 출구 농도가 수렴할 때까지 반복 (수십~수백회)
- Under-relaxation으로 안정적 수렴 유도

**Step 5. 결과 출력**
- 후액(Raffinate) 농도, 로딩 유기(Loaded Organic) 농도
- Stage별 pH 프로파일, NaOH 소비량
- 전체 추출률, 분리 계수 등 성능 지표
"""
    )

    st.markdown("---")
    st.markdown("### 9️⃣ 성능 지표")
    col_perf1, col_perf2 = st.columns(2)
    with col_perf1:
        st.markdown("**추출률 (Percent Extraction):**")
        st.latex(r"E_M = \frac{C_{M,\text{feed}} - C_{M,\text{raffinate}}}{C_{M,\text{feed}}} \times 100\%")
        st.markdown("**분리 계수 (Separation Factor):**")
        st.latex(r"\alpha_{M_1, M_2} = \frac{D_{M_1}}{D_{M_2}}")
    with col_perf2:
        if target_pH:
            st.markdown(f"**pH = {target_pH}에서 분리 계수:**")
            sep_data = []
            D_vals = {
                metal: distribution_coefficient(
                    target_pH,
                    metal,
                    extractant,
                    C_ext,
                    temperature=temperature,
                    extractant_params=active_extractant_params,
                )
                for metal in metals
            }
            for i, metal_1 in enumerate(metals):
                for metal_2 in metals[i + 1 :]:
                    if D_vals[metal_1] <= 1e-10 and D_vals[metal_2] <= 1e-10:
                        alpha_str = "N/A"
                    elif D_vals[metal_2] <= 1e-10:
                        alpha_str = "> 10⁶"
                    else:
                        alpha_str = f"{D_vals[metal_1] / D_vals[metal_2]:.2f}"
                    sep_data.append({"M₁/M₂": f"{metal_1}/{metal_2}", "α": alpha_str})
            st.dataframe(pd.DataFrame(sep_data), use_container_width=True, hide_index=True)

    st.markdown("---")
    st.markdown("### 📋 현재 파라미터 요약")
    param_rows = []
    for metal in metals:
        params = active_extractant_params[extractant][metal]
        ph50_eff = get_effective_pH50(
            metal,
            extractant,
            C_ext,
            temperature,
            extractant_params=active_extractant_params,
        )
        k_eff = get_effective_k(
            metal,
            extractant,
            temperature,
            extractant_params=active_extractant_params,
        )
        param_rows.append(
            {
                "금속": metal,
                "pH₅₀ (기준)": params["pH50"],
                "α": params["alpha"],
                "C_ref (M)": params["C_ref"],
                "β": params.get("beta", 0.0),
                "γ": params.get("gamma", 0.0),
                f"pH₅₀ ({temperature:.0f}°C)": round(ph50_eff, 3),
                f"k ({temperature:.0f}°C)": round(k_eff, 3),
                "E_max (%)": params["E_max"],
                "n_H": params["n_H"],
            }
        )
    st.dataframe(pd.DataFrame(param_rows), use_container_width=True, hide_index=True)


def render_fitting_tab(extractant: str, default_metals: list) -> None:
    st.subheader("📝 실험 데이터 피팅")
    st.markdown(
        """
    실험 데이터(pH vs 추출률)를 업로드하여 시그모이드 모델 파라미터를 자동으로 피팅합니다.
    피팅된 파라미터는 시뮬레이션에 적용할 수 있습니다.
    """
    )

    with st.expander("📥 샘플 데이터 다운로드"):
        st.markdown(
            """
        아래 버튼을 클릭하여 샘플 CSV 파일을 다운로드하세요. 이 파일을 참고하여 데이터를 준비하시면 됩니다.

        **필수 열**: `pH`, `E_pct` (추출률 %)
        **선택 열**: `metal` (금속명 - 복수 금속 동시 피팅 시)
        """
        )
        sample_csv = (
            "pH,E_pct,metal\n2.0,0.5,Co\n2.5,2.0,Co\n3.0,12.0,Co\n3.5,45.0,Co\n"
            "4.0,82.0,Co\n4.5,96.0,Co\n5.0,99.5,Co\n2.0,0.1,Ni\n3.0,0.5,Ni\n"
            "4.0,3.0,Ni\n5.0,25.0,Ni\n5.5,60.0,Ni\n6.0,90.0,Ni\n6.5,98.0,Ni"
        )
        st.download_button(
            "📥 샘플 CSV 다운로드",
            sample_csv,
            "sample_extraction_data.csv",
            "text/csv",
        )

    uploaded = st.file_uploader("📂 CSV 파일 업로드 (pH, E_pct 열 필수)", type=["csv"])

    if uploaded is None:
        st.info("👆 CSV 파일을 업로드하여 피팅을 시작하세요.")
        return

    try:
        df_raw = pd.read_csv(uploaded)
    except Exception as exc:
        st.error(f"CSV 파일 읽기 오류: {exc}")
        return

    required_cols = {"pH", "E_pct"}
    if not required_cols.issubset(set(df_raw.columns)):
        st.error(f"CSV에 필수 열이 없습니다: {required_cols}. 현재 열: {list(df_raw.columns)}")
        return

    st.markdown("#### 📊 업로드된 데이터")
    st.dataframe(df_raw.head(20), use_container_width=True, hide_index=True)

    if "metal" in df_raw.columns:
        metal_list = sorted(df_raw["metal"].unique().tolist())
    else:
        metal_list = ["Unknown"]
        df_raw["metal"] = "Unknown"

    st.markdown(f"감지된 금속: **{', '.join(metal_list)}** ({len(df_raw)}개 데이터 포인트)")

    st.markdown("---")
    st.markdown("#### ⚙️ 피팅 설정")
    fit_model = "Sigmoid"

    if not st.button("🚀 피팅 실행", type="primary", use_container_width=True):
        return

    st.markdown("---")
    st.markdown("#### 📈 피팅 결과")

    for metal_name in metal_list:
        df_m = df_raw[df_raw["metal"] == metal_name]
        pH_vals = df_m["pH"].values
        E_vals = df_m["E_pct"].values

        if len(pH_vals) < 3:
            st.warning(f"{metal_name}: 데이터 포인트가 3개 미만입니다.")
            continue

        st.markdown(f"##### 🔹 {metal_name}")

        fit_result = fit_sigmoid(pH_vals, E_vals)
        if not fit_result["success"]:
            st.error(f"피팅 실패: {fit_result.get('error', 'Unknown error')}")
            continue

        res_col1, res_col2 = st.columns([1, 1])
        params = fit_result["params"]
        errors = fit_result["errors"]

        with res_col1:
            st.markdown("**피팅 파라미터:**")
            param_df = pd.DataFrame(
                [
                    {"파라미터": "pH₅₀", "값": params["pH50"], "±95%CI": errors["pH50_err"]},
                    {"파라미터": "k", "값": params["k"], "±95%CI": errors["k_err"]},
                    {"파라미터": "E_max (%)", "값": params["E_max"], "±95%CI": errors["E_max_err"]},
                ]
            )
            st.dataframe(param_df, use_container_width=True, hide_index=True)
            st.metric("R²", f"{fit_result['r_squared']:.4f}")

        with res_col2:
            fig_fit = go.Figure()
            fig_fit.add_trace(
                go.Scatter(
                    x=pH_vals,
                    y=E_vals,
                    mode="markers",
                    name="실험 데이터",
                    marker=dict(size=10, color="#e94560"),
                )
            )

            pH_fine = np.linspace(float(min(pH_vals)) - 0.5, float(max(pH_vals)) + 0.5, 200)
            E_fit = sigmoid_model(pH_fine, params["pH50"], params["k"], params["E_max"])
            fig_fit.add_trace(
                go.Scatter(
                    x=pH_fine,
                    y=E_fit,
                    mode="lines",
                    name="피팅 곡선",
                    line=dict(color="#0f3460", width=2),
                )
            )
            fig_fit.update_yaxes(title="추출률 (%)")
            fig_fit.update_layout(
                title=f"{metal_name} 피팅 결과 ({fit_model})",
                xaxis_title="pH",
                height=350,
                template="plotly_white",
            )
            st.plotly_chart(fig_fit, use_container_width=True)

        if metal_name in default_metals:
            if st.button(f"✅ {metal_name} 피팅 파라미터 적용", key=f"apply_{metal_name}"):
                st.session_state.custom_params[extractant][metal_name]["pH50"] = params["pH50"]
                st.session_state.custom_params[extractant][metal_name]["k"] = params["k"]
                st.session_state.custom_params[extractant][metal_name]["E_max"] = params["E_max"]
                st.success(f"{metal_name} 파라미터가 적용되었습니다! 페이지를 새로고침하면 시뮬레이션에 반영됩니다.")
                st.rerun()

        st.markdown("---")
