"""
Dashboard Tab Renderers
=======================
Streamlit 대시보드의 주요 동적 탭 렌더링 로직을 분리합니다.

목표:
- sx_dashboard.py의 길이를 줄이고 책임을 분리
- 시뮬레이션 계산과 UI 렌더링을 느슨하게 연결
- 탭 단위 수정 시 충돌 범위 축소
"""

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
    loading_damping_factor,
)


def _metal_color_list(metals: list, metal_colors: dict) -> list:
    return [metal_colors.get(metal, "#636e72") for metal in metals]


def render_results_tab(
    result: dict,
    metals: list,
    loading_pct: float,
    scope_assessment: dict,
    profile_label: str,
    selected_preset_note: str | None,
    C_aq_feed: dict,
    pH_feed: float,
    n_stages: int,
    metal_colors: dict,
    scope_renderer,
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

    col_info1, col_info2, col_info3, col_info4 = st.columns(4)
    with col_info1:
        st.metric("후액 최종 pH", f"{result['pH_profile'][-1]:.2f}")
    with col_info2:
        st.metric("총 NaOH 소비", f"{result['total_NaOH_mol_hr']:.1f} mol/hr")
    with col_info3:
        if loading_pct > 85:
            st.metric(
                "최대 로딩률",
                f"{loading_pct:.1f}%",
                delta="⚠️ 포화 근접",
                delta_color="inverse",
            )
        else:
            st.metric("최대 로딩률", f"{loading_pct:.1f}%")
    with col_info4:
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
    target_pH: float | None,
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
        D = D * loading_damping_factor(loading)
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
