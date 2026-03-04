#!/usr/bin/env python3
"""
SX Simulator Web Dashboard
===========================
Streamlit 기반 Li/Ni/Co/Mn Mixer-Settler 용매추출 시뮬레이션 대시보드.

실행: streamlit run sx_dashboard.py
"""

import sys, os, math, copy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

APP_VERSION = "1.3.0"

# CHANGELOG 읽기
_changelog_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "CHANGELOG.md")
try:
    with open(_changelog_path, "r", encoding="utf-8") as _f:
        CHANGELOG_CONTENT = _f.read()
except FileNotFoundError:
    CHANGELOG_CONTENT = "CHANGELOG.md 파일을 찾을 수 없습니다."

import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

from sx_simulator.extraction_isotherm import (
    extraction_efficiency, compute_all_extractions, get_effective_pH50,
    get_effective_k, distribution_coefficient, calc_loading_fraction,
)
from sx_simulator.multistage_sx import solve_multistage_countercurrent
from sx_simulator.config import EXTRACTANT_PARAMS, MOLAR_MASS, DEFAULT_METALS, T_REF, MAX_LOADING_FRACTION

# =============================================================================
# 페이지 설정
# =============================================================================
st.set_page_config(
    page_title="SX Simulator",
    page_icon="⚗️",
    layout="wide",
    initial_sidebar_state="expanded",
)

# 커스텀 CSS
st.markdown("""
<style>
    .stMetric { background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
                padding: 15px; border-radius: 10px; border: 1px solid #0f3460; }
    .stMetric label { color: #e94560 !important; font-weight: 600; }
    .stMetric [data-testid="stMetricValue"] { color: #ffffff !important; }
    div[data-testid="stSidebar"] { background: linear-gradient(180deg, #0f0c29 0%, #302b63 50%, #24243e 100%); }
    div[data-testid="stSidebar"] label { color: #e0e0e0 !important; }
    h1 { color: #e94560; }
    h2, h3 { color: #0f3460; }
</style>
""", unsafe_allow_html=True)

# =============================================================================
# 사이드바: 입력 파라미터
# =============================================================================
st.sidebar.title("⚗️ SX 시뮬레이터")
st.sidebar.markdown("---")

# --- Feed 조건 ---
st.sidebar.header("📥 Feed (수계) 조건")
col_f1, col_f2 = st.sidebar.columns(2)
with col_f1:
    C_Li = st.number_input("Li (g/L)", 0.0, 50.0, 1.5, 0.1, key="li")
    C_Co = st.number_input("Co (g/L)", 0.0, 50.0, 3.0, 0.1, key="co")
with col_f2:
    C_Ni = st.number_input("Ni (g/L)", 0.0, 50.0, 5.0, 0.1, key="ni")
    C_Mn = st.number_input("Mn (g/L)", 0.0, 50.0, 2.0, 0.1, key="mn")

pH_feed = st.sidebar.slider("피드 pH", 0.5, 7.0, 3.0, 0.1)
Q_aq = st.sidebar.number_input("수계 유량 (L/hr)", 1.0, 10000.0, 100.0, 10.0)

st.sidebar.markdown("---")

# --- NaOH 조건 ---
st.sidebar.header("🧪 pH 제어 (NaOH)")
pH_mode = st.sidebar.radio("pH 제어 모드", ["목표 pH (자동 NaOH)", "고정 NaOH"], index=0)

if pH_mode == "목표 pH (자동 NaOH)":
    target_pH = st.sidebar.slider("목표 pH", 2.0, 8.0, 5.0, 0.1)
    use_staged_pH = st.sidebar.checkbox("Stage별 차등 pH", value=False)
else:
    target_pH = None
    C_NaOH = st.sidebar.number_input("NaOH 농도 (M)", 0.1, 20.0, 5.0, 0.5)
    Q_NaOH = st.sidebar.number_input("NaOH 유량 (L/hr)", 0.0, 500.0, 12.0, 1.0)

st.sidebar.markdown("---")

# --- 용매 조건 ---
st.sidebar.header("🫗 유기계 (용매) 조건")
extractant = st.sidebar.selectbox("추출제", ["Cyanex 272", "D2EHPA"])
C_ext = st.sidebar.number_input("추출제 농도 (M)", 0.05, 5.0, 0.5, 0.05)
Q_org = st.sidebar.number_input("유기계 유량 (L/hr)", 1.0, 10000.0, 100.0, 10.0)

# --- 온도 설정 ---
st.sidebar.markdown("---")
st.sidebar.header("🌡️ 온도 설정")
temperature = st.sidebar.slider("운전 온도 (°C)", 10.0, 60.0, 25.0, 1.0)
if abs(temperature - T_REF) > 0.5:
    st.sidebar.info(f"온도 보정 활성: ΔT = {temperature - T_REF:+.0f}°C")

st.sidebar.markdown("---")

# --- Stage 설정 ---
st.sidebar.header("🔄 Mixer-Settler 설정")
n_stages = st.sidebar.slider("Stage 수", 1, 10, 4)

# Stage별 pH 입력 (차등 모드)
staged_pHs = None
if pH_mode == "목표 pH (자동 NaOH)" and use_staged_pH:
    st.sidebar.markdown("**Stage별 목표 pH:**")
    staged_pHs = []
    for i in range(n_stages):
        pH_val = st.sidebar.slider(f"Stage {i+1} pH", 2.0, 8.0, 3.0 + i * (5.0 / n_stages), 0.1, key=f"staged_ph_{i}")
        staged_pHs.append(pH_val)

st.sidebar.markdown("---")

# --- 모델 파라미터 편집 ---
st.sidebar.header("📐 모델 파라미터 편집")
edit_params = st.sidebar.checkbox("Isotherm 파라미터 수정", value=False)

# 파라미터를 session_state에 저장
if "custom_params" not in st.session_state:
    st.session_state.custom_params = copy.deepcopy(EXTRACTANT_PARAMS)

if edit_params:
    st.sidebar.markdown(f"**{extractant} 파라미터:**")
    for metal in DEFAULT_METALS:
        with st.sidebar.expander(f"🔹 {metal}"):
            p = st.session_state.custom_params[extractant][metal]
            p["pH50"] = st.number_input(f"pH₅₀ ({metal})", 0.0, 14.0, p["pH50"], 0.1, key=f"pH50_{metal}")
            p["k"] = st.number_input(f"k ({metal})", 0.1, 10.0, p["k"], 0.1, key=f"k_{metal}")
            p["E_max"] = st.number_input(f"E_max ({metal})", 50.0, 100.0, p["E_max"], 0.5, key=f"Emax_{metal}")
            p["n_H"] = st.number_input(f"n_H ({metal})", 1, 4, p["n_H"], 1, key=f"nH_{metal}")
            p["alpha"] = st.number_input(f"α ({metal})", 0.0, 3.0, p["alpha"], 0.1, key=f"alpha_{metal}")
            p["beta"] = st.number_input(f"β 온도계수 ({metal})", -0.2, 0.1, float(p.get("beta", 0.0)), 0.005, key=f"beta_{metal}", format="%.3f")
            p["gamma"] = st.number_input(f"γ k-온도계수 ({metal})", -0.05, 0.05, float(p.get("gamma", 0.0)), 0.001, key=f"gamma_{metal}", format="%.3f")

    # 파라미터를 전역에 반영
    for ext_name in st.session_state.custom_params:
        for metal in st.session_state.custom_params[ext_name]:
            for k, v in st.session_state.custom_params[ext_name][metal].items():
                EXTRACTANT_PARAMS[ext_name][metal][k] = v

    if st.sidebar.button("🔄 기본값으로 초기화"):
        from sx_simulator import config as _cfg
        import importlib
        importlib.reload(_cfg)
        st.session_state.custom_params = copy.deepcopy(_cfg.EXTRACTANT_PARAMS)
        st.rerun()

# =============================================================================
# 메인 영역
# =============================================================================
st.title("⚗️ Li/Ni/Co/Mn Mixer-Settler SX 시뮬레이터")
st.caption("MSE Thermodynamic Framework 기반 (ALTA 2024 / Wang et al.)")

metals = DEFAULT_METALS
C_aq_feed = {"Li": C_Li, "Ni": C_Ni, "Co": C_Co, "Mn": C_Mn}

# =============================================================================
# 탭 구성
# =============================================================================
tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "📊 시뮬레이션 결과", "📐 모델 수식", "📈 pH Isotherm",
    "🔬 추출제 비교", "📋 상세 데이터", "📜 변경 이력"
])

# =============================================================================
# 시뮬레이션 실행
# =============================================================================
with st.spinner("시뮬레이션 계산 중..."):
    sim_kwargs = dict(
        C_aq_feed=C_aq_feed, pH_feed=pH_feed,
        Q_aq=Q_aq, Q_org=Q_org,
        extractant=extractant, C_ext=C_ext,
        n_stages=n_stages, metals=metals,
        temperature=temperature,
    )

    if pH_mode == "목표 pH (자동 NaOH)":
        if staged_pHs:
            sim_kwargs["target_pH_per_stage"] = staged_pHs
        else:
            sim_kwargs["target_pH"] = target_pH
    else:
        sim_kwargs["C_NaOH"] = C_NaOH
        sim_kwargs["Q_NaOH"] = Q_NaOH

    result = solve_multistage_countercurrent(**sim_kwargs)

# =============================================================================
# TAB 1: 시뮬레이션 결과
# =============================================================================
with tab1:
    # --- 요약 메트릭 ---
    st.subheader("🎯 추출 결과 요약")
    met_cols = st.columns(4)
    colors = {"Li": "#ff6b6b", "Ni": "#4ecdc4", "Co": "#45b7d1", "Mn": "#f7b731"}
    for i, metal in enumerate(metals):
        ext_pct = result["overall_extraction"][metal]
        raff = result["raffinate"][metal]
        met_cols[i].metric(
            label=f"{metal} 추출률",
            value=f"{ext_pct:.1f}%",
            delta=f"후액: {raff:.3f} g/L",
        )

    col_info1, col_info2, col_info3 = st.columns(3)
    with col_info1:
        st.metric("후액 최종 pH", f"{result['pH_profile'][-1]:.2f}")
    with col_info2:
        st.metric("총 NaOH 소비", f"{result['total_NaOH_mol_hr']:.1f} mol/hr")
    with col_info3:
        # Stage 1 (가장 높은 로딩) 로딩률 표시
        max_loading = max(sr.get('loading_fraction', 0.0) for sr in result['stages'])
        loading_pct = max_loading * 100
        if loading_pct > 85:
            st.metric("최대 로딩률", f"{loading_pct:.1f}%", delta="⚠️ 포화 근접", delta_color="inverse")
        else:
            st.metric("최대 로딩률", f"{loading_pct:.1f}%")

    st.markdown("---")

    # --- Stage별 추출률 바 차트 ---
    st.subheader("📊 금속별 추출률")
    ext_data = pd.DataFrame({
        "금속": metals,
        "추출률 (%)": [result["overall_extraction"][m] for m in metals],
        "피드 (g/L)": [C_aq_feed[m] for m in metals],
        "후액 (g/L)": [result["raffinate"][m] for m in metals],
    })

    fig_ext = px.bar(ext_data, x="금속", y="추출률 (%)",
                     color="금속", color_discrete_sequence=["#ff6b6b", "#4ecdc4", "#45b7d1", "#f7b731"],
                     text="추출률 (%)")
    fig_ext.update_traces(texttemplate='%{text:.1f}%', textposition='outside')
    fig_ext.update_layout(height=400, showlegend=False,
                          plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
                          yaxis=dict(range=[0, 110]))
    st.plotly_chart(fig_ext, use_container_width=True)

    st.markdown("---")

    # --- pH 프로파일 + NaOH 프로파일 ---
    col_ph, col_naoh = st.columns(2)

    with col_ph:
        st.subheader("📉 Stage별 pH 프로파일")
        stages_x = list(range(1, n_stages + 1))
        fig_ph = go.Figure()
        fig_ph.add_trace(go.Scatter(
            x=stages_x, y=result["pH_profile"],
            mode='lines+markers', name='pH',
            line=dict(color='#e94560', width=3),
            marker=dict(size=12, symbol='circle'),
        ))
        fig_ph.update_layout(
            xaxis_title="Stage", yaxis_title="pH",
            height=350, plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
            xaxis=dict(dtick=1),
        )
        st.plotly_chart(fig_ph, use_container_width=True)

    with col_naoh:
        st.subheader("🧪 Stage별 NaOH 소비량")
        fig_naoh = go.Figure()
        fig_naoh.add_trace(go.Bar(
            x=stages_x, y=result["NaOH_profile"],
            marker_color='#0f3460', text=[f"{v:.1f}" for v in result["NaOH_profile"]],
            textposition='outside',
        ))
        fig_naoh.update_layout(
            xaxis_title="Stage", yaxis_title="NaOH (mol/hr)",
            height=350, plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
            xaxis=dict(dtick=1),
        )
        st.plotly_chart(fig_naoh, use_container_width=True)

    # --- Stage별 금속 농도 변화 (수계) ---
    st.subheader("🔬 Stage별 수계 금속 농도 변화")
    conc_data = {"Stage": ["Feed"] + [f"Stage {i+1}" for i in range(n_stages)]}
    for metal in metals:
        conc_data[metal] = [C_aq_feed[metal]]
        for sr in result["stages"]:
            conc_data[metal].append(sr["C_aq_out"][metal])

    df_conc = pd.DataFrame(conc_data)
    fig_conc = go.Figure()
    color_list = ["#ff6b6b", "#4ecdc4", "#45b7d1", "#f7b731"]
    for i, metal in enumerate(metals):
        fig_conc.add_trace(go.Scatter(
            x=df_conc["Stage"], y=df_conc[metal],
            mode='lines+markers', name=metal,
            line=dict(color=color_list[i], width=2),
            marker=dict(size=8),
        ))
    fig_conc.update_layout(
        xaxis_title="Stage", yaxis_title="농도 (g/L)",
        height=400, plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    st.plotly_chart(fig_conc, use_container_width=True)

# =============================================================================
# TAB 2: pH Isotherm 곡선
# =============================================================================
with tab3:
    st.subheader(f"📈 pH-추출률 Isotherm ({extractant}, {C_ext}M, {temperature:.0f}°C)")

    fig_iso = go.Figure()
    pH_range = [x * 0.1 for x in range(10, 101)]
    for i, metal in enumerate(metals):
        E_values = [extraction_efficiency(pH, metal, extractant, C_ext, temperature) for pH in pH_range]
        fig_iso.add_trace(go.Scatter(
            x=pH_range, y=E_values,
            mode='lines', name=metal,
            line=dict(color=color_list[i], width=3),
        ))
    fig_iso.update_layout(
        xaxis_title="pH", yaxis_title="추출률 (%)",
        height=500, plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        yaxis=dict(range=[0, 105]),
    )
    # 현재 pH 표시
    if target_pH:
        fig_iso.add_vline(x=target_pH, line_dash="dash", line_color="red",
                          annotation_text=f"목표 pH = {target_pH}")

    st.plotly_chart(fig_iso, use_container_width=True)

    # pH50 테이블
    st.subheader("📋 금속별 pH₅₀ (보정값)")
    ph50_data = []
    for metal in metals:
        ph50_eff = get_effective_pH50(metal, extractant, C_ext, temperature)
        k_eff = get_effective_k(metal, extractant, temperature)
        params = EXTRACTANT_PARAMS[extractant][metal]
        ph50_data.append({
            "금속": metal, "pH₅₀ (기준)": params["pH50"],
            "pH₅₀ (보정)": round(ph50_eff, 2), "k (기준)": params["k"],
            "k (보정)": round(k_eff, 2),
            "E_max (%)": params["E_max"], "n_H": params["n_H"],
        })
    st.dataframe(pd.DataFrame(ph50_data), use_container_width=True, hide_index=True)

# =============================================================================
# TAB 3: 추출제 비교
# =============================================================================
with tab4:
    st.subheader("🔬 Cyanex 272 vs D2EHPA 비교")
    compare_pH = st.slider("비교 목표 pH", 2.0, 8.0, 5.0, 0.1, key="compare_ph")

    results_compare = {}
    for ext in ["Cyanex 272", "D2EHPA"]:
        ext_conc = 0.5 if ext == "Cyanex 272" else 0.64
        r = solve_multistage_countercurrent(
            C_aq_feed=C_aq_feed, pH_feed=pH_feed,
            Q_aq=Q_aq, Q_org=Q_org,
            extractant=ext, C_ext=ext_conc, n_stages=n_stages,
            target_pH=compare_pH, metals=metals,
            temperature=temperature,
        )
        results_compare[ext] = r

    # 비교 바 차트
    compare_data = []
    for ext in ["Cyanex 272", "D2EHPA"]:
        for m in metals:
            compare_data.append({
                "추출제": ext, "금속": m,
                "추출률 (%)": results_compare[ext]["overall_extraction"][m],
            })

    df_compare = pd.DataFrame(compare_data)
    fig_compare = px.bar(df_compare, x="금속", y="추출률 (%)", color="추출제",
                         barmode="group", color_discrete_sequence=["#e94560", "#0f3460"],
                         text="추출률 (%)")
    fig_compare.update_traces(texttemplate='%{text:.1f}%', textposition='outside')
    fig_compare.update_layout(height=450, yaxis=dict(range=[0, 115]),
                              plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
    st.plotly_chart(fig_compare, use_container_width=True)

    # NaOH 소비 비교
    col_c1, col_c2 = st.columns(2)
    with col_c1:
        st.metric("Cyanex 272 NaOH", f"{results_compare['Cyanex 272']['total_NaOH_mol_hr']:.1f} mol/hr")
    with col_c2:
        st.metric("D2EHPA NaOH", f"{results_compare['D2EHPA']['total_NaOH_mol_hr']:.1f} mol/hr")

    # Isotherm 비교
    st.markdown("---")
    st.subheader("📈 Isotherm 비교")
    selected_metal = st.selectbox("금속 선택", metals, index=2)

    fig_iso_cmp = go.Figure()
    for ext, color in [("Cyanex 272", "#e94560"), ("D2EHPA", "#0f3460")]:
        ext_c = 0.5 if ext == "Cyanex 272" else 0.64
        E_vals = [extraction_efficiency(pH, selected_metal, ext, ext_c, temperature) for pH in pH_range]
        fig_iso_cmp.add_trace(go.Scatter(
            x=pH_range, y=E_vals, mode='lines', name=ext,
            line=dict(color=color, width=3),
        ))
    fig_iso_cmp.add_vline(x=compare_pH, line_dash="dash", line_color="gray",
                           annotation_text=f"비교 pH = {compare_pH}")
    fig_iso_cmp.update_layout(
        title=f"{selected_metal} 추출 Isotherm 비교",
        xaxis_title="pH", yaxis_title="추출률 (%)",
        height=400, plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
    )
    st.plotly_chart(fig_iso_cmp, use_container_width=True)

# =============================================================================
# TAB 4: 상세 데이터
# =============================================================================
with tab5:
    st.subheader("📋 시뮬레이션 상세 결과")

    # 수렴 정보
    st.info(f"수렴: {'✅ Yes' if result['converged'] else '❌ No'} | "
            f"반복: {result['iterations']} | "
            f"후액 pH: {result['pH_profile'][-1]:.2f} | "
            f"NaOH 총: {result['total_NaOH_mol_hr']:.1f} mol/hr")

    # Stage별 상세 테이블
    for i, sr in enumerate(result["stages"]):
        with st.expander(f"Stage {i+1} (pH = {sr['pH_out']:.2f}, NaOH = {sr.get('NaOH_consumed_mol_hr', 0):.2f} mol/hr)"):
            rows = []
            for m in metals:
                rows.append({
                    "금속": m,
                    "수계 입구 (g/L)": C_aq_feed[m] if i == 0 else result["stages"][i-1]["C_aq_out"][m],
                    "수계 출구 (g/L)": sr["C_aq_out"][m],
                    "유기계 출구 (g/L)": sr["C_org_out"][m],
                    "Stage 추출률 (%)": sr["extraction"][m],
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    # 전체 물질수지 체크
    st.markdown("---")
    st.subheader("⚖️ 물질수지 검증")
    balance_data = []
    for m in metals:
        feed_flow = C_aq_feed[m] * Q_aq
        raff_flow = result["raffinate"][m] * Q_aq
        org_flow = result["loaded_organic"][m] * Q_org
        balance = feed_flow - raff_flow - org_flow
        balance_data.append({
            "금속": m,
            "피드 (g/hr)": round(feed_flow, 3),
            "후액 (g/hr)": round(raff_flow, 3),
            "유기계 (g/hr)": round(org_flow, 3),
            "잔차 (g/hr)": round(balance, 3),
        })
    st.dataframe(pd.DataFrame(balance_data), use_container_width=True, hide_index=True)


# =============================================================================
# TAB 5 (tab5): 모델 수식 확인
# =============================================================================
with tab2:
    st.subheader(f"📐 현재 모델 수식 ({extractant}, {C_ext}M, {temperature:.0f}°C)")
    st.markdown("현재 사이드바에 설정된 파라미터 값이 적용된 수식입니다. "
                "파라미터를 변경하면 수식도 실시간으로 업데이트됩니다.")

    # ---------------------------------------------------------------
    # 1. 추출률 Isotherm 수식
    # ---------------------------------------------------------------
    st.markdown("---")
    st.markdown("### 1️⃣ 추출률 (Extraction Efficiency)")
    st.markdown("각 금속의 pH에 따른 추출률을 **시그모이드(Sigmoid) 함수**로 근사합니다:")
    st.latex(r"E_M(\text{pH}, T) = \frac{E_{\max}}{1 + \exp\!\left(-k(T) \cdot (\text{pH} - \text{pH}_{50,\text{eff}}(T))\right)}")

    st.markdown("여기서 $\\text{pH}_{50,\\text{eff}}$는 추출제 농도와 **온도**에 따라 보정됩니다:")
    st.latex(r"\text{pH}_{50,\text{eff}}(T) = \text{pH}_{50,\text{ref}} - \alpha \cdot \log_{10}\!\left(\frac{C_{\text{ext}}}{C_{\text{ref}}}\right) + \beta \cdot (T - T_{\text{ref}})")

    st.markdown("시그모이드 기울기 $k$도 온도에 따라 보정됩니다:")
    st.latex(r"k(T) = k_{\text{ref}} \cdot \exp\!\left(\gamma \cdot (T - T_{\text{ref}})\right)")

    if abs(temperature - T_REF) > 0.5:
        st.info(f"🌡️ 온도 보정 활성: T = {temperature:.0f}°C (T_ref = {T_REF:.0f}°C, ΔT = {temperature - T_REF:+.0f}°C)")

    # 금속별 실제 수식 전개
    st.markdown("#### 금속별 적용 수식")
    for metal in metals:
        params = EXTRACTANT_PARAMS[extractant][metal]
        ph50_eff = get_effective_pH50(metal, extractant, C_ext, temperature)
        k_eff = get_effective_k(metal, extractant, temperature)
        col_eq, col_plot = st.columns([3, 2])

        with col_eq:
            st.markdown(f"**{metal}** ({extractant}, {temperature:.0f}°C):")
            # 보정 pH50 계산 과정
            if abs(C_ext - params['C_ref']) > 0.001 or abs(temperature - T_REF) > 0.5:
                parts = [f"{params['pH50']:.1f}"]
                if abs(C_ext - params['C_ref']) > 0.001:
                    log_ratio = math.log10(C_ext / params['C_ref'])
                    parts.append(rf"- {params['alpha']:.1f} \times \log_{{10}}\!\left(\frac{{{C_ext}}}{{{params['C_ref']}}}\right)")
                if abs(temperature - T_REF) > 0.5:
                    beta = params.get('beta', 0.0)
                    parts.append(rf"+ ({beta:.3f}) \times ({temperature - T_REF:+.0f})")
                st.latex(
                    rf"\text{{pH}}_{{50,\text{{eff}}}}^{{\text{{{metal}}}}} = "
                    + " ".join(parts) + rf" = \mathbf{{{ph50_eff:.2f}}}"
                )
            else:
                st.latex(
                    rf"\text{{pH}}_{{50,\text{{eff}}}}^{{\text{{{metal}}}}} = "
                    rf"{params['pH50']:.1f} \quad (C_{{\text{{ext}}}} = C_{{\text{{ref}}}}, T = T_{{\text{{ref}}}})"
                )

            # k 온도 보정
            if abs(temperature - T_REF) > 0.5:
                gamma = params.get('gamma', 0.0)
                st.latex(
                    rf"k^{{\text{{{metal}}}}}(T) = {params['k']:.1f} \times "
                    rf"\exp({gamma:.3f} \times {temperature - T_REF:+.0f}) = \mathbf{{{k_eff:.2f}}}"
                )

            # 최종 추출률 수식
            st.latex(
                rf"E_{{\text{{{metal}}}}}(\text{{pH}}) = "
                rf"\frac{{{params['E_max']:.1f}}}"
                rf"{{1 + \exp\!\left(-{k_eff:.2f} \cdot "
                rf"(\text{{pH}} - {ph50_eff:.2f})\right)}}"
            )

        with col_plot:
            # 해당 금속의 미니 Isotherm 차트
            pHs = [x * 0.1 for x in range(10, 101)]
            Es = [extraction_efficiency(pH, metal, extractant, C_ext, temperature) for pH in pHs]
            fig_mini = go.Figure()
            fig_mini.add_trace(go.Scatter(x=pHs, y=Es, mode='lines',
                line=dict(color=color_list[metals.index(metal)], width=2)))
            if target_pH:
                E_at_target = extraction_efficiency(target_pH, metal, extractant, C_ext, temperature)
                fig_mini.add_trace(go.Scatter(x=[target_pH], y=[E_at_target],
                    mode='markers', marker=dict(size=12, color='red', symbol='x'),
                    name=f'pH={target_pH}'))
                fig_mini.add_annotation(x=target_pH, y=E_at_target,
                    text=f'{E_at_target:.1f}%', showarrow=True, arrowhead=2, yshift=15)
            fig_mini.update_layout(height=180, margin=dict(l=30,r=10,t=10,b=30),
                xaxis_title='pH', yaxis_title='E(%)', showlegend=False,
                plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
                yaxis=dict(range=[0, 105]))
            st.plotly_chart(fig_mini, use_container_width=True)

        st.markdown("")

    # ---------------------------------------------------------------
    # 2. 분배 계수
    # ---------------------------------------------------------------
    st.markdown("---")
    st.markdown("### 2️⃣ 분배 계수 (Distribution Coefficient)")
    st.latex(r"D_M = \frac{C_{M,\text{org}}}{C_{M,\text{aq}}} = \frac{E_M / 100}{1 - E_M / 100}")

    if target_pH:
        st.markdown(f"**목표 pH = {target_pH}에서의 분배 계수:**")
        d_data = []
        for metal in metals:
            E_val = extraction_efficiency(target_pH, metal, extractant, C_ext, temperature)
            D_val = distribution_coefficient(target_pH, metal, extractant, C_ext, temperature=temperature)
            d_data.append({"금속": metal, f"E(pH={target_pH})": f"{E_val:.2f}%",
                           f"D(pH={target_pH})": f"{D_val:.4f}"})
        st.dataframe(pd.DataFrame(d_data), use_container_width=True, hide_index=True)

    # ---------------------------------------------------------------
    # 3. 물질수지
    # ---------------------------------------------------------------
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
    st.latex(
        r"C_{M,\text{org}}^{\text{out}} = D_M \cdot C_{M,\text{aq}}^{\text{out}}"
    )

    st.markdown(f"현재 설정: $Q_{{\\text{{aq}}}}$ = {Q_aq} L/hr, $Q_{{\\text{{org}}}}$ = {Q_org} L/hr, "
                f"A/O = {Q_aq/Q_org:.2f}")

    # ---------------------------------------------------------------
    # 4. pH 수지
    # ---------------------------------------------------------------
    st.markdown("---")
    st.markdown("### 4️⃣ pH 수지 (Proton Balance)")
    st.markdown("양이온 추출형 추출제의 금속-양성자 교환 반응:")
    st.latex(r"M^{n+}_{(\text{aq})} + n(\text{HA})_{2,(\text{org})} \rightarrow M(\text{A}_2)_{n,(\text{org})} + 2n\text{H}^+_{(\text{aq})}")

    st.markdown("**금속별 H⁺ 방출:**")
    for metal in metals:
        n_H = EXTRACTANT_PARAMS[extractant][metal]["n_H"]
        charge = "+" if n_H == 1 else f"{n_H}+"
        st.latex(
            rf"\text{{{metal}}}^{{{charge}}} \rightarrow {n_H} \text{{H}}^+ \text{{ 방출}}"
        )

    st.markdown("**수계 H⁺ 수지:**")
    st.latex(
        r"[\text{H}^+]_{\text{out}} \cdot Q_{\text{aq}} = "
        r"[\text{H}^+]_{\text{in}} \cdot Q_{\text{aq}} + "
        r"\sum_M n_{\text{H},M} \cdot \Delta C_M \cdot Q_{\text{aq}} - "
        r"C_{\text{NaOH}} \cdot Q_{\text{NaOH}}"
    )
    st.latex(
        r"\text{pH}_{\text{out}} = -\log_{10}\!\left([\text{H}^+]_{\text{out}}\right)"
    )

    if pH_mode == "목표 pH (자동 NaOH)":
        st.info(f"**목표 pH 모드**: pH = {target_pH if not staged_pHs else staged_pHs} 를 유지하기 위해 필요한 NaOH를 자동 역산합니다.\n\n"
                f"NaOH 필요량 = $[H^+]_{{in}} \\cdot Q_{{aq}} + \\Sigma(n_H \\cdot \\Delta C_M \\cdot Q_{{aq}}) - [H^+]_{{target}} \\cdot Q_{{aq}}$")

    # ---------------------------------------------------------------
    # 5. 추출률 공식
    # ---------------------------------------------------------------
    st.markdown("---")
    st.markdown("### 5️⃣ 성능 지표")
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
            D_vals = {m: distribution_coefficient(target_pH, m, extractant, C_ext, temperature=temperature) for m in metals}
            for i, m1 in enumerate(metals):
                for m2 in metals[i+1:]:
                    alpha = D_vals[m1] / D_vals[m2] if D_vals[m2] > 0 else float('inf')
                    sep_data.append({"M₁/M₂": f"{m1}/{m2}", "α": f"{alpha:.2f}"})
            st.dataframe(pd.DataFrame(sep_data), use_container_width=True, hide_index=True)

    # ---------------------------------------------------------------
    # 6. 파라미터 요약 테이블
    # ---------------------------------------------------------------
    st.markdown("---")
    st.markdown("### 📋 현재 파라미터 요약")
    param_rows = []
    for metal in metals:
        p = EXTRACTANT_PARAMS[extractant][metal]
        ph50_eff = get_effective_pH50(metal, extractant, C_ext, temperature)
        k_eff = get_effective_k(metal, extractant, temperature)
        param_rows.append({
            "금속": metal,
            "pH₅₀ (기준)": p["pH50"],
            "α": p["alpha"],
            "C_ref (M)": p["C_ref"],
            "β": p.get("beta", 0.0),
            "γ": p.get("gamma", 0.0),
            f"pH₅₀ ({temperature:.0f}°C)": round(ph50_eff, 3),
            f"k ({temperature:.0f}°C)": round(k_eff, 3),
            "E_max (%)": p["E_max"],
            "n_H": p["n_H"],
        })
    st.dataframe(pd.DataFrame(param_rows), use_container_width=True, hide_index=True)


# =============================================================================
# TAB 6: 변경 이력 (Changelog)
# =============================================================================
with tab6:
    st.subheader(f"📜 변경 이력 (현재 v{APP_VERSION})")
    st.markdown(CHANGELOG_CONTENT)


# =============================================================================
# 푸터
# =============================================================================
st.markdown("---")
st.caption(f"⚗️ SX Simulator v{APP_VERSION} | MSE Framework (ALTA 2024 / Wang et al. 2002, 2004, 2006)")
