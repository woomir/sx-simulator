#!/usr/bin/env python3
"""
SX Simulator Web Dashboard
===========================
Streamlit 기반 Li/Ni/Co/Mn Mixer-Settler 용매추출 시뮬레이션 대시보드.

실행: streamlit run sx_dashboard.py
"""

import sys, os, math, copy
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

APP_VERSION = "1.7.0"

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
from sx_simulator.config import (EXTRACTANT_PARAMS, MOLAR_MASS, DEFAULT_METALS,
                                  T_REF, MAX_LOADING_FRACTION)
from sx_simulator.fitting import (
    fit_sigmoid, sigmoid_model
)

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
st.sidebar.header("🚰 Feed (수계) 조건")
C_Li = st.sidebar.number_input("Li 농도 (g/L)", 0.0, 50.0, 1.5, 0.1)
C_Ni = st.sidebar.number_input("Ni 농도 (g/L)", 0.0, 150.0, 5.0, 0.5)
C_Co = st.sidebar.number_input("Co 농도 (g/L)", 0.0, 150.0, 3.0, 0.5)
C_Mn = st.sidebar.number_input("Mn 농도 (g/L)", 0.0, 150.0, 2.0, 0.5)

pH_feed = st.sidebar.number_input("Feed pH", 0.0, 14.0, 3.0, 0.1)
Q_aq = st.sidebar.number_input("Feed 수계 유량 (L/hr)", 1.0, 10000.0, 100.0, 10.0)

C_sulfate = st.sidebar.number_input("총 황산염 농도 (M)", 0.0, 5.0, 0.5, 0.1)
st.sidebar.caption("👉 Buffer Capacity 계산용 (0=순수물)")

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
    Q_NaOH = st.sidebar.number_input("총 NaOH 유량 (L/hr)", 0.0, 500.0, 12.0, 1.0)

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
    st.sidebar.caption("Stage별 목표 pH 입력 (1단~N단):")
    cols = st.sidebar.columns(min(n_stages, 4))
    staged_pHs = []
    for i in range(n_stages):
        with cols[i % len(cols)]:
            val = st.number_input(f"St.{i+1}", 2.0, 8.0, float(target_pH), 0.1, key=f"stage_ph_{i}")
            staged_pHs.append(val)

# NaOH 분배 전략 (고정 NaOH 모드)
naoh_strategy = "uniform"
naoh_weights = None
if pH_mode == "고정 NaOH":
    st.sidebar.markdown("---")
    st.sidebar.header("💧 NaOH 분배 전략")
    naoh_strategy_label = st.sidebar.selectbox("분배 방식", ["균등 (Uniform)", "전단집중 (Front-loaded)", "커스텀 (Custom)"])
    naoh_strategy_map = {"균등 (Uniform)": "uniform", "전단집중 (Front-loaded)": "front_loaded", "커스텀 (Custom)": "custom"}
    naoh_strategy = naoh_strategy_map[naoh_strategy_label]

    if naoh_strategy == "custom":
        st.sidebar.caption("각 Stage별 가중치 (비율) 입력:")
        cols = st.sidebar.columns(min(n_stages, 4))
        naoh_weights = []
        for i in range(n_stages):
            with cols[i % len(cols)]:
                w = st.number_input(f"St.{i+1}", 0.0, 10.0, 1.0, 0.1, key=f"naoh_w_{i}")
                naoh_weights.append(w)


st.sidebar.markdown("---")

# --- 고급 옵션 (Phase 3) ---
st.sidebar.header("⚙️ 고급 옵션 (Phase 3)")
use_competition = st.sidebar.checkbox("추출제 경쟁 분배 활성", value=False, help="Vasilyev et al. (2019) 모델 기반: 수계 금속 이온 간의 공유 유기 추출제 풀 경쟁 방식을 모델에 반영하여 고로딩 조건 정합성을 개선합니다.")
use_speciation = st.sidebar.checkbox("수계 종분화 효과 활성", value=False, help="금속 수산화물 착물(MOH⁺) 생성에 의한 추가적인 pH 완충 효과를 시뮬레이션에 반영합니다.")

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
tab1, tab2, tab3, tab8, tab4, tab5, tab6, tab9, tab7 = st.tabs([
    "📊 시뮬레이션 결과", "📐 모델 수식", "📈 pH Isotherm", "📉 McCabe-Thiele",
    "🔬 추출제 비교", "📋 상세 데이터", "📝 데이터 피팅", "📖 용어 해설", "📜 변경 이력"
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
        C_sulfate=C_sulfate,
        use_competition=use_competition,
        use_speciation=use_speciation,
    )

    if pH_mode == "목표 pH (자동 NaOH)":
        if staged_pHs:
            sim_kwargs["target_pH_per_stage"] = staged_pHs
        else:
            sim_kwargs["target_pH"] = target_pH
    else:
        sim_kwargs["C_NaOH"] = C_NaOH
        sim_kwargs["Q_NaOH"] = Q_NaOH
        sim_kwargs["naoh_strategy"] = naoh_strategy
        if naoh_weights is not None:
            sim_kwargs["naoh_weights"] = naoh_weights

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
        stages_x = [0] + list(range(1, n_stages + 1))
        pH_y = [float(pH_feed)] + result["pH_profile"]
        fig_ph = go.Figure()
        fig_ph.add_trace(go.Scatter(
            x=stages_x, y=pH_y,
            mode='lines+markers', name='pH',
            line=dict(color='#e94560', width=3),
            marker=dict(size=12, symbol='circle'),
        ))
        fig_ph.update_layout(
            xaxis_title="Stage", yaxis_title="pH",
            height=350, plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
            xaxis=dict(
                tickmode='array',
                tickvals=stages_x,
                ticktext=["Feed"] + [f"St.{i}" for i in range(1, n_stages + 1)]
            ),
        )
        st.plotly_chart(fig_ph, use_container_width=True)

    with col_naoh:
        st.subheader("🧪 Stage별 NaOH 소비량")
        fig_naoh = go.Figure()
        naoh_y = [0.0] + result["NaOH_profile"]
        fig_naoh.add_trace(go.Bar(
            x=stages_x, y=naoh_y,
            marker_color='#0f3460', text=[f"{v:.1f}" if v > 0 else "" for v in naoh_y],
            textposition='outside',
        ))
        fig_naoh.update_layout(
            xaxis_title="Stage", yaxis_title="NaOH (mol/hr)",
            height=350, plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
            xaxis=dict(
                tickmode='array',
                tickvals=stages_x,
                ticktext=["Feed"] + [f"St.{i}" for i in range(1, n_stages + 1)]
            ),
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
# TAB 8: McCabe-Thiele 다이어그램
# =============================================================================
with tab8:
    st.subheader("📉 McCabe-Thiele 다이어그램")
    st.markdown("선택한 금속에 대한 조작선(Operating Line)과 평형 곡선(Equilibrium Curve)을 표시합니다.")
    
    mt_metal = st.selectbox("다이어그램을 그릴 금속 선택", metals, index= metals.index("Co") if "Co" in metals else 0)
    
    # 평형 곡선 생성을 위한 기준 pH (최종 후액 pH 기준)
    ref_pH = result["pH_profile"][-1]
    st.caption(f"💡 평형 곡선은 최종 후액 기준 pH ({ref_pH:.2f}) 로 계산되었습니다.")
    
    # X축(수계 농도) 범위 설정 (후액 ~ 피드)
    x_max = C_aq_feed[mt_metal] * 1.2
    x_vals = np.linspace(0, x_max, 100)
    
    # 평형 곡선 계산
    eq_y_vals = []
    for x in x_vals:
        # 단일 금속만 존재한다고 단순 가정하여 Isotherm 곡선 생성
        D = distribution_coefficient(ref_pH, mt_metal, extractant, C_ext, temperature)
        # 로딩 감쇠 적용 (단일 금속 근사)
        from sx_simulator.extraction_isotherm import loading_damping_factor
        sim_org = {mt_metal: D * x}
        from sx_simulator.extraction_isotherm import calc_loading_fraction
        loading = calc_loading_fraction(sim_org, extractant, C_ext, [mt_metal])
        D = D * loading_damping_factor(loading)
        
        # y = D * x
        y = D * x
        eq_y_vals.append(y)
        
    fig_mt = go.Figure()
    
    # 1. 평형 곡선
    fig_mt.add_trace(go.Scatter(x=x_vals, y=eq_y_vals, mode='lines', name='Equilibrium Curve (평형선)', line=dict(color='blue', width=2)))
    
    # 2. 조작선 (Operating Line)
    # y = (Q_aq / Q_org) * x + (y_in - (Q_aq / Q_org) * x_in)
    # 여기서 x_in 은 feed, y_in 은 result['loaded_organic'][metal]
    slope = Q_aq / Q_org
    x_in_feed = C_aq_feed[mt_metal]
    y_out_org = result['loaded_organic'][mt_metal]
    x_out_raff = result['raffinate'][mt_metal]
    y_in_org = 0.0 # fresh solvent
    
    fig_mt.add_trace(go.Scatter(
        x=[x_out_raff, x_in_feed], 
        y=[y_in_org, y_out_org], 
        mode='lines', name='Operating Line (조작선)', 
        line=dict(color='red', width=2, dash='dash')
    ))
    
    # 3. Stage 단계 그리기
    stage_x = []
    stage_y = []
    
    # 역류 추출이므로 (x_out, y_in) 에서 시작.
    # Stage N(마지막 단, x_raff, y_in) -> 평형선 위 (x_raff, y_stage_N) 
    # -> 조작선 표면 -> ...
    current_x = x_out_raff
    current_y = y_in_org
    
    stage_x.append(current_x)
    stage_y.append(current_y)
    
        # 실제 시뮬레이션된 단수 포인트들을 기반으로 단계 그리기
    for s_idx in range(n_stages-1, -1, -1):
        # 수직선: 수계 농도는 유지, 유기상 농도가 평형(또는 실제 stage 출구)으로 이동
        actual_y = result['stages'][s_idx]['C_org_out'][mt_metal]
        stage_x.append(current_x)
        stage_y.append(actual_y)
        
        # 수평선: 유기상 농도는 유지, 수계 농도가 다음 stage(또는 feed) 입구로 이동
        # 조작선 상의 해당 x점 찾기
        current_y = actual_y
        if s_idx == 0:
            current_x = C_aq_feed[mt_metal]
        else:
            current_x = result['stages'][s_idx - 1]['C_aq_out'][mt_metal]
            
        stage_x.append(current_x)
        stage_y.append(current_y)
        
    fig_mt.add_trace(go.Scatter(
        x=stage_x, y=stage_y,
        mode='lines+markers', name='Stages (실제 단수)',
        line=dict(color='green', width=2, shape='hv'),
        marker=dict(size=8, symbol='circle')
    ))
    
    fig_mt.update_layout(
        title=f"McCabe-Thiele for {mt_metal}",
        xaxis_title=f"Aqueous {mt_metal} (g/L)",
        yaxis_title=f"Organic {mt_metal} (g/L)",
        height=600,
        plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
        legend=dict(x=0.01, y=0.99),
    )
    st.plotly_chart(fig_mt, use_container_width=True)

# =============================================================================
# TAB 4: 추출제 비교
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

            if use_competition:
                st.markdown("**Phase 3: 금속 추출제 경합 보정 활성**")
                st.latex(r"D_{M}^{\text{adj}} = D_{M}^{\text{sig}} \times \left(\frac{[\overline{\text{HA}}]_{\text{free}}}{C_{\text{ext}}}\right)^{n_{\text{ext}}}")
                st.caption("유기상 잔여 액체 추출제량에 비례하여 높은 로딩(Loading) 시 추출 효율을 동적으로 감소시킵니다 (Vasilyev 모델).")
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

    if use_speciation:
        st.markdown("**Phase 3: 수계 금속 종분화(Speciation) 효과 활성**")
        st.latex(r"[\text{H}^+]_{\text{hydrolysis}} \approx \sum_M K_{MOH} \cdot [\text{M}^{n+}] \cdot [\text{OH}^-]")
        st.caption("고 pH 역외 접근 시 금속 수산화물 착물(MOH⁺ 등)이 생성되며 양성자를 방출(OH⁻ 소비)해 직접적인 완충제로 작용함을 수지에 추가 반영했습니다.")

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
# TAB 6: 데이터 피팅 (Experimental Data Fitting)
# =============================================================================
with tab6:
    st.subheader("📝 실험 데이터 피팅")
    st.markdown("""
    실험 데이터(pH vs 추출률)를 업로드하여 시그모이드 모델 파라미터를 자동으로 피팅합니다.
    피팅된 파라미터는 시뮬레이션에 적용할 수 있습니다.
    """)

    # --- 샘플 데이터 다운로드 ---
    with st.expander("📥 샘플 데이터 다운로드"):
        st.markdown("""
        아래 버튼을 클릭하여 샘플 CSV 파일을 다운로드하세요. 이 파일을 참고하여 데이터를 준비하시면 됩니다.

        **필수 열**: `pH`, `E_pct` (추출률 %)
        **선택 열**: `metal` (금속명 - 복수 금속 동시 피팅 시)
        """)
        sample_csv = "pH,E_pct,metal\n2.0,0.5,Co\n2.5,2.0,Co\n3.0,12.0,Co\n3.5,45.0,Co\n4.0,82.0,Co\n4.5,96.0,Co\n5.0,99.5,Co\n2.0,0.1,Ni\n3.0,0.5,Ni\n4.0,3.0,Ni\n5.0,25.0,Ni\n5.5,60.0,Ni\n6.0,90.0,Ni\n6.5,98.0,Ni"
        st.download_button("📥 샘플 CSV 다운로드", sample_csv, "sample_extraction_data.csv", "text/csv")

    # --- 데이터 업로드 ---
    uploaded = st.file_uploader("📂 CSV 파일 업로드 (pH, E_pct 열 필수)", type=["csv"])

    if uploaded is not None:
        try:
            df_raw = pd.read_csv(uploaded)
        except Exception as e:
            st.error(f"CSV 파일 읽기 오류: {e}")
            df_raw = None

        if df_raw is not None:
            # 열 이름 검증
            required_cols = {"pH", "E_pct"}
            if not required_cols.issubset(set(df_raw.columns)):
                st.error(f"CSV에 필수 열이 없습니다: {required_cols}. 현재 열: {list(df_raw.columns)}")
            else:
                st.markdown("#### 📊 업로드된 데이터")
                st.dataframe(df_raw.head(20), use_container_width=True, hide_index=True)

                # 금속 분류 처리
                if "metal" in df_raw.columns:
                    metal_list = sorted(df_raw["metal"].unique().tolist())
                else:
                    metal_list = ["Unknown"]
                    df_raw["metal"] = "Unknown"

                st.markdown(f"감지된 금속: **{', '.join(metal_list)}** ({len(df_raw)}개 데이터 포인트)")

                # --- 피팅 설정 ---
                st.markdown("---")
                st.markdown("#### ⚙️ 피팅 설정")
                # Only Sigmoid model is supported for fitting.
                fit_model = "Sigmoid" # Implicitly set to Sigmoid

                # --- 피팅 실행 ---
                if st.button("🚀 피팅 실행", type="primary", use_container_width=True):
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

                        result = fit_sigmoid(pH_vals, E_vals)

                        if not result["success"]:
                            st.error(f"피팅 실패: {result.get('error', 'Unknown error')}")
                            continue

                        # --- 피팅 결과 표시 ---
                        res_col1, res_col2 = st.columns([1, 1])

                        with res_col1:
                            st.markdown("**피팅 파라미터:**")
                            params = result["params"]
                            errors = result["errors"]

                            param_df = pd.DataFrame([
                                {"파라미터": "pH\u2085\u2080", "값": params["pH50"],
                                 "±95%CI": errors["pH50_err"]},
                                {"파라미터": "k", "값": params["k"],
                                 "±95%CI": errors["k_err"]},
                                {"파라미터": "E_max (%)", "값": params["E_max"],
                                 "±95%CI": errors["E_max_err"]},
                            ])
                            st.dataframe(param_df, use_container_width=True, hide_index=True)
                            st.metric("R\u00b2", f"{result['r_squared']:.4f}")

                        with res_col2:
                            # 피팅 곡선 차트
                            import plotly.graph_objects as go
                            fig_fit = go.Figure()

                            # 실험 데이터 (scatter)
                            fig_fit.add_trace(go.Scatter(
                                x=pH_vals, y=E_vals,
                                mode="markers", name="실험 데이터",
                                marker=dict(size=10, color="#e94560")
                            ))

                            # 피팅 곡선
                            pH_fine = np.linspace(float(min(pH_vals)) - 0.5, float(max(pH_vals)) + 0.5, 200)
                            E_fit = sigmoid_model(pH_fine, params["pH50"], params["k"], params["E_max"])
                            fig_fit.add_trace(go.Scatter(
                                x=pH_fine, y=E_fit,
                                mode="lines", name="피팅 곡선",
                                line=dict(color="#0f3460", width=2)
                            ))
                            fig_fit.update_yaxes(title="추출률 (%)")

                            fig_fit.update_layout(
                                title=f"{metal_name} 피팅 결과 ({fit_model})",
                                xaxis_title="pH",
                                height=350,
                                template="plotly_white",
                            )
                            st.plotly_chart(fig_fit, use_container_width=True)

                        # --- 파라미터 적용 버튼 (Sigmoid만) ---
                        if metal_name in DEFAULT_METALS:
                            if st.button(f"✅ {metal_name} 피팅 파라미터 적용", key=f"apply_{metal_name}"):
                                EXTRACTANT_PARAMS[extractant][metal_name]["pH50"] = params["pH50"]
                                EXTRACTANT_PARAMS[extractant][metal_name]["k"] = params["k"]
                                EXTRACTANT_PARAMS[extractant][metal_name]["E_max"] = params["E_max"]
                                st.success(f"{metal_name} 파라미터가 적용되었습니다! 페이지를 새로고침하면 시뮬레이션에 반영됩니다.")
                                st.rerun()

                        st.markdown("---")
    else:
        st.info("👆 CSV 파일을 업로드하여 피팅을 시작하세요.")


# =============================================================================
# TAB 9: 용어 해설 (Glossary)
# =============================================================================
with tab9:
    st.subheader("📖 용어 해설")
    st.markdown("이 시뮬레이터에서 사용되는 전문 용어, 기호, 약어를 쉽게 풀어서 설명합니다.")

    # --- 공정 용어 ---
    st.markdown("### ⚗️ 공정 용어")
    st.markdown("""
| 용어 | 영어 | 뜻 |
|---|---|---|
| **역류** | Counter-current | 수계(물)와 유기계(기름)가 서로 반대 방향으로 흐르는 방식. 추출 효율을 극대화하는 운전법 |
| **수계** | Aqueous Phase | 금속이 녹아 있는 물(수용액) 쪽. 약자로 'Aq'라고 씀 |
| **유기계** | Organic Phase | 추출제가 녹아 있는 기름(유기용매) 쪽. 약자로 'Org'라고 씀 |
| **피드** | Feed | 공정에 처음 투입되는 원료 용액 |
| **후액** | Raffinate | 금속이 추출된 뒤 남은 수계 용액. 금속 농도가 낮아진 상태 |
| **로딩** | Loading | 추출제가 금속을 얼마나 담고 있는지의 비율. 100%에 가까우면 '포화'되어 더 이상 추출 못함 |
| **사포닌화** | Saponification | 추출제(HL)를 NaOH로 미리 처리하여 Na형(NaL)으로 바꿔주는 것. 추출 효율을 높이는 전처리 |
| **O/A 비** | Organic/Aqueous Ratio | 유기계 유량 ÷ 수계 유량. 이 비율이 크면 추출제가 상대적으로 많아져 추출률이 올라감 |
    """)

    # --- 화학/열역학 용어 ---
    st.markdown("### 🧪 화학·열역학 용어")
    st.markdown("""
| 용어 | 영어 | 뜻 |
|---|---|---|
| **추출률 (%)** | Extraction Efficiency | 피드 속 금속이 유기계로 얼마나 옮겨갔는지의 백분율 |
| **분배계수 (D)** | Distribution Coefficient | (유기계 금속 농도) ÷ (수계 금속 농도). 값이 클수록 금속이 유기계를 '좋아함' |
| **등온선** | Isotherm | 온도를 고정한 상태에서 pH에 따라 추출률이 어떻게 변하는지 그린 S자 곡선 |
| **완충 용량** | Buffer Capacity | 용액이 pH 변화에 저항하는 능력. 황산염이 많으면 pH가 쉽게 변하지 않음 |
| **Proton Balance** | Proton Balance | 반응 전후로 H⁺(수소이온)가 얼마나 생기고 사라지는지 계산하는 것 |
| **물질수지** | Mass Balance | '들어간 양 = 나온 양'이라는 보존 법칙. 시뮬레이터의 핵심 원리 |
    """)

    # --- 수학 기호 ---
    st.markdown("### 📐 수학 기호·파라미터")
    st.markdown("""
| 기호 | 이름 | 뜻 |
|---|---|---|
| **pH₅₀** | 반추출 pH | 추출률이 50%가 되는 pH 값. 금속마다 다름 |
| **k** | 기울기 계수 | S자 곡선이 얼마나 가파른지를 나타냄. 클수록 pH 변화에 민감 |
| **E_max** | 최대 추출률 | 아무리 pH를 올려도 넘지 못하는 추출률의 상한(보통 99~100%) |
| **n_H** | 수소이온 화학양론수 | 금속 1개가 추출될 때 방출되는 H⁺의 개수 (예: Ni²⁺→2, Co²⁺→2) |
| **α (알파)** | 추출제 농도 보정 계수 | 추출제 농도가 pH₅₀에 미치는 영향의 크기 |
| **β (베타)** | 온도 보정 계수 | 온도가 pH₅₀을 얼마나 변화시키는지의 크기 |
| **γ (감마)** | k-온도 보정 계수 | 온도가 기울기(k)를 얼마나 변화시키는지의 크기 |
| **C_ext** | 추출제 농도 (M) | 유기계에 녹인 추출제의 몰 농도 |
| **Q** | 유량 (L/hr) | 액체가 시간당 흐르는 양. Q_aq=수계, Q_org=유기계 |
    """)

    # --- 차트 용어 ---
    st.markdown("### 📈 차트·다이어그램 용어")
    st.markdown("""
| 용어 | 영어 | 뜻 |
|---|---|---|
| **McCabe-Thiele** | McCabe-Thiele Diagram | 다단 추출 공정의 이론 단수를 시각적으로 구하는 그래프. 화학공학 교과서의 대표적인 분석법 |
| **평형 곡선** | Equilibrium Curve | 수계 농도와 유기계 농도가 평형(더 이상 변하지 않는 상태)을 이룰 때의 관계를 나타낸 곡선 |
| **조작선** | Operating Line | 물질수지(들어간 양=나온 양)로부터 그려지는 직선. 실제 공정의 운전 조건을 나타냄 |
| **Pinch Point** | Pinch Point | 평형 곡선과 조작선이 거의 만나는 지점. 여기서는 단수를 아무리 늘려도 추출이 더 이상 안 됨 |
| **Sigmoid** | Sigmoid Curve | 'S'자 모양의 수학 함수. pH에 따른 추출률을 모사하는 데 사용 |
| **R²** | 결정계수 (Coefficient of Determination) | 피팅(회귀)이 얼마나 잘 맞는지의 지표. 1.0에 가까울수록 완벽한 적합 |
    """)

    # --- 추출제 이름 ---
    st.markdown("### 🫗 추출제 이름")
    st.markdown("""
| 이름 | 정식 명칭 | 설명 |
|---|---|---|
| **Cyanex 272** | Bis(2,4,4-trimethylpentyl)phosphinic acid | 코발트(Co)/니켈(Ni) 분리에 탁월한 인산 계열 추출제 |
| **D2EHPA** | Di-(2-ethylhexyl)phosphoric acid | 범용적인 인산 계열 추출제. 가격↓, 선택성은 Cyanex 272보다 낮음 |
    """)

    # --- 시뮬레이터 특수 기능 ---
    st.markdown("### ⚙️ 시뮬레이터 고급 기능 용어")
    st.markdown("""
| 용어 | 뜻 |
|---|---|
| **Bisection Solver** | 이분법 탐색. '정답'이 있는 범위를 반씩 좁혀가며 빠르고 안정적으로 해(solution)를 찾는 수치 알고리즘 |
| **NaOH 분배 전략** | 다단 공정에서 총 NaOH를 각 Stage에 어떻게 나눠줄지 결정하는 방식 (균등/전단집중/커스텀) |
| **균등 (Uniform)** | 모든 Stage에 NaOH를 동일하게 나눠주는 방식 |
| **전단집중 (Front-loaded)** | 앞쪽 Stage에 NaOH를 더 많이 주는 방식. pH를 빨리 올려 초기 추출을 강화 |
| **수렴 (Convergence)** | 다단 계산을 반복할 때, 결과가 더 이상 변하지 않고 안정되는 상태 |
| **피팅 (Fitting)** | 실험 데이터에 수학 모델을 맞추는(최적화하는) 과정 |
| **MSE Framework** | Mixed-Solvent Electrolyte. 이 시뮬레이터가 기반으로 하는 열역학 프레임워크 이름 |
    """)

# =============================================================================
# TAB 7: 변경 이력 (Changelog)
# =============================================================================
with tab7:
    st.subheader(f"📜 변경 이력 (현재 v{APP_VERSION})")
    st.markdown(CHANGELOG_CONTENT)


# =============================================================================
# 푸터
# =============================================================================
st.markdown("---")
st.caption(f"⚗️ SX Simulator v{APP_VERSION} | MSE Framework (ALTA 2024 / Wang et al. 2002, 2004, 2006)")
