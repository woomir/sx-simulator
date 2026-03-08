#!/usr/bin/env python3
"""
SX Simulator Web Dashboard
===========================
Streamlit 기반 Li/Ni/Co/Mn Mixer-Settler 용매추출 시뮬레이션 대시보드.

실행: streamlit run sx_dashboard.py
"""

import sys, os, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# 버전 정보 명시 (문서 및 캐시 관리에 활용 가능)
APP_VERSION = "2.1.1"
LAST_UPDATED = "2026-03-07"
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
    extraction_efficiency, get_effective_pH50,
    get_effective_k, distribution_coefficient, calc_loading_fraction,
)
from sx_simulator.config import (DEFAULT_METALS,
                                  T_REF,
                                  get_parameter_profile,
                                  SPECIATION_CONSTANTS)
from sx_simulator.dashboard_service import (
    SimulationInputs,
    build_scope_assessment,
    compute_loading_pct,
    run_compare_simulations,
    run_simulation,
)
from sx_simulator.datasets import (
    DASHBOARD_PRESETS as PRESETS,
    KNOWN_PRESET_NOTES,
    calc_sulfate_from_feed,
)
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

PROFILE_LABELS = {
    "현장 보정 (Isd 108)": "field_calibrated",
    "문헌 기본값": "literature_default",
}
METAL_COLORS = {
    "Li": "#ff6b6b",
    "Ni": "#4ecdc4",
    "Co": "#45b7d1",
    "Mn": "#f7b731",
    "Ca": "#a29bfe",
    "Mg": "#fd79a8",
    "Zn": "#6c5ce7",
}
PH_PLOT_RANGE = [x * 0.1 for x in range(10, 101)]


def get_active_extractant_params() -> dict:
    """현재 세션에서 사용하는 파라미터 세트를 반환합니다."""
    return st.session_state.custom_params


def apply_profile() -> None:
    """선택한 보정 프로필을 세션과 엔진에 반영합니다."""
    profile_key = PROFILE_LABELS[st.session_state.ui_param_profile]
    st.session_state.custom_params = get_parameter_profile(profile_key)


def load_markdown_doc(relative_path: str, missing_message: str) -> str:
    """docs 또는 루트의 markdown 파일을 읽어 반환합니다."""
    doc_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), relative_path)
    try:
        with open(doc_path, "r", encoding="utf-8") as file:
            return file.read()
    except FileNotFoundError:
        return missing_message


def render_scope_assessment(container, assessment: dict) -> None:
    """검증 커버리지 평가를 Streamlit 컨테이너에 렌더링합니다."""
    if assessment["level"] == "high":
        container.success("🧭 현장 검증 범위와 잘 겹치는 조건입니다. 내부 비교와 추세 해석에 유리합니다.")
    elif assessment["level"] == "medium":
        container.warning("🧭 일부만 현장 검증 범위와 겹칩니다. 방향성 해석은 가능하지만 정량값은 주의가 필요합니다.")
    else:
        container.error("🧭 현장 검증 범위를 벗어난 조건입니다. 정량 예측보다 경향 확인 용도로 해석하는 편이 안전합니다.")

    if assessment["highlights"]:
        container.markdown("**해석 포인트**")
        for item in assessment["highlights"]:
            container.markdown(f"- {item}")

    if assessment["cautions"]:
        container.markdown("**주의사항**")
        for item in assessment["cautions"]:
            container.markdown(f"- {item}")

def apply_preset():
    sel = st.session_state.preset_sel
    if sel == "사용자 직접 입력": return
    d = PRESETS[sel]
    # 희석 효과 반영하여 Feed 농도 입력 (BM Feed와 NaOH 수용액의 합산 유량 기준)
    dilution = d['feed_flow'] / (d['feed_flow'] + d['naoh_flow'])
    st.session_state.ui_C_Li = float(d['input']['Li'] * dilution)
    st.session_state.ui_C_Ni = float(d['input']['Ni'] * dilution)
    st.session_state.ui_C_Co = float(d['input']['Co'] * dilution)
    st.session_state.ui_C_Mn = float(d['input']['Mn'] * dilution)
    st.session_state.ui_C_Ca = float(d['input']['Ca'] * dilution)
    st.session_state.ui_C_Mg = float(d['input']['Mg'] * dilution)
    st.session_state.ui_C_Zn = float(d['input']['Zn'] * dilution)
    
    st.session_state.ui_Q_aq = float(d['feed_flow'] + d['naoh_flow'])
    st.session_state.ui_Q_org = float(d['org_flow'])
    st.session_state.ui_extractant = d['ext']
    st.session_state.ui_C_ext = d['C_ext']
    st.session_state.ui_target_pH = d['pH']
    st.session_state.ui_n_stages = d['n_stages']
    st.session_state.ui_temperature = d['T']
    st.session_state.ui_pH_mode = "목표 pH (자동 NaOH)"

# =============================================================================
# 사이드바: 입력 파라미터
# =============================================================================
st.sidebar.title("⚗️ SX 시뮬레이터")

if "ui_param_profile" not in st.session_state:
    st.session_state.ui_param_profile = "현장 보정 (Isd 108)"
if "custom_params" not in st.session_state:
    st.session_state.custom_params = get_parameter_profile(
        PROFILE_LABELS[st.session_state.ui_param_profile]
    )

st.sidebar.selectbox(
    "🧪 현장 실험 Data 프리셋", 
    ["사용자 직접 입력"] + list(PRESETS.keys()),
    key="preset_sel",
    on_change=apply_preset,
    help="선택 시 실험 조건 값(희석된 초기 농도, 유량, 추출제, 단수, 목표 pH 등)이 즉시 동기화됩니다."
)
selected_preset_note = KNOWN_PRESET_NOTES.get(st.session_state.get("preset_sel"))
if selected_preset_note:
    st.sidebar.warning(selected_preset_note)

st.sidebar.markdown("---")

st.sidebar.header("🧭 모델 보정 프로필")
st.sidebar.selectbox(
    "보정 기준",
    list(PROFILE_LABELS.keys()),
    key="ui_param_profile",
    on_change=apply_profile,
    help="현장 보정 프로필은 Data1~6과 같은 현장 조건에 더 가깝고, 문헌 기본값은 범용 비교용 기준선입니다.",
)
if st.session_state.ui_param_profile == "현장 보정 (Isd 108)":
    st.sidebar.info("현재 프로필은 Isd 108 현장 보정 파라미터를 사용합니다.")
else:
    st.sidebar.warning("현재 프로필은 문헌 기본값을 사용합니다. 현장 프리셋과의 직접 정합도는 낮아질 수 있습니다.")

st.sidebar.markdown("---")

# --- Feed 조건 ---
st.sidebar.header("🚰 Feed (수계) 조건")
C_Li = st.sidebar.number_input("Li 농도 (g/L)", 0.0, 50.0, 1.5, 0.1, key="ui_C_Li")
C_Ni = st.sidebar.number_input("Ni 농도 (g/L)", 0.0, 150.0, 5.0, 0.5, key="ui_C_Ni")
C_Co = st.sidebar.number_input("Co 농도 (g/L)", 0.0, 150.0, 3.0, 0.5, key="ui_C_Co")
C_Mn = st.sidebar.number_input("Mn 농도 (g/L)", 0.0, 150.0, 2.0, 0.5, key="ui_C_Mn")
C_Ca = st.sidebar.number_input("Ca 농도 (g/L)", 0.0, 50.0, 0.0, 0.1, key="ui_C_Ca")
C_Mg = st.sidebar.number_input("Mg 농도 (g/L)", 0.0, 50.0, 0.0, 0.1, key="ui_C_Mg")
C_Zn = st.sidebar.number_input("Zn 농도 (g/L)", 0.0, 50.0, 0.0, 0.1, key="ui_C_Zn")

pH_feed = st.sidebar.number_input("Feed pH", 0.0, 14.0, 3.0, 0.1)
Q_aq = st.sidebar.number_input("Feed 수계 유량 (L/hr)", 1.0, 10000.0, 100.0, 10.0, key="ui_Q_aq")

# 총 황산염 농도 — Feed 금속 황산염 조성으로부터 자동 계산
C_sulfate = calc_sulfate_from_feed({
    "Li": C_Li,
    "Ni": C_Ni,
    "Co": C_Co,
    "Mn": C_Mn,
    "Ca": C_Ca,
    "Mg": C_Mg,
    "Zn": C_Zn,
})
st.sidebar.info(f"💎 계산된 SO₄²⁻: **{C_sulfate:.3f} M**")

st.sidebar.markdown("---")

# --- NaOH 조건 ---
# --- NaOH 조건 ---
st.sidebar.header("🧪 pH 제어 (NaOH)")
pH_mode = st.sidebar.radio("pH 제어 모드", ["목표 pH (자동 NaOH)", "고정 NaOH"], index=0, key="ui_pH_mode", help="'목표 pH' 모드를 적극 권장합니다. 시뮬레이터가 원하는 산도에 도달하기 위한 최적 알칼리 투입량을 역산합니다.")

if pH_mode == "목표 pH (자동 NaOH)":
    target_pH = st.sidebar.slider("목표 pH", 2.0, 8.0, 5.0, 0.1, key="ui_target_pH")
    use_staged_pH = st.sidebar.checkbox("Stage별 차등 pH", value=False)
else:
    target_pH = None
    C_NaOH = st.sidebar.number_input("NaOH 농도 (M)", 0.1, 20.0, 5.0, 0.5)
    Q_NaOH = st.sidebar.number_input("총 NaOH 유량 (L/hr)", 0.0, 500.0, 12.0, 1.0)

st.sidebar.markdown("---")

# --- 용매 조건 ---
st.sidebar.header("🫗 유기계 (용매) 조건")
extractant = st.sidebar.selectbox("추출제", ["Cyanex 272", "D2EHPA"], key="ui_extractant", help="Cyanex 272는 Co를 Ni보다 먼저 추출하고, D2EHPA는 불순물(Ca, Mg, Zn, Mn)을 주로 선추출합니다.")
C_ext = st.sidebar.number_input("추출제 농도 (M)", 0.05, 5.0, 0.5, 0.05, key="ui_C_ext", help="사용할 용매의 실제 유효 농도입니다. 부족할 경우 추출 효율이 급감합니다.")
Q_org = st.sidebar.number_input("유기계 유량 (L/hr)", 1.0, 10000.0, 100.0, 10.0, key="ui_Q_org", help="유기상 투입 유량입니다. 수계 유량과 비례하여 O/A 비를 결정합니다.")

# --- 온도 설정 ---
st.sidebar.markdown("---")
st.sidebar.header("🌡️ 온도 설정")
temperature = st.sidebar.slider("운전 온도 (°C)", 10.0, 60.0, 25.0, 1.0, key="ui_temperature")
if abs(temperature - T_REF) > 0.5:
    st.sidebar.info(f"온도 보정 활성: ΔT = {temperature - T_REF:+.0f}°C")

st.sidebar.markdown("---")

# --- Stage 설정 ---
st.sidebar.header("🔄 Mixer-Settler 설정")
n_stages = st.sidebar.slider("Stage 수", 1, 10, 4, key="ui_n_stages", help="실제 현장에 구성될 연속 추출 다단 반응조의 갯수입니다. 많을수록 추출 한계치에 수렴합니다.")

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

# 고급 옵션 (Phase 3)은 성능 향상을 위해 기본 활성화 됨 (무조건 True)

# --- 모델 파라미터 편집 ---
st.sidebar.header("📐 모델 파라미터 편집")
edit_params = st.sidebar.checkbox("Isotherm 파라미터 수정", value=False)

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

    if st.sidebar.button("🔄 선택 프로필로 초기화"):
        st.session_state.custom_params = get_parameter_profile(
            PROFILE_LABELS[st.session_state.ui_param_profile]
        )
        st.rerun()

active_extractant_params = get_active_extractant_params()

# =============================================================================
# 메인 영역
# =============================================================================
st.title("⚗️ Li/Ni/Co/Mn Mixer-Settler SX 시뮬레이터")
st.caption("MSE Thermodynamic Framework 기반 (ALTA 2024 / Wang et al.)")

metals = DEFAULT_METALS
C_aq_feed = {"Li": C_Li, "Ni": C_Ni, "Co": C_Co, "Mn": C_Mn,
             "Ca": C_Ca, "Mg": C_Mg, "Zn": C_Zn}
simulation_inputs = SimulationInputs(
    C_aq_feed=C_aq_feed,
    pH_feed=pH_feed,
    Q_aq=Q_aq,
    Q_org=Q_org,
    extractant=extractant,
    C_ext=C_ext,
    n_stages=n_stages,
    temperature=temperature,
    C_sulfate=C_sulfate,
    pH_mode=pH_mode,
    target_pH=target_pH,
    staged_pHs=staged_pHs,
    C_NaOH=C_NaOH if pH_mode == "고정 NaOH" else 0.0,
    Q_NaOH=Q_NaOH if pH_mode == "고정 NaOH" else 0.0,
    naoh_strategy=naoh_strategy,
    naoh_weights=naoh_weights,
    metals=tuple(metals),
)

# =============================================================================
# 탭 구성
# =============================================================================
tab1, tab2, tab3, tab8, tab4, tab5, tab6, tab9, tab10, tab11, tab7, tab0 = st.tabs([
    "📊 시뮬레이션 결과", "📐 수식 및 메커니즘 알고리즘", "📈 pH Isotherm", "📉 McCabe-Thiele",
    "🔬 추출제 비교", "📋 상세 데이터", "📝 데이터 피팅", "📖 용어 해설", "📚 파라미터 및 문헌", "🧪 데이터 피팅 검증 이력", "📜 변경 이력", "📘 사용자 매뉴얼"
])

# =============================================================================
# TAB 0: 사용자 매뉴얼
# =============================================================================
with tab0:
    st.markdown(
        load_markdown_doc(
            os.path.join("docs", "user_manual.md"),
            "매뉴얼 파일을 찾을 수 없습니다.",
        )
    )

# =============================================================================
# 시뮬레이션 실행
# =============================================================================
with st.spinner("시뮬레이션 계산 중..."):
    result = run_simulation(simulation_inputs, active_extractant_params)

profile_key = PROFILE_LABELS[st.session_state.ui_param_profile]
loading_pct = compute_loading_pct(result)
scope_assessment = build_scope_assessment(profile_key, simulation_inputs, result)

st.sidebar.markdown("---")
st.sidebar.header("📎 검증 커버리지")
st.sidebar.caption(f"활성 프로필: {st.session_state.ui_param_profile}")
if result["converged"]:
    st.sidebar.success(f"계산 수렴: {result['iterations']} iter")
else:
    st.sidebar.error(f"계산 미수렴: {result['iterations']} iter")
render_scope_assessment(st.sidebar, scope_assessment)

# =============================================================================
# TAB 1: 시뮬레이션 결과
# =============================================================================
with tab1:
    if not result["converged"]:
        st.error("현재 조건에서 다단 계산이 수렴하지 않았습니다. 결과 수치를 의사결정에 직접 사용하기 전에 입력 조건과 전략을 다시 확인하세요.")

    # --- 요약 메트릭 ---
    st.subheader("🎯 추출 결과 요약")
    colors = METAL_COLORS
    # 7개 금속 카드: 상단 4개 + 하단 3개
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
            st.metric("최대 로딩률", f"{loading_pct:.1f}%", delta="⚠️ 포화 근접", delta_color="inverse")
        else:
            st.metric("최대 로딩률", f"{loading_pct:.1f}%")
    with col_info4:
        if result["converged"]:
            st.metric("수렴 상태", "수렴", delta=f"{result['iterations']} iter")
        else:
            st.metric("수렴 상태", "미수렴", delta=f"{result['iterations']} iter", delta_color="inverse")

    st.markdown("---")
    st.subheader("🧭 현재 시뮬레이션 해석 가이드")
    st.caption(f"활성 프로필: **{st.session_state.ui_param_profile}**")
    render_scope_assessment(st, scope_assessment)
    if selected_preset_note:
        st.warning(f"프리셋 주의: {selected_preset_note}")

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
                     color="금속", color_discrete_sequence=[colors.get(m, "#636e72") for m in metals],
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
    color_list = [colors.get(m, "#636e72") for m in metals]
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
    pH_range = PH_PLOT_RANGE
    for i, metal in enumerate(metals):
        E_values = [
            extraction_efficiency(
                pH, metal, extractant, C_ext, temperature,
                extractant_params=active_extractant_params,
            )
            for pH in pH_range
        ]
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
        ph50_eff = get_effective_pH50(
            metal, extractant, C_ext, temperature,
            extractant_params=active_extractant_params,
        )
        k_eff = get_effective_k(
            metal, extractant, temperature,
            extractant_params=active_extractant_params,
        )
        params = active_extractant_params[extractant][metal]
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
        D = distribution_coefficient(
            ref_pH, mt_metal, extractant, C_ext, temperature,
            extractant_params=active_extractant_params,
        )
        # 로딩 감쇠 적용 (단일 금속 근사)
        from sx_simulator.extraction_isotherm import loading_damping_factor
        sim_org = {mt_metal: D * x}
        from sx_simulator.extraction_isotherm import calc_loading_fraction
        loading = calc_loading_fraction(
            sim_org, extractant, C_ext, [mt_metal],
            extractant_params=active_extractant_params,
        )
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
    st.info(
        f"비교 조건: 공통 C_ext {C_ext:.4f} M, {n_stages}단, {temperature:.0f}°C, SO₄²⁻ {C_sulfate:.3f} M"
    )

    results_compare = run_compare_simulations(
        simulation_inputs,
        compare_pH,
        active_extractant_params,
    )

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
        E_vals = [
            extraction_efficiency(
                pH, selected_metal, ext, ext_c, temperature,
                extractant_params=active_extractant_params,
            )
            for pH in pH_range
        ]
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
    st.subheader(f"📐 시스템 수식 및 메커니즘 알고리즘 ({extractant}, {C_ext}M, {temperature:.0f}°C)")
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
        params = active_extractant_params[extractant][metal]
        ph50_eff = get_effective_pH50(
            metal, extractant, C_ext, temperature,
            extractant_params=active_extractant_params,
        )
        k_eff = get_effective_k(
            metal, extractant, temperature,
            extractant_params=active_extractant_params,
        )
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

            st.markdown("**Phase 3: 금속 추출제 경합 보정 (기본 활성)**")
            st.latex(r"D_{M}^{\text{adj}} = D_{M}^{\text{sig}} \times \left(\frac{[\overline{\text{HA}}]_{\text{free}}}{C_{\text{ext}}}\right)^{n_{\text{eff}}}")
            st.caption("유기상 잔여 액체 추출제량에 비례하여 높은 로딩(Loading) 시 다핵 착물 형성에 따른 추출 효율을 동적 감소시킵니다.")
            # 해당 금속의 미니 Isotherm 차트
            pHs = [x * 0.1 for x in range(10, 101)]
            Es = [
                extraction_efficiency(
                    pH, metal, extractant, C_ext, temperature,
                    extractant_params=active_extractant_params,
                )
                for pH in pHs
            ]
            fig_mini = go.Figure()
            fig_mini.add_trace(go.Scatter(x=pHs, y=Es, mode='lines',
                line=dict(color=color_list[metals.index(metal)], width=2)))
            if target_pH:
                E_at_target = extraction_efficiency(
                    target_pH, metal, extractant, C_ext, temperature,
                    extractant_params=active_extractant_params,
                )
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
            E_val = extraction_efficiency(
                target_pH, metal, extractant, C_ext, temperature,
                extractant_params=active_extractant_params,
            )
            D_val = distribution_coefficient(
                target_pH, metal, extractant, C_ext,
                temperature=temperature,
                extractant_params=active_extractant_params,
            )
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
        n_H = active_extractant_params[extractant][metal]["n_H"]
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

    st.markdown("**Phase 3: 수계 금속 종분화(Speciation) 효과 (기본 활성)**")
    st.latex(r"[\text{H}^+]_{\text{hydrolysis}} \approx \sum_M K_{MOH} \cdot [\text{M}^{n+}] \cdot [\text{OH}^-]")
    st.caption("고 pH 역외 접근 시 금속 수산화물 착물(MOH⁺ 등)이 생성되며 양성자를 방출(OH⁻ 소비)해 직접적인 완충제로 작용함을 수지에 기본 반영했습니다.")

    if pH_mode == "목표 pH (자동 NaOH)":
        st.info(f"**목표 pH 모드**: pH = {target_pH if not staged_pHs else staged_pHs} 를 유지하기 위해 필요한 NaOH를 자동 역산합니다.\n\n"
                f"NaOH 필요량 = $[H^+]_{{in}} \\cdot Q_{{aq}} + \\Sigma(n_H \\cdot \\Delta C_M \\cdot Q_{{aq}}) - [H^+]_{{target}} \\cdot Q_{{aq}}$")

    # ---------------------------------------------------------------
    # 5. 시스템 입력 보정 로직 (황산염)
    # ---------------------------------------------------------------
    st.markdown("---")
    st.markdown("### 5️⃣ 시스템 입력 보정 로직 (황산염 농도 자동 연산)")
    st.markdown("피드 용액에 포함된 양이온 금속들이 형태상 금속 황산염($M\\text{SO}_4$, 또는 $M_2(\\text{SO}_4)_3$)으로 투입된다고 가정하고 화학양론비에 따라 **초기 총 음이온 농도**를 연산합니다.")
    st.latex(r"C_{\text{sulfate}} (\text{M}) = \sum_{M} \left( \frac{C_{M,\text{feed}}}{\text{MW}_M} \times \frac{\text{valency}_M}{2} \right)")
    st.info(f"계산된 피드 용액 내 황산염 음이온($\\text{{SO}}_4^{{2-}}$) 농도: **{C_sulfate:.4f} M**")

    # ---------------------------------------------------------------
    # 6. 다단 교대 수렴 (McCabe-Thiele Simulator)
    # ---------------------------------------------------------------
    st.markdown("---")
    st.markdown("### 6️⃣ 다단 연속 추출 순환 해법 (McCabe-Thiele Iterator)")
    st.markdown("**교류 방식 (Counter-current) 반복 알고리즘:** 다단(Multi-stage) Mixer-Settler 시스템의 작동 원리에 따라 수계(Aqueous)는 $\\text{Stage } 1 \\rightarrow N$ 방향으로, 유기계(Organic)는 $\\text{Stage } N \\rightarrow 1$ 방향으로 엇갈리며 흐릅니다.")
    st.markdown(f"시뮬레이터는 이 동적 흐름 모델에서 {n_stages}개의 각 Stage가 평형에 도달할 때까지 **교대 수렴(Alternating Convergence)** 반복 루프를 구동합니다. 이전 Stage와 다음 Stage에서 넘어오는 물질들을 교차 연산하여 물질수지가 모두 일치할 때까지 수십 차례 왕복 연산합니다.")
    st.latex(r"C_{M,\text{org}}^{\text{stage } i} = f_{\text{eq}}\!\left( C_{M,\text{aq}}^{\text{stage } i}, \text{pH}_i \right)")
    st.latex(r"Q_{\text{aq}} \cdot C_{M,\text{aq}}^{\text{stage } i} + Q_{\text{org}} \cdot C_{M,\text{org}}^{\text{stage } i} = Q_{\text{aq}} \cdot C_{M,\text{aq}}^{\text{stage } i-1} + Q_{\text{org}} \cdot C_{M,\text{org}}^{\text{stage } i+1}")
    
    # ---------------------------------------------------------------
    # ---------------------------------------------------------------
    # 7. 추출제 한계 (Phase 3)
    # ---------------------------------------------------------------
    st.markdown("---")
    st.markdown("### 7️⃣ 추출 농도 물리적 한계 (Phase 3: 다핵 착물 모델)")
    st.markdown("금속 로딩량이 치솟아 남은 자유 추출제($\\overline{\\text{HA}}$)가 모자랄 때, $D_M$ 지표가 급격하게 저하되도록 억제하는 고도화 수식입니다. (v2.0 다핵 착물 $n_{\\text{eff}}$ 완화 적용)")
    st.latex(r"[\overline{\text{HA}}]_{\text{free}} = C_{\text{ext}} - \sum_{M} n_{\text{eff},M} \cdot C_{M,\text{org}}")
    st.latex(r"D_{M}^{\text{adj}} = D_{M}^{\text{sig}} \times \left( \max\left(10^{-4}, \frac{[\overline{\text{HA}}]_{\text{free}}}{C_{\text{ext}}}\right) \right)^{n_{\text{eff}}}")
    st.caption("고농도 피드 유입 시 모든 금속이 100% 추출되는 비현실적 결과를 방지하고 가용 유기제 내에서 실제적인 상호 조율점(Crowding Out)을 찾습니다.")

    # ---------------------------------------------------------------
    # 8. 추출률 공식
    # ---------------------------------------------------------------
    st.markdown("---")
    st.markdown("### 8️⃣ 성능 지표")
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
                m: distribution_coefficient(
                    target_pH, m, extractant, C_ext,
                    temperature=temperature,
                    extractant_params=active_extractant_params,
                )
                for m in metals
            }
            for i, m1 in enumerate(metals):
                for m2 in metals[i+1:]:
                    if D_vals[m1] <= 1e-10 and D_vals[m2] <= 1e-10:
                        alpha_str = "N/A"
                    elif D_vals[m2] <= 1e-10:
                        alpha_str = "> 10⁶"
                    else:
                        alpha_str = f"{D_vals[m1] / D_vals[m2]:.2f}"
                    sep_data.append({"M₁/M₂": f"{m1}/{m2}", "α": alpha_str})
            st.dataframe(pd.DataFrame(sep_data), use_container_width=True, hide_index=True)

    # ---------------------------------------------------------------
    # 6. 파라미터 요약 테이블
    # ---------------------------------------------------------------
    st.markdown("---")
    st.markdown("### 📋 현재 파라미터 요약")
    param_rows = []
    for metal in metals:
        p = active_extractant_params[extractant][metal]
        ph50_eff = get_effective_pH50(
            metal, extractant, C_ext, temperature,
            extractant_params=active_extractant_params,
        )
        k_eff = get_effective_k(
            metal, extractant, temperature,
            extractant_params=active_extractant_params,
        )
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

                        fit_result = fit_sigmoid(pH_vals, E_vals)

                        if not fit_result["success"]:
                            st.error(f"피팅 실패: {fit_result.get('error', 'Unknown error')}")
                            continue

                        # --- 피팅 결과 표시 ---
                        res_col1, res_col2 = st.columns([1, 1])

                        with res_col1:
                            st.markdown("**피팅 파라미터:**")
                            params = fit_result["params"]
                            errors = fit_result["errors"]

                            param_df = pd.DataFrame([
                                {"파라미터": "pH\u2085\u2080", "값": params["pH50"],
                                 "±95%CI": errors["pH50_err"]},
                                {"파라미터": "k", "값": params["k"],
                                 "±95%CI": errors["k_err"]},
                                {"파라미터": "E_max (%)", "값": params["E_max"],
                                 "±95%CI": errors["E_max_err"]},
                            ])
                            st.dataframe(param_df, use_container_width=True, hide_index=True)
                            st.metric("R\u00b2", f"{fit_result['r_squared']:.4f}")

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
                                st.session_state.custom_params[extractant][metal_name]["pH50"] = params["pH50"]
                                st.session_state.custom_params[extractant][metal_name]["k"] = params["k"]
                                st.session_state.custom_params[extractant][metal_name]["E_max"] = params["E_max"]
                                st.success(f"{metal_name} 파라미터가 적용되었습니다! 페이지를 새로고침하면 시뮬레이션에 반영됩니다.")
                                st.rerun()

                        st.markdown("---")
    else:
        st.info("👆 CSV 파일을 업로드하여 피팅을 시작하세요.")


# =============================================================================
# TAB 9: 용어 해설 (Glossary)
# =============================================================================
with tab9:
    st.markdown(
        load_markdown_doc(
            os.path.join("docs", "glossary.md"),
            "용어 해설 문서를 찾을 수 없습니다.",
        )
    )

# =============================================================================
# TAB 10: 파라미터 및 문헌 확인
# =============================================================================
with tab10:
    st.subheader("📚 시스템 파라미터 및 참고 문헌 (References)")
    st.markdown("이 시뮬레이터에 적용된 열역학 모델 및 파라미터 값들의 출처와 세부 데이터베이스 표를 제공합니다.")
    
    # 1. 참고 문헌
    st.markdown(
        load_markdown_doc(
            os.path.join("docs", "references.md"),
            "참고 문헌 문서를 찾을 수 없습니다.",
        )
    )

    # 2. 파라미터 데이터
    st.markdown("---")
    st.markdown("### 🔬 추출제 초기 파라미터 (Extractant Database)")
    
    t10_col1, t10_col2 = st.columns(2)
    with t10_col1:
        st.markdown("**Cyanex 272 기초 파라미터**")
        df_cyanex = pd.DataFrame(active_extractant_params["Cyanex 272"]).T
        st.dataframe(df_cyanex, use_container_width=True)
    with t10_col2:
        st.markdown("**D2EHPA 기초 파라미터**")
        df_d2ehpa = pd.DataFrame(active_extractant_params["D2EHPA"]).T
        st.dataframe(df_d2ehpa, use_container_width=True)
    st.caption("위 값들은 `$T_{\\text{ref}} = 25^{\\circ}\\text{C}$`, 기준 농도 단위에서 도출된 이분법/시그모이드 파라미터의 초기 기준점이며, 시뮬레이션 환경에 맞춰 좌측 사이드바 설정대로 수동/동적 보정이 이루어집니다.")

    st.markdown("---")
    st.markdown("### 🔬 화학 종분화 반응 상수 (Speciation Constants)")
    st.markdown("고-pH 환경에서 금속 수산화물($M\\text{OH}^+$, $M(\\text{OH})_2$) 생성 및 점유를 결정하는 $\\beta_1/K_w$ 또는 $K_{MOH}$, 그리고 황산염 결합 상수들입니다.")
    df_spec = pd.DataFrame(SPECIATION_CONSTANTS).T
    st.dataframe(df_spec, use_container_width=True)

# =============================================================================
# TAB 11: 데이터 피팅 검증 이력
# =============================================================================
with tab11:
    st.markdown(
        load_markdown_doc(
            os.path.join("docs", "validation_history.md"),
            "데이터 피팅 검증 이력 문서를 찾을 수 없습니다.",
        )
    )

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
