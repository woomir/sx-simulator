#!/usr/bin/env python3
"""
SX Simulator Web Dashboard
===========================
Streamlit 기반 Li/Ni/Co/Mn Mixer-Settler 용매추출 시뮬레이션 대시보드.

실행: streamlit run sx_dashboard.py
"""

import sys, os
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
import pandas as pd

from sx_simulator.config import (DEFAULT_METALS,
                                  T_REF,
                                  get_parameter_profile,
                                  SPECIATION_CONSTANTS)
from sx_simulator.dashboard_service import (
    SimulationInputs,
    build_scope_assessment,
    compute_loading_pct,
    run_simulation,
)
from sx_simulator.dashboard_tabs import (
    render_compare_tab,
    render_detail_tab,
    render_fitting_tab,
    render_formula_tab,
    render_isotherm_tab,
    render_mccabe_thiele_tab,
    render_results_tab,
)
from sx_simulator.datasets import (
    DASHBOARD_PRESETS as PRESETS,
    KNOWN_PRESET_NOTES,
    calc_sulfate_from_feed,
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
    render_results_tab(
        result=result,
        metals=metals,
        loading_pct=loading_pct,
        scope_assessment=scope_assessment,
        profile_label=st.session_state.ui_param_profile,
        selected_preset_note=selected_preset_note,
        C_aq_feed=C_aq_feed,
        pH_feed=pH_feed,
        n_stages=n_stages,
        metal_colors=METAL_COLORS,
        scope_renderer=render_scope_assessment,
    )

# =============================================================================
# TAB 2: pH Isotherm 곡선
# =============================================================================
with tab3:
    render_isotherm_tab(
        metals=metals,
        extractant=extractant,
        C_ext=C_ext,
        temperature=temperature,
        target_pH=target_pH,
        active_extractant_params=active_extractant_params,
        metal_colors=METAL_COLORS,
        ph_plot_range=PH_PLOT_RANGE,
    )

# =============================================================================
# TAB 8: McCabe-Thiele 다이어그램
# =============================================================================
with tab8:
    render_mccabe_thiele_tab(
        metals=metals,
        result=result,
        extractant=extractant,
        C_ext=C_ext,
        temperature=temperature,
        C_aq_feed=C_aq_feed,
        n_stages=n_stages,
        active_extractant_params=active_extractant_params,
    )

# =============================================================================
# TAB 4: 추출제 비교
# =============================================================================
with tab4:
    render_compare_tab(
        metals=metals,
        simulation_inputs=simulation_inputs,
        active_extractant_params=active_extractant_params,
        temperature=temperature,
        C_ext=C_ext,
        n_stages=n_stages,
        C_sulfate=C_sulfate,
        ph_plot_range=PH_PLOT_RANGE,
    )

# =============================================================================
# TAB 4: 상세 데이터
# =============================================================================
with tab5:
    render_detail_tab(
        result=result,
        metals=metals,
        C_aq_feed=C_aq_feed,
        Q_aq=Q_aq,
        Q_org=Q_org,
    )


# =============================================================================
# TAB 5 (tab5): 모델 수식 확인
# =============================================================================
with tab2:
    render_formula_tab(
        extractant=extractant,
        C_ext=C_ext,
        temperature=temperature,
        T_REF=T_REF,
        metals=metals,
        active_extractant_params=active_extractant_params,
        target_pH=target_pH,
        staged_pHs=staged_pHs,
        pH_mode=pH_mode,
        Q_aq=Q_aq,
        Q_org=Q_org,
        C_sulfate=C_sulfate,
        n_stages=n_stages,
        metal_colors=METAL_COLORS,
    )


# =============================================================================
# TAB 6: 데이터 피팅 (Experimental Data Fitting)
# =============================================================================
with tab6:
    render_fitting_tab(
        extractant=extractant,
        default_metals=DEFAULT_METALS,
    )


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
