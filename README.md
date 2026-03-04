# ⚗️ SX Simulator — Li/Ni/Co/Mn Mixer-Settler

MSE Thermodynamic Framework 기반 Li/Ni/Co/Mn 용매추출(SX) 시뮬레이터.

## Features

- **pH-기반 추출 Isotherm**: 시그모이드 함수를 이용한 금속별 추출률 예측
- **다단 역류 시뮬레이션**: Counter-current mixer-settler 공정 모사
- **목표 pH 모드**: 원하는 pH 유지에 필요한 NaOH 자동 계산
- **추출제 비교**: Cyanex 272 vs D2EHPA 성능 비교
- **모델 수식 확인**: 현재 파라미터가 대입된 수학 모델을 LaTeX로 표시
- **대화형 대시보드**: Streamlit + Plotly 기반 실시간 시각화

## Quick Start

```bash
pip install -r requirements.txt
streamlit run sx_dashboard.py
```

## Project Structure

```
sx_simulator/          # 시뮬레이션 엔진
├── config.py          # 금속/추출제 파라미터
├── extraction_isotherm.py  # Isotherm 모델
├── single_stage.py    # 단일 stage 솔버
└── multistage_sx.py   # 다단 역류 솔버
sx_dashboard.py        # Streamlit 웹 대시보드
sx_simulator_app.py    # CLI 실행 스크립트
```

## References

- ALTA 2024 / Wang et al. (2002, 2004, 2006)
- Mixed-Solvent Electrolyte (MSE) Model Framework
