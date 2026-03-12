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

```text
sx_dashboard.py                 # Streamlit 웹 대시보드 (호환 경로 유지)
sx_simulator_app.py             # CLI 실행 스크립트
sx_simulator/                   # 시뮬레이션 엔진
tests/                          # 자동 테스트
scripts/analysis/               # 수동 분석/비교 스크립트
docs/                           # 사용자 문서, 검토 기록, literature intake
docs/architecture/              # 저장소 구조 가이드
docs/analysis/                  # 기술 메모와 분석 요약
deliverables/                   # 발표자료, 보고 산출물, 초안 작업공간
tmp/                            # 로컬 임시 산출물 (gitignore)
```

## Organization Rules

- 웹앱 호환성을 위해 `sx_dashboard.py`와 `docs/user_manual.md` 등 웹앱이 직접 읽는 경로는 유지합니다.
- 문헌 PDF와 handoff 기록은 `docs/literature/` 아래에서 계속 관리합니다.
- 보고자료와 발표 초안은 루트가 아니라 `deliverables/` 아래에 둡니다.
- 일회성 분석 스크립트는 루트가 아니라 `scripts/analysis/` 아래에 둡니다.

자세한 규칙은 `AGENTS.md`와 `docs/architecture/repository_organization.md`를 참고합니다.

## References

- ALTA 2024 / Wang et al. (2002, 2004, 2006)
- Mixed-Solvent Electrolyte (MSE) Model Framework
