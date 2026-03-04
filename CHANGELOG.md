# Changelog

이 프로젝트의 모든 주요 변경 사항을 기록합니다.

형식: [Semantic Versioning](https://semver.org/lang/ko/)

## [v1.2.0] - 2026-03-04

### ✨ 추가

- 🌡️ **온도 의존성 모델** 추가
  - 반-경험적 온도 보정: pH₅₀(T) = pH₅₀(T_ref) + β·(T - T_ref)
  - 시그모이드 기울기 온도 보정: k(T) = k_ref · exp(γ·(T - T_ref))
  - 금속별 β, γ 온도 계수(사용자 조정 가능)
  - 대시보드 사이드바에 온도 슬라이더(10~60°C) 추가
  - 모델 수식 탭에 온도 보정 공식 표시

## [v1.1.0] - 2026-03-04

### ✨ 추가

- 📐 **모델 수식 탭** 추가
  - LaTeX 렌더링으로 시뮬레이션 수식 표시
  - 파라미터 변경 시 수식 실시간 업데이트
  - 금속별 미니 Isotherm 차트 표시
- 🚀 **GitHub 배포** 준비
  - `requirements.txt`, `.gitignore`, `README.md` 생성
  - Streamlit Community Cloud 배포 지원

## [v1.0.0] - 2026-03-03

### ✨ 추가

- 🧪 **핵심 시뮬레이션 엔진**
  - 시그모이드 pH-추출률 Isotherm 모델
  - 추출제 농도 보정 (α·log₁₀(C_ext/C_ref))
  - 단일 Stage 물질수지 + pH 수지
  - 역류 다단 Mixer-Settler 시뮬레이션
  - 목표 pH 모드 (NaOH 자동 역산)
- 📊 **Streamlit 대시보드**
  - 추출 결과 요약 (메트릭 카드)
  - Stage별 pH/NaOH 프로파일
  - pH-추출률 Isotherm 곡선
  - Cyanex 272 vs D2EHPA 비교
  - 모델 파라미터 편집 UI
- 🔬 **지원 금속/추출제**
  - 금속: Li, Ni, Co, Mn
  - 추출제: Cyanex 272, D2EHPA
