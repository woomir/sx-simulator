# PROJECT_STRUCTURE

이 문서는 이 저장소에 처음 들어온 개발자가 프로젝트의 목적, 구조, 실행 흐름, 주요 파일 역할을 빠르게 파악하기 위한 안내서다.

## 1. 프로젝트 목적

이 프로젝트는 Li/Ni/Co/Mn 및 일부 불순물(Ca/Mg/Zn)을 포함한 황산염계 수용액에 대해 다단 용매추출(SX, Solvent Extraction) 공정을 시뮬레이션하는 Streamlit 웹앱이다.

핵심 목적:

- 공정 조건에 따른 금속별 추출률 예측
- 다단 Mixer-Settler 역류 공정 거동 계산
- 현장 데이터와의 비교/보정
- pH isotherm, McCabe-Thiele, 추출제 비교 등 공정 해석 지원

## 2. 실행 진입점

### 웹앱

- `sx_dashboard.py`
  - Streamlit 대시보드 진입점
  - 사용자 입력 수집, 탭 렌더링, 시뮬레이션 실행, 결과 표시 담당

실행 예:

```bash
streamlit run sx_dashboard.py
```

### 검증 스크립트

- `validation_test.py`
  - 구조 검증
  - 물질수지, 단일단/다단 일관성 확인
- `test_verification.py`
  - 현장 데이터(Data1~6)와 시뮬레이션 결과 비교
  - 현재는 보고서형 검증에 가깝고, 허용오차 기반 자동 실패 테스트는 제한적

## 3. 런타임 구조

실행 흐름은 아래와 같다.

```text
sx_dashboard.py
  -> sx_simulator/dashboard_service.py
      -> sx_simulator/dashboard_tabs.py
      -> sx_simulator/multistage_sx.py
          -> sx_simulator/single_stage.py
              -> sx_simulator/extraction_isotherm.py
                  -> sx_simulator/config.py
                  -> sx_simulator/datasets.py
```

### 각 계층의 역할

#### `sx_dashboard.py`

- Streamlit 위젯과 탭 UI
- 세션 파라미터 관리
- 그래프/표 렌더링

#### `sx_simulator/dashboard_service.py`

- 대시보드 입력을 `SimulationInputs`로 묶음
- 메인 계산/비교 계산 kwargs 조립
- 검증 커버리지 평가
- UI와 엔진 사이의 얇은 서비스 계층

#### `sx_simulator/dashboard_tabs.py`

- 주요 동적 탭 렌더링 함수 모음
- 시뮬레이션 결과, Isotherm, McCabe-Thiele, 비교, 상세 데이터 탭 담당
- 메인 대시보드 파일의 UI 블록을 모듈 단위로 분리

#### `sx_simulator/multistage_sx.py`

- 다단 역류 시뮬레이션 엔진
- stage 간 유기상/수계 연결
- 수렴 반복과 relaxation 처리

#### `sx_simulator/single_stage.py`

- 단일 stage 평형 계산
- 목표 pH 모드 / 고정 NaOH 모드 처리
- 경쟁 추출 및 양성자 수지 반영

#### `sx_simulator/extraction_isotherm.py`

- pH 기반 시그모이드 추출률 모델
- 분배계수 계산
- 로딩 감쇠
- 경쟁 추출 보정
- 금속별 H+ 방출량 계산

#### `sx_simulator/config.py`

- 기본 금속 목록
- 추출제 파라미터
- 상수
- 검증 범위 설정
- 프로필(`field_calibrated`, `literature_default`) 제공

#### `sx_simulator/datasets.py`

- 대시보드 프리셋 데이터
- 검증 데이터와 공용으로 쓰는 조성 정보
- 피드 조성 기반 황산염 자동 계산 helper

#### `sx_simulator/fitting.py`

- 실험 데이터 기반 시그모이드 파라미터 피팅
- 피팅 탭에서 사용

## 4. 주요 데이터 흐름

### 1. 사용자 입력

`sx_dashboard.py`의 사이드바에서 다음을 입력한다.

- 피드 금속 농도
- 피드 pH
- 수계/유기계 유량
- 추출제와 추출제 농도
- 운전 온도
- stage 수
- 목표 pH 또는 고정 NaOH 조건
- 파라미터 프로필 및 수동 파라미터 편집값

### 2. 입력 조립

입력은 `dashboard_service.py`의 `SimulationInputs` 객체로 묶인다.

### 3. 메인 계산

`run_simulation()`이 `solve_multistage_countercurrent()` 호출용 kwargs를 조립하고, 다단 계산을 실행한다.

### 4. 결과 해석

대시보드는 결과를 바탕으로 다음을 렌더링한다.

- 금속별 추출률
- 후액 농도
- pH profile
- NaOH profile
- 물질수지 표
- Isotherm
- McCabe-Thiele
- 비교 차트

## 5. 사용자 관점 주요 탭

- `📊 시뮬레이션 결과`
  - 핵심 결과 요약
- `📐 수식 및 메커니즘 알고리즘`
  - 현재 설정값이 반영된 식 설명
- `📈 pH Isotherm`
  - 금속별 pH-추출률 곡선
- `📉 McCabe-Thiele`
  - 조작선/평형선 시각화
- `🔬 추출제 비교`
  - Cyanex 272 vs D2EHPA 비교
- `📋 상세 데이터`
  - stage별 결과와 물질수지
- `📝 데이터 피팅`
  - 실험 데이터 기반 파라미터 추정
- `📖 용어 해설`
  - 용어집
- `📚 파라미터 및 문헌`
  - 파라미터 표 및 참고 문헌
- `🧪 데이터 피팅 검증 이력`
  - 현장 데이터 기반 보정 기록
- `📜 변경 이력`
  - 사용자 기능 변경 로그
- `📘 사용자 매뉴얼`
  - 사용법 문서

## 6. docs 폴더 역할

`docs/`에는 크게 세 종류 문서가 있다.

- 사용자 문서
  - 예: `docs/user_manual.md`
  - 예: `docs/glossary.md`
  - 예: `docs/references.md`
  - 예: `docs/validation_history.md`
- 검토/분석 문서
  - 예: 일반화 검토, 검증 보고서
- 개발자 리팩토링 기록
  - 예: `docs/v2.1.1_리팩토링_기록.md`
  - 예: `docs/v2.1.1_2차_리팩토링_기록.md`
  - 예: `docs/v2.1.1_3차_리팩토링_기록.md`
  - 예: `docs/v2.1.1_4차_리팩토링_기록.md`
  - 예: `docs/v2.1.1_5차_리팩토링_기록.md`

원칙:

- 사용자 기능 변경은 `CHANGELOG.md`
- 내부 구조 변경은 refactoring 기록 문서

## 7. 현재 모델의 성격

이 프로젝트는 완전한 1원리 열역학 해석기라기보다 다음 성격에 가깝다.

- 시그모이드 기반 준경험 추출 모델
- 경쟁 추출과 로딩 효과를 반영한 경량 메커니즘 보정
- 종분화와 양성자 수지를 포함한 공정 시뮬레이터
- 일부 현장 데이터에 맞춘 보정 프로필 포함

따라서 "조건 비교와 경향 해석"에는 강하고, "모든 현장에 대한 범용 정량 예측"은 아직 제한이 있다.

## 8. 현재 알려진 한계

- 일부 현장 데이터에서 Ni/Co/Li 오차가 남아 있다.
- 활동도/비이상성은 제한적으로만 반영된다.
- `test_verification.py`는 자동 실패 기준보다 보고서 출력 비중이 크다.
- `sx_dashboard.py`는 여전히 큰 파일이며, 탭 렌더러 분리는 앞으로의 과제다.

## 9. 앞으로 코드를 볼 때 추천 순서

처음 보는 개발자라면 아래 순서로 읽는 것이 가장 빠르다.

1. `PROJECT_STRUCTURE.md`
2. `README.md`
3. `sx_dashboard.py`
4. `sx_simulator/dashboard_service.py`
5. `sx_simulator/dashboard_tabs.py`
6. `sx_simulator/multistage_sx.py`
7. `sx_simulator/single_stage.py`
8. `sx_simulator/extraction_isotherm.py`
9. `sx_simulator/config.py`
10. `validation_test.py`
11. `test_verification.py`

## 10. 최근 내부 구조 변경 기록

최근 리팩토링 관련 문서:

- `docs/v2.1.1_리팩토링_기록.md`
- `docs/v2.1.1_2차_리팩토링_기록.md`
- `docs/v2.1.1_3차_리팩토링_기록.md`
- `docs/v2.1.1_4차_리팩토링_기록.md`
- `docs/v2.1.1_5차_리팩토링_기록.md`

이 문서들을 함께 보면 왜 구조가 바뀌었는지 추적할 수 있다.
