# 🧪 현장 실험 데이터 피팅 및 검증 이력

연구원 제공 실제 현장 실험 데이터(Isd 108 희석제, 사포닌화 등)를 바탕으로 시뮬레이터 파라미터를 보정(Data Fitting)한 주요 연혁을 기록합니다.

---

### 📌 [post-v2.1.2] 2026-03-12 : 전문가 리뷰 후 검증/claim 정렬

#### 1. 이번 변경의 성격

- 이번 단계는 파라미터 재피팅이 아니라, 현재 로컬 코드와 사용자-facing claim 사이의 간격을 줄이기 위한 안전장치 보강이다.
- 따라서 legacy release gate는 유지하되, 어떤 경로가 `validated`, 어떤 경로가 `diagnostic`인지 더 명시적으로 드러내는 쪽에 초점을 맞췄다.

#### 2. 로컬 재실행 기준 확인한 사실

- `python3 test_verification.py`
  - `legacy_premixed_target_pH` 기준은 여전히 PASS
  - raw-feed 계열 basis는 계속 `diagnostic only`
- `python3 validation_test.py`
  - Step 1, 2, 4는 PASS
  - Step 3 multistage counter-current logic는 현재 Co balance diff 약 `0.000446 g`로 FAIL

#### 3. 이번에 반영한 안전장치

- 웹앱 설명 문구에서 `full MSE`처럼 읽힐 수 있는 표현을 축소
- `McCabe-Thiele` 탭을 formal 설계 도구가 아니라 해석용 보조 시각화로 명시
- `고정 NaOH + 사포니피케이션` 경로는 raw-feed replay/진단용이라는 경고 강화
- `D2EHPA + 목표 pH + 고황산염 + 다단` 조합은 검증 취약영역 경고 추가
- v2.1.1 종합검증 보고서 상단에 live 재실행 불일치 경고 추가

#### 4. 현재 해석 원칙

- 공식 회귀 보호 기준은 계속 `legacy_premixed_target_pH`
- raw-feed 결과는 NaOH 농도/주입 메타데이터 확보 전까지 정식 validation으로 승격하지 않음
- `physical_v2` 사포니피케이션 경로는 구현은 존재하지만 현재 field replay 정확도 기준으로는 experimental/shadow 해석이 더 적절함

---

### 📌 [v2.1.2] 2026-03-11 : 알칼리 인터페이스 정리 + Li/Co 잔여오차 미세 보정

#### 1. 알칼리 경로 명시화

- `aqueous_direct`와 `fresh_organic_saponification`을 코드 레벨에서 분리
- 사포니피케이션 경로를
  - `physical_v2`
  - `legacy_equivalent_target`
  로 명시적으로 구분
- 다단 사포니피케이션에 대해 `fresh sap inventory`를 실제 counter-current state로 전달하는 physical path를 추가

#### 2. field-calibrated 파라미터 미세 조정

- `Cyanex 272 / Li`
  - `pH50: 8.00 -> 8.10`
  - `k: 2.0 -> 1.8`
- `D2EHPA / Li`
  - `pH50: 6.50 -> 6.30`
- `D2EHPA / Co`
  - `E_max: 99.5 -> 99.0`

#### 3. 검증 결과

- legacy release gate 유지
  - Li MAE `0.816 -> 0.788 g/L`
  - Co MAE `0.060 -> 0.047 g/L`
- relative error 기준 개선
  - Cyanex 272 legacy aggregate: `14.2% -> 13.8%`
  - D2EHPA legacy aggregate: `27.1% -> 25.3%`
  - Cyanex 272 raw-feed fixed sap: `23.7% -> 23.4%`
  - D2EHPA raw-feed fixed sap: `51.0% -> 50.8%`

#### 4. 해석 주의사항

- `raw_feed_physical_saponification_v2`는 구현은 완료됐지만 아직 field 정확도가 충분하지 않다.
- 따라서 현재는 **experimental / shadow only**로 유지하고, 공식 release gate는 계속 `legacy_premixed_target_pH`를 사용한다.

---

### 📌 [v1.9.0] 2026-03-06 : CoSX 및 IMSX 현장 데이터(Data1~6) 기반 피팅

#### 1. Cyanex 272 (Ni 과추출 방지)

- **배경**: Data1~3 (CoSX) 검증 시 Ni이 이론 모델 대비 덜 추출되는 현상 관찰 (Isd 108 희석제 영향으로 추정)
- **조정 내역**: `Ni` 추출 시작점 `pH50` 파라미터를 **5.8 -> 6.3**으로 상향
- **결과**: Ni의 시뮬레이션 오차율 대폭 완화

#### 2. D2EHPA (Mg 저추출 방지)

- **배경**: Data4~6 (IMSX) 검증 시 D2EHPA에서 Mg이 조기 추출되는 현상 관찰
- **조정 내역**: `Mg` 추출 시작점 `pH50` 파라미터를 **4.5 -> 3.2**로 하향
- **결과**: Mg 피팅 오차율 0%에 근접하는 정확도 확보

---

### 📊 검증 데이터 프리셋 로드

좌측 사이드바의 **[🧪 현장 실험 Data 프리셋]** 메뉴를 통해, v1.9.0 피팅에 사용되었던 Data1~6 조건을 원클릭으로 호출하여 즉각 연산해 볼 수 있습니다.
