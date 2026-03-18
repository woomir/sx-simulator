# SX 시뮬레이터 종합 전문가 검토 보고서

**작성일:** 2026-03-16
**대상 버전:** main 브랜치 (commit c86c30e)
**검토 범위:** 평형 / 플로우시트 / 공정 / 시뮬레이션 4개 전문 영역 통합

---

## 1. 요약 (Executive Summary)

### 핵심 발견사항 5가지

1. **추출제 농도 범위 미스매치**: Ni-Cy 현장 데이터는 C_ext = 0.473 M이나, `VALIDATED_FIELD_WINDOW`는 0.605–0.631 M 범위만 커버한다. 현장 데이터가 검증 범위의 25% 아래에 있어 파라미터 외삽이 불가피하다.

2. **프리셋 과적합 (Overfitting)**: Data1~6 (5단, C_ext ≈ 0.63 M)에 맞춘 `SITE_PARAMETER_OVERRIDES`가 Ni-Cy (3단, C_ext = 0.473 M, O/A 6~20:1)에 적용될 때 Ni 후액 MAE 26 g/L의 과추출 오차를 생성한다. 프리셋 MAE < 1 g/L과의 극적인 격차는 파라미터 전이 실패를 의미한다.

3. **사포니피케이션 모델 미성숙**: 현장의 NaOH pre-saponified organic 조건을 재현하는 `physical_v2` 경로가 핵심 병목이다. `SAPONIFICATION_P_H50_SHIFT`의 Ni 계수 0.65는 pH50을 유효 5.90까지 낮추지만, 극고농도 Ni(62–75 g/L)에서의 로딩 포화와 결합된 효과를 과대평가한다.

4. **극고농도 Ni의 이온강도/activity 효과 미반영**: Ni 62–75 g/L (1.06–1.28 M)에서 이온강도가 매우 높아 (I > 3 M 추정) 자유 금속 분율이 급감하나, 현재 종분화 모델은 1:1 착물만 고려하며 activity coefficient 보정이 없다 (`extraction_isotherm.py:343–383`).

5. **비정상상태 데이터 혼재**: Ni-Cy 초기 운전일(25.03.11–03.13)의 Feed pH 변동(4.12→6.28), 극단적 1단 Ni 잔류(54.8 g/L), 그리고 NaOH 유량 결측은 비정상 상태를 시사한다. 안정화 후 데이터(04월 이후)와 혼재되어 단일 파라미터 세트로의 적합이 불가능하다.

### 최우선 권장사항

| 우선순위 | 권장사항 | 기대 효과 |
|---------|---------|----------|
| P0 | Ni-Cy 농도 범위에 맞는 별도 보정 프로필 구축 | Ni MAE 26→5 g/L 이하 목표 |
| P0 | 비정상상태 운전일 식별 및 분리 | 보정 데이터셋 품질 확보 |
| P1 | 이온강도 기반 activity 보정 도입 | 극고농도 정확도 향상 |
| P1 | 사포니피케이션 계면화학 모델 고도화 | pH50 shift 물리적 근거 강화 |
| P2 | 3단 역류 수렴 특성 전용 검증 | 단수별 solver 안정성 확보 |

---

## 2. 평형 전문가 검토 (sx-equilibrium-expert)

### 2.1 Ni pH50 = 6.55 SITE_OVERRIDE의 물리적 근거

**코드 위치:** `config.py:252–255`

```python
SITE_PARAMETER_OVERRIDES = {
    ("Cyanex 272", "Ni"): {
        "pH50": 6.55,  # 현장 Data1~3 기준 저 pH Ni 과추출 완화
        "k": 2.0,
    },
```

**평가:**
- 문헌 기본값 pH50 = 5.8에서 +0.75 상향 조정. 이는 Cyanex 272-Ni의 일반적 문헌 범위(pH50 = 5.5–6.5, Devi 1998, Santanilla 2021) 내에 있어 **방향적으로는 방어 가능**하다.
- 그러나 이 값은 Data1~3 (C_ext = 0.631 M, 5단, feed Ni ≈ 28–34 g/L)에만 최적화된 field-calibration 결과이다.
- Ni-Cy 현장(C_ext = 0.473 M, feed Ni ≈ 62–78 g/L)에서는 추출제 농도 보정 `alpha * log10(C_ext/C_ref)`에 의해 pH50_eff가 추가로 상승하지만, alpha = 0.7과 C_ref = 0.5 M 기준으로 Δ = 0.7 × log10(0.473/0.5) = +0.017에 불과하다. 이 미세한 보정으로는 농도 범위 외삽의 비선형 효과를 커버할 수 없다.

**위험:** Ni-Cy 현장의 과추출은 pH50 자체가 아니라, 극고농도에서의 이온강도 효과와 로딩 포화 상호작용이 pH50 보정만으로는 포착되지 않기 때문일 가능성이 높다.

**권장:** C_ext = 0.473 M 데이터에 대한 독립적 pH50 보정값 도출. 단, working-rules에 따라 "단일 캘리브레이션 슬라이스에서 개선된 파라미터를 '물리적으로 검증됨'이라 주장하지 않을 것."

### 2.2 n_eff_multiplier 고로딩 방어의 물리적 타당성

**코드 위치:** `extraction_isotherm.py:549–572`

```python
# n_eff_multiplier: 1.0 (낮은 로딩) -> 0.35 (극한 로딩)
n_eff_multiplier = max(0.35, 1.0 - (loading_ratio_clamped ** 2))
```

**평가:**
- 고로딩 시 다핵 착물 형성(Co₂A₄, NiA₂(HA)₂ 등)으로 인해 유효 화학양론이 감소한다는 개념은 물리적으로 **합리적**이다 (Santanilla 2021: Ni는 더 많은 octahedral hydration 유지, Co는 더 적은 hydration).
- 그러나 `max(0.35, 1 - loading²)` 공식의 바닥값 0.35와 제곱 감쇠 형태는 **순수 수치적 편의**이다. 이를 뒷받침하는 직접적 실험 증거나 문헌 근거는 확인되지 않는다.
- Data5 (Ni 71 g/L, D2EHPA)에서는 이 방어 알고리즘이 극단적 과로딩을 방지하는 데 효과적이었으나, 이것이 Cyanex 272 Ni-Cy 조건(다른 추출제, 다른 농도)에서도 동일하게 작동하는지는 별도 검증이 필요하다.

**위험:** n_eff_multiplier가 경쟁추출 D 보정과 loading_damping_factor의 이중 감쇠와 결합될 때, 효과의 중복이나 상쇄가 발생할 수 있다.

### 2.3 loading_damping_factor 파라미터의 물리적 근거

**코드 위치:** `extraction_isotherm.py:290–325`

```python
# 시그모이드 감쇠: 로딩률이 85%를 넘으면 급격히 감소
ratio = loading_fraction / L_max
exponent = 12.0 * (ratio - 0.85)
```

**평가:**
- 중심점 0.85와 steepness 12는 수치적 안정성을 위한 경험적 선택이다.
- 물리적으로, 추출제 로딩이 85%에 근접하면 잔여 자유 추출제 접근성이 steric hindrance로 급감하는 것은 합리적이나, 정확한 임계점과 전이 속도는 추출제/금속/용매 조합에 따라 다르다.
- Vasilyev 2019의 공유 추출제 풀 개념에서는 로딩 감쇠가 `(HA_free/C_ext)^n_ext` 형태로 내재적으로 표현되므로, 경쟁추출 모드 ON 시 이 별도 damping과의 관계가 모호하다.

**권장:** 경쟁추출 모드와 loading_damping_factor의 역할을 명확히 분리하거나, 경쟁 모드 ON 시 damping을 비활성화하는 옵션 검토.

### 2.4 종분화 모델의 적절성과 선택적 적용

**코드 위치:** `config.py:435–460`, `extraction_isotherm.py:343–437`

**현재 상태:**
- 7개 금속 모두에 K_MOH, K_MSO4 상수가 정의되어 있으나, sulfate D correction은 D2EHPA의 Co/Ni/Li 3개 조합에만 적용 (`SULFATE_D_CORRECTION_RULES`, `config.py:456–460`).
- Cyanex 272에 대한 sulfate correction rule이 없다.

**평가:**
- 선택적 적용은 "검증된 효과만 적용" 원칙에서는 보수적이고 안전하다.
- 그러나 Ni-Cy 현장에서 Ni 1.06–1.28 M, 이온강도 I > 3 M일 때, NiSO₄⁰ 착물 비율이 상당하다. 현재 모델의 `denominator = max(1.0, 1.0 + hydroxo_term + sulfate_term)` (`extraction_isotherm.py:372`)에서 `max(1.0, ...)` 클램핑은 low-sulfate 조건에서의 불안정 방지용이나, 고농도 Ni 조건에서는 자유 금속 분율이 실제보다 높게 계산될 수 있다.
- deep-dive-notes: "황산염 효과는 기울기 왜곡으로 나타남, pH 윈도우 이동만이 아님" — 현재 모델은 이를 반영하지 못한다.

**권장:** Cyanex 272-Ni 조합에 대한 sulfate D correction rule 추가를 실험 데이터 기반으로 검토.

### 2.5 사포니피케이션 pH50 shift: 계면화학 vs bulk pH offset

**코드 위치:** `config.py:392–415`

```python
SAPONIFICATION_P_H50_SHIFT = {
    "Cyanex 272": {
        "Ni": 0.65,
        ...
    },
```

**평가:**
- Keller 2022 (deep-dive-notes): "비누화는 계면 재고 문제이지 단순 벌크 pH 오프셋이 아니다. 부분 중화 = 추출제 가용성, 양성자 방출, 공추출 개시 변경."
- 현재 구현은 pH50 shift 0.65 × sap_fraction이라는 단순 선형 모델을 사용한다. 이는 사포니피케이션의 복잡한 계면화학을 bulk pH50 이동으로 근사한 것으로, 방향성은 맞지만 메커니즘적 정확성은 부족하다.
- 특히 sap_fraction이 높을 때(>0.5), 이 선형 근사가 실제 계면 평형과 괴리될 가능성이 있다.
- Devi 1998, Swain 2006: Na-form Cyanex 272가 황산염 시스템에서 선택도를 실질적으로 변경한다는 증거가 있으나, 이를 pH50 shift 계수 하나로 포착하는 것은 과도한 단순화이다.

### 2.6 경쟁추출 공유풀 모델의 한계

**코드 위치:** `extraction_isotherm.py:483–651`

**현재 구현:**
```
D_adjusted(M) = D_sigmoid(M) × (HA_free / C_ext)^n_eff
```

**평가:**
- Vasilyev 2019의 공유 추출제 풀 개념을 시그모이드와 결합한 것은 합리적인 Level 2 근사이다.
- 그러나 Lu et al. 2024의 다성분 ESI 모델과 비교하면:
  - Lu 모델은 금속 간 직접 상호작용 항을 포함
  - 현재 모델은 순차적 추출제 차감만으로 간접 경쟁을 표현
  - 고로딩 조건에서 금속 간 교환 반응(예: Co가 Ni를 displaced)을 표현할 수 없음
- BM-CY Li MAE 32.4% (동료 코드 기준)는 Li 추출에 대한 Mn/Co/Ni 경쟁 효과가 과소평가되고 있음을 시사한다.

---

## 3. 플로우시트 전문가 검토 (sx-flowsheet-simulation-expert)

### 3.1 현재 기능 단계 (Capability Level)

**판정: Level 2–3** (decision-rubric 기준)

- Level 2: 광범위한 회로 블록 없는 역류 캐스케이드 ✓
- Level 3 부분: 단계적 pH 제어 있음 (`target_pH_per_stage`) ✓
- Level 4 미달: scrub/strip/recycle 블록 미포함 ✗

**함의:** 현재 시뮬레이터는 단일 추출 블록만 모델링한다. Strip return 조성이 신선 유기상에 영향을 미치는 효과(잔류 금속, pH 변화)가 누락되어 있다. `C_org_fresh`가 기본 0으로 설정된다 (`multistage_sx.py:211`).

### 3.2 Tear-stream 누락의 영향

**현재 코드:**
```python
if C_org_fresh is None:
    C_org_fresh = {m: 0.0 for m in metals}  # multistage_sx.py:211
```

**평가:**
- 실제 공정에서 strip 후 유기상에는 소량의 잔류 금속이 남아있다. 특히 Ni-Cy 공정에서 Cyanex 272의 Ni 선택적 스트리핑이 완전하지 않으면, 유기 재순환에 의한 Ni 축적이 발생한다.
- 이 효과를 무시하면 시뮬레이터가 추출 용량을 과대평가한다 — 이는 관찰된 Ni 과추출 편향과 일치하는 방향이다.
- decision-rubric의 "Tear-stream 필요성 체크리스트"에서 strip return 조성이 핵심 항목으로 지정되어 있다.

**권장:** 단기적으로 `C_org_fresh`에 strip 후 잔류 Ni 농도를 수동 입력하는 옵션 추가 (예: 0.1–0.5 g/L). 장기적으로 strip 블록 연결.

### 3.3 3단 vs 5단 수렴 특성

**코드 위치:** `multistage_sx.py:62–105` (relaxation alpha)

**관찰:**
- 프리셋 Data1~6은 모두 5단이며 정상적으로 수렴한다.
- Ni-Cy 현장은 3단이다. 3단 시스템에서는:
  - 역류 수렴에 필요한 반복 횟수가 적어 수렴 자체는 빠르다
  - 그러나 각 stage의 추출 부하가 5단 대비 더 크므로 로딩 감쇠의 비선형 효과가 증폭된다
  - relaxation alpha의 3단계 warmup (`iteration < 40: 0.10, < 120: 0.18, < 250: 0.28`)이 3단에서 최적인지 검증되지 않았다

**평가:** 3단 시스템 전용 수렴 테스트가 부재하다. 기존 검증은 5단에만 집중되어 있다.

### 3.4 NaOH 투입 위치/방식

**코드 위치:** `multistage_sx.py:42–59` (`_build_naoh_distribution`)

**현재 모델:**
- `uniform` (기본): 모든 stage에 균등 분배
- `front_loaded`: 앞 stage에 가중 분배 (0.5^i)
- 사포니피케이션 모드: NaOH가 유기상 전처리에 사용, 수계 직접 투입 아님

**Ni-Cy 현장 실제:**
- NaOH 반응 유형: "Saponification" (모든 날짜)
- NaOH가 유기상에 사전 반응되어 각 stage로 들어감
- 투입 위치(어느 stage 앞에서 유기상과 혼합하는지)는 CSV에 기록되지 않음

**평가:** 사포니피케이션 NaOH의 stage별 분배가 시뮬레이터 출력에 큰 영향을 미치는데, 실제 분배 방식이 불명확하다. 이는 공정 전문가 검토의 메타데이터 불완전성과 직결된다 (§4.1 참조).

### 3.5 비정상상태 데이터 가능성

**Ni-Cy 현장 데이터 분석:**

| 날짜 | Feed Ni (g/L) | Raff 1단 Ni (g/L) | O/A | NaOH (mL/min) | 비고 |
|------|-------------|-----------------|-----|--------------|------|
| 03.11 | <1 | 54.85 | 6.0 | 2.6 | Feed Ni < 1 → 비정상 |
| 03.12 | 77.58 | 34.21 | 6.8 | 5.0 | 초기 안정화 |
| 03.13 | 68.11 | 11.54 | 7.9 | (결측) | NaOH 결측 |
| ... | ... | ... | ... | ... | ... |
| 04.16 | 66.99 | 0.73 | 12.3 | 12.95 | 안정 운전 |
| 04.18 | 68.49 | 0.63 | 8.9 | 15.5 | 안정 운전 |

- 03.11은 Feed Ni가 검출 한계 미만(<1 mg/L)이면서 1단 Raff'에 54,849 mg/L가 나오는 명백한 비정상 상태 (이전 운전의 잔류).
- 03.12~03.24: O/A 비 6~9:1, NaOH 25 wt% — 초기 운전 조건.
- 04.04 이후: O/A 비 8~16:1, NaOH 10 wt%로 전환 — 운전 조건 변경.
- field_validation_main.py에서 이미 03.11을 제외하고 있으나 (`field_validation_main.py:362`), 03.12~03.24의 초기 안정화 데이터도 분리가 필요하다.

---

## 4. 공정 전문가 검토 (sx-process-expert)

### 4.1 메타데이터 완성도 평가

decision-rubric의 메타데이터 완성도 점수표에 따른 평가:

**Ni-Cy 현장 데이터:**

| 메타데이터 항목 | 상태 | 비고 |
|---------------|------|------|
| Feed 조성 | ✓ 알려짐 | 8종 금속 + pH |
| 추출제/농도 | ✓ 알려짐 | Cyanex 272, 0.473 M |
| O/A 비 | ✓ 알려짐 | 계산 가능 (Org/Aq mL/min) |
| 단수 | ✓ 알려짐 | 3단 |
| NaOH 농도 | ✓ 알려짐 | 25 wt% → 10 wt% 전환 |
| NaOH 투입 위치/분배 | ✗ 불명확 | "Saponification" 표기만 |
| pH 측정점 | ✗ 불명확 | Feed pH / Raff pH만 기록 |
| 온도 | ✗ 불명확 | CSV에 온도 정보 없음 |
| 유기상 strip 조건 | ✗ 불명확 | strip 후 유기 조성 미기록 |
| Stage efficiency | ✗ 불명확 | 평형 가정 |

**점수: 5/10 → 중간 신뢰도** — 진단과 방향성 분석은 가능하나, 정밀 보정은 메타데이터 보완 필요.

**BM-CY 현장 데이터:**

- 가동조건이 일별로 기록되어 있어 Ni-Cy보다 양호
- 5단 각 stage의 금속 농도와 pH가 기록되어 있어 stage별 검증 가능
- 그러나 NaOH 반응 방식("Saponification"), 유기상 재순환 상태 등은 동일하게 불명확
- **점수: 6/10 → 중간-상 신뢰도**

### 4.2 O/A 비 해석

**Ni-Cy 현장:**
- O/A 범위: 6.0–20.4 (Org 150–250 mL/min, Aq 11.95–31 mL/min)
- 특히 04.09: O/A = 206/11.95 = 17.2 — 극단적 비율

**물리적 의미:**
- O/A >> 1은 유기상이 수계 대비 대과잉이며, 각 stage에서 유기상의 로딩 변화가 미미하다는 뜻이다.
- 이는 추출 driving force를 유지하는 데 유리하지만, 시뮬레이터의 역류 수렴에서는 유기상 프로파일 변화가 작아 수렴 판정이 민감해질 수 있다.
- 극단적 O/A에서는 수계 유량이 매우 적어 NaOH 희석 효과가 증폭되므로, `naoh_mode = "saponification"` (수계 희석 없음)과 `naoh_mode = "aqueous_direct"` (수계 희석 포함)의 차이가 극대화된다.

### 4.3 극고농도 Ni에서의 Stage Efficiency

**문제:**
- Feed Ni 62–78 g/L (1.06–1.33 M)은 극히 높은 농도이다.
- 이 농도에서 mixer-settler의 stage efficiency가 100%에 도달하지 못할 가능성이 있다:
  - 고점도로 인한 혼합 효율 저하
  - 고이온강도로 인한 계면 장력 변화
  - 유기상 로딩 포화 근접으로 인한 평형 접근 시간 증가

**현재 모델:** Stage efficiency = 100% 가정 (평형 모델). 이는 프리셋 조건(Ni ≈ 28–34 g/L)에서는 합리적이나, 극고농도에서는 과대평가 요인이다.

**권장:** Stage efficiency factor (0.8–0.95) 도입 검토. 이는 D값에 직접 곱하는 단순한 방식으로 구현 가능하다.

### 4.4 BM-CY의 Mn 63–91 mg/L 역할

**BM-CY 데이터에서:**
- Feed Mn: 63–92 mg/L (0.063–0.092 g/L)
- Cyanex 272 선택성: Mn > Co >> Ni → Mn이 가장 먼저 추출됨

**평가:**
- Mn이 0.06–0.09 g/L로 양은 적지만, pH50 = 3.5 (Cyanex 272)로 feed pH 3.0에서 이미 상당량 추출된다.
- 경쟁추출 모델에서 Mn이 추출제를 먼저 소비하여 Co/Ni 추출에 미치는 간접 영향은 미미하다 (Mn의 추출제 소비량 = 0.09/54.938 × 2 ≈ 0.003 M << C_ext 0.631 M).
- 그러나 Co 3–4 g/L이 추출제를 소비하는 양은 상당하다: 4.0/58.933 × 2 ≈ 0.136 M (C_ext의 21.5%).

### 4.5 단위 문제: Ni-Cy 1단 후액 값

**CSV 원본 확인:**
- Ni_Cy.csv의 Raff 1단 Ni: "635", "634", "1,434", "34,207" 등
- 단위는 mg/L (ppm) — 이는 field_validation_main.py에서 `/1000.0` 변환으로 g/L로 처리된다 (`field_validation_main.py:227`).
- 최솟값 634 mg/L = 0.634 g/L (04.18), 최댓값 54,849 mg/L = 54.85 g/L (03.11)
- 단위는 **mg/L로 확인됨**. "g/L"이 아닌 mg/L이므로 혼동 주의.

---

## 5. 시뮬레이션 전문가 검토 (sx-simulation-expert)

### 5.1 프리셋 양호 + 현장 오차 = Overfitting

**증거:**

| 데이터셋 | 단수 | C_ext (M) | Feed Ni (g/L) | O/A | MAE 목표 | 실제 MAE |
|---------|------|-----------|-------------|-----|---------|---------|
| Data1 (CoSX) | 5 | 0.631 | 33.9 | 4.7 | < 1 g/L | ~0.3 g/L |
| Data2 (CoSX) | 5 | 0.631 | 33.9 | 4.7 | < 1 g/L | ~0.2 g/L |
| Data3 (CoSX) | 5 | 0.631 | 28.5 | 5.9 | < 1 g/L | ~0.15 g/L |
| Ni-Cy | 3 | 0.473 | 62–78 | 6–20 | - | ~26 g/L |
| BM-CY | 5 | 0.631 | 28–38 | 4.4–5.2 | - | ~6.5 g/L |

**분석:**
- `SITE_PARAMETER_OVERRIDES`의 5개 규칙 중 2개가 Cyanex 272-Ni/Li에 적용되며, 이들은 Data1~3에 맞춰져 있다.
- BM-CY는 C_ext와 단수가 프리셋과 동일하지만 MAE 6.5 g/L → 다금속 경쟁 효과(특히 Li)에서 차이 발생.
- Ni-Cy는 C_ext, 단수, 농도 범위가 모두 다르며 MAE 26 g/L → 파라미터 전이 완전 실패.

**진단:** 전형적인 narrow-window overfitting. working-rules: "기존 기반만 보호되는데 일반적 검증이라고 말하지 말 것."

### 5.2 Bisection pH 범위의 물리적 부적절성

**코드 위치:** `single_stage.py:885–886`

```python
pH_low = -1.0
pH_high = 14.0
```

**평가:**
- pH < 0과 pH > 14는 물리적으로 의미 없다 (수용액 기준).
- 극단적 범위가 수렴에 직접 해를 끼치지는 않으나, 수렴 속도를 저하시킨다: log₂(15/1e-4) ≈ 17 iterations vs log₂(10/1e-4) ≈ 17 → 차이는 미미하다.
- 그러나 음수 pH에서 시그모이드 함수가 수치적으로 극단적 D값을 생성할 수 있으며, 이는 물질수지 계산에서 overflow/underflow 위험을 야기한다.

**권장:** `pH_low = 0.0`, `pH_high = 12.0`으로 제한 (또는 feed pH 기반 동적 범위 설정).

### 5.3 Relaxation Factor 3단계 Warmup의 효과성

**코드 위치:** `multistage_sx.py:62–105`

```python
def _get_relaxation_alpha(iteration, last_max_diff, target_mode, relaxation_scale):
    if target_mode:
        if iteration < 40:   alpha = 0.10
        elif iteration < 120: alpha = 0.18
        elif iteration < 250: alpha = 0.28
        else:                 alpha = 0.40
```

**평가:**
- 3단계 warmup은 비선형 시스템의 수렴 안정성을 위한 합리적 전략이다.
- 그러나 `last_max_diff`에 기반한 적응적 감쇠(`> 10.0 → min 0.05`, `> 3.0 → min 0.08` 등)가 warmup과 상호작용할 때, 초기 단계에서 alpha가 극도로 낮아져(0.03) 수렴이 매우 느려질 수 있다.
- MAX_ITERATIONS = 500에서 수렴하지 않는 경우, best_residual 기반 fallback이 작동하나, 이 때의 결과는 fully-converged가 아닌 best-effort이다.

**Ni-Cy 3단에서의 관찰:** 3단 시스템에서는 stage 수가 적어 유기상 프로파일 자유도가 낮으므로, warmup 단계를 줄이는 것이 효율적일 수 있다.

### 5.4 VALIDATED_FIELD_WINDOW의 C_ext 범위 문제

**코드 위치:** `config.py:465–472`

```python
VALIDATED_FIELD_WINDOW = {
    "pH": (3.9, 7.0),
    "C_ext_m": (0.6053, 0.6308),
    "temperature_c": 25.0,
    "temperature_margin_c": 5.0,
    "n_stages": 5,
    "total_metals_warning_gL": 80.0,
}
```

**Ni-Cy 현장 실제:**
- C_ext = 0.473 M → VALIDATED_FIELD_WINDOW 범위(0.605–0.631)의 **75%** 수준
- n_stages = 3 → 검증된 5가 아닌 3
- Feed Ni = 62–78 g/L → total_metals_warning_gL = 80 g/L에 근접

**평가:**
- 이 VALIDATED_FIELD_WINDOW는 현재 경고용으로만 사용되고 있으며, 시뮬레이션 로직 자체를 제한하지는 않는다.
- 그러나 사용자가 이 범위 밖의 조건을 입력할 때 명시적 경고가 필요하다.
- 현재 대시보드(`sx_dashboard.py`)에서 이 범위 밖 경고가 표시되는지 확인이 필요하다.

### 5.5 물질수지 정합성 검증

**현재 메커니즘:**
- 각 stage에서 `total_metal_flow = C_aq_mol_in × Q_aq + C_org_mol_in × Q_org` 보존 (`single_stage.py:72`, `extraction_isotherm.py:621`)
- 역류 수렴 후 final refinement pass 실행 (`multistage_sx.py:580–618`)
- `CONVERGENCE_TRACE_ORGANIC_ABS_TOL_G_L = 2e-5` 이하의 미세 변동은 무시

**평가:**
- 단일 stage 내 물질수지는 수학적으로 정확하다 (D값 기반 분배).
- 역류 반복에서 relaxation으로 인한 일시적 물질수지 불균형이 발생하나, final exact pass에서 해소된다.
- 그러나 사포니피케이션 모드에서 OH⁻ 수지와 Na⁺ 수지의 정합성은 명시적으로 검증되지 않는다. `saponification_inventory_error_mol_hr` 필드가 존재하나 (`single_stage.py:1011–1015`), 이 값이 큰 경우의 경고 메커니즘이 없다.

### 5.6 NaOH 사포니피케이션 모드 역산 정확도

**현재 경로:**
1. `physical_v2`: stage별 유기상 sap inventory를 역류 반복으로 전파 (`multistage_sx.py:306–324`)
2. `legacy_equivalent_target`: NaOH 조건을 등가 target pH로 환산하는 회귀식 (`single_stage.py:394–438`)

**평가:**
- `legacy_equivalent_target`의 회귀 계수 (`config.py:375–388`)는 Data1~6에 맞춰진 것으로, Ni-Cy 조건에서의 정확도가 검증되지 않았다.
- `physical_v2`의 stage별 sap 전파는 더 물리적이나, `SAPONIFICATION_DIRECT_NEUTRALIZATION_FACTOR`가 Cyanex 272에서 0.25로 설정되어 있어 대부분의 sap capacity가 금속 교환에 사용된다. 이 비율이 극고농도 조건에서도 유지되는지 검증이 필요하다.

---

## 6. 종합 개선 로드맵

### 6.1 단기 (1–2주)

| # | 개선안 | 대상 파일 | 기대 효과 | 위험 수준 |
|---|-------|---------|----------|----------|
| S1 | Ni-Cy 운전 데이터 전처리: 비정상상태 일자 식별 및 제외 | `field_validation_main.py` | 보정 데이터 품질 확보 | structure_safe |
| S2 | C_ext=0.473 M 전용 SITE_PARAMETER_OVERRIDES 프로필 추가 | `config.py` | Ni MAE 대폭 감소 | recalibration_risk |
| S3 | VALIDATED_FIELD_WINDOW에 Ni-Cy 범위 추가 및 경고 메시지 | `config.py`, `sx_dashboard.py` | 사용자 인식 개선 | structure_safe |
| S4 | Bisection pH 범위를 [0.0, 12.0]으로 제한 | `single_stage.py:885–886` | 수치 안정성 향상 | structure_safe |
| S5 | Stage efficiency factor 입력 옵션 추가 (기본값 1.0) | `single_stage.py` | 극고농도 조건 대응 | recalibration_risk |

### 6.2 중기 (1–2개월)

| # | 개선안 | 대상 | 기대 효과 | 위험 수준 |
|---|-------|------|----------|----------|
| M1 | 이온강도 기반 activity coefficient 보정 (Davies 식 또는 SIT 모델) | `extraction_isotherm.py` | 극고농도 자유금속 분율 정확도 | validation_breaking |
| M2 | Cyanex 272-Ni/Li sulfate D correction rule 추가 | `config.py`, `extraction_isotherm.py` | Ni-Cy 정확도 향상 | recalibration_risk |
| M3 | 사포니피케이션 pH50 shift의 비선형 모델 (포화 곡선) | `config.py`, `extraction_isotherm.py` | sap 효과 정확도 향상 | validation_breaking |
| M4 | 3단 역류 전용 수렴 벤치마크 구축 | 테스트 코드 | solver 안정성 검증 | structure_safe |
| M5 | C_org_fresh에 strip 후 잔류 농도 입력 옵션 | `multistage_sx.py`, `sx_dashboard.py` | 유기 재순환 효과 반영 | recalibration_risk |

### 6.3 장기 (3–6개월)

| # | 개선안 | 대상 | 기대 효과 | 위험 수준 |
|---|-------|------|----------|----------|
| L1 | Strip 블록 추가 → Level 4 기능 달성 | 신규 모듈 | 완전한 회로 시뮬레이션 | flowsheet-required |
| L2 | Lu et al. 2024 ESI 모델 기반 경쟁추출 고도화 | `extraction_isotherm.py` | 다금속 정확도 대폭 향상 | validation_breaking |
| L3 | 동적 상태 추적 (transient 모드) | 신규 모듈 | startup/shutdown 시뮬레이션 | dynamic-required |
| L4 | McCabe-Thiele 다이어그램 자동 생성 | `sx_dashboard.py` | 공정 해석 도구 | structure_safe |

---

## 부록 A: 핵심 코드-이슈 매핑

| 이슈 | 코드 위치 | 섹션 |
|------|----------|------|
| Ni pH50 override | `config.py:252–255` | §2.1 |
| n_eff_multiplier | `extraction_isotherm.py:549–572` | §2.2 |
| loading_damping_factor | `extraction_isotherm.py:290–325` | §2.3 |
| 종분화 1:1 한계 | `extraction_isotherm.py:343–383` | §2.4 |
| Sap pH50 shift | `config.py:392–415` | §2.5 |
| 경쟁추출 공유풀 | `extraction_isotherm.py:483–651` | §2.6 |
| C_org_fresh 기본값 | `multistage_sx.py:211` | §3.2 |
| Relaxation warmup | `multistage_sx.py:62–105` | §3.3, §5.3 |
| NaOH 분배 | `multistage_sx.py:42–59` | §3.4 |
| Bisection pH 범위 | `single_stage.py:885–886` | §5.2 |
| VALIDATED_FIELD_WINDOW | `config.py:465–472` | §5.4 |
| Sap inventory error | `single_stage.py:1011–1015` | §5.5 |
| 등가 target pH 회귀 | `single_stage.py:394–438`, `config.py:375–388` | §5.6 |

## 부록 B: 데이터 요약

### Ni-Cy (Ni_Cy.csv)
- 22개 날짜 (25.03.11 – 25.05.13)
- 추출제: Cyanex 272, C_ext = 0.473 M
- 3단 Mixer-Settler
- Feed Ni: 62–78 g/L (03.11 제외)
- O/A: 6.0–20.4
- NaOH: Saponification, 25 wt% → 10 wt% 전환

### BM-CY (BM-CY.csv)
- 11일차 (가동일수 1–11)
- 추출제: Cyanex 272, C_ext = 0.631 M
- 5단 Mixer-Settler
- Feed Ni: 28–38 g/L, Co: 3.0–4.1 g/L, Li: 5.3–6.3 g/L, Mn: 0.063–0.092 g/L
- O/A: 4.4–5.2
- NaOH: Saponification, 10 wt%, 8.2–12.0 mL/min

## 부록 C: Guardrails 준수 확인

| 규칙 | 준수 여부 | 비고 |
|------|---------|------|
| 현재 모델을 '완전한 MSE'로 부르지 않을 것 | ✓ | "반경험적 시그모이드 모델"로 일관 표기 |
| 단일 캘리브레이션에서 '물리적 검증' 주장 금지 | ✓ | §2.1에서 명시적으로 경고 |
| 황산염→염화물 직접 전이 금지 | ✓ | 황산염 시스템 내로 한정 |
| Cyanex 272↔D2EHPA 파라미터 직접 전이 금지 | ✓ | 추출제별 독립 논의 |
| 배치 평형→역류 회로 성능 주장 불가 | ✓ | 회로 수준 주장 시 Level 명시 |
| 정상상태 결과를 제어 준비로 표현 금지 | ✓ | §3.1에서 Level 2–3 명시 |
| 보고서만 수정으로 일반 검증 주장 금지 | ✓ | §5.1에서 narrow-window overfitting 진단 |
| overclaim 금지 | ✓ | 모든 권장사항에 위험 라벨 부착 |

---

*본 보고서는 sx-equilibrium-expert, sx-flowsheet-simulation-expert, sx-process-expert, sx-simulation-expert 4개 스킬의 working-rules, decision-rubric, deep-dive-notes, paper-summaries를 참조하여 작성되었습니다.*
