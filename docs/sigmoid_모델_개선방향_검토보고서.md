# SX Simulator Sigmoid 모델 개선방향 검토 보고서

**문서 버전:** v1.0
**작성일:** 2026-03-18
**대상 시스템:** sx-simulator (다금속 역류 용매추출 시뮬레이터)

---

## 목차

1. [검토 배경 및 목적](#1-검토-배경-및-목적)
2. [현행 모델 분석](#2-현행-모델-분석)
3. [문제 식별: E_max < 100% 금속의 모델링 한계](#3-문제-식별-e_max--100-금속의-모델링-한계)
4. [대안 1: Richards 함수 (5PL) 도입](#4-대안-1-richards-함수-5pl-도입)
5. [대안 2: McCabe-Thiele 등온선 기반 모델](#5-대안-2-mccabe-thiele-등온선-기반-모델)
6. [대안 비교 및 종합 평가](#6-대안-비교-및-종합-평가)
7. [권고 방향](#7-권고-방향)
8. [참고문헌](#8-참고문헌)

---

## 1. 검토 배경 및 목적

### 1.1 배경

현행 SX Simulator는 sigmoid(시그모이드) 함수 기반의 반경험적 모델로 pH-추출률 관계를 근사하고 있다. 이 모델은 Vasilyev 경쟁 추출 메커니즘과 결합하여 다금속 역류 mixer-settler 공정을 시뮬레이션한다.

그러나 **pH를 증가시켜도 추출률이 100%에 도달하지 않는 금속**(Li, Ca, Mg)의 경우, 현행 sigmoid 모델의 구조적 한계로 인해 예측 정확도가 저하될 수 있다는 문제가 제기되었다.

### 1.2 목적

본 보고서는 다음을 목적으로 한다:

1. 현행 sigmoid 모델에서 E_max < 100% 금속에 대한 구체적 오류 가능성을 분석
2. 대안적 모델링 접근법(Richards/5PL, McCabe-Thiele)의 이론적 타당성과 적용 가능성을 평가
3. 시뮬레이터 개선을 위한 권고 방향을 제시

---

## 2. 현행 모델 분석

### 2.1 핵심 수식

현행 시뮬레이터의 추출 효율 계산은 4-Parameter Logistic(4PL) 시그모이드 함수에 기반한다:

```
E(pH, T) = E_max / (1 + exp(-k(T) × (pH - pH50(T))))
```

각 파라미터의 정의:

| 파라미터 | 의미 | 보정 방법 |
|----------|------|-----------|
| `pH50` | 반추출 pH (E = E_max/2 지점) | 추출제 농도(α) 및 온도(β) 보정 |
| `k` | 시그모이드 기울기 (전이 급격도) | 온도(γ) 보정 |
| `E_max` | 최대 추출률 (%) | 고정값 (금속/추출제별) |

온도 및 농도 보정 수식:

```
pH50_eff = pH50_ref - α × log₁₀(C_ext / C_ref) + β × (T - T_ref)
k_eff    = k_ref × exp(γ × (T - T_ref))
```

### 2.2 분배계수 계산

시그모이드에서 산출된 추출률(E)을 물질수지에 사용하기 위해 분배계수(D)로 변환한다:

```
D = E / (100 - E)    (E는 %, D는 무차원)
```

수치 안정성을 위한 클리핑:
- E >= 99.99% → D = 1,000,000
- E <= 0.01% → D = 0.000001

### 2.3 경쟁 추출 모델 (Vasilyev)

다금속 시스템에서 추출제 풀을 공유하는 경쟁 효과를 반영한다:

```
D_adjusted(M) = D_sigmoid(M) × (HA_free / C_ext)^n_eff(M)
```

- pH50가 낮은 금속(Zn, Mn)이 우선 추출되어 추출제를 소비
- 후순위 금속(Ni, Li)은 잔여 추출제 비율에 의해 D가 감소

추출 우선순위:
- **Cyanex 272:** Zn → Mn → Co → Ni → Ca → Mg → Li
- **D2EHPA:** Zn → Mn → Ca → Co → Mg → Ni → Li

### 2.4 금속별 파라미터 현황 (Field-Calibrated 기준)

| 금속 | 추출제 | pH50 | k | E_max (%) | n_H | n_ext |
|------|--------|------|-----|-----------|-----|-------|
| Mn | Cyanex 272 | 3.5 | 3.0 | 99.5 | 2 | 2 |
| Co | Cyanex 272 | 4.0 | 3.5 | 99.5 | 2 | 2 |
| Ni | Cyanex 272 | 6.55 | 2.0 | 99.0 | 2 | 2 |
| Li | Cyanex 272 | 6.5 | 2.0 | **95.0** | 1 | 1 |
| Zn | Cyanex 272 | 1.85 | 3.5 | 99.5 | 2 | 2 |
| Ca | Cyanex 272 | 5.5 | 2.0 | **90.0** | 2 | 2 |
| Mg | Cyanex 272 | 5.7 | 2.0 | **85.0** | 2 | 2 |

| 금속 | 추출제 | pH50 | k | E_max (%) | n_H | n_ext |
|------|--------|------|-----|-----------|-----|-------|
| Mn | D2EHPA | 2.5 | 3.5 | 99.5 | 2 | 2 |
| Co | D2EHPA | 3.5 | 3.0 | 99.5 | 2 | 2 |
| Ni | D2EHPA | 4.5 | 2.5 | 99.0 | 2 | 2 |
| Li | D2EHPA | 4.3 | 2.5 | **95.0** | 1 | 1 |
| Zn | D2EHPA | 1.5 | 4.0 | 99.5 | 2 | 2 |
| Ca | D2EHPA | 3.5 | 2.0 | **90.0** | 2 | 2 |
| Mg | D2EHPA | 3.2 | 2.0 | **85.0** | 2 | 2 |

**E_max < 100% 금속: Li(95%), Ca(90%), Mg(85%)**

---

## 3. 문제 식별: E_max < 100% 금속의 모델링 한계

### 3.1 곡선 대칭성 문제 (가장 핵심)

표준 시그모이드(4PL)는 **pH50을 중심으로 대칭인 S자 곡선**을 생성한다. 그러나 E_max < 100%인 금속의 실제 추출 곡선은 비대칭인 경우가 많다:

- **상승 구간** (0% → E_max/2): 비교적 급격한 전이
- **포화 구간** (E_max/2 → E_max): 열역학적 한계에 의한 완만한 수렴

예를 들어 Mg(E_max = 85%)의 경우:

```
실제 곡선:     ___________──────── 85% (완만하게 수렴)
              /
             /
            /
───────────/
           pH →

sigmoid 곡선:        ┌────────── 85% (급격하게 수렴)
                    /
                   /
                  /
─────────────────/
                 pH →
```

단일 기울기 파라미터(k)로는 상승부와 포화부의 **서로 다른 곡률**을 독립적으로 제어할 수 없다.

### 3.2 분배계수(D) 상한 제약

E_max가 D의 상한을 구조적으로 결정한다:

| 금속 | E_max (%) | D_max = E_max/(100-E_max) |
|------|-----------|---------------------------|
| Co | 99.5 | 199.0 |
| Ni | 99.0 | 99.0 |
| Li | 95.0 | 19.0 |
| Ca | 90.0 | 9.0 |
| Mg | 85.0 | **5.67** |

Vasilyev 경쟁 모델에서 `D_adjusted = D_sigmoid × competition_factor`이므로, E_max가 낮은 금속은 경쟁 상황에서 **구조적으로 불리**하다. 이것이 물리적으로 정확한 측면도 있으나, 곡선 형태 오차와 결합되면 실제보다 추출률을 **과소 예측**할 위험이 있다.

### 3.3 pH50 해석의 모호성

E_max < 100%일 때, 모델의 pH50은 **E = E_max/2** 지점을 의미한다:

| 금속 | E_max | pH50에서의 실제 E값 | 문헌상 pH50 정의(E=50%) |
|------|-------|---------------------|------------------------|
| Co | 99.5% | 49.75% | ≈ 50% (거의 일치) |
| Li | 95.0% | **47.5%** | 50% (2.5%p 차이) |
| Mg | 85.0% | **42.5%** | 50% (7.5%p 차이) |

문헌에서 정의하는 pH50(추출률 50% 지점)과 모델의 pH50이 서로 다른 값을 가리키므로, 문헌 파라미터를 직접 사용할 때 체계적 오차가 발생할 수 있다.

### 3.4 피팅 안정성

현행 `fitting.py`의 bounds 설정:

```python
bounds = ([0, 0.1, 10], [14, 20, 105])
#          pH50_min, k_min, E_max_min → pH50_max, k_max, E_max_max
```

k와 E_max 사이에 상관관계가 존재하여, 실험 데이터가 포화 구간을 충분히 포함하지 않으면 **비물리적 파라미터 조합으로 수렴**할 수 있다 (높은 k + 낮은 E_max ↔ 낮은 k + 높은 E_max).

### 3.5 열역학적 배경

Li, Ca, Mg가 100% 추출에 도달하지 못하는 것은 공정 조건의 문제가 아닌 **열역학적 한계**이다:

- **Li-Cyanex 272:** Li⁺와 인산계 추출제(Cyanex 272) 간의 결합 친화력(affinity)이 본질적으로 약하여, 추출 평형상수(Kex)가 작다. pH를 무한히 올려도 열역학적 한계를 초과할 수 없음 (Nguyen & Lee, 2021).
- **Ca, Mg:** 알칼리토류 금속은 전이금속 대비 유기상 복합체 안정도가 낮아 추출제 선택성이 제한됨.

---

## 4. 대안 1: Richards 함수 (5PL) 도입

### 4.1 수학적 정의

Richards 함수(1959)는 Generalized Logistic Function, 또는 5-Parameter Logistic(5PL)이라고도 불리며, 기존 4PL에 **비대칭 파라미터(ν)**를 추가한 형태이다:

```
4PL (현행):   E(pH) = E_max / (1 + exp(-k × (pH - pH50)))
5PL (제안):   E(pH) = E_max / (1 + exp(-k × (pH - pH50)))^(1/ν)
```

| 파라미터 | 의미 |
|----------|------|
| ν = 1 | 기존 4PL과 동일 (대칭 시그모이드) |
| ν > 1 | 포화 구간이 더 완만 (E_max에 천천히 접근) |
| ν < 1 | 상승 구간이 더 완만 |

### 4.2 이론적 근거

#### 4.2.1 타 분야 검증 현황

5PL은 **생의학/약학 분야에서 가장 널리 검증**된 비대칭 곡선 모델이다:

- **FDA 권장 모델:** 면역분석(immunoassay) 및 생물분석(bioassay)의 calibration curve fitting에 FDA가 5PL 사용을 권고
- **대규모 실증 연구:** Gottschalk & Dunn(2005)이 **60,000개 이상의 면역분석 곡선**을 분석한 결과:
  - 14%만이 대칭 (ν = 0.9~1.1)
  - **86%가 비대칭** (ν 값이 다양하게 분포)
- **정확도 개선:** 비대칭 데이터에 5PL 적용 시 4PL 대비 농도 추정 정확도가 크게 향상되며, 특히 곡선 양 극단(저농도/고농도)에서 개선이 현저 (PMC, 2014)

#### 4.2.2 SX 분야 적용 선례

**용매추출 pH-추출률 곡선에 Richards/5PL을 직접 적용한 문헌은 확인되지 않았다.**

그러나 pH-추출률 곡선은 본질적으로 dose-response 곡선과 동일한 S자 형태를 가지며, 동일한 수학적 고려사항이 적용된다. 따라서 타 분야의 검증 결과를 유추 적용하는 것은 수학적으로 합리적이다.

#### 4.2.3 물리적 타당성

E_max < 100% 금속에서 포화 구간이 완만해지는 현상은 물리적으로 설명 가능하다:

- 추출 평형상수(Kex)가 작은 금속은 **평형에 도달하는 구동력(driving force)이 약함**
- 고추출률 영역에서 수상 금속 농도가 감소하면 **정반응 속도가 떨어지고 역반응이 상대적으로 강화**
- 이로 인해 E_max에 접근하는 곡선이 대칭 sigmoid보다 **완만한 형태**를 보일 수 있음

### 4.3 장점

| 항목 | 설명 |
|------|------|
| 하위 호환성 | ν = 1이면 기존 4PL과 동일. Co, Mn 등 E_max ≈ 100% 금속은 영향 없음 |
| 최소 코드 변경 | 기존 `extraction_efficiency()` 함수에 지수 `^(1/ν)` 한 줄 추가로 구현 가능 |
| 표현력 향상 | 단일 파라미터(ν) 추가로 상승부/포화부 곡률을 독립 제어 |
| 통계적 근거 | 60,000+ 곡선 분석에서 86%가 비대칭이라는 실증 결과 |

### 4.4 한계 및 위험

| 항목 | 설명 |
|------|------|
| SX 분야 선례 부족 | 용매추출에서 5PL을 적용한 문헌이 없어 동료 심사(peer review) 근거 미비 |
| 과적합 위험 | 파라미터 증가(ν)로 제한된 실험 데이터에서 과적합 가능성 |
| 피팅 난이도 | 5PL 피팅은 4PL보다 수렴이 어려우며 초기값 의존성이 높음 |
| 물리적 해석 한계 | ν 값 자체에 명확한 물리화학적 의미를 부여하기 어려움 |

---

## 5. 대안 2: McCabe-Thiele 등온선 기반 모델

### 5.1 McCabe-Thiele 방법 개요

McCabe-Thiele 도표법은 용매추출 플랜트 설계에서 수십 년간 사용된 **업계 표준 방법**이다.

핵심 구성 요소:

```
1. 평형 등온선 (Equilibrium Isotherm)
   - 실험적으로 측정된 수상/유기상 금속 농도 평형 관계
   - x축: 수상 금속 농도 (C_aq)
   - y축: 유기상 금속 농도 (C_org)

2. 작동선 (Operating Line)
   - 물질수지에서 도출, 기울기 = Q_aq/Q_org (유량비의 역수)
   - 입구/출구 조건을 연결

3. 단수 결정 (Stage Stepping)
   - 등온선과 작동선 사이를 계단식으로 이동하며 필요 단수 결정
```

### 5.2 SX 분야 적용 사례

McCabe-Thiele 방법은 다수의 SX 연구에서 사용되었다:

- **Co/Cyanex 272 추출:** McCabe-Thiele 도표로 2단계 추출 시 99.7% Co 추출 달성 확인 (ResearchGate, 2023)
- **Mixed-metals isotherm:** Ritcey(2006)는 2개 금속이 그룹으로 추출되는 시스템에서 복합 금속 등온선과 log-log McCabe-Thiele 도표를 적용
- **Cu SX 최적화:** 구리 용매추출 회로 데이터 평가에 McCabe-Thiele 도표 활용 (Du Preez, 2015)

### 5.3 현행 시뮬레이터와의 관계

현행 시뮬레이터의 다단 역류 솔버(`multistage_sx.py`)는 **McCabe-Thiele와 본질적으로 동일한 계산을 수치적으로 수행**하고 있다:

| McCabe-Thiele 구성 요소 | 현행 시뮬레이터 대응 |
|-------------------------|---------------------|
| 평형 등온선 (isotherm) | sigmoid 함수로 계산된 D값 |
| 작동선 (operating line) | 물질수지 (Q_aq, Q_org, 유량비) |
| 단수 결정 (stage stepping) | fixed-point 반복 수렴 (최대 500회) |

**핵심 차이점은 오직 하나:** 평형 등온선을 sigmoid 함수로 **계산**하느냐, 실험 데이터에서 **직접 읽어오느냐**의 차이이다.

### 5.4 장점

| 항목 | 설명 |
|------|------|
| 모델 오차 제거 | 실험 데이터를 직접 사용하므로 sigmoid 형태 가정 불필요 |
| E_max 문제 해소 | 데이터 자체가 실제 최대 추출률을 내재적으로 반영 |
| 비대칭 곡선 처리 | 어떤 형태의 곡선이든 그대로 반영 |
| 업계 검증 | SX 플랜트 설계에서 수십 년간 사용된 표준 방법 |

### 5.5 치명적 한계

#### 5.5.1 pH 의존성 문제

현 시스템은 **각 stage마다 pH가 변화**한다 (H⁺ 방출 + NaOH 첨가). McCabe-Thiele는 고정된 하나의 등온선을 전제로 하나, 실제로는 **pH마다 별도의 등온선이 필요**하다:

```
Stage 1 (pH 5.0) → Stage 2 (pH 5.5) → Stage 3 (pH 6.0) → ...
  등온선 A 적용        등온선 B 적용        등온선 C 적용
  Co 추출 시작         Co 거의 완료         Ni 추출 시작
```

#### 5.5.2 다금속 경쟁 반영 불가

단일 금속 등온선으로는 경쟁 효과를 반영할 수 없다. Ritcey(2006)의 mixed-metals isotherm도 **2개 금속이 그룹으로 추출되는 경우**에만 적용 가능하며, 7개 금속이 서로 다른 pH에서 경쟁하는 현 시스템에는 직접 적용이 어렵다.

#### 5.5.3 실험 데이터 요구량

필요한 등온선 데이터 조합:

```
7개 금속 × N개 pH 포인트 × M개 온도 × L개 추출제 농도 = 수백~수천 개 데이터
```

이는 현실적으로 확보가 극히 어렵다.

#### 5.5.4 보간 및 외삽 한계

- 데이터 포인트 사이의 **보간 오차** 존재
- 실험 범위 밖의 조건에 대한 **외삽이 원천적으로 불가능**
- 새로운 금속/추출제 조합 추가 시 추가 실험 필수

#### 5.5.5 업계 인식

McCabe-Thiele의 정확한 최적화를 위해서는 조건 변경 시마다 새로운 distribution curve를 측정해야 하며, 다금속 동시추출, stage별 pH 차이, 화학종 변화, 부피 비율 변화가 있는 경우 **예측 정확도가 떨어진다**는 것이 알려져 있다. 이러한 이유로 최근 연구에서는 McCabe-Thiele를 단독으로 사용하기보다 시뮬레이션 프로그램(RSP 등)과 **병행**하여 사용하는 추세이다.

---

## 6. 대안 비교 및 종합 평가

### 6.1 평가 기준별 비교

| 평가 기준 | 현행 4PL Sigmoid | 5PL (Richards) | McCabe-Thiele |
|-----------|-----------------|----------------|---------------|
| E_max<100% 곡선 정확도 | **낮음** (대칭 가정) | **높음** (비대칭 제어) | **높음** (데이터 직접 사용) |
| 다금속 경쟁 반영 | **우수** (Vasilyev) | **우수** (Vasilyev 유지) | **매우 제한적** |
| pH 가변 stage 대응 | **우수** (함수 기반) | **우수** (함수 기반) | **곤란** (등온선 재측정 필요) |
| 실험 데이터 요구량 | 낮음 (3개 파라미터/금속) | 낮음 (4개 파라미터/금속) | **매우 높음** (수백 데이터) |
| 외삽 능력 | 우수 | 우수 | **불가** |
| 구현 복잡도 | 기존 | 최소 변경 | **대규모 재설계** 필요 |
| 이론적 근거 (SX) | 관행적 사용 | 타 분야 검증 | **SX 업계 표준** |
| 물리적 해석 가능성 | 중간 | 낮음 (ν의 의미) | 높음 (실측 데이터) |
| 하위 호환성 | - | **완전 호환** (ν=1) | 비호환 |

### 6.2 시스템 적합도 종합 평가

| 접근법 | 적합도 | 핵심 근거 |
|--------|--------|-----------|
| 순수 McCabe-Thiele | **낮음** | 다금속 경쟁 + pH 변화 시스템에 구조적 한계 |
| 현행 4PL Sigmoid | **중간** | E_max < 100% 금속의 곡선 형태 오차 가능 |
| 5PL (Richards) 도입 | **중상** | 최소 변경으로 곡선 표현력 향상, 단 SX 선례 미비 |
| 실험 등온선 → sigmoid 피팅 보정 | **높음** | 기존 구조 유지 + 데이터 기반 정확도 확보 |
| 하이브리드 (등온선 + sigmoid) | **높음** | 가장 유연하나 구현 복잡도 증가 |

---

## 7. 권고 방향

### 7.1 단기 (즉시 적용 가능)

#### 7.1.1 현행 모델의 데이터 기반 검증

기존 4PL sigmoid 모델이 실제 데이터와 어떤 패턴의 잔차를 보이는지 검증한다:

- 현행 6개 field data(Data1~Data6)에 대한 잔차 분석 수행
- E_max < 100% 금속(Li, Ca, Mg)에 대해 **체계적 편향(systematic bias)** 존재 여부 확인
- 잔차가 랜덤 분포이면 현행 모델 유지, 체계적 편향이 있으면 개선 진행

#### 7.1.2 pH50 정의 명확화

문헌 호환성을 위해, pH50의 정의를 문서화하고 필요 시 보정 계수를 도입한다:

```
pH50_literature = pH at E = 50%
pH50_model      = pH at E = E_max/2

보정: pH50_model = pH50_literature + (1/k) × ln(E_max/50 - 1) - (1/k) × ln(1)
```

### 7.2 중기 (검증 결과에 따라)

#### 7.2.1 5PL 도입 (잔차 편향이 확인된 경우)

E_max < 100% 금속에 한해 ν 파라미터를 도입한다:

```python
# 변경 전 (extraction_isotherm.py:236)
E = E_max / (1.0 + math.exp(exponent))

# 변경 후
nu = params.get("nu", 1.0)  # 기본값 1.0 = 기존 동작 유지
E = E_max / (1.0 + math.exp(exponent)) ** (1.0 / nu)
```

적용 대상:
- Li: ν 추정 필요 (실험 데이터 기반 피팅)
- Ca: ν 추정 필요
- Mg: ν 추정 필요
- 나머지 금속: ν = 1.0 (변경 없음)

#### 7.2.2 실험 등온선 데이터 활용 피팅

McCabe-Thiele용 실험 등온선 데이터가 확보되면, 이를 sigmoid(4PL 또는 5PL) 파라미터 피팅에 활용하여 데이터 기반 모델 보정을 수행한다:

```
실험 등온선 데이터 → fit_sigmoid() → 보정된 (pH50, k, E_max, ν) → 시뮬레이션
```

### 7.3 장기 (연구 과제)

#### 7.3.1 열역학적 평형 모델 검토

sigmoid 근사를 넘어서 질량작용법칙(mass-action law) 기반의 열역학적 모델 도입을 검토한다:

```
Kex = [ML_n]_org × [H⁺]^n / ([M^n⁺]_aq × [(HA)₂]^m_org)
```

- 장점: 물리적 의미 명확, 온도/농도 변화에 대한 근본적 예측 가능
- 단점: 종 분화(speciation) 정보 필요, 구현 복잡도 높음
- 참고: Nguyen & Lee(2021)가 Li/Cyanex 272 시스템에 대해 유사한 접근을 시도한 바 있음

---

## 8. 참고문헌

1. **Richards, F.J.** (1959). "A Flexible Growth Function for Empirical Use." *Journal of Experimental Botany*, 10(2), 290-301.
   - Richards 함수(Generalized Logistic) 원저

2. **Gottschalk, P.G. & Dunn, J.R.** (2005). "The five-parameter logistic: A characterization and comparison with the four-parameter logistic." *Analytical Biochemistry*, 343(1), 54-65.
   - 60,000+ 곡선 분석, 86% 비대칭 실증 ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0003269705003313))

3. **Cumberland, W.G. et al.** (2014). "Nonlinear Calibration Model Choice between the Four and Five-Parameter Logistic Models." *PMC*.
   - 4PL vs 5PL 비교, 비대칭 데이터에서 5PL 우위 ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4263697/))

4. **Ritcey, G.M.** (2006). "Use of mixed-metals isotherm and log-log McCabe Thiele's diagram in solvent extraction — A case study." *Hydrometallurgy*, 81(1), 1-12.
   - 복합 금속 등온선 및 McCabe-Thiele 적용 ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0304386X05002239))

5. **Nguyen, T.H. & Lee, M.S.** (2021). "Development of heterogeneous equilibrium model for lithium solvent extraction using organophosphinic acid." *Separation and Purification Technology*, 274, 118970.
   - Li/Cyanex 272 열역학적 평형 모델 ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S1383586621010170))

6. **Sole, K.C.** (2023). "Design of Cobalt Solvent Extraction using Cyanex 272 as Extractant by Combining McCabe Thiele Procedure and Simulation Program."
   - McCabe-Thiele + RSP 병행 설계 ([ResearchGate](https://www.researchgate.net/publication/374695961))

7. **GraphPad** (2024). "Asymmetrical (five parameter) logistic dose-response curves."
   - 5PL 실용 가이드 ([GraphPad](https://www.graphpad.com/guides/prism/latest/curve-fitting/reg_asymmetric_dose_response_ec.htm))

8. **911 Metallurgist.** "Solvent Extraction Plants: Thiele Diagram & Theoretical Design Aspects."
   - SX McCabe-Thiele 설계 개론 ([911Metallurgist](https://www.911metallurgist.com/blog/solvent-extraction-plant-design/))

---

*본 보고서는 현행 시뮬레이터의 코드 분석, 문헌 조사 및 이론적 검토를 기반으로 작성되었습니다.*
