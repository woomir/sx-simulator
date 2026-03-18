# SX 시뮬레이터 종합 전문가 검토 보고서 (하이브리드 병렬 생성)

**작성일:** 2026-03-16 (하이브리드 재생성: 2026-03-17)
**대상 버전:** main 브랜치 (commit c86c30e)
**검토 범위:** 평형 / 플로우시트 / 공정 / 시뮬레이션 4개 전문 영역 독립 병렬 분석 후 통합
**생성 방식:** 4개 전문가 Agent 컨텍스트 격리 병렬 실행 → 메인 통합

---

## 1. 요약 (Executive Summary)

### 생성 방식 설명

본 보고서는 4개 SX 전문가(평형, 플로우시트, 공정, 시뮬레이션)를 **독립된 컨텍스트에서 병렬 실행**하여 생성하였다. 각 전문가는 자신의 스킬 도메인 지식(SKILL.md + references/)을 직접 읽고, 다른 전문가의 분석 결과 없이 독립적으로 코드를 검토하였다. 이를 통해:
- 교차 오염 없는 독립적 발견사항 도출
- 다수 전문가가 독립적으로 지적한 항목의 교차 검증
- 후반 전문가의 컨텍스트 압박 해소

### 핵심 발견사항 Top-5

다수 전문가가 독립적으로 지적한 항목을 우선 배치하였다.

1. **비정상상태 데이터 혼재 및 필터링 부재** [플로우시트 §3.5 `validation_breaking` + 공정 §4.5]
   - 플로우시트 전문가: 03.11 단일 날짜만 하드코딩 제외, 자동 판별 로직 부재
   - 공정 전문가: NaOH 결측 패턴이 25 wt% 구간에서 50%(4/8건)로 집중
   - **독립 교차 확인**: 두 전문가가 서로 다른 관점(데이터 필터링 vs 결측 패턴)에서 동일 문제를 식별

2. **극고농도 Ni의 이온강도/activity 효과 미반영** [평형 §2.4 `validation_breaking` + 시뮬레이션 §5.1]
   - 평형 전문가: Ni 1.06-1.28 M에서 I > 3 M, 종분화 1:1 착물만 고려, activity coefficient 부재
   - 시뮬레이션 전문가: 프리셋 vs 현장 MAE 격차의 구조적 원인 중 하나로 독립 식별
   - **독립 교차 확인**: 평형(화학적 근거) + 시뮬레이션(수치적 결과)에서 같은 결론

3. **사포니피케이션 모델의 다층적 한계** [평형 §2.5 + 플로우시트 §3.4 + 시뮬레이션 §5.5 `validation_breaking`, §5.6]
   - 평형 전문가: pH50 shift 선형 근사 한계, proton release와의 이중 계산 위험
   - 플로우시트 전문가: stage별 분배 불명확, sap inventory 역류 전파 수렴 미검증
   - 시뮬레이션 전문가: OH⁻/Na⁺ 수지 검증 부재, NaOH 전달 경로에서 sap truncation 미기록
   - **독립 교차 확인**: 3개 전문가가 사포니피케이션의 서로 다른 층위(화학/구조/수치)에서 문제 식별

4. **Stage efficiency 100% 가정의 극고농도 비검증** [공정 §4.4 `validation_breaking`]
   - 공정 전문가: Ni 62-78 g/L에서 수용액 점도 증가, 유기상 로딩 포화, phase disengagement 지연
   - loading_damping_factor와 stage efficiency는 물리적으로 다른 메커니즘이나 현재 코드에서 구분 없음

5. **3단 시스템 solver 수렴 특성 미최적화** [플로우시트 §3.3 + 시뮬레이션 §5.3, §5.7]
   - 플로우시트 전문가: 5단 기준 warmup 구간(40/120/250)이 3단에서 과보수적
   - 시뮬레이션 전문가: alpha 0.03 하한에서 500회 내 수렴 불가 위험, 중간 상태 해 수렴 가능성
   - **독립 교차 확인**: 두 전문가가 동일한 `_get_relaxation_alpha` 코드를 독립 분석하여 같은 결론

### 최우선 권장사항

| 우선순위 | 권장사항 | 기대 효과 | 교차 확인 |
|---------|---------|----------|----------|
| P0 | 비정상상태 운전일 식별 기준 구축 (NaOH 농도 전환점 + O/A 이상치) | 보정 데이터셋 품질 확보 | 2개 전문가 |
| P0 | C_ext=0.473 M 전용 파라미터 프로필 구축 | Ni MAE 26→5 g/L 이하 목표 | 2개 전문가 |
| P1 | 이온강도 기반 activity 보정 도입 (SIT 모델) | 극고농도 정확도 향상 | 2개 전문가 |
| P1 | 사포니피케이션 수지 검증 + sap truncation 기록 | 모델 투명성 확보 | 3개 전문가 |
| P1 | 3단 전용 warmup 스케일링 + alpha 하한 0.05 상향 | solver 안정성 확보 | 2개 전문가 |
| P2 | Stage efficiency factor 도입 (기본 1.0) | 극고농도 Ni 대응 | 1개 전문가 |
| P2 | C_org_fresh 수동 입력 → tear-stream 부분 보정 | MAE 3-8 g/L 감소 | 1개 전문가 |

### 위험 라벨 통계

| 라벨 | 건수 | 섹션 |
|------|------|------|
| `validation_breaking` | 5 | §2.4, §2.6, §3.5, §4.4, §5.5 |
| `accuracy_limiting` | 18 | §2.1-2.3, §2.5, §2.7, §3.1-3.4, §3.6, §4.1-4.3, §4.5-4.6, §5.1, §5.3-5.4, §5.6-5.7 |
| `cosmetic` | 2 | §4.1(BM-CY), §5.2 |

---

## 2. 평형 전문가 검토 (sx-equilibrium-expert)

> **모델 정체성 선언**: 본 시뮬레이터는 pH-시그모이드 기반 반경험적 보정 모델이다. `distribution_coefficient()`는 질량수지 공식 내부의 pseudo-D이며, 완전한 열역학적 엄밀성을 가진 분배 계수가 아니다 (working-rules 준수). 아래 분석은 코드 증거와 문헌 근거를 분리하여 제시한다.

---

### 2.1 Ni pH50=6.55 override의 C_ext=0.473M 적용성

**코드 위치:** `config.py:252-255` (SITE_PARAMETER_OVERRIDES), `extraction_isotherm.py:84` (pH50 보정식)

**현재 상태:**
```python
SITE_PARAMETER_OVERRIDES = {
    ("Cyanex 272", "Ni"): {"pH50": 6.55, "k": 2.0},
}
```
pH50 보정식 (`extraction_isotherm.py:84`):
```python
pH50_eff = pH50_ref - alpha * math.log10(C_ext / C_ref) + beta * (T - T_REF)
```
Cyanex 272 Ni의 보정 파라미터: alpha=0.7, C_ref=0.5 M, beta=-0.03.

**정량적 분석:**

Data1~3 보정 조건 (C_ext=0.631 M):
- alpha × log10(0.631/0.5) = 0.7 × 0.1013 = -0.0709
- pH50_eff = 6.55 + 0.071 = 6.479

Ni-Cy 현장 (C_ext=0.473 M):
- alpha × log10(0.473/0.5) = 0.7 × (-0.0241) = +0.0169
- pH50_eff = 6.55 - 0.017 = 6.567

**두 조건 간 pH50_eff 차이: 0.088 pH 단위.** C_ext가 25% 감소했는데 pH50 이동은 고작 0.09 pH 단위이다. 시그모이드 k=2.0에서 0.09 pH 차이는 추출률 변화 약 4.5%p에 해당하므로, 이 보정은 C_ext 변화를 사실상 흡수하지 못한다.

**이론적 alpha 추정:** 고전 평형식 logD = logK_ext + n×log[HA]₂ + 2×pH에서, pH50 (D=1)을 C_ext로 미분하면 dpH50/dlog(C_ext) = -n_ext/2. Ni-Cyanex 272의 n_ext=2이므로 이론적 alpha = 1.0이다. 현재 alpha=0.7은 이론값 대비 30% 과소하여 넓은 C_ext 범위에서의 보정력이 구조적으로 부족하다.

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. C_ext=0.473 M 데이터 기반의 독립적 Ni pH50 override 프로필 구축
2. alpha 값을 이론적 근거(n_ext/2=1.0)에 가깝게 상향 검토 (기존 Data1~3 적합도와의 균형 필요)
3. 또는 C_ext 의존성을 pH50 보정이 아닌 D 보정 계수로 분리 처리

### 2.2 n_eff_multiplier 바닥값 0.35의 물리적 근거

**코드 위치:** `extraction_isotherm.py:571-572`

**현재 상태:**
```python
n_eff_multiplier = max(0.35, 1.0 - (loading_ratio_clamped ** 2))
```

`loading_ratio_clamped`는 0~1 범위 유기상 로딩비. 로딩비 약 81% 이상에서 바닥값 0.35에 고정된다.

**수학적 형태 분석:**
- loading=0.0: multiplier=1.0 | loading=0.5: 0.75 | loading=0.7: 0.51 | loading≥0.81: 0.35

바닥값 0.35는 n_ext=2일 때 n_eff=0.7을 의미하며, 극한 로딩에서 금속 1몰당 추출제 dimer 소비가 2.0에서 0.7으로 65% 감소한다는 가정이다.

**물리적 근거 평가:**
- **방향적 타당성**: Santanilla 2021 — Ni/Co 추출 착물 구조 차이(Ni: octahedral hydrated), 고로딩 시 다핵 착물 형성으로 유효 화학양론 감소 가능
- **정량적 근거 부재**: 0.35 바닥값과 제곱 감쇠 형태에 대한 직접적 실험/문헌 근거 없음. Lu et al. 2024에서도 이 형태의 n_eff 감쇠 미사용
- **수치적 편의성**: Data5 (Ni 71 g/L, D2EHPA) 과로딩 방지를 위한 실용적 선택으로 추정

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. 바닥값 0.2~0.5 범위 민감도 분석 (Data5, Ni-Cy)
2. 금속별 차별적 바닥값 도입 검토 (Ni vs Co 착물 구조 차이 반영)
3. 문서에 "수치적 정규화 파라미터"임을 명시

### 2.3 loading_damping_factor와 경쟁추출 모드의 이중 감쇠

**코드 위치:** `extraction_isotherm.py:290-325` (loading_damping_factor), `extraction_isotherm.py:606-615` (competition_factor)

**현재 상태:**

메커니즘 1 — loading_damping_factor: 로딩 85%에서 factor=0.5, 95%에서 ~0.0 (시그모이드 cliff)
메커니즘 2 — competition_factor: `(HA_free/C_ext)^n_eff`, n_eff = n_ext × n_eff_multiplier

**이중 감쇠 분석:**

경쟁추출 모드 ON 시, n_eff_multiplier가 **이중으로 작용**한다:
1. HA_free를 높이고 (추출제 소비 계산 완화)
2. 지수를 낮춘다 (감쇠 효과 약화)

두 효과가 같은 방향(감쇠 완화)으로 작용하므로, 고로딩에서 감쇠가 예상보다 약해져 과추출을 야기할 수 있다. 반면 경쟁 OFF에서는 loading_damping_factor의 sharp sigmoid cliff가 작동한다. 두 경로의 감쇠 특성이 질적으로 다르다.

Vasilyev 2019의 공유 추출제 풀 개념에서는 `(HA_free/C_ext)^n_ext`가 로딩 감쇠를 내재적으로 표현하므로, 별도의 loading_damping_factor와 역할 중복이다.

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. 경쟁추출 모드 ON 시 loading_damping_factor 적용 여부를 명확히 문서화
2. n_eff_multiplier의 이중 완화 효과 sensitivity 분석 수행
3. 장기적으로 단일 통합 프레임워크 (Vasilyev 스타일 공유풀만으로 감쇠 표현) 검토

### 2.4 종분화 1:1 착물 한계: 고농도 Ni의 activity coefficient 미반영

**코드 위치:** `extraction_isotherm.py:343-383` (get_aqueous_speciation_state), `config.py:435-444` (SPECIATION_CONSTANTS)

**현재 상태:**
```python
# extraction_isotherm.py:370-376
hydroxo_term = K_MOH * free_OH
sulfate_term = K_MSO4 * free_sulfate
denominator = max(1.0, 1.0 + hydroxo_term + sulfate_term)
free_fraction = 1.0 / denominator
```

Ni 상수: K_MOH = 10^4.1, K_MSO4 = 195.0

**정량적 평가 (Ni-Cy 조건):**

Ni 70 g/L = 1.193 M → SO₄²⁻ ~ 1.19 M → 이온강도 I ≈ 4.77 M (완전 해리 가정). 실제 NiSO₄⁰ 이온쌍으로 유효 I는 낮지만 여전히 I >> 1 M.

**문제점:**
1. **1:1 착물만 고려**: `Ni(SO₄)₂²⁻` 등 고차 착물 누락
2. **Activity coefficient 부재**: 모든 평형 계산이 농도 기반. I > 3 M에서 Davies 식조차 적용 불가 (I < 0.5 M 범위), SIT 또는 Pitzer 모델 필요
3. **Free sulfate 피드백 부재**: `get_free_sulfate_concentration` (`extraction_isotherm.py:328-340`)이 금속-황산염 착물 형성에 의한 free sulfate 감소를 미반영
4. **Cyanex 272-Ni sulfate D correction 미정의**: `SULFATE_D_CORRECTION_RULES` (`config.py:456-460`)에 D2EHPA Co/Ni/Li만 존재

**위험 라벨:** `validation_breaking`

**권장사항:**
1. 단기: Cyanex 272-Ni SULFATE_D_CORRECTION_RULES 추가
2. 중기: SIT 모델 기반 이온강도 보정 도입 (I < 4 M 적용 가능, Pitzer보다 구현 간단)
3. Free sulfate 계산에 금속-황산염 착물 피드백 루프 추가

### 2.5 사포니피케이션 pH50 shift 선형 근사의 한계

**코드 위치:** `config.py:392-415` (SAPONIFICATION_P_H50_SHIFT), `extraction_isotherm.py:141-156` (get_saponification_pH50_shift)

**현재 상태:** pH50 shift = coefficient × sap_fraction (순수 선형). Ni-Cyanex 272: shift = 0.65 × sap_fraction.

**비선형 한계:**
1. **계면 화학 포화**: Keller 2022 — Na-form 분율이 높아질수록 marginal 효과 감소 (추출제 site 유한). 선형 모델은 포화 무시
2. **양성자 수지 이중 계산**: `get_proton_release` (`extraction_isotherm.py:670-672`)에서 `n_H × (1 - sap_fraction)` 감쇠와 pH50 shift가 독립 적용되어 이중 계산 위험
3. **D2EHPA 비교**: D2EHPA Ni/Co sap shift 계수 1.40 (Cyanex 272의 2배) — 같은 선형 근사 적용 시 고사포니피케이션에서 오차 증폭

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. 포화형 함수 교체 검토: `shift = shift_max × sap_fraction / (K_half + sap_fraction)`
2. proton release 감쇠와 pH50 shift의 이중 효과 분리 정량화

### 2.6 BM-CY Li MAE 32.4% 원인 분석

**코드 위치:** `extraction_isotherm.py:483-651` (compute_competitive_extractions), `config.py:257-259` (Li override)

**현재 상태:** Li 시뮬레이션 0.1~7.8% vs 실측 27.2~46.2%, MAE 32.4%p.

Li pH50 = 8.10 (override), BM-CY stage pH ≈ 5.6-6.3 → pH 6.0에서 시그모이드 E = 2.1%, pH 6.3에서 3.6%. 사포니피케이션 shift 0.40도 부족하여 유효 pH50 = 7.90 → E ≈ 2.7%.

**과소예측 메커니즘:**
1. 경쟁추출에서 Li는 최저 우선순위 (`config.py:348` EXTRACTION_PRIORITY) — Co의 추출제 소비(C_ext의 ~21.5%)로 Li 잔여 추출제 감소
2. **실측 Li 추출의 물리적 원인**: Swain 2006 — "과잉 Na-Cyanex 272가 Co 추출 완료 후 Li를 끌어들임". Na-form 추출제의 Li 추출력이 HL form보다 훨씬 강함
3. 실측 35-46% 재현에는 유효 pH50 ≈ 6.0-6.5 필요 → sap shift 1.6-2.1 필요 (현재 0.40의 4-5배)

**위험 라벨:** `validation_breaking`

**권장사항:**
1. Li sap shift 대폭 상향 (0.40 → 1.5-2.0) 또는
2. Na-form 추출 경로를 별도 모델링: `D_Li = D_sigmoid_HL + D_NaL`

### 2.7 (신규) C_ext alpha 보정의 비선형 한계 정량 분석

**코드 위치:** `extraction_isotherm.py:84`, `config.py:95-106`

**현재 상태:** `pH50_eff = pH50_ref - alpha × log10(C_ext/C_ref)`

**정량적 범위별 괴리:**

| C_ext (M) | alpha=0.7 보정 | 이론적 (n/2=1.0) | 괴리 (pH) | 추출률 오차 |
|-----------|--------------|-----------------|----------|-----------|
| 0.30 | +0.155 | +0.222 | 0.067 | ~3.3%p |
| 0.473 | +0.017 | +0.024 | 0.007 | ~0.3%p |
| 0.631 | -0.071 | -0.101 | 0.030 | ~1.5%p |
| 1.00 | -0.211 | -0.301 | 0.090 | ~4.5%p |

alpha=0.7 vs 이론값 1.0의 괴리는 C_ref(0.5M)에서 멀어질수록 증가. Ni alpha=0.7만이 이론값 대비 유의미하게 낮음 (Li 0.5=이론, Co 0.8=80%).

추출제 자기 회합(dimerization) 효과: C_ext 변화 시 단량체/이합체 비율 변화로 유효 추출제 농도가 C_ext와 비선형 관계. 현재 모델은 미반영.

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. Ni alpha를 0.85-1.0 범위로 상향 검토
2. C_ref를 Data1~3와 Ni-Cy 중간값(~0.55 M)으로 이동 검토

### 2.8 평형 전문가 종합 위험 매트릭스

| # | 이슈 | 위험 라벨 | 코드 위치 | 영향 범위 |
|---|------|----------|----------|----------|
| 2.1 | Ni pH50 C_ext 전이 실패 | `accuracy_limiting` | config.py:252, extraction_isotherm.py:84 | Ni-Cy 전체 |
| 2.2 | n_eff_multiplier 0.35 근거 부재 | `accuracy_limiting` | extraction_isotherm.py:571-572 | 고로딩 |
| 2.3 | 이중 감쇠 상호작용 미검증 | `accuracy_limiting` | extraction_isotherm.py:290-325, 606-615 | 경쟁+고로딩 |
| 2.4 | 종분화 1:1 + activity 미반영 | `validation_breaking` | extraction_isotherm.py:343-383 | Ni>1M, I>3M |
| 2.5 | Sap pH50 shift 선형 한계 | `accuracy_limiting` | config.py:392-415 | sap>0.5 |
| 2.6 | Li 과소예측 32.4%p | `validation_breaking` | extraction_isotherm.py:483-651 | Li 전체 |
| 2.7 | C_ext alpha 비선형 괴리 | `accuracy_limiting` | extraction_isotherm.py:84 | C_ext<0.4 or >0.8M |

---

## 3. 플로우시트 전문가 검토 (sx-flowsheet-simulation-expert)

> 본 섹션은 decision-rubric, working-rules, deep-dive-notes, paper-summaries에 기반하여 현재 시뮬레이터의 플로우시트 수준 기능을 독립적으로 평가한다.

---

### 3.1 Capability Level 판정

**코드 위치:** `multistage_sx.py:169-812` (역류 solver), `alkali_contract.py:34-101`

**현재 상태:**
- 단일 추출 블록 내 N-stage 역류 counter-current 계산
- Stage별 pH 제어(`target_pH_per_stage`), 균일 target pH, 고정 NaOH 세 가지 모드
- AlkaliContract dataclass로 알칼리 적용 경로 명시적 분류
- Scrub, strip, organic recycle, bleed 블록 부재. `C_org_fresh`는 `{m: 0.0 for m in metals}` 기본값

**판정: Level 2, 부분적 Level 3 요소 보유**

- **Level 2** 완전 충족: 역류 다단 추출 solver 작동, relaxation 기반 수렴
- **Level 3 부분**: stage별 pH 목표 + NaOH 분배 전략 존재. 그러나 stage efficiency, entrainment, bypass 미구현
- **Level 4 미달**: scrub/strip 블록, organic recycle, tear-stream convergence 전무

Galvez et al. 2004: extraction-stripping coupled solvent-loop 필요. Vasilyev et al. 2019: loading-scrubbing-stripping 전체 포괄 모델. 현재는 loading 부분만 구현.

**위험 라벨:** `accuracy_limiting`

**권장사항:** "단일 추출 블록의 역류 cascade solver"로 명확히 포지셔닝. Circuit-level 주장은 strip 블록 + organic recycle 구현 후에만 가능.

### 3.2 (신규) Tear-stream 부재의 과추출 기여도 정량 추정

**코드 위치:** `multistage_sx.py:210-211` (C_org_fresh 기본값), `multistage_sx.py:278, 467`

**현재 상태:**
```python
C_org_fresh = {m: 0.0 for m in metals}  # line 211
```

마지막 stage에서 항상 fresh organic(금속 0) 투입 — strip 완벽 가정.

**Ni-Cy 정량 추정 (04.16일 대표 조건):**

- Feed Ni 67 g/L, Q_aq 1.2 L/hr, Q_org 14.7 L/hr, O/A=12.3
- 유기상 총 Ni 용량 ≈ 204 g/hr, Ni 투입량 = 80.4 g/hr
- Strip 후 잔류 Ni 0.5 g/L 가정 → 7.35 g/hr (투입량의 9.1%)
- **Raffinate Ni 변화 추정: 3-8 g/L (MAE 26 g/L의 12-31%)**

Tear-stream 부재 과추출 편향은 유의미하지만 지배적이지 않은 기여 요인. deep-dive-notes #4: "recycle과 tear stream은 아키텍처 문제이지 단순 bookkeeping이 아니다."

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. 단기: `C_org_fresh`를 대시보드에 노출, strip 후 잔류 Ni 수동 입력 (MAE 3-8 g/L 감소 기대)
2. 장기: Strip 블록 + organic recycle + tear-stream convergence (Wegstein/direct substitution)

### 3.3 3단 vs 5단 relaxation warmup 적합성

**코드 위치:** `multistage_sx.py:62-105` (`_get_relaxation_alpha`)

**현재 상태:**

Target pH 모드 warmup: it<40: α=0.10 | it<120: α=0.18 | it<250: α=0.28 | else: α=0.40

적응 감쇠: diff>10.0→cap 0.05 | diff>3.0→cap 0.08 | diff>1.0→cap 0.12

**평가:**
- 3단: 자유도 21개(3×7 금속), 5단: 35개. 3단이 본질적으로 빠르게 수렴
- 그러나 warmup(40/120/250)이 n_stages에 무관하게 고정 → 3단에서 과보수적
- Ni-Cy 3단 초기: seed=0 → diff>10 → alpha=min(0.10, 0.05)=0.05 → 극보수적 수렴
- continuation_levels = [0.35, 0.65, 0.85, 1.0] (`multistage_sx.py:765`)도 5단 튜닝으로 추정

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. Warmup을 n_stages 비례로: `warmup_1 = 8×n_stages` → 3단(24, 72, 150), 5단(40, 120, 250)
2. 3단 전용 수렴 벤치마크 구축

### 3.4 NaOH 사포니피케이션 stage별 분배 불명확성

**코드 위치:** `multistage_sx.py:42-59` (`_build_naoh_distribution`), `multistage_sx.py:155-166` (`_resolve_stage_saponified_input_mol_hr`)

**현재 상태:**

**Physical_v2 경로** (현장 기본): 총 sap mol flow를 마지막 stage(N)에서만 유기상에 투입. 앞 stage들은 sap_profile_seed를 통해 간접 전파.

**Aqueous_direct 경로**: uniform/front_loaded/custom 분배.

**평가:**
- Physical_v2의 단일 투입점은 물리적으로 합리적 (saponification mixer → extraction bank 마지막 stage)
- 그러나 sap inventory 역류 전파의 3단 수렴 특성에 대한 독립 테스트 부재
- Ni_Cy.csv에 NaOH 물리적 반응 위치(단일 vs 분산) 미기록

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. 현장에서 NaOH-유기상 반응 위치 메타데이터 확인
2. 다지점 사포니피케이션 지원 확장 검토

### 3.5 비정상상태 데이터 분리 기준

**코드 위치:** `scripts/field_validation_main.py:361-366, 428-433`

**현재 상태:**

비정상 필터링:
1. `if entry["date"] == "25.03.11": continue` — 하드코딩 단일 날짜
2. `if entry["Q_naoh_lhr"] is None: continue` — NaOH 결측 제외
3. `if entry["raff3_ph"] is None: continue` — pH 결측 제외

**평가:**
- 03.12~03.24 구간에서도 비정상 징후 뚜렷: 1단 Raff Ni 34,207 mg/L(안정기 대비 10배), NaOH 25→10 wt% 전환, O/A 6.0~9.1(안정기 8.9~20.4)
- "startup", "transient" 등 키워드로 코드 검색 결과: 자동 판별 로직 전무
- NaOH 농도 전환점, O/A 비 변화, Feed pH 변동이 자연스러운 데이터 분류 경계이나 미활용

**위험 라벨:** `validation_breaking`

**권장사항:**
1. NaOH 농도 전환점(25→10 wt%)으로 Group A(03.12~03.24) vs Group B(04.04 이후) 분리
2. 자동 판별: Raff 1st Ni > Feed Ni×50% → 비정상 플래그 등
3. Aspen/MetSIM에서 표준인 통계적 outlier detection 전처리 파이프라인 도입

### 3.6 (신규) O/A 6~20 범위가 역류 수렴에 미치는 영향

**코드 위치:** `multistage_sx.py:450-576`, `multistage_sx.py:108-138`

**현재 상태:**

Ni-Cy O/A = 6.0~20.4 (3.4배 범위). 시뮬레이터에서 O/A 적응 없음 — relaxation alpha가 `last_max_diff`에만 반응.

**수치적 특성:**
- O/A=6: 유기 로딩 집중 → loading_damping 빨리 활성화 → D 급변 → 수렴 진동 가능
- O/A=20: 유기 로딩 낮음 → damping 거의 비활성 → 수계 농도 변화 큼 → pH 변동
- 넓은 O/A를 단일 파라미터 세트로 피팅하는 것은 물리적으로 비현실적 (O/A 6과 20에서 로딩 프로파일 근본적 차이)

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. O/A 범위 그룹 분류 검증: 저(6~10), 고(10~20)
2. O/A>15 시 initial alpha 상향으로 수렴 가속
3. O/A 의존 사포니피케이션 효과 민감도 분석

### 3.7 플로우시트 종합

| 비교 축 | 현재 시뮬레이터 | Aspen/gPROMS급 | MetSIM/SysCAD급 |
|---------|---------------|---------------|----------------|
| 공정 블록 | 추출 1개 | 추출+scrub+strip+결정화 | 추출+strip+EW |
| Recycle closure | 없음 (C_org_fresh=0) | Wegstein/Broyden | Direct substitution |
| 정상/동적 | 정상 only | 정상+동적 | 정상 (일부 동적) |
| Capability Level | Level 2 (부분 3) | Level 5-6 | Level 4-5 |

---

## 4. 공정 전문가 검토 (sx-process-expert)

---

### 4.1 메타데이터 완성도 점수

#### 4.1.1 Ni-Cy

| # | decision-rubric 항목 | 상태 | 근거 |
|---|---------------------|------|------|
| 1 | NaOH 농도 | **부분** | 25→10 wt% 전환 기록, 실측이 아닌 설정값 |
| 2 | NaOH injection split | **불명** | "Saponification" 표기만, 물리적 위치 불명 |
| 3 | Feed basis | **알려짐** | 8종 금속+pH, mg/L, raw feed |
| 4 | 보고 pH 기준 | **부분** | stage 출구 평형값 vs mixer setpoint 불명 |
| 5 | 실제 O/A | **알려짐** | Org/Aq mL/min 기록 |
| 6 | 단수/블록 duty | **알려짐** | 3단 고정 |
| 7 | 증거 유형 | **알려짐** | 파일럿 mixer-settler 연속 운전 |

**점수: 4.5/7 → "medium confidence, caution required"**

추가 항목: 주요 불순물 세트(알려짐), 중화 위치(불명), 하류 목표(불명) → **종합: 5.5/10**

#### 4.1.2 BM-CY

**점수: 6.5/7 → "high process interpretability"**, 종합 8.5/10

- **위험 라벨:** Ni-Cy `accuracy_limiting` / BM-CY `cosmetic`
- **권장사항:** NaOH injection 위치, pH 측정 위치, 제품 규격을 현장 확인

### 4.2 O/A 비 해석 및 시뮬레이터 매핑

**코드 위치:** `field_validation_main.py:222-223, 230-233`

**현재 상태:**

Ni-Cy O/A 범위:
- 초기(03.11~03.21): 6.0~9.5 | 전환기(04.04~04.09): 9.4~17.2 | 안정기(04.10~05.13): 7.4~13.2

시뮬레이터 매핑: `mL/min × 0.06 = L/hr` → `Q_aq`, `Q_org` 직접 전달. Stage별 O/A 변동 없음.

**공정적 해석:**
- O/A=17.2 (04.09): 유기상 극과잉, 로딩률 ~31%, 추출 driving force 충분
- **산업적 타당성**: Cyanex 272 CoSX 일반 O/A = 4-8 (Swain 2008, Zhang 2025). O/A>10은 stripping 부하 과도 → 과도기 운전으로 추정

**위험 라벨:** `accuracy_limiting`

**권장사항:** O/A>12 구간 데이터는 과도기 분류, 시뮬레이터에 O/A 범위 경고 추가

### 4.3 NaOH 25wt%→10wt% 전환의 공정적 의미

**코드 위치:** `Ni_Cy.csv:10`, `field_validation_main.py:179`

**현재 상태:**

| 구간 | NaOH wt% | 유량 (mL/min) | mol/min | 추정 전환 이유 |
|------|----------|-------------|---------|--------------|
| 03.11~03.24 | 25 | 2.6-5.2 | 0.82-1.64 | — |
| 04.04~05.13 | 10 | 9.47-19.12 | 0.66-1.34 | 유량 제어 정밀도 향상 + 추출제 분해 위험 완화 |

25 wt% NaOH: 점도 ~6 cP, 소량 유량(2.6~5.2 mL/min)에서 pump 정확도 저하. 유기상 직접 접촉 시 국소 과포화 → 추출제 degradation 위험.

시뮬레이터에서 밀도 보간/몰농도 환산은 정확함 (25 wt%→7.96 M, 10 wt%→2.77 M, 문헌 대비 <0.2% 오차).

**위험 라벨:** `accuracy_limiting`

**권장사항:** 25 wt% / 10 wt% 구간을 별도 코호트로 분석

### 4.4 극고농도 Ni에서의 Stage Efficiency 가정

**코드 위치:** `single_stage.py:258-332`, `single_stage.py:684-718`

**현재 상태:** 모든 stage에서 100% efficiency (완전 평형) 가정. `C_aq_mol_out = total_metal_flow / (Q_aq_out + D × Q_org)`.

**평가:**
1. NiSO₄ 1M 수용액 점도 ~1.5-2.0 cP (순수 대비 50-100% 높음) → settler emulsion band 두꺼워짐
2. C_ext=0.473 M에서 Ni 최대 로딩 ~13.9 g/L, 3단 역류에서 한계 접근 시 D 급감
3. 프리셋 Data1~3 Ni 28~34 g/L에서 보정 → 62~78 g/L은 **검증되지 않은 외삽**
4. 현장 증거: 03.12(O/A=6.8, Ni=77.6 g/L)에서 1단 추출률 55.9% — NaOH 부족 + stage efficiency<100% 복합 효과

**위험 라벨:** `validation_breaking`

**권장사항:**
1. Stage efficiency factor `eta` 도입: `D = D_raw × damping × eta`, 기본값 1.0
2. Ni-Cy 조건 eta = 0.85-0.95 보정 시도
3. 또는 농도 의존 함수: `eta = 1 - beta × (C_total/C_threshold)^gamma`

### 4.5 (신규) NaOH 결측 패턴 분석

**데이터 위치:** `Ni_Cy.csv:9`, **코드 위치:** `field_validation_main.py:178, 213, 224`

**현재 상태:** 22일 중 7일 NaOH 결측.

| 결측 날짜 | 1단 Raff pH | 1단 Raff Ni (mg/L) | 패턴 |
|----------|------------|-------------------|------|
| 03.13 | 6.25 | 11,541 | 높은 추출률 |
| 03.19 | 6.21 | 11,651 | 높은 추출률 |
| 03.24 | 6.17 | 10,146 | 높은 추출률 |
| 04.04 | 6.07 | 10,260 | 높은 추출률 |
| 04.10 | 6.12 | 4,498 | 높은 추출률 |
| 04.15 | 6.28 | 1,692 | 매우 높은 추출률 |
| 04.17 | 6.37 | 634 | 극히 높은 추출률 |

**분석:**
1. **비정상 운전과 무관**: 결측 날짜의 pH 5.3~6.4로 NaOH 분명 투입됨 → 유량계 기록 누락
2. **25 wt% 구간 집중**: 결측 7건 중 4건이 25 wt% 구간 (8일 중 4일=50% 결측률). 저유량(2.6~5.2 mL/min) → 유량계 읽기 어려움
3. **시뮬레이터 처리**: NaOH 모드(Report 1)에서 완전 제외 (올바름). pH 평형 모드(Report 2)에서는 포함 (NaOH 불필요)

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. 결측 NaOH를 인접일 운전 조건+실측 pH로부터 역산하는 보간 로직 추가
2. 25 wt% 구간(50% 결측) 전체를 NaOH 모드 검증에서 제외, 진단용으로만 사용

### 4.6 (신규) field_validation_main.py NaOH 환산 방식 검토

**코드 위치:** `field_validation_main.py:370`, `datasets.py:215-238`

**NaOH 환산 경로:**
```
CSV "NaOH, mL/min" → ×0.06 → Q_naoh_lhr (L/hr)
→ estimate_naoh_molarity_from_wt_pct(wt%) → C_NaOH (M)
→ multistage_sx: C_NaOH, Q_NaOH 전달
→ single_stage: C_NaOH × Q_NaOH = NaOH mol/hr
```

**정확성:** 밀도 보간/몰농도 환산 오차 <0.2% — **정확하다**.

**문제점:**
1. **naoh_wt_pct 기본값**: `field_validation_main.py:213` — 결측 시 10.0 기본값 → 25 wt% 구간 결측 시 mol/hr 2.87배 과소평가 위험
2. **사포니피케이션 반응 수율 100% 가정**: 접촉 시간, 혼합 강도, NaOH 농도에 무관. 25 wt%의 국소 과농축은 10 wt% 대비 반응 수율 다를 수 있음
3. **BM-CY NaOH wt% 파싱**: "10wt%" 문자열 → `_parse_numeric` 파싱 불가 → None. BM-CY는 pH 평형 모드 사용으로 현재 무해하나 사포니피케이션 모드 확장 시 문제

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. naoh_wt_pct 기본값을 None으로 변경, 결측 시 경고 출력
2. BM-CY "10wt%" 문자열 파싱 지원 추가

### 4.7 공정 전문가 종합

| 항목 | 위험 라벨 | 우선순위 |
|------|----------|---------|
| 4.1 메타데이터 Ni-Cy | `accuracy_limiting` | 현장 확인 |
| 4.2 O/A 해석 | `accuracy_limiting` | 중 |
| 4.3 NaOH 전환 | `accuracy_limiting` | 중 |
| 4.4 Stage efficiency | `validation_breaking` | **고** |
| 4.5 NaOH 결측 | `accuracy_limiting` | 중 |
| 4.6 NaOH 환산 | `accuracy_limiting` | 저-중 |

**Process realism verdict:**
- Ni-Cy: `directionally usable` — 추세 해석 가능, 정밀 공정 권장 불가
- BM-CY: `well-specified` — 의미 있는 매핑 가능

---

## 5. 시뮬레이션 전문가 검토 (sx-simulation-expert)

> 본 섹션은 수치 시뮬레이션/solver 관점에서 구조적 건전성, 수렴 안정성, 검증 범위를 독립 분석한다.

---

### 5.1 프리셋 vs 현장 MAE 격차의 구조적 원인 진단

**코드 위치:** `config.py:252-255`, `config.py:465-472`, `field_validation_main.py:37-49`

**현재 상태:** Data1~6 MAE < 1 g/L vs Ni-Cy MAE 26.3 g/L.

**세 가지 구조적 원인:**
1. **파라미터 외삽 실패**: pH50=6.55는 C_ext=0.631 M 보정. C_ext=0.473 M에서 alpha 보정 -0.017 pH 단위 (무의미)
2. **로딩 체제 전환**: Ni-Cy stage당 유기상 로딩이 프리셋 대비 2-3배. loading_damping_factor 임계점도 프리셋 calibrate
3. **O/A 불일치**: Ni-Cy O/A 6-20 vs 프리셋 ~4.4-5.2. 사포니피케이션 상호작용 미포착

**결론:** overfitting + 구조적 모델 한계의 복합. 파라미터 재보정만으로는 MAE < 5 g/L 불가, 모델 구조 개선 선행 필요.

**위험 라벨:** `accuracy_limiting`

**권장사항:** C_ext=0.473 M + 3단 독립 파라미터 프로필 (기존 프리셋 보호 병렬 구성)

### 5.2 Bisection pH 범위 [-1, 14] 수치적 위험

**코드 위치:** `single_stage.py:885-886`

**현재 상태:** `pH_low = -1.0, pH_high = 14.0`

**평가:**
- 수렴 속도: 범위 15에서 tol=1e-4까지 bisection 17.2회 (범위 [0,12]면 16.9회 — 차이 미미)
- pH=-1에서 free_H=10 M → 첫 1-2회 비물리적 D값 계산 (D→0, overflow 아님)
- pH=14에서 OH⁻ 보정 = Q_aq_eff × 1.0 → H_out_actual 음수 가능 → error 부호 불안정
- 실제 발산 보고 없음

**위험 라벨:** `cosmetic`

**권장사항:** `pH_low = 0.0, pH_high = 12.0` 정적 제한, 또는 `max(0, pH_in-2) ~ min(12, pH_in+5)` 동적 범위

### 5.3 Relaxation Warmup의 3단 시스템 적합성

**코드 위치:** `multistage_sx.py:62-105`

**현재 상태:**

| 구간 | iteration | alpha (target_mode) |
|------|-----------|-----|
| 1단계 | < 40 | 0.10 |
| 2단계 | 40-119 | 0.18 |
| 3단계 | 120-249 | 0.28 |
| 4단계 | ≥ 250 | 0.40 |

적응 감쇠 cap: diff>10→0.05, diff>3→0.08, diff>1→0.12. 최종 `max(0.03, min(0.85, alpha))`.

**3단 분석:**
- 자유도 21개(3×7) vs 5단 35개 → 본질적으로 빠른 수렴
- Ni-Cy 초기: seed=0 → diff>10 → alpha=min(0.10, 0.05)=0.05 → 극보수적
- alpha=0.05에서 수렴: n > log(1e-6)/log(0.95) = 269회 → warmup 3단계 범위에서도 수렴 미달 가능

**위험 라벨:** `accuracy_limiting`

**권장사항:**
- Warmup 스케일링: `threshold = f(n_stages)` 또는 초기 seed 개선 (단일 stage 결과로 추정)
- 3단 전용 수렴 테스트 `test_verification.py`에 추가

### 5.4 VALIDATED_FIELD_WINDOW C_ext 범위 경고 메커니즘

**코드 위치:** `config.py:465-472`, `dashboard_service.py:136-295`

**현재 상태:**

`VALIDATED_FIELD_WINDOW.C_ext_m = (0.6053, 0.6308)` — 매우 좁은 범위.

`build_scope_assessment` (line 190-195)에서 범위 이탈 경고 → **대시보드 UI에서만 표시**. 엔진 수준(`solve_multistage_countercurrent`)에는 검증/경고 없음. `field_validation_main.py`에서 C_ext=0.473 M이 무경고 처리.

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. `solve_multistage_countercurrent` 반환에 `warnings` 리스트 추가
2. 외삽 정도 정량화: `distance_ratio = (C_ext - cext_min) / (cext_max - cext_min)`
3. n_stages 등급화 경고 (1단 차이 vs 2단 이상 차이)

### 5.5 사포니피케이션 OH⁻/Na⁺ 수지 검증 부재

**코드 위치:** `single_stage.py:501-584`, `single_stage.py:1011-1015`

**현재 상태:**

`saponification_inventory_error_mol_hr`가 결과 dict에 포함되나 **어디서도 검사/경고하지 않음**.
- `_build_result_payload`에서 미집계
- `test_verification.py`, `validation_test.py`에서 미검사
- NaOH 총 소비 = Na-extractant 전환 + 금속 교환 + 직접 중화 합산 검증 없음

`SAPONIFICATION_DIRECT_NEUTRALIZATION_FACTOR` = 0.25 (Cyanex 272) — pH/금속농도/로딩률에 무관하게 고정.

**위험 라벨:** `validation_breaking`

**권장사항:**
1. `_build_result_payload`에 총 sap inventory error 합계 추가, 5% 초과 시 경고
2. `test_verification.py`에 사포니피케이션 모드 전용 Na 수지 테스트 추가
3. DIRECT_NEUTRALIZATION_FACTOR를 pH 의존적으로 검토

### 5.6 (신규) NaOH 파라미터 전달 경로 전체 추적

**코드 위치:** `dashboard_service.py:45-87` → `multistage_sx.py:169-258` → `single_stage.py:831-876` → `extraction_isotherm.py:141-156`

**전달 경로:**

1. dashboard_service: `naoh_mode`, `C_NaOH`, `Q_NaOH`, `saponification_model` 조립 → **정보 손실 없음**
2. multistage_sx: `build_alkali_contract()` → AlkaliContract 결정 → **정보 손실 없음**
3. multistage_sx → single_stage (physical_v2): `fresh_sap_input_mol_hr` 계산 후 seed profile 분배 → **C_NaOH/Q_NaOH 맥락 소실** (의도적이나 부작용 있음)
4. single_stage: `_resolve_stage_saponification_state` (line 492) → `sap_mol_hr = min(requested, free_site_capacity)` → **sap truncation 미기록**

**핵심 정보 손실:**
- Step 3에서 원래 NaOH 조건의 맥락 소실 → stage 수준에서 sap inventory 현실성 판단 불가
- Step 4에서 sap truncation 미기록 → 투입 NaOH 중 실제 활용률 불가시

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. `sap_truncation_mol_hr` 기록 추가
2. 총 truncation > 입력 20% 시 경고
3. `saponification_utilization_pct` 결과에 포함

### 5.7 (신규) 3단 극저감쇠(alpha 0.03) 시나리오 분석

**코드 위치:** `multistage_sx.py:62-105`

**현재 상태:** alpha 하한 0.03 도달 조건: `relaxation_scale < 1.0`(fallback 0.70) × diff>10 cap 0.05 = 0.035 → max(0.03, 0.035) = 0.035.

**수렴 분석:**
- alpha=0.03: n > log(1e-6)/log(0.97) = 454회. MAX_ITERATIONS=500 → 한도 근접
- **중간 상태 해 위험**: loading_damping sigmoid + alpha 0.03 조합 → 유기상이 "낮은 로딩" 상태에서 탈출 불가 → 500회 후 best_org_out이 완전 수렴 해가 아닌 중간 상태
- Fallback (`multistage_sx.py:796-809`)에서 `relaxation_scale=0.70` → alpha = 0.05×0.70 = 0.035 → 하한 근처

**위험 라벨:** `accuracy_limiting`

**권장사항:**
1. alpha 하한 0.03 → 0.05 상향 (필요 반복 454→269)
2. 또는 MAX_ITERATIONS 1000으로 증가
3. `_get_relaxation_alpha`에 n_stages 의존 추가
4. 수렴 실패 시 `convergence_risk_flag = True` + best_residual > tolerance×100 경고

### 5.8 시뮬레이션 전문가 종합 위험 매트릭스

| 항목 | 위험 라벨 | 코드 위치 | 핵심 이슈 |
|------|----------|----------|----------|
| 5.1 프리셋 vs 현장 MAE | `accuracy_limiting` | config.py:252, 465 | 파라미터 외삽+로딩 전환+O/A 불일치 |
| 5.2 Bisection pH | `cosmetic` | single_stage.py:885 | 비물리적 범위, 발산 위험 낮음 |
| 5.3 Warmup 3단 | `accuracy_limiting` | multistage_sx.py:62-105 | 5단 기준 과보수적 |
| 5.4 C_ext 경고 | `accuracy_limiting` | dashboard_service.py:190 | 엔진 수준 경고 부재 |
| 5.5 Sap 수지 미검증 | `validation_breaking` | single_stage.py:1011 | inventory error 미검사 |
| 5.6 NaOH 전달 경로 | `accuracy_limiting` | multistage→single_stage | sap truncation 미기록 |
| 5.7 극저감쇠 수렴 | `accuracy_limiting` | multistage_sx.py:104 | alpha 0.03에서 수렴 불가 위험 |

---

## 6. 종합 개선 로드맵

### 6.1 단기 (1-2주)

| # | 개선안 | 대상 파일 | 기대 효과 | 근거 섹션 | 교차 확인 |
|---|-------|---------|----------|----------|----------|
| S1 | 비정상상태 데이터 전처리: NaOH 전환점+O/A 이상치 기준 코호트 분리 | `field_validation_main.py` | 보정 데이터 품질 확보 | §3.5, §4.5 | 2개 전문가 |
| S2 | C_ext=0.473 M 전용 SITE_PARAMETER_OVERRIDES 프로필 | `config.py` | Ni MAE 대폭 감소 | §2.1, §5.1 | 2개 전문가 |
| S3 | VALIDATED_FIELD_WINDOW 범위 확장 + 엔진 수준 경고 | `config.py`, `multistage_sx.py` | 외삽 인식 개선 | §5.4 | 1개 전문가 |
| S4 | Bisection pH 범위 [0.0, 12.0] 제한 | `single_stage.py:885-886` | 수치 안정성 | §5.2 | 1개 전문가 |
| S5 | C_org_fresh 대시보드 노출 (strip 후 잔류 수동 입력) | `multistage_sx.py`, `sx_dashboard.py` | MAE 3-8 g/L 감소 | §3.2 | 1개 전문가 |
| S6 | Sap inventory error 집계 + 5% 초과 경고 | `single_stage.py`, `multistage_sx.py` | 모델 투명성 | §5.5, §5.6 | 2개 전문가 |

### 6.2 중기 (1-2개월)

| # | 개선안 | 대상 | 기대 효과 | 근거 섹션 | 교차 확인 |
|---|-------|------|----------|----------|----------|
| M1 | 이온강도 SIT 모델 기반 activity 보정 | `extraction_isotherm.py` | 극고농도 정확도 | §2.4, §5.1 | 2개 전문가 |
| M2 | Cyanex 272-Ni sulfate D correction rule 추가 | `config.py`, `extraction_isotherm.py` | Ni-Cy 정확도 | §2.4 | 1개 전문가 |
| M3 | 사포니피케이션 pH50 shift 비선형 모델 (Langmuir형) | `config.py`, `extraction_isotherm.py` | sap 정확도 | §2.5 | 1개 전문가 |
| M4 | 3단 warmup n_stages 비례 스케일링 + alpha 하한 0.05 | `multistage_sx.py` | solver 안정성 | §3.3, §5.3, §5.7 | 2개 전문가 |
| M5 | Stage efficiency factor eta 도입 (기본 1.0) | `single_stage.py` | 극고농도 대응 | §4.4 | 1개 전문가 |
| M6 | NaOH 결측 역산 보간 로직 | `field_validation_main.py` | 데이터 활용도 | §4.5 | 1개 전문가 |
| M7 | Li sap shift 상향 또는 Na-form 별도 경로 | `config.py`, `extraction_isotherm.py` | Li MAE 32%→10% | §2.6 | 1개 전문가 |

### 6.3 장기 (3-6개월)

| # | 개선안 | 대상 | 기대 효과 | 근거 섹션 |
|---|-------|------|----------|----------|
| L1 | Strip 블록 + organic recycle + tear-stream convergence | 신규 모듈 | Level 4 달성 | §3.1, §3.2 |
| L2 | Lu et al. 2024 ESI 모델 기반 경쟁추출 고도화 | `extraction_isotherm.py` | 다금속 정확도 | §2.3 |
| L3 | 동적 상태 추적 (transient 모드) | 신규 모듈 | startup/shutdown | §3.5 |
| L4 | McCabe-Thiele 자동 생성 | `sx_dashboard.py` | 공정 해석 도구 | — |
| L5 | n_eff_multiplier 바닥값 금속별 차별화 | `extraction_isotherm.py` | 고로딩 정확도 | §2.2 |

---

## 부록 A: 핵심 코드-이슈 매핑

| 이슈 | 코드 위치 | 섹션 | 위험 라벨 |
|------|----------|------|----------|
| Ni pH50 override | `config.py:252-255` | §2.1 | accuracy_limiting |
| n_eff_multiplier 바닥값 | `extraction_isotherm.py:571-572` | §2.2 | accuracy_limiting |
| loading_damping + 경쟁추출 이중 감쇠 | `extraction_isotherm.py:290-325, 606-615` | §2.3 | accuracy_limiting |
| 종분화 1:1 한계 | `extraction_isotherm.py:343-383` | §2.4 | validation_breaking |
| Sap pH50 shift 선형 | `config.py:392-415, extraction_isotherm.py:141-156` | §2.5 | accuracy_limiting |
| Li 과소예측 | `extraction_isotherm.py:483-651, config.py:257-259` | §2.6 | validation_breaking |
| C_ext alpha 비선형 | `extraction_isotherm.py:84, config.py:95-106` | §2.7 | accuracy_limiting |
| C_org_fresh 기본값 | `multistage_sx.py:210-211` | §3.2 | accuracy_limiting |
| Relaxation warmup | `multistage_sx.py:62-105` | §3.3, §5.3, §5.7 | accuracy_limiting |
| NaOH stage별 분배 | `multistage_sx.py:42-59, 155-166` | §3.4 | accuracy_limiting |
| 비정상상태 필터링 | `field_validation_main.py:361-366` | §3.5 | validation_breaking |
| O/A 범위 영향 | `multistage_sx.py:450-576` | §3.6 | accuracy_limiting |
| Stage efficiency | `single_stage.py:258-332` | §4.4 | validation_breaking |
| NaOH 환산 | `field_validation_main.py:370, datasets.py:215-238` | §4.6 | accuracy_limiting |
| Bisection pH 범위 | `single_stage.py:885-886` | §5.2 | cosmetic |
| VALIDATED_FIELD_WINDOW | `config.py:465-472, dashboard_service.py:190` | §5.4 | accuracy_limiting |
| Sap inventory error | `single_stage.py:1011-1015` | §5.5 | validation_breaking |
| NaOH 전달 경로 | `multistage_sx.py:306-324 → single_stage.py:857-876` | §5.6 | accuracy_limiting |

## 부록 B: 데이터 요약

### Ni-Cy (Ni_Cy.csv)
- 22개 날짜 (25.03.11 – 25.05.13)
- 추출제: Cyanex 272, C_ext = 0.473 M
- 3단 Mixer-Settler
- Feed Ni: 62–78 g/L (03.11 제외)
- O/A: 6.0–20.4 (산업 일반 4-8 대비 과도기 포함)
- NaOH: Saponification, 25 wt%(8일) → 10 wt%(14일) 전환
- NaOH 결측: 7/22일 (25 wt% 구간 4/8 = 50%)
- 비정상 판별 대상: 03.11-03.24 (NaOH 전환+O/A 저범위)

### BM-CY (BM-CY.csv)
- 11일차 (가동일수 1–11)
- 추출제: Cyanex 272, C_ext = 0.631 M
- 5단 Mixer-Settler
- Feed Ni: 28–38 g/L, Co: 3.0–4.1 g/L, Li: 5.3–6.3 g/L, Mn: 0.063–0.092 g/L
- O/A: 4.4–5.2 (산업 일반 범위 이내)
- NaOH: Saponification, 10 wt%, 8.2–12.0 mL/min (결측 없음)
- 메타데이터 완성도: 8.5/10 (Ni-Cy 5.5/10 대비 현저히 높음)

## 부록 C: Guardrails 준수 확인

| 규칙 | 준수 여부 | 비고 |
|------|---------|------|
| 현재 모델을 '완전한 MSE'로 부르지 않을 것 | ✓ | §2 서두에서 "반경험적 시그모이드 모델"로 선언 |
| 단일 캘리브레이션에서 '물리적 검증' 주장 금지 | ✓ | §2.1에서 명시적 경고, working-rules 인용 |
| 황산염→염화물 직접 전이 금지 | ✓ | 황산염 시스템 내 한정 |
| Cyanex 272↔D2EHPA 파라미터 직접 전이 금지 | ✓ | 추출제별 독립 논의 |
| 배치 평형→역류 회로 성능 주장 불가 | ✓ | §3.1에서 Level 2-3 명시, circuit 미주장 |
| 정상상태 결과를 제어 준비로 표현 금지 | ✓ | §3.7에서 동적 모드 부재 명시 |
| 보고서만 수정으로 일반 검증 주장 금지 | ✓ | §5.1에서 narrow-window + 구조적 한계 진단 |
| Overclaim 금지 | ✓ | 모든 하위 섹션에 위험 라벨 부착, 25건 전수 라벨링 |

## 부록 D: 기존 보고서 대비 변경/심화/신규 발견 분류

### 신규 발견 (기존 보고서에 없는 항목)

| 항목 | 섹션 | 발견 내용 |
|------|------|----------|
| C_ext alpha 비선형 한계 | §2.7 | alpha=0.7 vs 이론 1.0 괴리의 정량적 범위별 분석, 추출제 dimerization 효과 |
| Tear-stream 부재 정량 추정 | §3.2 | Ni-Cy 조건에서 MAE 3-8 g/L (12-31%) 기여 정량화 |
| O/A 범위 solver 영향 | §3.6 | O/A 6-20의 수렴 특성/파라미터 적합성 분석 |
| NaOH 결측 패턴 | §4.5 | 25 wt% 구간 50% 결측, 비정상 운전과 무관 확인 |
| NaOH 환산 경로 검토 | §4.6 | naoh_wt_pct 기본값 문제, BM-CY 파싱 한계 |
| NaOH 전달 경로 추적 | §5.6 | sap truncation 미기록, C_NaOH 맥락 소실 식별 |
| 극저감쇠 수렴 분석 | §5.7 | alpha 0.03에서 454회 필요, 중간 상태 해 위험 |

### 심화 분석 (기존 보고서에 있으나 더 깊게 분석)

| 항목 | 섹션 | 심화 내용 |
|------|------|----------|
| Ni pH50 적용성 | §2.1 | 이론적 alpha=n_ext/2=1.0 도출, 두 조건 간 0.088 pH 괴리 정량화 |
| n_eff_multiplier | §2.2 | loading ratio별 수학적 형태 완전 분석, 81% 바닥 도달점 식별 |
| 이중 감쇠 | §2.3 | n_eff_multiplier의 이중 완화 효과(HA_free↑ + 지수↓) 경로 분석 |
| 종분화 한계 | §2.4 | 이온강도 I 정량 계산, Davies 식 적용 불가 확인, free sulfate 피드백 부재 |
| Li MAE | §2.6 | pH별 시그모이드 E 정량 계산, Swain 2006 Na-form 메커니즘으로 sap shift 4-5배 필요 도출 |
| 비정상 필터링 | §3.5 | 03.12-03.24 구간별 정량 근거, 운전 조건 변경점 분류 기준 |
| Stage efficiency | §4.4 | 현장 증거 기반(03.12: 추출률 55.9%) 비현실성 실증 |
| Warmup 적합성 | §5.3 | alpha=0.05 → 수렴 269회 정량 추정, 3단 자유도 분석 |

### 유지 항목 (기존과 동일 결론, 독립 재확인)

| 항목 | 섹션 | 비고 |
|------|------|------|
| Capability Level 2-3 | §3.1 | 독립 판정에서 동일 결론 |
| NaOH stage별 분배 | §3.4 | Physical_v2 경로 추적 재확인 |
| O/A 매핑 방식 | §4.2 | mL/min→L/hr 환산 정확성 재확인 |
| NaOH wt% 전환 | §4.3 | 밀도 보간 정확성 재확인 |
| Bisection pH 범위 | §5.2 | cosmetic 판정 재확인 |
| VALIDATED_FIELD_WINDOW | §5.4 | 엔진 수준 경고 부재 재확인 |
| Sap 수지 미검증 | §5.5 | validation_breaking 재확인 |

### 정량 비교

| 지표 | 기존 보고서 | 하이브리드 보고서 | 변화 |
|------|-----------|----------------|------|
| validation_breaking 건수 | 3 | 5 | +2 (§4.4 Stage efficiency, §2.6 Li) |
| accuracy_limiting 건수 | 10 | 18 | +8 |
| 총 코드 참조 수 | 13 | 18 | +5 |
| 신규 발견 수 | 0 | 7 | +7 |
| 심화 분석 수 | 0 | 8 | +8 |
| 독립 교차 확인 항목 | 0 | 5 | +5 |

---

*본 보고서는 4개 전문가 Agent(sx-equilibrium-expert, sx-flowsheet-simulation-expert, sx-process-expert, sx-simulation-expert)를 독립된 컨텍스트에서 병렬 실행하여 생성하였다. 각 Agent는 해당 스킬의 SKILL.md, working-rules, decision-rubric, paper-summaries, deep-dive-notes를 직접 읽고, 다른 전문가의 분석 결과 없이 독립적으로 코드를 검토하였다. 이후 메인 대화에서 4개 섹션을 통합하고, Executive Summary와 교차 검증 분석을 추가하였다.*
