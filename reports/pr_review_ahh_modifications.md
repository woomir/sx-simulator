# PR Review: `colleague/ahh-modifications` — pKa 해리 모델 추가 및 Vasilyev 모델 기본화

## 변경 요약

동료 브랜치(`colleague/ahh-modifications`)의 핵심 변경 사항:

1. **Vasilyev 경쟁 추출 모델을 유일한 엔진으로 강제화** — `use_competition` 토글 제거, 경쟁 OFF 경로(`_partition_stage_with_damping`) 삭제
2. **pKa 해리 모델 추가** — `EXTRACTANT_PKA`, `pka_dissociation_factor()` 함수 추가 (단, 실제 계산에 미사용)
3. **대시보드 UI 대폭 단순화** — 고정 NaOH 모드, 수계 직접 투입, Stage별 차등 pH 등 제거
4. **브랜딩 변경** — "Wang/ALTA/MSE Framework" → "Vasilyev 경쟁 추출 모델 기반"으로 전면 교체
5. **Scope assessment 경고 다수 삭제**

---

## 1. 현장 데이터 성능 비교 — 핵심 우려 사항

### Report 1: NaOH 사포니피케이션 모드 (Ni-Cy 3단, 14개 데이터)

| 지표 | Main OFF | Main ON | Colleague |
|------|----------|---------|-----------|
| 최종 제품 Ni MAE (g/L) | 26.0 | 26.2 | **26.3** |
| 첫 접촉 Ni MAPE (%) | 10.7 | 10.5 | 10.5 |
| 최종 pH MAE | 1.3 | 1.3 | **0.3** |

- 세 버전 모두 1단 Ni 농도 예측에서 큰 오차 (실제 ~1-15 g/L → 예측 ~27-43 g/L)
- Colleague가 Main ON과 거의 동일 (차이 0.1-0.3 g/L) — 경쟁 모델 강제화의 자연스러운 결과
- pH 예측은 Colleague가 소폭 개선 (MAE 0.3 vs 1.3)
- **판정: NaOH sap 모드에서는 유의미한 차이 없음. 단, 근본적인 모델 오차는 어느 버전에서도 미해결**

### Report 2: pH 평형 모드 통합 비교 (Ni-Cy 20개 + BM-CY 11개 = 31개)

| 지표 | Main OFF | Main ON | Colleague |
|------|----------|---------|-----------|
| 전체 Raff Ni MAE (g/L) | **4.7** | 5.0 | 6.0 |
| 전체 추출률 MAE (%p) | **8.9** | **8.9** | 12.1 |
| 평균 추출률 편차 (%p) | +2.1 | -7.5 | **-11.6** |

데이터셋별 상세:

| Dataset | 지표 | Main OFF | Main ON | Colleague |
|---------|------|----------|---------|-----------|
| Ni-Cy | Raff MAE | **4.31** | 5.10 | 5.76 |
| Ni-Cy | Ex% MAE | **6.1** | 7.3 | 8.2 |
| BM-CY | Raff MAE | 5.33 | **4.72** | 6.54 |
| BM-CY | Ex% MAE | **14.0** | 11.9 | 19.1 |

- **Colleague 버전이 모든 통합 지표에서 최하위**
- 평균 추출률 편차 -11.6%p → **체계적으로 추출률을 과소추정** (실제보다 추출이 덜 되는 것으로 예측)
- BM-CY에서 특히 열화 심각 (Ex% MAE: 19.1 vs Main ON 11.9)

### Report 3: BM-CY 다금속 (5단, 11개) — 가장 심각한 문제

| 금속 | Main OFF MAE | Main ON MAE | Colleague MAE |
|------|-------------|-------------|---------------|
| Ni | 14.8 | 14.2 | **19.1** |
| Co | 0.2 | 0.2 | 0.1 |
| Li | 8.5 | 10.4 | **32.4** |
| Mn | 0.1 | 0.1 | 0.1 |
| Mg | 9.0 | 9.3 | 8.9 |

**Li 추출률 — 치명적 회귀:**
- 실제 Li 추출률: 26~49%
- Main OFF 예측: 28~49% (양호)
- Main ON 예측: 25~46% (양호)
- **Colleague 예측: 4~8%** (거의 0에 가까움)

이것은 **명백한 모델 회귀(regression)**입니다. Colleague 버전에서 Li가 사실상 추출되지 않는 것으로 예측되는데, 실제 현장에서는 26-49%가 추출됩니다.

경쟁 간섭 분석:
- Colleague의 Ni 경쟁 간섭: 3.5%p (Main OFF 0.5%p, Main ON 2.6%p)
- **경쟁 모델이 과도하게 작동하여 약한 금속(Li)의 추출을 억압**하는 것으로 보임

---

## 2. Colleague ≠ Main ON인 이유 분석

이론적으로 Colleague 코드는 Main + `use_competition=True`와 동일해야 합니다. 그러나 Report 3에서 Li MAE가 Main ON 10.4 vs Colleague 32.4로 대폭 차이가 납니다. 이 차이의 원인을 추적해야 합니다:

**가능한 원인:**
1. 검증 스크립트(`compare_phase3.py`)의 변경 — 동료가 이 파일도 수정했으므로, 입력 파라미터나 호출 방식이 달라졌을 가능성
2. 대시보드 UI 변경으로 인한 입력값 차이 — 사포니피케이션 전용 경로로 바뀌면서 이전과 다른 NaOH/sap 파라미터가 전달되었을 가능성
3. `build_simulation_kwargs` 변경 — 고정 NaOH 모드 경로 삭제로 kwargs 구성이 달라짐

**리뷰어 질문:** Colleague 데이터를 생성할 때 사용한 정확한 입력 파라미터와 스크립트 호출 방식을 확인해야 합니다. Main ON과 동일한 조건이었는지 검증 필요.

---

## 3. 코드 변경 상세 리뷰

### 3-1. 경쟁 모델 강제화 (single_stage.py)

**삭제된 코드:**
- `_partition_stage_with_damping()` (80줄) — 경쟁 OFF의 로딩 댐핑 경로
- `_solve_target_pH_stage_equilibrium()`의 이중 경로 (if use_competition / else)
- `_solve_with_fixed_NaOH()`의 이중 경로

**리네이밍:**
- `_solve_competitive_stage_state` → `_solve_stage_state`
- `compute_competitive_extractions` → `compute_stage_extractions`

**판정:** 코드 단순화 자체는 긍정적이나, **성능 저하가 확인된 상태에서 기존 경로를 완전 삭제하는 것은 위험합니다.** 최소한 실험적으로 경쟁 OFF 결과가 필요한 비교/검증 시나리오를 위해 유지하는 것이 바람직합니다.

### 3-2. pKa 해리 모델 (extraction_isotherm.py, config.py)

```python
# config.py — 새로 추가
EXTRACTANT_PKA = {
    "D2EHPA": 3.24,
    "Cyanex 272": 6.37,
}

# extraction_isotherm.py — 새 함수
def pka_dissociation_factor(pH, extractant):
    """비해리 분율: f = 1 / (1 + 10^(pH - pKa))"""
```

그러나 `compute_stage_extractions()` 내에서:
```python
# pKa 참고: 추출제 해리 상태는 시그모이드 파라미터(pH50, k)에
# 이미 내재되어 있으므로 경쟁 계수에 별도 적용하지 않습니다.
# pka_dissociation_factor(pH, extractant) 함수는 대시보드 표시용으로 제공됩니다.
```

**문제점:**
- 함수를 만들고 사용하지 않음 — 코드 정리 필요
- "이미 내재되어 있다"는 주장이 맞다면, 왜 별도 함수와 상수를 추가했는지 목적이 불명확
- 대시보드 표시용이라면 `dashboard_tabs.py`에 국한시키는 것이 적절

### 3-3. 대시보드 UI 변경 (sx_dashboard.py)

**삭제된 기능:**
- "고정 NaOH" pH 제어 모드 — `target_pH = None` 경로
- "수계 직접 투입" NaOH 적용 방식 — `naoh_mode = "aqueous_direct"` 옵션
- Stage별 차등 pH 입력 (`staged_pHs`)
- NaOH 분배 전략 (균등/전단집중/커스텀)
- NaOH 유량 직접 입력 (Q_NaOH → 항상 0.0)

**우려:**
- 고정 NaOH 모드는 현장에서 NaOH 유량을 직접 측정/입력하는 시나리오에 필요
- 사용자 옵션 대폭 축소로 유연성 상실
- 기존 프리셋이 "고정 NaOH"로 설정되어 있었는데, 강제로 "목표 pH"로 전환

### 3-4. 브랜딩/문서 변경

**변경 전:** "Wang/ALTA 계열 문헌을 참조한 준경험 황산염계 SX 모델"
**변경 후:** "Vasilyev 경쟁 추출 모델 기반 (Vasilyev et al. 2019)"

**문제:**
- 시뮬레이터의 근간은 **시그모이드(OLI/MSE 파라미터) + Vasilyev 경쟁 보정**의 하이브리드 모델
- 시그모이드 파라미터(pH50, k, E_max)는 Wang/ALTA/MSE 프레임워크에서 도출
- Vasilyev는 경쟁 보정 레이어일 뿐, 전체 모델의 기반이 아님
- 원 모델 출처를 삭제하는 것은 학술적으로 부적절하며, 외부 리뷰어에게 오해를 줄 수 있음

### 3-5. Scope Assessment 경고 삭제 (dashboard_service.py)

삭제된 주의/경고 사항들:
- "physical v2 사포니피케이션 경로는 현재 공식 회귀 보호 기준이 아님" 경고
- "고정 NaOH + 사포니피케이션 해석은 raw-feed replay/진단용 경로" 경고
- "D2EHPA + 목표 pH + 고황산염 + 다단 조합은 대표 취약 구간" 경고
- McCabe-Thiele 한계 설명 경고

**문제:** 이러한 경고는 사용자가 모델의 한계를 인지하게 하는 안전장치입니다. 제거하면 사용자가 모델 결과를 과신할 위험이 있습니다.

### 3-6. 코드 품질 이슈

- `from typing import Optional`이 `config.py`에서 두 번 import됨 (Line 1, Line 272)
- `from typing import Optional`이 모듈 docstring 위에 위치 (관례 위반)
- `multistage_sx.py`에서 `use_competition` 제거 후 trailing whitespace 남음
- McCabe-Thiele의 로딩 댐핑이 `loading_damping_factor()` 대신 인라인 `max(0.0, 1.0 - loading) if loading > 0.85 else 1.0`으로 대체 — 원래 시그모이드 감쇠와 다른 선형 감쇠로 동작이 달라짐

---

## 4. 종합 판정

### Merge 불가 (Request Changes)

**차단 사유:**
1. **Li 추출률 치명적 회귀** — BM-CY 다금속 시나리오에서 Li MAE 32.4%p (Main ON 10.4%p 대비 3배 악화). Colleague 코드가 Main ON과 이론적으로 동일해야 하는데 큰 차이가 발생하므로, 숨겨진 버그 또는 검증 조건 차이가 있음
2. **전 데이터셋에서 일관된 성능 하락** — 단일 지표도 Main 대비 개선된 것이 없음 (Report 1 pH MAE 제외)
3. **안전장치(경고/주의) 무단 삭제** — 모델 한계 인지를 위한 scope assessment 경고 제거

**개선 필요 사항:**
1. Li 추출률 회귀 원인 규명 및 수정
2. 경쟁 OFF 경로 완전 삭제 대신, deprecated 처리 또는 내부 비교용 유지
3. 브랜딩: "Vasilyev 경쟁 추출 모델"만 표기 → "MSE/시그모이드 + Vasilyev 경쟁 보정" 하이브리드 표기
4. Scope assessment 경고 복원
5. pKa 함수 — 사용하지 않을 거면 제거하거나, 향후 적용 계획을 명시
6. 코드 품질 (중복 import, trailing whitespace, docstring 위치)

**긍정적 측면:**
- 이중 코드 경로 제거로 유지보수 단순화 방향은 좋음
- 수식 탭에 계산 흐름 설명(Step-by-Step) 추가는 사용자 이해에 도움
- Saponification % 메트릭 추가는 유용
- `Optional[]` 타입 힌트 통일 (Python 3.9 호환)
