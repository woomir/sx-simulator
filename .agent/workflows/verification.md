---
description: SX 시뮬레이터 검증 절차 워크플로우
---

# SX 시뮬레이터 검증 워크플로우

새 기능을 구현하거나 모델을 변경한 후에는 반드시 아래 절차를 따라 검증합니다.

> ⚠️ **필수 규칙: UI 검증은 반드시 배포된 URL에서 수행합니다.**
>
> - 배포 URL: **https://sx-simulator.streamlit.app**
> - localhost에서의 검증은 허용되지 않습니다.
> - Git push 후 Streamlit Cloud가 자동 배포할 때까지 대기 후 확인합니다.

## 1. 검증 계획 작성

- `docs/verification/v{버전}_plan.md` 파일을 생성합니다.
- 아래 항목을 포함해야 합니다:
  - 검증 매트릭스 (조건 조합)
  - 탭별 입력 의존성 확인 항목
  - 물질수지/pH수지 정합성 체크
  - Phase 3 기능별 체크 (경쟁, 종분화 등)

## 2. 자동화 테스트 실행

// turbo

- `test_verification.py` 스크립트를 실행합니다:
  ```bash
  cd /Users/woomir/Project/회사/08.\ SX시뮬레인션\ -\ AHH
  .venv/bin/python test_verification.py
  ```
- 결과를 `docs/verification/v{버전}_results.md`에 기록합니다.

## 3. Git Push 및 배포 확인

- 모든 변경사항을 Git commit & push합니다.
- Streamlit Cloud 자동 배포 완료를 기다립니다 (약 1~2분).

## 4. 배포된 앱 UI 검증

- **반드시** https://sx-simulator.streamlit.app 에서 확인합니다.
- 핵심 시나리오를 브라우저에서 직접 확인합니다.
- 스크린샷을 캡처하여 결과 문서에 첨부합니다.
- 확인 항목:
  - 버전 번호가 올바른지
  - 사이드바 설정 변경 시 모든 탭이 업데이트되는지
  - pH Isotherm과 시뮬레이션 결과의 정합성

## 5. 결과 기록

- 발견된 이슈는 `v{버전}_results.md`의 "발견 이슈" 섹션에 기록합니다.
- 이슈 수정 후 재검증하여 체크리스트를 업데이트합니다.

## 6. 검증 문서 위치

- 검증 계획/결과: `docs/verification/`
- 테스트 스크립트: 프로젝트 루트 `test_verification.py`
