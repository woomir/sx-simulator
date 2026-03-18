#!/usr/bin/env python3
"""
Jantunen et al. (2022) Figure 3(b) 자동 디지타이징
====================================================
PDF 내장 이미지에서 색상 필터링 + 연결 컴포넌트 분석으로
pH-추출률(E%) 데이터 포인트를 자동 추출합니다.

파이프라인:
  1. PyMuPDF로 PDF 내장 JPEG 추출 (xref=158, page 7)
  2. Panel (b) 크롭 (col >= 1100)
  3. PIL HSV 변환 → 색상 필터링 (Ni/Li/Co/Mn)
  4. scipy.ndimage.label() → 컴포넌트 중심점 → (row, col)
  5. 보정 상수로 (row, col) → (pH, E%) 변환
  6. Na는 자동 감지 불가 → 기존 CSV 값 유지

Usage:
    python scripts/digitize_fig3b.py
"""

import sys
from pathlib import Path

import fitz
import numpy as np
from PIL import Image
from scipy import ndimage

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PDF_FILE = PROJECT_ROOT / "docs" / "literature" / "inbox" / "metals-12-01445.pdf"
CSV_FILE = PROJECT_ROOT / "data" / "jantunen2022_fig3b_cyanex272_OA2.5.csv"

# ─────────────────────────────────────────────────────────────────────
# 이미지 추출 상수
# ─────────────────────────────────────────────────────────────────────
XREF = 158                  # 내장 JPEG xref (2182×773px, page 7)
PANEL_B_CROP_LEFT = 1100    # Panel (b)는 우측 절반

# ─────────────────────────────────────────────────────────────────────
# 좌표계 보정 상수 (Plan agent 실험 + 틱마크 검증)
# ─────────────────────────────────────────────────────────────────────
# Panel (b) 좌표계 (panel_b = full_image[:, 1100:])
# X축: pH = (col_panel_b - PH_ORIGIN) / PH_SCALE
# Y축: E% = E_SLOPE * row + E_INTERCEPT
PH_ORIGIN = 124       # panel_b 내 pH=0 에 해당하는 col
PH_SCALE = 85.8       # pixels per pH unit
E_SLOPE = -0.18143    # E% per pixel (음수: 위=100%, 아래=0%)
E_INTERCEPT = 116.12

# 유효 플롯 영역 (E% 범위)
E_MIN_CLAMP = 0.0
E_MAX_CLAMP = 100.0

# ─────────────────────────────────────────────────────────────────────
# 색상 필터 정의 (PIL HSV: H=0-255, S=0-255, V=0-255)
# ─────────────────────────────────────────────────────────────────────
# H 범위는 PIL HSV 기준 (0-255 스케일)
COLOR_FILTERS = {
    "Ni": {
        "h_range": (60, 80),    # 녹색 다이아몬드
        "s_min": 100, "v_min": 100,
        "min_size": 80, "max_size": 300,
    },
    "Li": {
        "h_range": (20, 30),    # 주황색 사각형
        "s_min": 100, "v_min": 100,
        "min_size": 50, "max_size": 200,
    },
    "Co": {
        "h_range": (195, 210),  # 보라색(자주) 원형
        "s_min": 80, "v_min": 70,
        "min_size": 15, "max_size": 300,
    },
    "Mn": {
        "h_range": (230, 245),  # 핑크색 삼각형
        "s_min": 50, "v_min": 80,
        "min_size": 10, "max_size": 200,
    },
}

# Na 회색 마커: 자동 감지 불가 → 기존 CSV 값 사용
NA_FALLBACK = [
    (1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0), (5.0, 0.0),
    (5.5, 2.0), (6.0, 5.0), (6.5, 10.0), (7.0, 18.0),
    (7.5, 23.0), (8.0, 28.0),
]

# ─────────────────────────────────────────────────────────────────────
# 논문 텍스트 기반 앵커 포인트 (p.7, O/A=2.5)
# Figure 3(b) 마커 감지로는 정확한 pH=5.27 데이터가 없음
# 논문 본문에서 명시된 값을 앵커로 추가
# ─────────────────────────────────────────────────────────────────────
TEXT_ANCHORS_PH527 = {
    "Mn": 99.8,
    "Co": 99.0,
    "Ni": 13.9,
    "Li": 3.6,
}

# Li E=0% 기저선: 자동 감지 불가 (주황 마커가 x축과 구분 안 됨)
LI_ZERO_BASELINE = [(1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0), (5.0, 0.0)]


# ─────────────────────────────────────────────────────────────────────
# Step 1: 이미지 추출
# ─────────────────────────────────────────────────────────────────────
def extract_panel_b(pdf_path: Path) -> np.ndarray:
    """PDF에서 내장 JPEG를 추출하고 Panel (b)를 크롭합니다."""
    doc = fitz.open(str(pdf_path))
    pix = fitz.Pixmap(doc, XREF)
    img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
    arr = np.array(img)
    doc.close()

    panel_b = arr[:, PANEL_B_CROP_LEFT:]
    print(f"  이미지 추출 완료: {arr.shape[1]}×{arr.shape[0]}px → "
          f"Panel (b) 크롭: {panel_b.shape[1]}×{panel_b.shape[0]}px")
    return panel_b


# ─────────────────────────────────────────────────────────────────────
# Step 2-3: 좌표 변환
# ─────────────────────────────────────────────────────────────────────
def pixel_to_ph(col: float) -> float:
    """Panel (b) 열 좌표 → pH 값."""
    return (col - PH_ORIGIN) / PH_SCALE


def pixel_to_e(row: float) -> float:
    """Panel (b) 행 좌표 → E% 값 (클램핑 적용)."""
    e = E_SLOPE * row + E_INTERCEPT
    return max(E_MIN_CLAMP, min(E_MAX_CLAMP, e))


# ─────────────────────────────────────────────────────────────────────
# Step 4: 색상 필터링 + 컴포넌트 분석
# ─────────────────────────────────────────────────────────────────────
def detect_markers(panel_b: np.ndarray, metal: str) -> list[tuple[float, float]]:
    """지정 금속의 색상 필터를 적용하여 마커 중심점을 (pH, E%) 쌍으로 반환."""
    filt = COLOR_FILTERS[metal]
    h_lo, h_hi = filt["h_range"]

    # PIL HSV 변환
    img = Image.fromarray(panel_b)
    hsv = np.array(img.convert("HSV"))
    h, s, v = hsv[:, :, 0], hsv[:, :, 1], hsv[:, :, 2]

    # 색상 마스크
    mask = (h >= h_lo) & (h <= h_hi) & (s >= filt["s_min"]) & (v >= filt["v_min"])

    # 연결 컴포넌트 라벨링
    labeled, n_components = ndimage.label(mask)

    points = []
    for i in range(1, n_components + 1):
        component = np.argwhere(labeled == i)
        size = len(component)

        if size < filt["min_size"] or size > filt["max_size"]:
            continue

        # 컴포넌트 중심점
        center_row, center_col = component.mean(axis=0)
        ph = pixel_to_ph(center_col)
        e_pct = pixel_to_e(center_row)

        # pH 유효 범위 필터 (0-9)
        if ph < -0.5 or ph > 9.5:
            continue

        points.append((round(ph, 2), round(e_pct, 1)))

    # pH 순으로 정렬
    points.sort(key=lambda p: p[0])
    return points


# ─────────────────────────────────────────────────────────────────────
# Step 5: Mn/Co 보정 (E=0% 영역 추가)
# ─────────────────────────────────────────────────────────────────────
def supplement_data(
    points: list[tuple[float, float]], metal: str
) -> list[tuple[float, float]]:
    """자동 감지 데이터에 추론/앵커 포인트를 보충합니다.

    추가 항목:
    - E=0% 기저선: Mn/Co는 낮은 pH에서 x축 마커 미감지
    - Li E=0% 기저선: 주황 마커가 x축과 혼동
    - pH 5.27 앵커: 논문 텍스트(p.7) 기반 검증 포인트
    """
    if not points:
        return points

    existing_phs = {round(ph, 2) for ph, _ in points}
    additions: list[tuple[float, float]] = []

    # E=0% 기저선 추가
    min_ph = points[0][0]
    if metal == "Mn":
        for ph in [1.0, 2.0]:
            if ph < min_ph - 0.3:
                additions.append((ph, 0.0))
    elif metal == "Co":
        for ph in [1.0, 2.0, 3.0]:
            if ph < min_ph - 0.3:
                additions.append((ph, 0.0))
    elif metal == "Li":
        for ph, e in LI_ZERO_BASELINE:
            if ph not in existing_phs:
                additions.append((ph, e))

    # pH 5.27 앵커 포인트 (논문 텍스트 기반)
    if metal in TEXT_ANCHORS_PH527 and 5.27 not in existing_phs:
        additions.append((5.27, TEXT_ANCHORS_PH527[metal]))

    all_points = points + additions
    all_points.sort(key=lambda p: p[0])
    return all_points


# ─────────────────────────────────────────────────────────────────────
# Step 6: CSV 출력 + 비교
# ─────────────────────────────────────────────────────────────────────
def load_existing_csv(path: Path) -> dict[str, list[tuple[float, float]]]:
    """기존 CSV를 {metal: [(pH, E%), ...]} 형태로 로드합니다."""
    data: dict[str, list[tuple[float, float]]] = {}
    if not path.exists():
        return data

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("metal"):
                continue
            parts = line.split(",")
            if len(parts) != 3:
                continue
            metal = parts[0].strip()
            ph = float(parts[1])
            e = float(parts[2])
            data.setdefault(metal, []).append((ph, e))

    return data


def write_csv(
    path: Path,
    all_data: dict[str, list[tuple[float, float]]],
    metal_order: list[str],
):
    """CSV 파일을 작성합니다."""
    header = (
        "# Source: Jantunen et al. (2022), Metals 12, 1445, Figure 3(b)\n"
        "# Extractant: 0.8 M Cyanex 272 (diluted in Exxsol D80)\n"
        "# O/A = 2.5, T = 25±1°C, t_eq = 15 min\n"
        "# Feed [g/L]: Li 2.52, Ni 2.05, Co 16.5, Mn 2.09\n"
        "# Na introduced via NaOH pH adjustment (no fixed initial conc.)\n"
        "# Extraction method: auto-digitized from PDF (scripts/digitize_fig3b.py)\n"
        "# Calibration: PH_SCALE=85.8 px/pH, E_SLOPE=-0.18143, E_INTERCEPT=116.12\n"
        "# Ni/Li/Co: color-filtered marker detection; Mn: pink-filtered; Na: manual (grey undetectable)\n"
        "# Supplemented: pH 5.27 anchors from paper text (p.7); Li/Mn/Co E=0% baselines inferred\n"
    )

    with open(path, "w") as f:
        f.write(header)
        f.write("metal,pH,E_pct\n")
        for metal in metal_order:
            for ph, e in all_data.get(metal, []):
                f.write(f"{metal},{ph:.2f},{e:.1f}\n")


def print_comparison(
    old_data: dict[str, list[tuple[float, float]]],
    new_data: dict[str, list[tuple[float, float]]],
    metal_order: list[str],
):
    """기존 CSV vs 자동 감지 비교를 출력합니다."""
    print("\n" + "=" * 72)
    print("  기존 CSV vs 자동 디지타이징 비교")
    print("=" * 72)

    for metal in metal_order:
        old_pts = sorted(old_data.get(metal, []))
        new_pts = sorted(new_data.get(metal, []))

        print(f"\n  [{metal}] 기존 {len(old_pts)}개 → 자동 {len(new_pts)}개")
        print(f"  {'소스':>4} {'pH':>6} | {'E%':>8}")
        print(f"  {'-'*4}-{'-'*6}-+-{'-'*8}")

        for ph, e in old_pts:
            print(f"  {'기존':>4} {ph:>6.2f} | {e:>8.1f}")
        print(f"  {'-'*4}-{'-'*6}-+-{'-'*8}")
        for ph, e in new_pts:
            print(f"  {'자동':>4} {ph:>6.2f} | {e:>8.1f}")

    print("\n" + "=" * 72)


# ─────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────
METAL_ORDER = ["Mn", "Co", "Ni", "Li", "Na"]


def main():
    print("\n" + "=" * 72)
    print("  Jantunen et al. (2022) Figure 3(b) 자동 디지타이징")
    print("  0.8 M Cyanex 272, O/A = 2.5")
    print("=" * 72)

    if not PDF_FILE.exists():
        print(f"[오류] PDF 파일을 찾을 수 없습니다: {PDF_FILE}")
        sys.exit(1)

    # Step 1: 이미지 추출
    panel_b = extract_panel_b(PDF_FILE)

    # Step 4: 색상별 마커 감지
    print("\n  마커 감지 중...")
    all_data: dict[str, list[tuple[float, float]]] = {}

    for metal in ["Ni", "Li", "Co", "Mn"]:
        points = detect_markers(panel_b, metal)
        n_raw = len(points)

        # Step 5: 보충 포인트 추가 (E=0% 기저선, pH 5.27 앵커)
        points = supplement_data(points, metal)

        all_data[metal] = points
        added = len(points) - n_raw
        extra = f" (+{added} supplemented)" if added > 0 else ""
        print(f"    {metal:>2}: {n_raw}개 감지{extra}")
        for ph, e in points:
            print(f"        pH {ph:5.2f}  E {e:5.1f}%")

    # Na: 기존 CSV 값 사용
    all_data["Na"] = list(NA_FALLBACK)
    print(f"    Na: {len(NA_FALLBACK)}개 (기존 CSV 값 유지, 회색 마커 감지 불가)")

    # Step 6: 기존 CSV 로드 & 비교
    old_data = load_existing_csv(CSV_FILE)
    if old_data:
        print_comparison(old_data, all_data, METAL_ORDER)

    # CSV 저장
    CSV_FILE.parent.mkdir(parents=True, exist_ok=True)
    write_csv(CSV_FILE, all_data, METAL_ORDER)

    total = sum(len(pts) for pts in all_data.values())
    print(f"\n  → CSV 저장 완료: {CSV_FILE}")
    print(f"  → 총 {total}개 데이터 포인트 ({len(all_data)}종 금속)")
    print(f"\n  검증 명령:")
    print(f"    python scripts/plot_jantunen_fig3b.py")


if __name__ == "__main__":
    main()
