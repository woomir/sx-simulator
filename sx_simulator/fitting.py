"""
Experimental Data Fitting Module
=================================
실험 데이터(pH vs 추출률)를 이용하여 시그모이드 모델
파라미터를 자동으로 피팅하는 모듈.

지원 모델:
1. Sigmoid: E(pH) = E_max / (1 + exp(-k × (pH - pH₅₀)))
"""

import math
import numpy as np
from scipy.optimize import curve_fit


# =========================================================================
# 시그모이드 모델 피팅
# =========================================================================

def sigmoid_model(pH, pH50, k, E_max):
    """시그모이드 추출률 함수 (피팅용)."""
    x = -k * (pH - pH50)
    # overflow 방지
    x = np.clip(x, -500, 500)
    return E_max / (1.0 + np.exp(x))


def fit_sigmoid(pH_data, E_data, p0=None):
    """
    실험 데이터에 시그모이드 모델을 피팅합니다.

    Parameters
    ----------
    pH_data : array-like — pH 값 배열
    E_data : array-like — 추출률(%) 배열
    p0 : tuple, optional — 초기 추정값 (pH50, k, E_max)

    Returns
    -------
    dict:
        params: {pH50, k, E_max}
        errors: {pH50_err, k_err, E_max_err} — 95% 신뢰구간
        r_squared: R² 결정계수
        residuals: 잔차 배열
    """
    pH_arr = np.array(pH_data, dtype=float)
    E_arr = np.array(E_data, dtype=float)

    # 초기값 자동 추정
    if p0 is None:
        # pH50: 추출률 50% 근처 pH 추정
        E_mid = (np.max(E_arr) + np.min(E_arr)) / 2
        idx_mid = np.argmin(np.abs(E_arr - E_mid))
        pH50_guess = pH_arr[idx_mid]
        k_guess = 3.0
        E_max_guess = np.max(E_arr) * 1.05
        p0 = (pH50_guess, k_guess, E_max_guess)

    try:
        popt, pcov = curve_fit(
            sigmoid_model, pH_arr, E_arr, p0=p0,
            bounds=([0, 0.1, 10], [14, 20, 105]),
            maxfev=10000,
        )
        perr = np.sqrt(np.diag(pcov)) * 1.96  # 95% CI

        # R² 계산
        E_pred = sigmoid_model(pH_arr, *popt)
        ss_res = np.sum((E_arr - E_pred) ** 2)
        ss_tot = np.sum((E_arr - np.mean(E_arr)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

        return {
            "success": True,
            "params": {
                "pH50": round(float(popt[0]), 3),
                "k": round(float(popt[1]), 3),
                "E_max": round(float(popt[2]), 2),
            },
            "errors": {
                "pH50_err": round(float(perr[0]), 3),
                "k_err": round(float(perr[1]), 3),
                "E_max_err": round(float(perr[2]), 2),
            },
            "r_squared": round(float(r_squared), 4),
            "residuals": (E_arr - E_pred).tolist(),
        }
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
        }

