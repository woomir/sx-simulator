"""
Experimental Data Fitting Module
=================================
실험 데이터(pH vs 추출률)를 이용하여 시그모이드 또는 ESI 모델
파라미터를 자동으로 피팅하는 모듈.

지원 모델:
1. Sigmoid: E(pH) = E_max / (1 + exp(-k × (pH - pH₅₀)))
2. ESI: log(D) = log(C) + a × log([NaL])
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


# =========================================================================
# ESI 모델 피팅
# =========================================================================

def esi_model_D(NaL, log_C, a):
    """ESI 분배계수 함수 (피팅용). log(D) = log(C) + a × log([NaL])"""
    # NaL이 0인 경우 방지
    NaL_safe = np.maximum(NaL, 1e-20)
    log_D = log_C + a * np.log10(NaL_safe)
    return log_D


def fit_esi(pH_data, D_data, C_ext, K_a=7.674e-8, p0=None):
    """
    실험 데이터에 ESI 모델을 피팅합니다.

    pH → [NaL] 변환 후, log(D) vs log([NaL]) 선형 회귀.

    Parameters
    ----------
    pH_data : array-like — pH 값 배열
    D_data : array-like — 분배계수(D) 배열
    C_ext : float — 추출제 농도 (mol/L)
    K_a : float — 계면 해리 상수
    p0 : tuple, optional — 초기 추정값 (log_C, a)

    Returns
    -------
    dict:
        params: {log_C, a}
        errors: {log_C_err, a_err}
        r_squared: R²
    """
    pH_arr = np.array(pH_data, dtype=float)
    D_arr = np.array(D_data, dtype=float)

    # pH → [NaL] 변환
    ratio = K_a * (10.0 ** pH_arr)
    NaL = C_ext * ratio / (1.0 + ratio)

    # D가 0인 데이터 제거
    valid = D_arr > 0
    if np.sum(valid) < 2:
        return {"success": False, "error": "유효한 데이터 포인트가 부족합니다 (D > 0 필요)."}

    NaL_valid = NaL[valid]
    logD_valid = np.log10(D_arr[valid])

    if p0 is None:
        p0 = (1.0, 1.5)

    try:
        popt, pcov = curve_fit(
            esi_model_D, NaL_valid, logD_valid, p0=p0,
            bounds=([-5, 0.1], [10, 5]),
            maxfev=10000,
        )
        perr = np.sqrt(np.diag(pcov)) * 1.96

        logD_pred = esi_model_D(NaL_valid, *popt)
        ss_res = np.sum((logD_valid - logD_pred) ** 2)
        ss_tot = np.sum((logD_valid - np.mean(logD_valid)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

        return {
            "success": True,
            "params": {
                "log_C": round(float(popt[0]), 3),
                "a": round(float(popt[1]), 3),
            },
            "errors": {
                "log_C_err": round(float(perr[0]), 3),
                "a_err": round(float(perr[1]), 3),
            },
            "r_squared": round(float(r_squared), 4),
        }
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
        }


def extraction_to_D(E_pct, ao_ratio=1.0):
    """
    추출률(%)을 분배계수(D)로 변환합니다.
    E = D × (Q_org/Q_aq) / (1 + D × (Q_org/Q_aq)) × 100
    → D = E / ((100 - E) × ao_ratio)  (여기서 ao_ratio = Q_aq/Q_org)

    Parameters
    ----------
    E_pct : float or array — 추출률 (%)
    ao_ratio : float — A/O ratio (Q_aq / Q_org)

    Returns
    -------
    float or array — 분배 계수 D
    """
    E = np.array(E_pct, dtype=float)
    E = np.clip(E, 0.01, 99.99)  # 0과 100 방지
    D = (E / 100.0) / ((1.0 - E / 100.0) * (1.0 / ao_ratio))
    return D
