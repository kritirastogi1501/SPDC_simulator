# phase_match.py

import sys
sys.path.append("C:/Users/Asus/OneDrive/Desktop/Sanya Personal/Projects/DRDO internship/spdc_simulator/src")

import numpy as np
from scipy.optimize import brentq
from optics.sellmeier import n_z, n_y, poling_period

def delta_k(lambda_s_um, T, lambda_p_um, poling_period_25C, interaction_type='type0'):
    """
    Compute Δk = kₚ - kₛ - kᵢ - k_QPM,
    where kⱼ = (2π/λⱼ)·nⱼ(λⱼ,T), and
    k_QPM = 2π / Λ(T).
    interaction_type: 'type0' or 'type2'
    """
    try:
        # 1/λₚ = 1/λₛ + 1/λᵢ  ⇒  λᵢ = 1 / (1/λₚ - 1/λₛ)
        lambda_i_um = 1.0 / (1.0 / lambda_p_um - 1.0 / lambda_s_um)
        if lambda_i_um <= 0 or np.isnan(lambda_i_um):
            return np.nan

        # pick refractive indices depending on interaction type
        if interaction_type == 'type0':
            n_p = n_z(lambda_p_um, T)
            n_s = n_z(lambda_s_um, T)
            n_i = n_z(lambda_i_um, T)
        else:  # 'type2'
            n_p = n_z(lambda_p_um, T)
            n_s = n_y(lambda_s_um, T)
            n_i = n_z(lambda_i_um, T)

        if np.isnan(n_p) or np.isnan(n_s) or np.isnan(n_i):
            return np.nan

        k_p = 2.0 * np.pi * n_p / lambda_p_um
        k_s = 2.0 * np.pi * n_s / lambda_s_um
        k_i = 2.0 * np.pi * n_i / lambda_i_um
        k_qpm = 2.0 * np.pi / poling_period(poling_period_25C, T)

        return k_p - k_s - k_i - k_qpm

    except Exception:
        return np.nan


def find_phase_match_lambda_s(
    T,
    lambda_p_um,
    poling_period_25C,
    interaction_type='type0',
    search=(0.5, 1.2)
):
    """
    Find λₛ ∈ [search[0], search[1]] such that Δk(λₛ,T)=0 using Brent’s method.
    Returns None if there is no sign‐change in that interval.
    """
    def f(λs):
        return delta_k(λs, T, lambda_p_um, poling_period_25C, interaction_type)

    a, b = search
    fa = f(a)
    fb = f(b)

    if np.isnan(fa) or np.isnan(fb) or (fa * fb > 0):
        return None

    return brentq(f, a, b)


def find_degenerate_temperature(
    lambda_s_target_um,
    lambda_p_um,
    poling_period_25C,
    interaction_type='type2',
    search=(30.0, 100.0)
):
    """
    Find T ∈ [search[0], search[1]] so that Δk(lambda_s_target_um, T)=0.
    Raises ValueError if no sign‐change is found in that T‐range.
    """
    def f(T):
        val = delta_k(lambda_s_target_um, T, lambda_p_um, poling_period_25C, interaction_type)
        if np.isnan(val):
            raise ValueError(f"Δk(T) returned NaN at T={T:.2f}")
        return val

    a, b = search
    fa = f(a)
    fb = f(b)
    if fa * fb > 0:
        raise ValueError(
            f"No root for Δk=0 in T ∈ [{a:.2f}, {b:.2f}]: f({a:.2f})={fa:.3e}, f({b:.2f})={fb:.3e}"
        )
    return brentq(f, a, b)

