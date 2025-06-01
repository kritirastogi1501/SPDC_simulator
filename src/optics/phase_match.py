import sys
sys.path.append("C:/Users/Asus/OneDrive/Desktop/Sanya Personal/Projects/DRDO internship/spdc_simulator/src")

import numpy as np
from scipy.optimize import brentq
from optics.sellmeier import n_z, n_y, poling_period

# ------------------------------------------------------------------------------
# Phase-matching module for collinear SPDC in PPKTP crystals
# Implements:
#   - Eq. (5): Δk = kp - ks - ki - 2π / Λ(T)
#   - Eq. (6): Energy conservation: 1/λp = 1/λs + 1/λi
#   - Eq. (4): Λ(T) = Λ₀ × [1 + A(T–25) + B(T–25)²]
# ------------------------------------------------------------------------------

def delta_k(lambda_s_um, T, lambda_p_um, poling_period_25C, interaction_type='type0'):
    """
    Computes phase mismatch Δk in µm⁻¹ for collinear SPDC in PPKTP.

    Args:
        lambda_s_um: signal wavelength in µm
        T: temperature in °C
        lambda_p_um: pump wavelength in µm
        poling_period_25C: poling period at 25°C (Λ₀)
        interaction_type: 'type0' = all extraordinary
                          'type2' = pump e, signal o, idler e

    Returns:
        Δk in µm⁻¹
    """

    # --- Eq. (6): Energy conservation → λi
    lambda_i_um = 1 / (1 / lambda_p_um - 1 / lambda_s_um)

    # --- Get indices based on polarization (Sellmeier)
    if interaction_type == 'type0':
        n_p = n_z(lambda_p_um, T)
        n_s = n_z(lambda_s_um, T)
        n_i = n_z(lambda_i_um, T)
    elif interaction_type == 'type2':
        n_p = n_z(lambda_p_um, T)
        n_s = n_y(lambda_s_um, T)
        n_i = n_z(lambda_i_um, T)
    else:
        raise ValueError("interaction_type must be 'type0' or 'type2'")

    # --- Compute wavevectors kj = 2π * n / λ
    k_p = 2 * np.pi * n_p / lambda_p_um
    k_s = 2 * np.pi * n_s / lambda_s_um
    k_i = 2 * np.pi * n_i / lambda_i_um

    # --- QPM grating at temperature T
    Λ_T = poling_period(poling_period_25C, T)
    k_qpm = 2 * np.pi / Λ_T

    return k_p - k_s - k_i - k_qpm


def find_phase_match_lambda_s(T, lambda_p_um, poling_period_25C, interaction_type='type0', search=(0.6, 1.1)):
    """
    Find λs such that Δk(λs, T) ≈ 0.

    Returns:
        λs (µm), or None if root not found
    """
    f = lambda λs: delta_k(λs, T, lambda_p_um, poling_period_25C, interaction_type)
    a, b = search

    # Check if root is bracketed
    if f(a) * f(b) > 0:
        print(f"Warning: No root at T={T}°C in range {a}-{b}")
        return None

    return brentq(f, a, b)
