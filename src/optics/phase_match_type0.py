# phase_match.py

# Add custom module path (adjust based on your system)
import sys
sys.path.append("C:/Users/Asus/OneDrive/Desktop/Sanya Personal/Projects/DRDO internship/spdc_simulator/src")

import math
from optics.sellmeier import n_z_T, poling_period

# ----------------------------
# Constants and conversions
# ----------------------------

# Pump wavelength fixed at 405 nm
λp_um = 0.405


# ----------------------------
# Eq. 9 – Wavelength relation:
# 1/λp = 1/λs + 1/λi
# ----------------------------

def λi_from_λs(λs_um):
    """
    Calculate idler wavelength λi [nm] from signal wavelength λs [nm] using Equation 9.
    """
    return 1.0 / (1.0/λp_um - 1.0/λs_um)

def λs_from_λi(λi_um):
    """
    Calculate signal wavelength λs [nm] from idler wavelength λi [nm] using Equation 9.
    """
    return 1.0 / (1.0/λp_um - 1.0/λi_um)

# ----------------------------
# Eq. 8 – Δk expression:
# Δk = (2π·n_y(λp)/λp) - (2π·n_z(λs)/λs) - (2π·n_y(λi)/λi) - (2π/Λ(T))
#
# Notes:
# - n_y(λp) and n_y(λi) from Equation 1 → temp-adjusted using Equation 5
# - n_z(λs) from Equation 2 → temp-adjusted using Equation 6
# - Λ(T) from Equation 7
# ----------------------------

def delta_k_s(λs_um, T):
    """
    Compute phase mismatch Δk [1/m] as a function of signal wavelength λs [nm] and temperature T [°C].
    λi is computed via Equation 9.
    """
    conv = 1e-6  # Convert μm to meters
    λi_um = λi_from_λs(λs_um)

    # Each term corresponds to one term in Equation 8
    kp = 2 * math.pi * n_z_T(λp_um, T) / (λp_um * conv)   
    ks = 2 * math.pi * n_z_T(λs_um, T) / (λs_um * conv)   
    ki = 2 * math.pi * n_z_T(λi_um, T) / (λi_um * conv)   
    kG = 2 * math.pi / (poling_period(T) * conv)          

    return kp - ks - ki - kG

def delta_k_i(λi_um, T):
    """
    Compute phase mismatch Δk [1/m] as a function of idler wavelength λi [nm] and temperature T [°C].
    λs is computed via Equation 9.
    """
    λs_um = λs_from_λi(λi_um)
    conv = 1e-6  # Convert μm to meters

    # Terms as per Equation 8
    kp = 2 * math.pi * n_z_T(λp_um, T) / (λp_um * conv)   # n_y(λp)
    ks = 2 * math.pi * n_z_T(λs_um, T) / (λs_um * conv)   # n_z(λs)
    ki = 2 * math.pi * n_z_T(λi_um, T) / (λi_um * conv)   # n_y(λi)
    kG = 2 * math.pi / (poling_period(T) * conv)          # Λ(T)

    return kp - ks - ki - kG


def delta_k(λi_um, λs_um, T):
    conv = 1e-6  # Convert μm to meters
    kp = 2 * math.pi * n_z_T(λp_um, T) / (λp_um * conv)   # n_y(λp)
    ks = 2 * math.pi * n_z_T(λs_um, T) / (λs_um * conv)   # n_z(λs)
    ki = 2 * math.pi * n_z_T(λi_um, T) / (λi_um * conv)   # n_y(λi)
    kG = 2 * math.pi / (poling_period(T) * conv)          # Λ(T)

    return kp - ks - ki - kG
