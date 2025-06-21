# -*- coding: utf-8 -*-
"""
Updated on: 2025-06-20
Updated by: Sanya Garg

This script computes and plots the wavelength vs. temperature tuning curves
for Type-II SPDC in PPKTP (405 nm pump → ~810 nm signal/idler) using:
  - Sellmeier equations for KTP
  - Thermo-optic corrections (Emanueli et al. form)
  - Thermal expansion of the poling period

All units: wavelength in micrometers, temperature in °C.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#------------------------ Sellmeier coefficients for KTP ------------------------
# From Hamel thesis Appendix B at 25°C
Y0 = [2.19229,  0.83547,  0.04970, 0.01621]  
Z0 = [2.12725, 1.18431, 0.0514852, 0.6603, 100.00507, 9.68956e-3]  

#--------------------------- Sellmeier equations --------------------------------
def nY(lam):
    return np.sqrt(
        Y0[0]
        + (Y0[1] / (lam**2 - Y0[2]) - Y0[3]) * lam**2
    )

def nZ(lam):
    return np.sqrt(
        Z0[0]
        + ((Z0[1] / (lam**2 - Z0[2])
            + Z0[3] / (lam**2 - Z0[4])
            - Z0[5]
           ) * lam**2)
    )

#----------------------- Thermo-optic corrections --------------------------------
# coefficients a_m for n1 and n2, y‑ and z‑polarizations
aY1 = [6.2897e-6, 6.3061e-6, -6.0629e-6, 2.6486e-6]
aY2 = [-1.4445e-9, 2.2244e-8, -3.5770e-8, 1.3470e-8]

aZ1 = [9.9587e-6, 9.9228e-6, -8.9603e-6, 4.1010e-6]
aZ2 = [-1.1882e-8, 1.0459e-7, -9.8136e-8, 3.1481e-8]

def poly_inv(lam, coeffs):
    """Σ_{m=0}^3 coeffs[m] / lam**m"""
    return sum(c / lam**m for m, c in enumerate(coeffs))

def delta_nY(lam, T):
    dT = T - 25.0
    n1 = poly_inv(lam, aY1)
    n2 = poly_inv(lam, aY2)
    return n1 * dT + n2 * dT**2

def delta_nZ(lam, T):
    dT = T - 25.0
    n1 = poly_inv(lam, aZ1)
    n2 = poly_inv(lam, aZ2)
    return n1 * dT + n2 * dT**2

def nY_T(lam, T):
    """y‑polarized refractive index at temperature T"""
    return nY(lam) + delta_nY(lam, T)

def nZ_T(lam, T):
    """z‑polarized refractive index at temperature T"""
    return nZ(lam) + delta_nZ(lam, T)

#----------------------- Poling period expansion --------------------------------
def Lambda_QPM(T, Λ0=9.925):
    alpha = 6.7e-6   # /°C
    beta  = 1.1e-8   # /°C^2
    dT = T - 25.0
    return Λ0 * (1 + alpha * dT + beta * dT**2)

#------------------ Energy conservation helper ----------------------------------
lambda_p = 0.405  # pump in micrometers

def lambda_signal(lambda_idler):
    return 1 / (1/lambda_p - 1/lambda_idler)

#------------------ Phase-matching equation solver -------------------------------
def phase_match_eq(lambda_idler, T):
    lam_s = lambda_signal(lambda_idler)
    return (
        nY_T(lambda_p, T)/lambda_p
        - nZ_T(lam_s,   T)/lam_s
        - nY_T(lambda_idler, T)/lambda_idler
        - 1/Lambda_QPM(T)
    )

def find_lambda_idler(T):
    sol = fsolve(phase_match_eq, x0=1.5, args=(T,), xtol=1e-6)
    return sol[0]

#--------------------------- Compute and Plot ------------------------------------
if __name__ == "__main__":
    T_vals    = np.linspace(0, 120, 400)
    idler_vals  = [find_lambda_idler(T) for T in T_vals]
    signal_vals = [lambda_signal(li)  for li in idler_vals]

    plt.plot(T_vals, signal_vals, label='Signal (λ_s)')
    plt.plot(T_vals, idler_vals,  label='Idler (λ_i)')
    plt.xlabel('Temperature (°C)')
    plt.ylabel('Wavelength (µm)')
    plt.title('Type-II SPDC Wavelength vs Temperature for PPKTP')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def compute_typeii_wavelengths(T_vals, lambda_p=0.405):
    idler_vals = [find_lambda_idler(T) for T in T_vals]
    signal_vals = [lambda_signal(li) for li in idler_vals]
    return signal_vals, idler_vals


def find_signal_idler_at_temp(T, lambda_p=0.405):
    """Compute signal and idler wavelengths at a specific temperature."""
    lambda_idler = find_lambda_idler(T)
    lambda_signal_val = lambda_signal(lambda_idler)
    return lambda_signal_val, lambda_idler

def find_degenerate_temperature(lambda_p=0.405):
    """Find the temperature where signal ≈ idler (degeneracy)."""
    def degenerate_condition(T):
        lam_i = find_lambda_idler(T)
        lam_s = lambda_signal(lam_i)
        return lam_s - lam_i

    T_guess = 50  # Good starting point for 405 nm pump
    T_deg = fsolve(degenerate_condition, x0=T_guess)[0]
    return T_deg
