# -*- coding: utf-8 -*-
"""
Updated on: 2025-06-20
Updated by: Sanya Garg

This script computes and plots the wavelength vs. temperature tuning curves
for Type-II SPDC in PPKTP (405 nm pump → ~810 nm signal/idler) using:
  - Sellmeier equations for KTP
  - Thermo-optic corrections
  - Thermal expansion of the poling period

All units: wavelength in micrometers, temperature in °C.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#------------------------ Sellmeier coefficients for KTP ------------------------
# From Hamel thesis Appendix B at 25°C
Y0 = [2.09930, 0.922683, 0.0467695, 0.0138404]  # Eq1 for n_y
Z0 = [2.12725, 1.18431, 0.0514852, 0.6603, 100.00507, 9.68956e-3]  # Eq2 for n_z

#--------------------------- Sellmeier equations --------------------------------
def nY(lam):
    return np.sqrt(Y0[0] + (Y0[1] / (lam**2 - Y0[2]) - Y0[3]) * lam**2)

def nZ(lam):
    return np.sqrt(Z0[0] + ((Z0[1] / (lam**2 - Z0[2]) + Z0[3] / (lam**2 - Z0[4]) - Z0[5]) * lam**2))

#----------------------- Thermo-optic corrections --------------------------------
def delta_nY(lam, T):
    dT = T - 25.0
    return (
        6.2897e-6 * dT - 1.4445e-9 * dT**2 +
        ((6.3061e-6 + 2.2244e-8 * dT) * dT) / lam +
        ((6.0629e-6 + 3.5770e-8 * dT) * dT) / lam**2 +
        ((2.6486e-6 + 1.3470e-8 * dT) * dT) / lam**3
    )

def delta_nZ(lam, T):
    dT = T - 25.0
    return (
        9.9587e-6 * dT - 1.1882e-8 * dT**2 +
        ((9.9228e-6 + 1.0459e-7 * dT) * dT) / lam +
        ((8.9603e-6 + 9.8136e-8 * dT) * dT) / lam**2 +
        ((4.1010e-6 + 3.1481e-8 * dT) * dT) / lam**3
    )

def nY_T(lam, T):
    return nY(lam) + delta_nY(lam, T)

def nZ_T(lam, T):
    return nZ(lam) + delta_nZ(lam, T)

#----------------------- Poling period expansion --------------------------------
def Lambda_QPM(T, Λ0=9.925):
    alpha = 6.7e-6  # /°C
    beta  = 1.1e-8  # /°C^2
    return Λ0 * (1 + alpha * (T - 25) + beta * (T - 25)**2)

#------------------ Energy conservation helper ----------------------------------
lambda_p = 0.405  # pump in micrometers

def lambda_signal(lambda_idler):
    return 1 / (1 / lambda_p - 1 / lambda_idler)

#------------------ Phase-matching equation solver --------------------------------
def phase_match_eq(lambda_idler, T):
    lam_s = lambda_signal(lambda_idler)
    return (
        nY_T(lambda_p, T) / lambda_p
        - nZ_T(lam_s, T) / lam_s
        - nY_T(lambda_idler, T) / lambda_idler
        - 1 / Lambda_QPM(T)
    )

def find_lambda_idler(T):
    sol = fsolve(phase_match_eq, x0=0.81, args=(T,), xtol=1e-6)
    return sol[0]

#--------------------------- Compute and Plot ------------------------------------
T_vals = np.linspace(0, 120, 400)
idler_vals = [find_lambda_idler(T) for T in T_vals]
signal_vals = [lambda_signal(li) for li in idler_vals]

plt.plot(T_vals, signal_vals, label='Signal (λ_s)')
plt.plot(T_vals, idler_vals, label='Idler (λ_i)')
plt.xlabel('Temperature (°C)')
plt.ylabel('Wavelength (µm)')
plt.title('Type-II SPDC Wavelength vs Temperature for PPKTP')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
