# tuning_curves.py

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from optics.phase_match import delta_k_s, delta_k_i

# Temperature grid and search bracket
T_vals  = np.linspace(0, 70, 141)      # 0–70 °C
λ_brack = (700.0, 900.0)               # nm

# Solve roots Δk_s=0 and Δk_i=0 for each T
λs_roots = []
λi_roots = []
for T in T_vals:
    # find λ_s such that Δk_s(λ_s, T)=0
    λs0 = brentq(lambda ls: delta_k_s(ls, T), *λ_brack)
    # find λ_i such that Δk_i(λ_i, T)=0
    λi0 = brentq(lambda li: delta_k_i(li, T), *λ_brack)
    λs_roots.append(λs0)
    λi_roots.append(λi0)

# Plot 1: λ_s vs T
plt.figure()
plt.plot(T_vals, λs_roots)
plt.xlabel('Temperature (°C)')
plt.ylabel('Signal λₛ (nm)')
plt.title('Signal Wavelength vs. Temperature')
plt.grid(True)

# Plot 2: λ_i vs T
plt.figure()
plt.plot(T_vals, λi_roots)
plt.xlabel('Temperature (°C)')
plt.ylabel('Idler λᵢ (nm)')
plt.title('Idler Wavelength vs. Temperature')
plt.grid(True)

# Plot 3: Both curves
plt.figure()
plt.plot(T_vals, λs_roots, label='λₛ (signal)')
plt.plot(T_vals, λi_roots, label='λᵢ (idler)')
plt.xlabel('Temperature (°C)')
plt.ylabel('Wavelength (nm)')
plt.title('Tuning Curves (Δk=0)')
plt.legend()
plt.grid(True)

plt.show()
