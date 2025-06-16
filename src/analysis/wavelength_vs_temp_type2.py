# tuning_curves.py

# -------------------------
# Path setup for local module imports
# -------------------------
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq  # Brent’s method for root-finding
from optics.phase_match_type2 import delta_k_s, delta_k_i  # Δk functions from Equation 8

# -------------------------
# Define temperature and wavelength search ranges
# -------------------------

# Temperature range [°C] over which to compute tuning curves
T_vals = np.linspace(50, 500, 10000)  # 141 points from 0°C to 500°C

# Wavelength bracket [nm] for root-finding (where Δk = 0)
λ_brack = (0.7, 1)

# -------------------------
# Root-finding: solve Δk = 0 (Equation 8)
# -------------------------
# For each temperature T:
# - Find λs where Δk_s(λs, T) = 0 → signal phase matching
# - Find λi where Δk_i(λi, T) = 0 → idler phase matching

λs_roots = []  # Signal wavelength roots
λi_roots = []  # Idler wavelength roots

for T in T_vals:
    # Use Brent’s method to find signal wavelength λs such that Δk_s(λs, T) = 0
    λs0 = brentq(lambda ls: delta_k_s(ls, T), *λ_brack)

    # Similarly, find idler wavelength λi such that Δk_i(λi, T) = 0
    λi0 = brentq(lambda li: delta_k_i(li, T), *λ_brack)

    # Store roots
    λs_roots.append(λs0)
    λi_roots.append(λi0)

# -------------------------
# Plotting: Signal and Idler Tuning Curves
# -------------------------

# Plot 1 – Signal wavelength vs temperature
plt.figure()
plt.plot(T_vals, λs_roots)
plt.xlabel('Temperature (°C)')
plt.ylabel('Signal λₛ ')
plt.title('Signal Wavelength vs. Temperature (Δk = 0)')
plt.grid(True)

# Plot 2 – Idler wavelength vs temperature
plt.figure()
plt.plot(T_vals, λi_roots)
plt.xlabel('Temperature (°C)')
plt.ylabel('Idler λᵢ')
plt.title('Idler Wavelength vs. Temperature (Δk = 0)')
plt.grid(True)

# Plot 3 – Combined tuning curves
plt.figure()
plt.plot(T_vals, λs_roots, label='λₛ (signal)')
plt.plot(T_vals, λi_roots, label='λᵢ (idler)')
plt.xlabel('Temperature (°C)')
plt.ylabel('Wavelength')
plt.title('Tuning Curves (Δk = 0)')
plt.legend()
plt.grid(True)

# Display all plots
plt.show()

#94.12