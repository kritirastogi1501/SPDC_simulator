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
from optics.phase_match_type0 import delta_k_s, delta_k_i  # Δk functions from Equation 8

# -------------------------
# Define temperature and wavelength search ranges
# -------------------------

# Temperature range [°C] over which to compute tuning curves
T_vals = np.linspace(0, 50, 50)  # 141 points from 0°C to 500°C

# Wavelength bracket [nm] for root-finding (where Δk = 0)
λ_brack = (0.7, 0.9)

# -------------------------
# Root-finding: solve Δk = 0 (Equation 8)
# -------------------------
# For each temperature T:
# - Find λs where Δk_s(λs, T) = 0 → signal phase matching
# - Find λi where Δk_i(λi, T) = 0 → idler phase matching

λs_roots = []  # Signal wavelength roots
λi_roots = []  # Idler wavelength roots

# Debugging: Check function behavior at the bracket bounds
for T in T_vals:
    f_a = delta_k_s(λ_brack[0], T)
    f_b = delta_k_s(λ_brack[1], T)

    print(f"Temperature {T}: f({λ_brack[0]}) = {f_a}, f({λ_brack[1]}) = {f_b}")

    if f_a * f_b > 0:
        print(f"Warning: No sign change for T = {T}, skipping...")
        continue  # Skip if no sign change

    # Perform root-finding only if sign change is detected
    λs0 = brentq(lambda ls: delta_k_s(ls, T), *λ_brack)
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
plt.ylabel('Signal λₛ (nm)')
plt.title('Signal Wavelength vs. Temperature (Δk = 0)')
plt.grid(True)

# Plot 2 – Idler wavelength vs temperature
plt.figure()
plt.plot(T_vals, λi_roots)
plt.xlabel('Temperature (°C)')
plt.ylabel('Idler λᵢ (nm)')
plt.title('Idler Wavelength vs. Temperature (Δk = 0)')
plt.grid(True)

# Plot 3 – Combined tuning curves
plt.figure()
plt.plot(T_vals, λs_roots, label='λₛ (signal)')
plt.plot(T_vals, λi_roots, label='λᵢ (idler)')
plt.xlabel('Temperature (°C)')
plt.ylabel('Wavelength (nm)')
plt.title('Tuning Curves (Δk = 0)')
plt.legend()
plt.grid(True)

# Display all plots
plt.show()

