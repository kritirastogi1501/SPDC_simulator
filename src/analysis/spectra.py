# spdc_emission_spectrum.py
"""
Compute and plot the SPDC emission spectra for both signal and idler using
existing phase-matching routines.

Outputs two curves on the same figure:
 - Signal spectrum: I_s(λ_s) ∝ sinc²[Δk_s(λ_s, T)·L/2]
 - Idler spectrum:  I_i(λ_i) ∝ sinc²[Δk_i(λ_i, T)·L/2]

Usage:
    python spdc_emission_spectrum.py

Dependencies:
    numpy, matplotlib, optics.phase_match (delta_k_s, delta_k_i)
"""
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import matplotlib.pyplot as plt
import math
from optics.phase_match import delta_k_s, delta_k_i

# ----------------------
# Parameters
# ----------------------
L = 25e-3           # Crystal length in meters (25 mm)
T = 132            # Temperature in °C (choose your operating T)

# Wavelength ranges (μm)
lambda_s_um = np.linspace(0.75, 0.9, 500)   # Signal: 750–900 nm
lambda_i_um = np.linspace(0.75, 0.9, 500)   # Idler:  750–900 nm

# ----------------------
# Compute intensities
# ----------------------
# sinc²(x) = (sin(x)/x)²

def sinc_sq(x):
    # handle x=0
    return np.where(x==0, 1.0, (np.sin(x)/x)**2)

I_s = []
for ls in lambda_s_um:
    dk_s = delta_k_s(ls, T)
    I_s.append(sinc_sq(dk_s * L / 2))

I_i = []
for li in lambda_i_um:
    dk_i = delta_k_i(li, T)
    I_i.append(sinc_sq(dk_i * L / 2))

# Convert to numpy arrays
I_s = np.array(I_s)
I_i = np.array(I_i)

# ----------------------
# Plotting
# ----------------------
plt.figure(figsize=(8, 5))
plt.plot(lambda_s_um * 1e3, I_s, label='Signal')
plt.plot(lambda_i_um * 1e3, I_i, label='Idler', linestyle='--')
plt.xticks(np.arange(750, 901, 50))
plt.xlabel('Wavelength (nm)')
plt.ylabel('Normalized Intensity')
plt.title(f'SPDC Emission Spectra at T={T}°C, L={L*1e3:.0f} mm')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
