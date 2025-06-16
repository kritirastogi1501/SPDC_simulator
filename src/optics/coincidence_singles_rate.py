# coincidence_single_rate.py
# Calculates R_pair_filtered, R_coincidence, and R_single rates

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from scipy.integrate import quad
from optics.phase_match_type2 import delta_k_s

# -------------------------------
# Constants and Parameters
# -------------------------------
c = 3e8  # Speed of light (m/s)
L = 25e-3  # Crystal length in meters
T = 75  # Temperature in Celsius

# Detection/collection efficiencies
eta_s = 0.6  # Signal arm efficiency
eta_i = 0.6  # Idler arm efficiency

# -------------------------------
# sinc^2 helper function
# -------------------------------
def sinc_squared(x):
    return np.sinc(x / np.pi)**2

# -------------------------------
# Compute R_pair_filtered
# -------------------------------
def R_pair_filtered(λ0, Δλ, T=94.2, L=0.025, resolution=500):
    """
    Computes the normalized filtered pair rate across the filter bandwidth centered at λ0.

    λ0: Central signal wavelength (μm)
    Δλ: Filter bandwidth (μm)
    T: Temperature (°C)
    L: Crystal length (m)
    resolution: Number of points to integrate over
    """

    λs_vals = np.linspace(λ0 - Δλ/2, λ0 + Δλ/2, resolution)
    dλ = λs_vals[1] - λs_vals[0]

    integrand_vals = []

    for λs in λs_vals:
        Δk = delta_k_s(λs, T)  # in 1/m
        x = Δk * L / 2
        sinc_squared = np.sinc(x / np.pi)**2
        integrand_vals.append(sinc_squared)

        if abs(λs - λ0) < 0.0001:  # around center
            print(f"λs = {λs:.6f} μm, Δk = {Δk:.3e}, sinc² = {sinc_squared:.3e}")

    R_pair = np.trapz(integrand_vals, λs_vals)
    return R_pair


# -------------------------------
# Coincidence rate
# -------------------------------
def R_coincidence(lambda0_um, delta_lambda_um):
    R_pair = R_pair_filtered(lambda0_um, delta_lambda_um)
    return eta_s * eta_i * R_pair

# -------------------------------
# Singles rate
# -------------------------------
def R_singles(lambda0_um, delta_lambda_um):
    R_pair = R_pair_filtered(lambda0_um, delta_lambda_um)
    R_s = eta_s * R_pair
    R_i = eta_i * R_pair
    return R_s, R_i

# -------------------------------
# Example/test usage (can be deleted/commented out)
# -------------------------------
if __name__ == "__main__":
    lambda0 = 0.81  # Center signal wavelength in um
    delta_lambda = 0.01  # Filter bandwidth in um (10 nm)

    Rpair = R_pair_filtered(lambda0, delta_lambda)
    Rcoinc = R_coincidence(lambda0, delta_lambda)
    Rsingle_s, Rsingle_i = R_singles(lambda0, delta_lambda)

    print(f"R_pair_filtered = {Rpair:.5f} (normalized)")
    print(f"R_coincidence = {Rcoinc:.5f}")
    print(f"R_single_s = {Rsingle_s:.5f}, R_single_i = {Rsingle_i:.5f}")
