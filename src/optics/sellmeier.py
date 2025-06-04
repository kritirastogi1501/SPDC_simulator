# src/optics/sellmeier.py

import numpy as np

# -----------------------------------------------------------------------------
# Sellmeier Model for PPKTP with Thermal Corrections
# Equations from ISRO paper
# Implements:
#   - Eq. (1): Sellmeier dispersion relation with temperature-independent terms
#   - Eq. (2) and (3): dn/dT thermo-optic relations (λ-dependent)
#   - Eq. (4): Temperature-dependent thermal expansion of the poling period
# -----------------------------------------------------------------------------

# === Sellmeier Coefficients from Table I ===
# Note: λ in µm
# For y-axis (ordinary ray)
A_y, B_y, C_y, D_y = 2.19229, 0.83547, 0.04970, 0.01621
# For z-axis (extraordinary ray)
A_z, B_z, C_z, D_z = 2.25411, 1.06543, 0.05486, 0.02140

# === Eq. (2): dn_y/dT — Temperature-dependent correction for n_y ===
def dn_dT_y(lambda_um):
    l = lambda_um
    return (1.997 / l**3 - 4.067 / l**2 + 5.154 / l - 5.425) * 1e-6

# === Eq. (3): dn_z/dT — Temperature-dependent correction for n_z ===
def dn_dT_z(lambda_um):
    l = lambda_um
    return (9.221 / l**3 - 29.220 / l**2 + 36.667 / l - 1.897) * 1e-6

# === Eq. (1) + (2): Full temperature-dependent Sellmeier for n_y(λ, T) ===
def n_y(lambda_um, T=25.0):
    """
    Ordinary refractive index n_y(λ,T) using Sellmeier + thermo-optic correction.
    Implements Eq. (1) and Eq. (2).
    """
    l2 = lambda_um**2
    n2 = A_y + (B_y*l2 )/ (l2 - C_y) - D_y * l2
    n0 = np.sqrt(n2)
    return n0 + dn_dT_y(lambda_um) * (T - 25.0)

# === Eq. (1) + (3): Full temperature-dependent Sellmeier for n_z(λ, T) ===
def n_z(lambda_um, T=25.0):
    """
    Extraordinary refractive index n_z(λ,T) using Sellmeier + thermo-optic correction.
    Implements Eq. (1) and Eq. (3).
    """
    l2 = lambda_um**2
    n2 = A_z + (B_z*l2)/ (l2 - C_z) - D_z * l2
    n0 = np.sqrt(n2)
    return n0 + dn_dT_z(lambda_um) * (T - 25.0)

# === Eq. (4): Poling period thermal expansion ===
def poling_period(Λ0_um, T, A_coeff=6.7e-6, B_coeff=11e-9):
    """
    Computes Λ(T) = Λ₀ * [1 + A(T - 25°C) + B(T - 25°C)^2]
    Implements Eq. (4).
    
    Parameters:
        Λ0_um: nominal poling period at 25°C (e.g. 3.42 or 10 µm)
        T: current temperature in °C
        A_coeff, B_coeff: expansion coefficients
    """
    ΔT = T - 25.0
    return Λ0_um * (1 + A_coeff * ΔT + B_coeff * ΔT**2)

# === Self-test and validation ===
if __name__ == "__main__":
    print("λ (µm)\tT (°C)\tn_y\t\tn_z")
    for lam in (0.405, 0.810):
        for T in (25, 50, 75):
            ny = n_y(lam, T)
            nz = n_z(lam, T)
            print(f"{lam:.3f}\t{T}\t{ny:.6f}\t{nz:.6f}")
    print("\nPoling period @ 60°C for Λ₀=3.42 µm:", poling_period(3.42, 60.0))