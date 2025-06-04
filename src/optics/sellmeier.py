import math

# --- Coefficients for Sellmeier Equations ---
A_y, B_y, C_y, D_y = 2.19229, 0.83547, 0.04970, 0.01621
A_z, B_z, C_z, D_z = 2.25411, 1.06543, 0.05486, 0.02140

# Grating‐period constants
Λ0 = 10.0         # initial poling period in microns
A_grate = 6.7e-6
B_grate = 11e-9

# Convert nm to microns
def nm_to_um(lam_nm):
    return lam_nm * 1e-3

# Eq.1 & 2: refractive index squared
def n2_y(λ_um):
    return A_y + (B_y * λ_um**2) / (λ_um**2 - C_y) - D_y * λ_um**2

def n2_z(λ_um):
    return A_z + (B_z * λ_um**2) / (λ_um**2 - C_z) - D_z * λ_um**2

# sqrt → refractive index
def n_y(λ_um):
    return math.sqrt(n2_y(λ_um))

def n_z(λ_um):
    return math.sqrt(n2_z(λ_um))

# Eq.3 & 4: temperature derivatives
def dn_dT_y(λ_um):
    return (1.997/λ_um**3 - 4.067/λ_um**2 + 5.154/λ_um - 5.425) * 1e-6

def dn_dT_z(λ_um):
    return (9.221/λ_um**3 - 29.220/λ_um**2 + 36.667/λ_um - 1.897) * 1e-6

# Eq.5 & 6: temperature-adjusted refractive indices
def n_y_T(λ_um, T):
    return n_y(λ_um) + dn_dT_y(λ_um) * (T - 25.0)

def n_z_T(λ_um, T):
    return n_z(λ_um) + dn_dT_z(λ_um) * (T - 25.0)

# Eq.7: temperature-dependent grating period
def grating_period(T):
    return Λ0 * (1 + A_grate*(T - 25.0) + B_grate*(T - 25.0)**2)

# Sanity-check print when run directly
if __name__ == "__main__":
    λp_um = nm_to_um(405)
    λi_um = nm_to_um(810)
    print(f"n_y @405 nm: {n_y(λp_um):.6f}")
    print(f"n_z @810 nm: {n_z(λi_um):.6f}")
    print(f"Grating period @60°C: {grating_period(60):.6f} μm")
