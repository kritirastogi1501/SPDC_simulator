import math  # provides sqrt and pi

# --- Coefficients for Sellmeier Equations ---
# These constants define the material dispersion along the Y and Z crystal axes
A_y, B_y, C_y, D_y = 2.19229, 0.83547, 0.04970, 0.01621
A_z, B_z, C_z, D_z = 2.25411, 1.06543, 0.05486, 0.02140

# Poling‐period thermal-expansion constants
# Λ0: initial period at 25 °C, A_pol/B_pol: linear/quadratic expansion coefficients
Λ0 = 9.925         # initial poling period in microns
A_pol = 6.7e-6    # linear thermal expansion coefficient
B_pol = 11e-9     # quadratic thermal expansion coefficient

# Utility conversion

def nm_to_um(lam_nm):
    """
    Convert wavelength from nanometers to microns.
    Parameters:
        lam_nm (float): Wavelength in nanometers.
    Returns:
        float: Wavelength in microns.
    """
    return lam_nm * 1e-3  # 1 nm = 1e-3 μm

# Sellmeier refractive-index equations (Eqs. 1 & 2)

def n2_y(lambda_um):
    """
    Compute n_y^2 according to Sellmeier equation for Y-axis.
    n_y^2 = A_y + (B_y * λ^2)/(λ^2 - C_y) - D_y * λ^2
    Parameters:
        lambda_um (float): Wavelength in microns.
    Returns:
        float: Square of refractive index along Y-axis.
    """
    return A_y + (B_y * lambda_um**2) / (lambda_um**2 - C_y) - D_y * lambda_um**2


def n2_z(lambda_um):
    """
    Compute n_z^2 according to Sellmeier equation for Z-axis.
    n_z^2 = A_z + (B_z * λ^2)/(λ^2 - C_z) - D_z * λ^2
    Parameters:
        lambda_um (float): Wavelength in microns.
    Returns:
        float: Square of refractive index along Z-axis.
    """
    return A_z + (B_z * lambda_um**2) / (lambda_um**2 - C_z) - D_z * lambda_um**2

# Taking square root to get refractive indices (n = sqrt(n^2))

def n_y(lambda_um):
    """
    Compute refractive index n_y by taking sqrt of n2_y.
    Parameters:
        lambda_um (float): Wavelength in microns.
    Returns:
        float: Refractive index along Y-axis.
    """
    return math.sqrt(n2_y(lambda_um))


def n_z(lambda_um):
    """
    Compute refractive index n_z by taking sqrt of n2_z.
    Parameters:
        lambda_um (float): Wavelength in microns.
    Returns:
        float: Refractive index along Z-axis.
    """
    return math.sqrt(n2_z(lambda_um))

# Thermo-optic coefficients (Eqs. 3 & 4)

def dn_dT_y(lambda_um):
    """
    Temperature derivative of n_y: dn_y/dT at a given wavelength.
    Uses polynomial fit: (1.997/λ^3 - 4.067/λ^2 + 5.154/λ - 5.425) × 10^-6.
    Parameters:
        lambda_um (float): Wavelength in microns.
    Returns:
        float: dn_y/dT (per °C).
    """
    return (1.997/lambda_um**3 - 4.067/lambda_um**2 + 5.154/lambda_um - 5.425) * 1e-6


def dn_dT_z(lambda_um):
    """
    Temperature derivative of n_z: dn_z/dT at a given wavelength.
    Uses polynomial fit: (9.221/λ^3 - 29.220/λ^2 + 36.667/λ - 1.897) × 10^-6.
    Parameters:
        lambda_um (float): Wavelength in microns.
    Returns:
        float: dn_z/dT (per °C).
    """
    return (9.221/lambda_um**3 - 29.220/lambda_um**2 + 36.667/lambda_um - 1.897) * 1e-6

# Temperature-adjusted refractive indices (Eqs. 5 & 6)

def n_y_T(lambda_um, T):
    """
    Compute temperature-adjusted refractive index n_y(T).
    n_y(T) = n_y(λ) + (dn_y/dT)(λ) × (T - 25°C).
    Parameters:
        lambda_um (float): Wavelength in microns.
        T (float): Temperature in °C.
    Returns:
        float: Temperature-adjusted n_y.
    """
    return n_y(lambda_um) + dn_dT_y(lambda_um) * (T - 25.0)


def n_z_T(lambda_um, T):
    """
    Compute temperature-adjusted refractive index n_z(T).
    n_z(T) = n_z(λ) + (dn_z/dT)(λ) × (T - 25°C).
    Parameters:
        lambda_um (float): Wavelength in microns.
        T (float): Temperature in °C.
    Returns:
        float: Temperature-adjusted n_z.
    """
    return n_z(lambda_um) + dn_dT_z(lambda_um) * (T - 25.0)

# Temperature-dependent poling period (Eq. 7)

def poling_period(T):
    """
    Compute the poling period Λ(T).
    Λ(T) = Λ0 × [1 + A_pol(T - 25) + B_pol(T - 25)^2]
    Parameters:
        T (float): Temperature in °C.
    Returns:
        float: Poling period in microns at temperature T.
    """
    return Λ0 * (1 + A_pol * (T - 25.0) + B_pol * (T - 25.0)**2)

# test
if __name__ == "__main__":
    # Convert pump and idler wavelengths to microns
    λp_um = nm_to_um(405)
    λi_um = nm_to_um(810)
    # Print baseline refractive indices and poling period
    print(f"n_y @405 nm: {n_y(λp_um):.6f}")
    print(f"n_z @810 nm: {n_z(λi_um):.6f}")
    print(f"Poling period @30°C: {poling_period(30):.6f} μm")