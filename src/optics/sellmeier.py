import math  # provides sqrt and pi

# --- Constants for Modified Sellmeier Equations (from user-specified equations) ---
def n2_y(lambda_um):
    """
    Compute n_y^2 using the modified Sellmeier equation.
    """
    term1 = 3.0333
    term2 = 0.04154 / ((lambda_um)**2 - 0.04547)
    term3 = -0.01408 * (lambda_um)**2
    return term1 + term2 + term3


def n2_z(lambda_um):
    """
    Compute n_z^2 using the modified Sellmeier equation.
    """
    term1 = 3.3134
    term2 = 0.05694 / ((lambda_um)**2 - 0.05658)
    term3 = -0.01682 * (lambda_um)**2
    return term1 + term2 + term3

# Taking square root to get refractive indices (n = sqrt(n^2))

def n_y(lambda_um):
    """
    Compute refractive index n_y by taking sqrt of n2_y.
    """
    return math.sqrt(n2_y(lambda_um))


def n_z(lambda_um):
    """
    Compute refractive index n_z by taking sqrt of n2_z.
    """
    return math.sqrt(n2_z(lambda_um))

# --- Thermo-optic Coefficients (from user-specified equations) ---

def dn_dT_y(lambda_um, T):
    """
    Compute dn_y/dT from user-provided formula.
    """
    L = lambda_um
    return (0.5014 * L**-3 - 2.0030 * L**-2 + 3.3016 * L**-1 + 0.7498) * 1e-5 * (T - 25.0)


def dn_dT_z(lambda_um, T):
    """
    Compute dn_z/dT from user-provided formula.
    """
    L = lambda_um
    return (0.3896 * L**-3 - 1.3332 * L**-2 + 2.2762 * L**-1 + 2.1151) * 1e-5 * (T - 25.0)

# --- Temperature-adjusted refractive indices ---

def n_y_T(lambda_um, T):
    """
    Compute temperature-adjusted refractive index n_y(T).
    """
    return n_y(lambda_um) + dn_dT_y(lambda_um, T)



def n_z_T(lambda_um, T):
    """
    Compute temperature-adjusted refractive index n_z(T).
    """
    return n_z(lambda_um) + dn_dT_z(lambda_um, T)

# --- Poling period stays unchanged ---

# Poling‐period thermal-expansion constants
Λ0 = 10.0         # initial poling period in microns
A_pol = 6.7e-6    # linear thermal expansion coefficient
B_pol = 11e-9     # quadratic thermal expansion coefficient

def poling_period(T):
    """
    Compute the poling period Λ(T).
    """
    return Λ0 * (1 + (A_pol * (T - 25.0)) + (B_pol * (T - 25.0)**2))

# --- Test section ---
if __name__ == "__main__":
    λp_um = 0.405
    λi_um = 0.810
    T_test = 30.0

    print(f"n_y @405 nm: {n_y(λp_um):.6f}")
    print(f"n_z @810 nm: {n_z(λi_um):.6f}")
    print(f"n_y(T=30°C) @405 nm: {n_y_T(λp_um, T_test):.6f}")
    print(f"n_z(T=30°C) @810 nm: {n_z_T(λi_um, T_test):.6f}")
    print(f"Poling period @30°C: {poling_period(T_test):.6f} μm")
