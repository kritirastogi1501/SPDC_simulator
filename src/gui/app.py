import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import io

# ----------------- Sellmeier & thermo-optic coefficients -----------------
# Sellmeier coefficients for refractive index of KTP crystal (Y and Z axes)
Y0 = [2.19229, 0.83547, 0.04970, 0.01621]
Z0 = [2.12725, 1.18431, 0.0514852, 0.6603, 100.00507, 9.68956e-3]

# Thermo-optic coefficients for polynomial temperature dependence corrections
aY1 = [6.2897e-6, 6.3061e-6, -6.0629e-6, 2.6486e-6]
aY2 = [-1.4445e-9, 2.2244e-8, -3.5770e-8, 1.3470e-8]
aZ1 = [9.9587e-6, 9.9228e-6, -8.9603e-6, 4.1010e-6]
aZ2 = [-1.1882e-8, 1.0459e-7, -9.8136e-8, 3.1481e-8]

# ---------------- Helper function for polynomial thermo-optic terms -------------
def poly_inv(lam, coeffs):
    """
    Computes the sum of coefficients divided by increasing powers of wavelength.
    This is used to calculate thermo-optic corrections.

    Args:
        lam (float): Wavelength in micrometers.
        coeffs (list): List of polynomial coefficients.

    Returns:
        float: Result of the polynomial sum.
    """
    return sum(c / lam**i for i, c in enumerate(coeffs))

# ---------------- Refractive index functions for Y and Z polarizations ---------
def nY(lam):
    """
    Sellmeier equation for Y-polarized refractive index at 25°C.

    Args:
        lam (float): Wavelength in micrometers.

    Returns:
        float: Refractive index nY.
    """
    return np.sqrt(Y0[0] + (Y0[1] / (lam**2 - Y0[2]) - Y0[3]) * lam**2)

def nZ(lam):
    """
    Sellmeier equation for Z-polarized refractive index at 25°C.

    Args:
        lam (float): Wavelength in micrometers.

    Returns:
        float: Refractive index nZ.
    """
    return np.sqrt(Z0[0] + ((Z0[1] / (lam**2 - Z0[2]) + Z0[3] / (lam**2 - Z0[4]) - Z0[5]) * lam**2))

# ------------- Thermo-optic corrections to refractive indices -----------------
def delta_nY(lam, T):
    """
    Temperature-dependent change in nY refractive index.

    Args:
        lam (float): Wavelength in micrometers.
        T (float): Temperature in °C.

    Returns:
        float: Change in refractive index due to temperature.
    """
    dT = T - 25.0  # reference temperature is 25°C
    return poly_inv(lam, aY1) * dT + poly_inv(lam, aY2) * dT**2

def delta_nZ(lam, T):
    """
    Temperature-dependent change in nZ refractive index.

    Args:
        lam (float): Wavelength in micrometers.
        T (float): Temperature in °C.

    Returns:
        float: Change in refractive index due to temperature.
    """
    dT = T - 25.0
    return poly_inv(lam, aZ1) * dT + poly_inv(lam, aZ2) * dT**2

def nY_T(lam, T):
    """
    Total refractive index nY at temperature T.

    Args:
        lam (float): Wavelength in micrometers.
        T (float): Temperature in °C.

    Returns:
        float: Temperature-corrected refractive index.
    """
    return nY(lam) + delta_nY(lam, T)

def nZ_T(lam, T):
    """
    Total refractive index nZ at temperature T.

    Args:
        lam (float): Wavelength in micrometers.
        T (float): Temperature in °C.

    Returns:
        float: Temperature-corrected refractive index.
    """
    return nZ(lam) + delta_nZ(lam, T)

# ------------- Poling period expansion with temperature -----------------------
def Lambda_QPM(T, L0):
    """
    Calculates the temperature-dependent poling period due to thermal expansion.

    Args:
        T (float): Temperature in °C.
        L0 (float): Poling period at reference temperature (25°C) in micrometers.

    Returns:
        float: Adjusted poling period at temperature T.
    """
    dT = T - 25.0
    alpha, beta = 6.7e-6, 1.1e-8  # linear and quadratic thermal expansion coefficients
    return L0 * (1 + alpha * dT + beta * dT**2)

# ----------------- Energy conservation for signal wavelength ------------------
def lambda_signal(lambda_idler, lambda_p):
    """
    Calculates signal wavelength from idler and pump wavelengths based on energy conservation.

    Args:
        lambda_idler (float): Idler wavelength in micrometers.
        lambda_p (float): Pump wavelength in micrometers.

    Returns:
        float: Signal wavelength in micrometers.
    """
    return 1 / (1/lambda_p - 1/lambda_idler)

# ---------------- Phase-matching equation for root finding ---------------------
def phase_match_eq(lambda_idler, T, lambda_p, L0):
    """
    Phase matching condition for Type-II SPDC, solved to find idler wavelength.

    Args:
        lambda_idler (float): Trial idler wavelength in micrometers.
        T (float): Temperature in °C.
        lambda_p (float): Pump wavelength in micrometers.
        L0 (float): Poling period at 25°C in micrometers.

    Returns:
        float: Difference from perfect phase matching (should be zero at solution).
    """
    lam_s = lambda_signal(lambda_idler, lambda_p)  # calculate signal wavelength from idler
    # phase mismatch calculation
    return (
        nY_T(lambda_p, T)/lambda_p
        - nZ_T(lam_s, T)/lam_s
        - nY_T(lambda_idler, T)/lambda_idler
        - 1/Lambda_QPM(T, L0)
    )

# -------------- Find idler wavelength by solving phase-match eq ----------------
def find_lambda_idler(T, lambda_p, L0):
    """
    Numerically finds the idler wavelength for a given temperature and pump wavelength
    by solving the phase-matching equation.

    Args:
        T (float): Temperature in °C.
        lambda_p (float): Pump wavelength in micrometers.
        L0 (float): Poling period at 25°C in micrometers.

    Returns:
        float: Idler wavelength in micrometers.
    """
    # Initial guess for root: twice the pump wavelength (arbitrary but reasonable)
    sol = fsolve(phase_match_eq, x0=2*lambda_p, args=(T, lambda_p, L0), xtol=1e-6)
    return sol[0]

# ----------------------- Streamlit User Interface ------------------------------
st.set_page_config(page_title="SPDC Tuning Curves", layout="wide")
st.title("Type-II SPDC Wavelength vs. Temperature")

# 1. Display settings sidebar
st.sidebar.header("Display Settings")
decimals = st.sidebar.slider("Decimal places:", min_value=0, max_value=10, value=4)

# 2. Pump, poling period and temperature sweep inputs
st.sidebar.header("Pump, Poling & Temp Sweep")
lambda_p = st.sidebar.number_input(
    "Pump wavelength (µm):", min_value=0.1, max_value=2.0,
    value=0.405, format=f"%.{decimals}f"
)
L0 = st.sidebar.number_input(
    "Poling period L0 (µm):", min_value=1.0, max_value=20.0,
    value=9.925, format=f"%.{decimals}f"
)
T_min = st.sidebar.number_input("Min temperature (°C):", value=0)
T_max = st.sidebar.number_input("Max temperature (°C):", value=120)
points = st.sidebar.slider("Resolution (# points):", min_value=10, max_value=1000, value=400)

# Validate temperature range order
if T_max <= T_min:
    st.sidebar.error("Max temperature must be greater than min temperature.")

# 3. Compute wavelength tuning curves over temperature range
temps = np.linspace(T_min, T_max, points)
idlers = [find_lambda_idler(T, lambda_p, L0) for T in temps]      # idler wavelengths at each temperature
signals = [lambda_signal(li, lambda_p) for li in idlers]          # corresponding signal wavelengths

# 4. Plot the tuning curves (signal and idler wavelengths vs temperature)
fig, ax = plt.subplots()
ax.plot(temps, signals, label="Signal (λs)", color='blue')
ax.plot(temps, idlers, label="Idler (λi)", color='green')
ax.set_xlabel("Temperature (°C)")
ax.set_ylabel("Wavelength (µm)")
ax.set_title(f"Tuning Curve (λp={lambda_p:.{decimals}f} µm, L0={L0:.{decimals}f} µm)")
ax.grid(True)
ax.legend()

# Save plot to buffer and display in Streamlit
buf = io.BytesIO()
fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
buf.seek(0)
st.image(buf, width=700)

# 5. Find degenerate temperature where λi = 2 * λp (signal and idler equal)
lambda_deg = 2 * lambda_p  # Degenerate wavelength for signal/idler (since λs=λi=λp/2)
f_deg = lambda T: find_lambda_idler(T, lambda_p, L0) - lambda_deg  # function whose root is degenerate temp

try:
    # Find root in the middle of temperature range
    T_deg = fsolve(f_deg, x0=(T_min + T_max) / 2)[0]
    # Display only if root within selected temperature range
    if T_min <= T_deg <= T_max:
        st.write(f"**Degenerate Temperature:** {T_deg:.{decimals}f} °C")
        st.write(f"At degeneracy λs = λi = {lambda_deg:.{decimals}f} µm")
    else:
        st.write("No degenerate temperature found within the specified range.")
except Exception:
    st.write("Could not compute degenerate temperature.")

# 6. Single temperature analysis section
st.sidebar.header("Single-Temperature Analysis")
T_set = st.sidebar.number_input(
    "Select temperature (°C):", min_value=float(T_min), max_value=float(T_max),
    value=float((T_min + T_max) / 2)
)

# Button to compute signal and idler wavelengths at chosen temperature
if st.sidebar.button("Compute λs & λi"):
    li = find_lambda_idler(T_set, lambda_p, L0)
    ls = lambda_signal(li, lambda_p)
    st.write(f"At T = {T_set:.{decimals}f} °C → λi = {li:.{decimals}f} µm, λs = {ls:.{decimals}f} µm")
