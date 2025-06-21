import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# ----------------- Sellmeier & thermo-optic, from user script -----------------
Y0 = [2.19229, 0.83547, 0.04970, 0.01621]
Z0 = [2.12725, 1.18431, 0.0514852, 0.6603, 100.00507, 9.68956e-3]

aY1 = [6.2897e-6, 6.3061e-6, -6.0629e-6, 2.6486e-6]
aY2 = [-1.4445e-9, 2.2244e-8, -3.5770e-8, 1.3470e-8]
aZ1 = [9.9587e-6, 9.9228e-6, -8.9603e-6, 4.1010e-6]
aZ2 = [-1.1882e-8, 1.0459e-7, -9.8136e-8, 3.1481e-8]

# helper for polynomial thermo-optic terms
def poly_inv(lam, coeffs):
    return sum(c / lam**i for i, c in enumerate(coeffs))

# refractive indices
def nY(lam):
    return np.sqrt(Y0[0] + (Y0[1] / (lam**2 - Y0[2]) - Y0[3]) * lam**2)

def nZ(lam):
    return np.sqrt(Z0[0] + ((Z0[1] / (lam**2 - Z0[2]) + Z0[3] / (lam**2 - Z0[4]) - Z0[5]) * lam**2))

# thermo-optic shifts
def delta_nY(lam, T):
    dT = T - 25.0
    return poly_inv(lam, aY1) * dT + poly_inv(lam, aY2) * dT**2

def delta_nZ(lam, T):
    dT = T - 25.0
    return poly_inv(lam, aZ1) * dT + poly_inv(lam, aZ2) * dT**2

def nY_T(lam, T):
    return nY(lam) + delta_nY(lam, T)

def nZ_T(lam, T):
    return nZ(lam) + delta_nZ(lam, T)

# poling period expansion
def Lambda_QPM(T, L0):
    dT = T - 25.0
    alpha, beta = 6.7e-6, 1.1e-8
    return L0 * (1 + alpha * dT + beta * dT**2)

# energy conservation
def lambda_signal(lambda_idler, lambda_p):
    return 1 / (1/lambda_p - 1/lambda_idler)

# phase-matching equation
def phase_match_eq(lambda_idler, T, lambda_p, L0):
    lam_s = lambda_signal(lambda_idler, lambda_p)
    return (
        nY_T(lambda_p, T)/lambda_p
        - nZ_T(lam_s, T)/lam_s
        - nY_T(lambda_idler, T)/lambda_idler
        - 1/Lambda_QPM(T, L0)
    )

# find idler root
def find_lambda_idler(T, lambda_p, L0):
    sol = fsolve(phase_match_eq, x0=2*lambda_p, args=(T, lambda_p, L0), xtol=1e-6)
    return sol[0]

# ----------------------- Streamlit UI ------------------------
st.set_page_config(page_title="SPDC Tuning Curves", layout="wide")
st.title("Type-II SPDC Wavelength vs. Temperature")

# 1. Display settings
st.sidebar.header("Display Settings")
decimals = st.sidebar.slider("Decimal places:", min_value=0, max_value=10, value=4)

# 2. Pump, poling & temperature sweep
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

# compute tuning curves
temps = np.linspace(T_min, T_max, points)
idlers = [find_lambda_idler(T, lambda_p, L0) for T in temps]
signals = [lambda_signal(li, lambda_p) for li in idlers]

# plot curves
fig, ax = plt.subplots()
ax.plot(temps, signals, label="Signal (λs)")
ax.plot(temps, idlers, label="Idler (λi)")
ax.set_xlabel("Temperature (°C)")
ax.set_ylabel("Wavelength (µm)")
ax.set_title(f"Tuning Curve (λp={lambda_p:.{decimals}f} µm, L0={L0:.{decimals}f} µm)")
ax.grid(True)
ax.legend()
st.pyplot(fig)

# degenerate temperature finder
lambda_deg = 2 * lambda_p
f_deg = lambda T: find_lambda_idler(T, lambda_p, L0) - lambda_deg
try:
    T_deg = fsolve(f_deg, x0=(T_min+T_max)/2)[0]
    if T_min <= T_deg <= T_max:
        st.write(f"**Degenerate Temperature:** {T_deg:.{decimals}f} °C")
        st.write(f"At degeneracy λs = λi = {lambda_deg:.{decimals}f} µm")
    else:
        st.write("No degenerate temperature within range.")
except:
    st.write("Could not compute degenerate temperature.")

# 3. Single-temperature cross-section
st.sidebar.header("Single-Temperature Analysis")
T_set = st.sidebar.number_input(
    "Select temperature (°C):", min_value=float(T_min),
    max_value=float(T_max), value=float((T_min+T_max)/2)
)
if st.sidebar.button("Compute λs & λi"):
    li = find_lambda_idler(T_set, lambda_p, L0)
    ls = lambda_signal(li, lambda_p)
    st.write(f"At T = {T_set:.{decimals}f} °C → λi = {li:.{decimals}f} µm, λs = {ls:.{decimals}f} µm")
