# app.py

import sys
sys.path.append("C:/Users/Asus/OneDrive/Desktop/Sanya Personal/Projects/DRDO internship/spdc_simulator/src")

# app.py

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

# functions
poly_inv = lambda lam, c: sum(cx / lam**i for i, cx in enumerate(c))

def nY(lam):
    return np.sqrt(Y0[0] + (Y0[1] / (lam**2 - Y0[2]) - Y0[3]) * lam**2)

def nZ(lam):
    return np.sqrt(Z0[0] + ((Z0[1] / (lam**2 - Z0[2]) + Z0[3] / (lam**2 - Z0[4]) - Z0[5]) * lam**2))

def delta_nY(lam, T):
    dT = T - 25
    return poly_inv(lam, aY1) * dT + poly_inv(lam, aY2) * dT**2

def delta_nZ(lam, T):
    dT = T - 25
    return poly_inv(lam, aZ1) * dT + poly_inv(lam, aZ2) * dT**2

def nY_T(lam, T): return nY(lam) + delta_nY(lam, T)
def nZ_T(lam, T): return nZ(lam) + delta_nZ(lam, T)

# poling period
def Lambda_QPM(T, L0=9.925):
    dT = T - 25
    return L0 * (1 + 6.7e-6 * dT + 1.1e-8 * dT**2)

# energy conservation
def lambda_signal(lambda_idler, lambda_p):
    return 1 / (1/lambda_p - 1/lambda_idler)

# phase matching
def phase_match_eq(lambda_idler, T, lambda_p):
    lam_s = lambda_signal(lambda_idler, lambda_p)
    return (nY_T(lambda_p, T)/lambda_p - nZ_T(lam_s, T)/lam_s - nY_T(lambda_idler, T)/lambda_idler - 1/Lambda_QPM(T))

def find_lambda_idler(T, lambda_p):
    return fsolve(phase_match_eq, x0=2*lambda_p, args=(T, lambda_p), xtol=1e-6)[0]

# ----------------------- Streamlit UI ------------------------
st.title("Type-II SPDC Wavelength vs. Temperature")

# 1. Display settings
decimals = st.sidebar.slider("Decimal places", min_value=0, max_value=10, value=4)

# 2. Pump & scan settings
st.sidebar.header("Pump & Temperature Sweep")
lambda_p = st.sidebar.number_input("Pump wavelength (µm)", min_value=0.1, max_value=1.5, value=0.405, format="%.5f")
T_min = st.sidebar.number_input("Min temperature (°C)", value=0)
T_max = st.sidebar.number_input("Max temperature (°C)", value=120)
points = st.sidebar.slider("Resolution (# points)", min_value=10, max_value=1000, value=400)

# compute tuning curve
T_vals = np.linspace(T_min, T_max, points)
idler_vals = [find_lambda_idler(T, lambda_p) for T in T_vals]
signal_vals = [lambda_signal(li, lambda_p) for li in idler_vals]

# plot
fig, ax = plt.subplots()
ax.plot(T_vals, signal_vals, label='Signal (λs)')
ax.plot(T_vals, idler_vals, label='Idler (λi)')
ax.set_xlabel('Temperature (°C)')
ax.set_ylabel('Wavelength (µm)')
ax.set_title(f'Tuning curve for λp={lambda_p:.{decimals}f} µm')
ax.grid(True)
ax.legend()
st.pyplot(fig)

# degenerate temperature
lambda_deg = 2 * lambda_p
f_deg = lambda T: find_lambda_idler(T, lambda_p) - lambda_deg
try:
    T_deg = fsolve(f_deg, x0=(T_min+T_max)/2)[0]
    if T_min <= T_deg <= T_max:
        st.write(f"**Degenerate temperature:** {T_deg:.{decimals}f} °C  \nAt degeneracy λs = λi = {lambda_deg:.{decimals}f} µm")
    else:
        st.write("No degenerate temperature found in the selected range.")
except Exception:
    st.write("Error finding degenerate temperature.")



# end of app
