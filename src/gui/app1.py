import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import io

# Import core index, thermo-optic, and QPM functions
from optics.SPDC_Phase_Matching import (
    poly_inv, nY, nZ, delta_nY, delta_nZ,
    nY_T, nZ_T, Lambda_QPM
)
# energy conservation
from optics.SPDC_Phase_Matching import lambda_signal  # uses global lambda_p ignored here

# ---------------------------------------------------------------------------- #
# Module: Wavelength vs Temperature
# ---------------------------------------------------------------------------- #

def run():
    """
    Streamlit UI for plotting SPDC signal/idler wavelengths
    versus temperature for a given pump λ and poling L0.
    """
    # Sidebar: display precision
    st.sidebar.header("Display Settings")
    decimals = st.sidebar.slider("Decimal places:", 0, 10, 4)

    # Sidebar: SPDC parameters
    st.sidebar.header("Simulation Parameters")
    lambda_p = st.sidebar.number_input(
        "Pump wavelength λp (µm):", 0.1, 2.0,
        value=0.405, format=f"%.{decimals}f"
    )
    L0 = st.sidebar.number_input(
        "Poling period L0 (µm):", 1.0, 20.0,
        value=9.925, format=f"%.{decimals}f"
    )
    T_min = st.sidebar.number_input("Min temperature (°C):", value=0)
    T_max = st.sidebar.number_input("Max temperature (°C):", value=120)
    points = st.sidebar.slider("Resolution (# points):", 10, 1000, 400)

    # Validate input range
    if T_max <= T_min:
        st.sidebar.error("Max temperature must exceed min temperature.")
        return

    # Local solver using passed λp and L0
    def find_idler(T):
        def phase_mismatch(li):
            ls = 1/(1/lambda_p - 1/li)
            return (
                nY_T(lambda_p, T)/lambda_p
                - nZ_T(ls, T)/ls
                - nY_T(li, T)/li
                - 1/Lambda_QPM(T, L0)
            )
        sol = fsolve(phase_mismatch, x0=2*lambda_p, xtol=1e-6)
        return sol[0]

    # Compute tuning curves
    temps = np.linspace(T_min, T_max, points)
    idlers = [find_idler(T) for T in temps]
    signals = [1/(1/lambda_p - 1/li) for li in idlers]

    # Plot
    fig, ax = plt.subplots()
    ax.plot(temps, signals, label="Signal (λs)")
    ax.plot(temps, idlers, label="Idler (λi)")
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Wavelength (µm)")
    ax.set_title(f"Tuning Curve: λp={lambda_p:.{decimals}f} µm, L0={L0:.{decimals}f} µm")
    ax.grid(True)
    ax.legend()

    # Display figure
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    buf.seek(0)
    st.image(buf, width = 700)

    # Degenerate temperature
    lambda_deg = 2 * lambda_p
    def f_deg(T): return find_idler(T) - lambda_deg
    try:
        T_deg = fsolve(f_deg, x0=(T_min + T_max)/2)[0]
        if T_min <= T_deg <= T_max:
            st.write(f"**Degenerate Temperature:** {T_deg:.{decimals}f} °C")
            st.write(f"λs = λi = {lambda_deg:.{decimals}f} µm")
        else:
            st.write("No degenerate temperature found within range.")
    except:
        st.write("Could not compute degenerate temperature.")

    # Single-temperature analysis
    st.sidebar.header("Single-Temperature Analysis")
    T_set = st.sidebar.number_input(
        "Select temperature (°C):",
        float(T_min), float(T_max), (T_min + T_max)/2
    )
    if st.sidebar.button("Compute λs & λi"):
        li = find_idler(T_set)
        ls = 1/(1/lambda_p - 1/li)
        st.write(f"At T={T_set:.{decimals}f}°C → λi={li:.{decimals}f} µm, λs={ls:.{decimals}f} µm")
