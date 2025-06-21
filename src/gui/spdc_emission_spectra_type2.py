# gui/app3.py
import streamlit as st
import numpy as np
from scipy.optimize import fsolve
import plotly.graph_objects as go

# Import core thermo-optic & QPM functions
from optics.SPDC_Phase_Matching import nY_T, nZ_T, Lambda_QPM

"""
Module: SPDC Signal & Idler Spectra vs Wavelength

Uses:
  I_s(λ_s) ∝ sinc²[Δk(λ_s,T)·L/2]
  I_i(λ_i) ∝ sinc²[Δk(λ_i,T)·L/2]

where Δk = k_p - k_s - k_i - 2π/ΛQPM(T)
and L is crystal length.

Features:
 1. Decimal precision slider
 2. Pump λ, poling L0, crystal length L
 3. Temperature slider
 4. Wavelength range and resolution sliders
 5. Interactive Plotly spectrum plot
 6. Peak intensity lookup at selected T
"""

def run():
    # --- Sidebar inputs ---
    st.sidebar.header("Display Settings")
    decimals = st.sidebar.slider("Decimal places:", 0, 6, 4)

    st.sidebar.header("SPDC Parameters")
    lambda_p = st.sidebar.number_input(
        "Pump wavelength λp (µm):", 0.1, 2.0,
        value=0.405, format=f"%.{decimals}f"
    )
    L0 = st.sidebar.number_input(
        "Poling period L0 (µm):", 1.0, 20.0,
        value=9.925, format=f"%.{decimals}f"
    )
    L = st.sidebar.number_input(
        "Crystal length L (mm):", 0.1, 100.0,
        value=10.0, format=f"%.{decimals}f"
    )
    # convert L to µm
    L_um = L * 1e3

    st.sidebar.header("Scan Settings")
    T_set = st.sidebar.number_input(
        "Temperature (°C):", min_value=0.0, max_value=200.0, value=66.8832,
        format=f"%.{decimals}f"
    )
    lam_min = st.sidebar.number_input(
        "Min wavelength (µm):", 0.1, 5.0,
        value=0.7, format=f"%.{decimals}f"
    )
    lam_max = st.sidebar.number_input(
        "Max wavelength (µm):", 0.1, 5.0,
        value=1.0, format=f"%.{decimals}f"
    )
    pts = st.sidebar.slider(
        "Resolution (# points):", 100, 2000, 500
    )

    if lam_max <= lam_min:
        st.sidebar.error("Max wavelength must exceed min wavelength.")
        return

    # --- Compute spectra ---
    lam_s = np.linspace(lam_min, lam_max, pts)
    lam_i = 1.0 / (1.0 / lambda_p - 1.0 / lam_s)

    def delta_k(ls, li, T):
        kp = 2 * np.pi * nY_T(lambda_p, T) / lambda_p
        ks = 2 * np.pi * nZ_T(ls, T) / ls
        ki = 2 * np.pi * nY_T(li, T) / li
        kq = 2 * np.pi / Lambda_QPM(T, L0)
        return kp - ks - ki - kq

    dks = delta_k(lam_s, lam_i, T_set)
    Is = (np.sinc(dks * L_um / (2 * np.pi)))**2

    lam_i_scan = np.linspace(lam_min, lam_max, pts)
    lam_s_from_i = 1.0 / (1.0 / lambda_p - 1.0 / lam_i_scan)
    dk_i = delta_k(lam_s_from_i, lam_i_scan, T_set)
    Ii = (np.sinc(dk_i * L_um / (2 * np.pi)))**2

    # --- Plot spectra ---
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=lam_s, y=Is,
        mode='lines', name='Signal I_s',
        hovertemplate=f'λs=%{{x:.{decimals}f}} µm<br>I_s=%{{y:.{decimals}f}}'
    ))
    fig.add_trace(go.Scatter(
        x=lam_i_scan, y=Ii,
        mode='lines', name='Idler I_i',
        hovertemplate=f'λi=%{{x:.{decimals}f}} µm<br>I_i=%{{y:.{decimals}f}}'
    ))
    fig.update_layout(
        title=f"SPDC Spectra at T={T_set:.{decimals}f}°C, L={L:.{decimals}f} mm",
        xaxis_title='Wavelength (µm)',
        yaxis_title='Normalized Intensity',
        hovermode='x unified'
    )
    st.plotly_chart(fig, use_container_width=True)

    # --- Peak lookup ---
    idx_s = np.argmax(Is)
    idx_i = np.argmax(Ii)

    # Display the temperature at which these peaks were computed
    st.write(f"**Analysis Temperature:** {T_set:.{decimals}f} °C")
    st.write(f"Signal peak at λs = {lam_s[idx_s]:.{decimals}f} µm with I_s = {Is[idx_s]:.{decimals}f}")
    st.write(f"Idler peak at λi = {lam_i_scan[idx_i]:.{decimals}f} µm with I_i = {Ii[idx_i]:.{decimals}f}")

    # --- Peak vs Temperature Scan ---
    T_scan = np.linspace(0, 200, 200)  # You can change the range/resolution
    peak_Is = []
    peak_Ii = []

    for T in T_scan:
        dks = delta_k(lam_s, lam_i, T)
        Is_temp = (np.sinc(dks * L_um / (2 * np.pi)))**2
        peak_Is.append(np.max(Is_temp))

        dk_i = delta_k(lam_s_from_i, lam_i_scan, T)
        Ii_temp = (np.sinc(dk_i * L_um / (2 * np.pi)))**2
        peak_Ii.append(np.max(Ii_temp))

    T_max_Is = T_scan[np.argmax(peak_Is)]
    T_max_Ii = T_scan[np.argmax(peak_Ii)]

    st.write("---")
    st.write("### Peak Intensity Temperature Scan")
    st.write(f"Signal intensity **maximum** at T = {T_max_Is:.{decimals}f} °C")
    st.write(f"Idler intensity **maximum** at T = {T_max_Ii:.{decimals}f} °C")

