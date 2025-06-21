# gui/app2.py
import streamlit as st
import numpy as np
import plotly.graph_objects as go

"""
Module: Poling Period vs Temperature

This Streamlit module computes and plots the temperature-dependent
poling period using the relation:

    Pol(T) = Pol_0 * (1 + A_pol*(T - 25.0) + B_pol*(T - 25.0)**2)

where:
  - Pol_0: poling period at 25°C (user-defined)
  - A_pol = 6.7e-6 (°C⁻¹)
  - B_pol = 11e-9  (°C⁻²)

Features:
  1. Slider for decimal precision
  2. Input for Pol_0 with formatting
  3. Temperature sweep range and resolution
  4. Interactive Plotly chart with hover tooltips
  5. Single-temperature lookup
"""

def run():
    """
    Launch the poling-period vs. temperature UI.
    """
    # --- Sidebar: Display settings ---
    st.sidebar.header("Display Settings")
    decimals = st.sidebar.slider("Decimal places:", 0, 10, 4)

    # --- Sidebar: Simulation parameters ---
    st.sidebar.header("Simulation Parameters")
    # Poling period at 25°C
    Pol_0 = st.sidebar.number_input(
        "Poling period Pol_0 (µm):",
        min_value=0.1, max_value=50.0,
        value=9.925,
        format=f"%.{decimals}f"
    )
    # Temperature sweep
    T_min = st.sidebar.number_input(
        "Min temperature (°C):", value=0.0,
        format=f"%.{decimals}f"
    )
    T_max = st.sidebar.number_input(
        "Max temperature (°C):", value=120.0,
        format=f"%.{decimals}f"
    )
    points = st.sidebar.slider(
        "Resolution (# points):",
        min_value=10, max_value=1000, value=10
    )

    # Validate range
    if T_max <= T_min:
        st.sidebar.error("Max temperature must exceed min temperature.")
        return

    # --- Constants for thermal expansion ---
    A_pol = 6.7e-6   # linear coefficient (°C⁻¹)
    B_pol = 11e-9    # quadratic coefficient (°C⁻²)

    # --- Compute poling period curve ---
    temps = np.linspace(T_min, T_max, points)
    # Pol(T) = Pol_0 * [1 + A_pol*(T-25) + B_pol*(T-25)^2]
    delT = temps - 25.0
    Pol_vals = Pol_0 * (1 + A_pol * delT + B_pol * delT**2)

    # --- Build interactive Plotly figure ---
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=temps, y=Pol_vals,
            mode='lines+markers',
            name='Poling Period',
            hovertemplate='T=%{x:.2f}°C<br>Pol= %{y:.{decimals}f} µm'
        )
    )
    fig.update_layout(
        title=f"Poling Period vs Temperature (Pol_0={Pol_0:.{decimals}f} µm)",
        xaxis_title='Temperature (°C)',
        yaxis_title='Poling Period (µm)',
        hovermode='x unified'
    )

    # --- Display chart ---
    st.plotly_chart(fig, use_container_width=True)

    # --- Single-temperature lookup ---
    st.sidebar.header("Single-Temperature Lookup")
    T_sel = st.sidebar.number_input(
        "Select temperature (°C):",
        min_value=float(T_min), max_value=float(T_max),
        value=float((T_min + T_max) / 2),
        format=f"%.{decimals}f"
    )
    if st.sidebar.button("Compute Pol(T)"):
        Pol_T = Pol_0 * (
            1 + A_pol * (T_sel - 25.0) + B_pol * (T_sel - 25.0)**2
        )
        st.write(
            f"At T = {T_sel:.{decimals}f} °C, Pol(T) = {Pol_T:.{decimals}f} µm"
        )
