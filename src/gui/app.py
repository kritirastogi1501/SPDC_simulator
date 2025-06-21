# app.py

import sys
sys.path.append("C:/Users/Asus/OneDrive/Desktop/Sanya Personal/Projects/DRDO internship/spdc_simulator/src")

# app.py

# app.py

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from optics.SPDC_Phase_Matching import (
    compute_typeii_wavelengths,
    find_signal_idler_at_temp,
    find_degenerate_temperature,
    lambda_signal,
    find_lambda_idler
)

# ---------------- Streamlit App UI ----------------

st.set_page_config(page_title="SPDC PPKTP Simulator", layout="wide")
st.title("🔬 SPDC in PPKTP — Wavelength vs Temperature")
st.markdown("""
Simulate **Type-II SPDC** in PPKTP:  
View how signal and idler wavelengths depend on temperature for a given pump wavelength.
""")

# Sidebar settings
st.sidebar.markdown("### ⚙️ Display Settings")
decimal_places = st.sidebar.slider("Decimal Places", min_value=2, max_value=10, value=4)
format_str = f"{{:.{decimal_places}f}}"

# User Inputs
lambda_p = st.number_input("Pump Wavelength (µm)", min_value=0.300, max_value=1.000, value=0.405, step=0.001, format="%.5f")
temp_min = st.slider("Minimum Temperature (°C)", min_value=0, max_value=150, value=20)
temp_max = st.slider("Maximum Temperature (°C)", min_value=temp_min+1, max_value=200, value=120)
num_points = st.slider("Resolution (points)", 50, 1000, 400)

# ---------------- Compute Main Tuning Curve ----------------

T_vals = np.linspace(temp_min, temp_max, num_points)
signal_vals, idler_vals = compute_typeii_wavelengths(T_vals, lambda_p=lambda_p)

fig, ax = plt.subplots()
ax.plot(T_vals, signal_vals, label="Signal (λₛ)", color="blue")
ax.plot(T_vals, idler_vals, label="Idler (λᵢ)", color="green")
ax.set_xlabel("Temperature (°C)")
ax.set_ylabel("Wavelength (µm)")
ax.set_title("Type-II SPDC Wavelength vs Temperature")
ax.grid(True)
ax.legend()
st.pyplot(fig)

# ---------------- Instant Value Lookup ----------------

st.markdown("### 📍 Signal & Idler at a Specific Temperature")
lookup_temp = st.number_input("Enter temperature (°C) to inspect", min_value=0.0, max_value=200.0, value=25.0, step=0.1)

sig, idi = find_signal_idler_at_temp(lookup_temp, lambda_p=lambda_p)

st.write(
    f"At **{format_str.format(lookup_temp)} °C** → "
    f"Signal: **{format_str.format(sig)} µm**, "
    f"Idler: **{format_str.format(idi)} µm**"
)

# ---------------- Degenerate Temperature ----------------

st.markdown("### 🎯 Find Degenerate Temperature (λₛ = λᵢ)")
if st.button("Find Degenerate Temperature"):
    T_deg = find_degenerate_temperature(lambda_p=lambda_p)
    lam_deg = lambda_signal(find_lambda_idler(T_deg))
    st.success(
        f"Degenerate temperature: **{format_str.format(T_deg)} °C**, "
        f"Degenerate λ: **{format_str.format(lam_deg)} µm**"
    )
