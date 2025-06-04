import sys
import numpy as np
import matplotlib.pyplot as plt

# Ensure Python can find `optics`:
sys.path.append("C:/Users/Asus/OneDrive/Desktop/Sanya Personal/Projects/DRDO internship/spdc_simulator/src")

from optics.phase_match import (
    find_phase_match_lambda_s,
    find_degenerate_temperature
)

# -----------------------------
# 1. Constants
# -----------------------------
lambda_p_um = 0.40533       # pump wavelength = 405.33 nm → 0.40533 µm
# ───► **Use 2.508 µm (≈ poling period that phase‐matches 405 nm → 810 nm at 25 °C)**
poling_period_um = 2.5081181150762486  
interaction = 'type2'       # 'type0' or 'type2'

# Degenerate signal = 2 × λₚ = 2 × 0.40533 µm = 0.81066 µm (≈ 810.66 nm)
degenerate_lambda_s_um = 2.0 * lambda_p_um
degenerate_lambda_nm = degenerate_lambda_s_um * 1e3  # in nm

# Temperature sweep (°C)
temps = np.linspace(10, 80, 300)  # you can widen or shift this as needed

# Preallocate lists for plotting (in nm)
lambda_s_list_nm = []
lambda_i_list_nm = []

# -----------------------------
# 2. Find phase‐matched λₛ(T) and λᵢ(T)
# -----------------------------
for T in temps:
    # Since we know λₛ will be around 0.8 µm for Type‐II PPKTP with Λ≈2.508 µm,
    # we narrow the search window to [0.7, 0.9] µm. Brent’s method will run faster,
    # and we will capture the small wavelength shift around 810 nm as T varies.
    λs_um = find_phase_match_lambda_s(
        T,
        lambda_p_um,
        poling_period_um,
        interaction_type=interaction,
        search=(0.70, 0.90)
    )

    if λs_um is not None:
        # Idler: 1/λₚ = 1/λₛ + 1/λᵢ  ⇒  λᵢ = 1 / (1/λₚ − 1/λₛ)
        λi_um = 1.0 / (1.0 / lambda_p_um - 1.0 / λs_um)
        lambda_s_list_nm.append(λs_um * 1e3)
        lambda_i_list_nm.append(λi_um * 1e3)
    else:
        # No real root at this T → plot NaN (Matplotlib will skip it)
        lambda_s_list_nm.append(np.nan)
        lambda_i_list_nm.append(np.nan)

# -----------------------------
# 3. Find degenerate temperature T₍deg₎ so that λₛ=λᵢ=0.81066 µm
# -----------------------------
try:
    T_deg = find_degenerate_temperature(
        degenerate_lambda_s_um,
        lambda_p_um,
        poling_period_um,
        interaction_type=interaction,
        search=(0.0, 100.0)  # widen if needed
    )
except ValueError:
    T_deg = None

# -----------------------------
# 4. Plotting
# -----------------------------
plt.figure(figsize=(8, 5))
plt.plot(temps, lambda_s_list_nm, 'o', label="Signal λₛ", markersize=3)
plt.plot(temps, lambda_i_list_nm, 'o', label="Idler λᵢ", markersize=3)

if T_deg is not None:
    # Draw the degenerate point (vertical + horizontal lines + black dot)
    plt.axhline(degenerate_lambda_nm, linestyle='--', color='gray', linewidth=1)
    plt.axvline(T_deg, linestyle=':', color='gray', linewidth=1)
    plt.plot(T_deg, degenerate_lambda_nm, 'ko', label="Degeneracy Point")

plt.xticks([10, 20, 30, 40, 50, 60, 70, 80])
plt.yticks([750,760,770,780, 790, 800, 810, 820, 830, 840])
plt.xlim(10, 80)
plt.ylim(780, 840)

plt.xlabel("Crystal Temperature (°C)")
plt.ylabel("Wavelength (nm)")
plt.title("SPDC Phase-Matching (Type-II) in PPKTP: Signal & Idler vs Temperature")
plt.grid(True)
plt.legend(loc="upper right")
plt.tight_layout()
plt.show()

if T_deg is not None:
    print(f"✅ Degenerate phase‐matching at T = {T_deg:.2f} °C "
          f"(λₛ = λᵢ = {degenerate_lambda_nm:.2f} nm)")
else:
    print("⚠️ No degenerate temperature found in [0 °C, 100 °C].")
