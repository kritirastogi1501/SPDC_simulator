import sys
sys.path.append("C:/Users/Asus/OneDrive/Desktop/Sanya Personal/Projects/DRDO internship/spdc_simulator/src")

import numpy as np
import matplotlib.pyplot as plt
from optics.phase_match import find_phase_match_lambda_s

# ------------------------------------------------------------------------------
# Plot tuning curve: λs vs T for Type-0 and Type-2 SPDC in PPKTP
# Each point is found by solving Δk = 0 numerically
# ------------------------------------------------------------------------------

lambda_p = 0.405  # pump wavelength in µm
temps = np.linspace(40, 80, 100)

results = {
    'type0': {
        'poling_period': 3.42,
        'label': 'Type-0 (Λ = 3.42 µm)',
        'color': 'blue',
        'values': []
    },
    'type2': {
        'poling_period': 10.0,
        'label': 'Type-2 (Λ = 10.0 µm)',
        'color': 'orange',
        'values': []
    }
}

for interaction in results:
    period = results[interaction]['poling_period']
    λs_list = []
    for T in temps:
        λs = find_phase_match_lambda_s(T, lambda_p, period, interaction)
        λs_list.append(λs * 1e3 if λs else np.nan)  # convert to nm
    results[interaction]['values'] = λs_list

# --- Plotting ---
plt.figure(figsize=(8, 5))
for interaction in results:
    plt.plot(temps, results[interaction]['values'], label=results[interaction]['label'], color=results[interaction]['color'])

plt.xlabel("Crystal Temperature (°C)", fontsize=12)
plt.ylabel("Signal Wavelength λs (nm)", fontsize=12)
plt.title("PPKTP SPDC Tuning Curves (Δk = 0)", fontsize=13)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
