import numpy as np
import matplotlib.pyplot as plt

# Constants
Lambda_0 = 10  # micrometers
A_pol = 6.7e-6  # °C⁻¹
B_pol = 11e-9   # °C⁻²

# Temperature range from -50°C to 150°C with bigger steps (e.g., 5°C steps)
T = np.arange(0, 500, 500)

# Calculate Λ for each temperature
Lambda = Lambda_0 * (1 + A_pol * (T - 25.0) + B_pol * (T - 25.0)**2)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(T, Lambda, marker='o', color='purple', label='Λ vs T')
plt.xlabel('Temperature (°C)')
plt.ylabel('Λ (μm)')
plt.title('Plot of Λ vs Temperature')
plt.grid(True)
plt.ylim(9, 11)
plt.xlim(0, 500)
plt.xticks(np.arange(0, 500, 500))  # Larger x-axis steps
plt.legend()
plt.tight_layout()
plt.show()
