#Sellmeier equation & thermo-optics

"""
sellmeier.py

Sellmeier equations for PPKTP crystal refractive indices with temperature dependence.
Supports calculation of n_x, n_y, n_z given wavelength (in microns) and temperature (in Celsius).
"""

import math

# Sellmeier coefficients from QDIT Labs Entangled Photon Source Simulation Suite
# Wavelength in microns, temperature in Celsius

# Coefficients for n_x, n_y, n_z at 25°C
# n^2(λ) = A + B/(λ^2 - C) - D*λ^2

sellmeier_coeffs = {
    'x': {  # n_x (ordinary)
        'A': 3.29100,
        'B': 0.04140,
        'C': 0.03978,
        'D': 9.35522,
        'E': 31.45571
    },
    'y': {  # n_y (ordinary)
        'A': 3.45018,
        'B': 0.04341,
        'C': 0.04597,
        'D': 9.35522,
        'E': 16.98825
    },
    'z': {  # n_z (extraordinary)
        'A': 4.59423,
        'B': 0.06206,
        'C': 0.04763,
        'D': 110.80672,
        'E': 86.12171
    }
}

# Thermo-optic coefficients (dn/dT) approximate for PPKTP (per axis)
thermo_optic_coeffs = {
    'x': 1.7e-5,  # per °C
    'y': 1.8e-5,
    'z': 1.2e-5
}


def sellmeier_n_squared(axis: str, wavelength_um: float) -> float:
    """
    Calculate n^2 using Sellmeier equation for given axis and wavelength at reference temp 25°C.
    """
    cfs = sellmeier_coeffs[axis] # Get coefficients for the given axis
    λ2 = wavelength_um ** 2
    n2 = cfs['A'] + cfs['B'] / (λ2 - cfs['C']) + cfs['D'] / (λ2 - cfs['E'])
    return n2


def refractive_index(axis: str, wavelength_um: float, temperature_C: float) -> float:
    """
    Calculate refractive index for given axis, wavelength (μm), and temperature (°C).
    Applies thermo-optic correction linearly.
    """
    # Step 1: Calculate n(λ) at 25°C
    n_ref = math.sqrt(sellmeier_n_squared(axis, wavelength_um))

    # Step 2: Apply linear temperature correction using dn/dT
    dndT = thermo_optic_coeffs[axis]         # Get the temperature coefficient
    delta_T = temperature_C - 25.0           # ΔT from reference temperature (25°C)
    n_temp = n_ref + dndT * delta_T          # Adjust n using: n(T) = n_ref + dn/dT * ΔT
    return n_temp


if __name__ == "__main__":
    # Simple test/demo
    wavelengths = [0.4, 0.8, 1.064, 1.55]  # μm
    temperature = 30.0  # °C

    print(f"Refractive indices at T={temperature} °C:")
    for wl in wavelengths:
        nx = refractive_index('x', wl, temperature)
        ny = refractive_index('y', wl, temperature)
        nz = refractive_index('z', wl, temperature)
        print(f"λ = {wl:.3f} μm: n_x={nx:.5f}, n_y={ny:.5f}, n_z={nz:.5f}")
