# dispersion.py

from sellmeier import refractive_index  # Import the refractive index function from sellmeier.py
import numpy as np

def dn_dlambda(axis: str, wavelength_um: float, temperature_C: float, delta_um: float = 1e-4) -> float:
    """
    Calculates the derivative dn/dλ numerically using the central difference method.

    Parameters:
    - axis: 'x', 'y', or 'z' — which crystal axis
    - wavelength_um: wavelength in microns (μm)
    - temperature_C: temperature in Celsius
    - delta_um: small change in wavelength for derivative (default: 0.0001 μm)

    Returns:
    - Approximate dn/dλ in units of 1/μm
    """
    wl1 = wavelength_um - delta_um  # Slightly lower wavelength
    wl2 = wavelength_um + delta_um  # Slightly higher wavelength
    n1 = refractive_index(axis, wl1, temperature_C)  # Refractive index at λ - Δλ
    n2 = refractive_index(axis, wl2, temperature_C)  # Refractive index at λ + Δλ
    return (n2 - n1) / (wl2 - wl1)  # Central difference approximation


def group_index(axis: str, wavelength_um: float, temperature_C: float) -> float:
    """
    Calculates the group index n_g = n - λ * dn/dλ.

    Group index is used to determine group velocity and dispersion effects.

    Parameters:
    - axis: 'x', 'y', or 'z'
    - wavelength_um: wavelength in microns
    - temperature_C: temperature in Celsius

    Returns:
    - Group index (unitless)
    """
    n = refractive_index(axis, wavelength_um, temperature_C)       # Refractive index at this λ and T
    dn_dl = dn_dlambda(axis, wavelength_um, temperature_C)         # Derivative of refractive index
    return n - wavelength_um * dn_dl                               # Group index formula


def group_velocity(axis: str, wavelength_um: float, temperature_C: float) -> float:
    """
    Calculates group velocity: v_g = c / n_g

    Parameters:
    - axis: 'x', 'y', or 'z'
    - wavelength_um: wavelength in microns
    - temperature_C: temperature in Celsius

    Returns:
    - Group velocity in m/s
    """
    c = 299792458  # Speed of light in vacuum (m/s)
    n_g = group_index(axis, wavelength_um, temperature_C)  # Compute group index
    return c / n_g  # Group velocity


if __name__ == "__main__":
    # Example usage and test output
    wavelength = 1.064  # in μm
    temperature = 30.0  # in °C

    for axis in ['x', 'y', 'z']:
        n = refractive_index(axis, wavelength, temperature)
        dn_dl = dn_dlambda(axis, wavelength, temperature)
        n_g = group_index(axis, wavelength, temperature)
        v_g = group_velocity(axis, wavelength, temperature)

        print(f"\nAxis: {axis}")
        print(f"Refractive index n = {n:.6f}")
        print(f"dn/dλ = {dn_dl:.6e} /μm")
        print(f"Group index n_g = {n_g:.6f}")
        print(f"Group velocity v_g = {v_g:.3e} m/s")
