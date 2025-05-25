# dispersion.py

from sellmeier import refractive_index  # Import the temperature-aware refractive index function
import numpy as np  # For numerical operations

def dn_dlambda(axis: str, wavelength_um: float, temperature_C: float, delta_um: float = 1e-4) -> float:
    """
    Calculates the numerical derivative dn/dλ using the central difference method.

    This tells how fast the refractive index changes with wavelength, which is crucial for group velocity
    and dispersion analysis in optics.

    Parameters:
    - axis: 'x', 'y', or 'z' (principal optical axis of the crystal)
    - wavelength_um: wavelength in microns (μm)
    - temperature_C: temperature in degrees Celsius
    - delta_um: small step size for λ used in central difference approximation

    Returns:
    - dn/dλ (unit: 1/μm)
    """
    wl1 = wavelength_um - delta_um  # λ - Δλ
    wl2 = wavelength_um + delta_um  # λ + Δλ

    # Compute n at both slightly shifted wavelengths
    n1 = refractive_index(axis, wl1, temperature_C)
    n2 = refractive_index(axis, wl2, temperature_C)

    # Use central difference to estimate derivative
    return (n2 - n1) / (wl2 - wl1)


def d2n_dlambda2(axis: str, wavelength_um: float, temperature_C: float, delta_um: float = 1e-4) -> float:
    """
    Calculates the second derivative d²n/dλ² numerically.

    This is used for computing Group Velocity Dispersion (GVD).

    Parameters:
    - axis: 'x', 'y', or 'z'
    - wavelength_um: wavelength in microns (μm)
    - temperature_C: temperature in Celsius
    - delta_um: step size for central difference (default: 0.0001 μm)

    Returns:
    - d²n/dλ² (unit: 1/μm²)
    """
    wl0 = wavelength_um
    wl1 = wl0 - delta_um
    wl2 = wl0 + delta_um

    n0 = refractive_index(axis, wl0, temperature_C)
    n1 = refractive_index(axis, wl1, temperature_C)
    n2 = refractive_index(axis, wl2, temperature_C)

    return (n2 - 2 * n0 + n1) / (delta_um ** 2)


def group_index(axis: str, wavelength_um: float, temperature_C: float) -> float:
    """
    Computes the group index: n_g = n - λ * (dn/dλ)

    The group index characterizes the speed of the pulse envelope (group velocity).
    It is used in analyzing ultrafast laser pulses and dispersion.

    Parameters:
    - axis: 'x', 'y', or 'z'
    - wavelength_um: wavelength in μm
    - temperature_C: temperature in Celsius

    Returns:
    - Group index (unitless)
    """
    n = refractive_index(axis, wavelength_um, temperature_C)         # Get base refractive index
    dn_dl = dn_dlambda(axis, wavelength_um, temperature_C)           # Compute how n varies with λ
    return n - wavelength_um * dn_dl                                 # Apply formula for group index


def group_velocity(axis: str, wavelength_um: float, temperature_C: float) -> float:
    """
    Computes group velocity: v_g = c / n_g

    Parameters:
    - axis: 'x', 'y', or 'z'
    - wavelength_um: wavelength in microns
    - temperature_C: temperature in Celsius

    Returns:
    - Group velocity in meters per second (m/s)
    """
    c = 299_792_458  # Speed of light in vacuum (m/s)
    n_g = group_index(axis, wavelength_um, temperature_C)
    return c / n_g


def group_velocity_dispersion(axis: str, wavelength_um: float, temperature_C: float) -> float:
    """
    Computes Group Velocity Dispersion (GVD) in fs²/mm.

    GVD = (λ³ / (2πc²)) * d²n/dλ²

    Parameters:
    - axis: 'x', 'y', or 'z'
    - wavelength_um: wavelength in microns
    - temperature_C: temperature in Celsius

    Returns:
    - GVD in femtoseconds squared per millimeter (fs²/mm)
    """
    c = 299792458  
    d2n_dl2 = d2n_dlambda2(axis, wavelength_um, temperature_C)  # 1/μm²
    λ = wavelength_um * 1e-6  # Convert to meters

    # GVD in s²/m
    gvd_s2_per_m = (λ ** 3 / (2 * np.pi * c ** 2)) * d2n_dl2 * 1e12  # 1e12 for μm² to m²

    # Convert s²/m to fs²/mm
    return gvd_s2_per_m * 1e30 / 1e3  # s²→fs² (1e30), m→mm (1e3)


if __name__ == "__main__":
    # Test inputs
    wavelength = 1.064  # μm
    temperature = 30.0  # °C

    for axis in ['x', 'y', 'z']:
        n = refractive_index(axis, wavelength, temperature)
        dn_dl = dn_dlambda(axis, wavelength, temperature)
        n_g = group_index(axis, wavelength, temperature)
        v_g = group_velocity(axis, wavelength, temperature)
        gvd = group_velocity_dispersion(axis, wavelength, temperature)

        print(f"\nAxis: {axis}")
        print(f"Refractive index n = {n:.6f}")
        print(f"dn/dλ = {dn_dl:.6e} /μm")
        print(f"Group index n_g = {n_g:.6f}")
        print(f"Group velocity v_g = {v_g:.3e} m/s")
        print(f"Group Velocity Dispersion GVD = {gvd:.3f} fs²/mm")
