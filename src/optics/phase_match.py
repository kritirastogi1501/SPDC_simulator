# phase_match.py

import sys
sys.path.append("C:/Users/Asus/OneDrive/Desktop/Sanya Personal/Projects/DRDO internship/spdc_simulator/src")

import math
from optics.sellmeier import n_y_T, n_z_T, nm_to_um, grating_period

# Pump wavelength fixed at 405 nm
λp_nm = 405.0
λp_um = nm_to_um(λp_nm)

# Eq.9 inversions
def λi_from_λs(λs_nm):
    return 1.0 / (1.0/λp_nm - 1.0/λs_nm)

def λs_from_λi(λi_nm):
    return 1.0 / (1.0/λp_nm - 1.0/λi_nm)

# Eq.8 → Δk as function of (λs, T), **using grating_period(T) from Eq7**
def delta_k_s(λs_nm, T):
    λs_um = nm_to_um(λs_nm)
    λi_nm = λi_from_λs(λs_nm)
    λi_um = nm_to_um(λi_nm)

    conv = 1e-6  # μm → m
    kp =  2*math.pi * n_y_T(λp_um, T) / (λp_um * conv)
    ks =  2*math.pi * n_z_T(λs_um, T) / (λs_um * conv)
    ki =  2*math.pi * n_y_T(λi_um, T) / (λi_um * conv)
    kG =  2*math.pi / (grating_period(T) * conv)  # <-- Eq7 here!

    return kp - ks - ki - kG

# Eq.8 → Δk as function of (λi, T)
def delta_k_i(λi_nm, T):
    λi_um = nm_to_um(λi_nm)
    λs_nm = λs_from_λi(λi_nm)
    λs_um = nm_to_um(λs_nm)

    conv = 1e-6
    kp =  2*math.pi * n_y_T(λp_um, T) / (λp_um * conv)
    ks =  2*math.pi * n_z_T(λs_um, T) / (λs_um * conv)
    ki =  2*math.pi * n_y_T(λi_um, T) / (λi_um * conv)
    kG =  2*math.pi / (grating_period(T) * conv)

    return kp - ks - ki - kG
