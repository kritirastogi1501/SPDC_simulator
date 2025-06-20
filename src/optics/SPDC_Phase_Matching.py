# -*- coding: utf-8 -*-
"""
Updated on: 2025-06-20
Updated by: ChatGPT
Modified to use Eq1 to Eq4 for refractive index calculations and 
thermo-optic corrections for KTP Type-II SPDC simulation.
"""

#=================================Library======================================
from qutip import *
import numpy as np
from numpy import linalg as LA
import math as ma
import cmath as ca
import scipy as sc
from scipy.optimize import fsolve
from scipy.integrate import odeint
import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex": False})
from mpl_toolkits import mplot3d
import sympy as sp
import random as ra
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from tqdm import tqdm
import time

I = complex(0,1)
I_2 = qeye(2)

#================================Contants_for_KTP==============================
Y0 = [2.09930, 0.922683, 0.0467695, 0.0138404]  # Eq1
Z0 = [2.12725, 1.18431, 0.0514852, 0.6603, 100.00507, 9.68956e-3]  # Eq2

#==============================================================================
#===============================Sellmeier Equations=============================
def nY(lam):
    # Eq1
    return np.sqrt(Y0[0] + ((Y0[1] / (lam**2 - Y0[2])) - Y0[3]) * lam**2)

def nZ(lam):
    # Eq2
    return np.sqrt(Z0[0] + ((Z0[1] / (lam**2 - Z0[2]) + Z0[3] / (lam**2 - Z0[4]) - Z0[5]) * lam**2))

#==============================Temp Corrections (Eq3 & Eq4)=====================
def delta_nY(lam, T):
    dT = T - 25.0
    return (
        6.2897e-6 * dT - 1.4445e-9 * dT**2 +
        ((6.3061e-6 + 2.2244e-8 * dT) * dT) / lam +
        ((6.0629e-6 + 3.5770e-8 * dT) * dT) / lam**2 +
        ((2.6486e-6 + 1.3470e-8 * dT) * dT) / lam**3
    )

def delta_nZ(lam, T):
    dT = T - 25.0
    return (
        9.9587e-6 * dT - 1.1882e-8 * dT**2 +
        ((9.9228e-6 + 1.0459e-7 * dT) * dT) / lam +
        ((8.9603e-6 + 9.8136e-8 * dT) * dT) / lam**2 +
        ((4.1010e-6 + 3.1481e-8 * dT) * dT) / lam**3
    )

#===============================Final Refractive Index==========================
def nY_T(lam, T):
    return nY(lam) + delta_nY(lam, T)

def nZ_T(lam, T):
    return nZ(lam) + delta_nZ(lam, T)

#===============================Poling Period==================================
def Lambda(T):
    p = 9.925  # micron
    alpha = 6.7e-6
    beta = 11e-9
    return p * (1 + alpha * (T - 25) + beta * (T - 25)**2)

#===============================Wavelength Calculations=========================
def Lambda_s(lam_p, lam_i):
    return 1 / (1 / lam_p - 1 / lam_i)

lambda_p = 0.405  # micron

def eqn1(lam_i, T):
    lam_s = Lambda_s(lambda_p, lam_i)
    return 2 * np.pi *(
        nY_T(lambda_p, T) / lambda_p -
        nZ_T(lam_s, T) / lam_s -
        nY_T(lam_i, T) / lam_i -
        1 / Lambda(T)
    )

def find_lambda_i(T):
    sol = fsolve(eqn1, x0=0.81, args=(T,), xtol=1e-4)
    return sol[0]

#===============================Plotting=======================================
Tmin, Tmax = 0, 120
T_vals = np.linspace(Tmin, Tmax, 500)
lambda_i_vals = [find_lambda_i(T) for T in T_vals]
Lambda_s_vec = np.vectorize(Lambda_s)

plt.plot(T_vals, lambda_i_vals, label=r"$\lambda_{i}$")
plt.plot(T_vals, Lambda_s_vec(lambda_p, lambda_i_vals), label=r"$\lambda_{s}$")
plt.xlabel("Temperature (°C)")
plt.ylabel("Wavelength (µm)")
plt.title("Signal and Idler Wavelength vs. Temperature for Type-II SPDC")
plt.grid(True)
plt.legend()
plt.show()
