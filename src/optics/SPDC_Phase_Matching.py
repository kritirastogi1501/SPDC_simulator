
#=================================Library======================================

from qutip import *
import numpy as np
from numpy import linalg as LA
import  math as ma ; import cmath as ca
import scipy as sc
from scipy.optimize import fsolve
from scipy.integrate import odeint
import matplotlib.pyplot as plt 
plt.rcParams.update({
    "text.usetex": False
})
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

# Y0 = [2.19229, 0.83547, 0.04970, 0.01621] # Paper
Y0 = [2.09930, 0.922683, 0.0467695, 0.0138404] # Thorlabs 
Z0 = [2.12725, 1.18431, 0.0514852, 0.6603, 100.00507, 9.68956e-3]
n1Y = [6.2897, 6.3061, -6.0629, 2.6486]
n2Y = [-0.14445, 2.2244, -3.5770, 1.3470]
n1Z = [9.9587, 9.9228, -8.9603, 4.1010]
n2Z = [-1.1882, 10.459, -9.8136, 3.1481]
#==============================================================================

#================================Sellmeier_Equation============================

def nY(lam):
    k1 = Y0[0] + Y0[1] / (1 - Y0[2] / (lam**2)) - Y0[3] * lam**2
    return np.sqrt(k1)

def nZ(lam):
    k1 = Z0[0] + Z0[1] / (1 - Z0[2] / (lam**2)) - Z0[3] / (1 - Z0[4] / (lam**2)) - Z0[5] * lam**2
    return np.sqrt(k1)
#==============================================================================

def Lambda(T):
     p = 10 ; α = 6.7e-6 ; β = 11e-9
     return p*(1 + α*(T - 21) + β*(T - 21)**2)

#===============================Temp_Denpen====================================


def n1(lam, pol):
    coeffs = n1Y if pol == 'y' else n1Z
    return sum(a / (lam**m) for m, a in enumerate(coeffs)) * 1e-6  # Scale n1 correctly

def n2(lam, pol):
    coeffs = n2Y if pol == 'y' else n2Z
    return sum(a / (lam**m) for m, a in enumerate(coeffs)) * 1e-8  # Scale n2 correctly

def Delta_n(lam, T, pol):
    return n1(lam, pol) * (T - 21) + n2(lam, pol) * (T - 21)**2

#==============================================================================

def Lambda_s(lam_p, lam_i):
    return 1 / (1 / lam_p - 1 / lam_i)


λ_p = 0.405

# Define the equation to solve
def eqn1(lam_i, T):
    return ((nY(λ_p) + Delta_n(λ_p, T, pol ='y')) / λ_p - \
           (nZ(Lambda_s(λ_p, lam_i)) + Delta_n(Lambda_s(λ_p, lam_i), T, pol ='z')) / Lambda_s(λ_p, lam_i) - \
           (nY(lam_i) + Delta_n(lam_i, T, pol ='y')) / lam_i - 1 / Lambda(T))


def find_lambda_i(T):
    sol = fsolve(eqn1, x0 = 1.5, args=(T,), xtol=1e-4)
    return sol[0]  # Return the first (and only) solution

# Define temperature range
Tmin, Tmax = 0, 150
T_vals = np.linspace(Tmin, Tmax, 500)
lambda_i_vals = [find_lambda_i(T) for T in T_vals]

Lambda_s_vec = np.vectorize(Lambda_s)

# Plotting
plt.plot(T_vals, lambda_i_vals,label=r"$\lambda_{i}$")
plt.plot(T_vals, Lambda_s_vec(λ_p, lambda_i_vals),label=r"$\lambda_{s}$")
plt.xlabel("Temperature (°C)")
plt.ylabel("Wavelength (µm)")
plt.title("Signal and Idler Wavelength vs. Temperature for Type-II SPDC")
plt.grid(True) ; plt.legend()
plt.show()



















