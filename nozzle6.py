# -*- coding: utf-8 -*-
"""
Created on Sun Sep 28 23:29:57 2025

@author: trefr
"""

import cantera as ct
import numpy as np
from scipy.integrate import solve_ivp

# Reactor parameters
A_in = 0.018
A_out = 0.003
L = 1.284 * 0.0254
n = 100
mdot = 1.125

# Gas object
gas = ct.Solution('gri30.yaml')

T_calc = np.zeros(n+1)
rho_calc = np.zeros(n+1)
nsp = gas.n_species
Y_calc = np.zeros((n+1, nsp))


# Example initial composition and T
gas.TPX = 1473, 4.47*101325, {'CH4':0.2899,'O2':2.0,'N2':7.52}
T_calc[0] = gas.T
rho_calc[0] = gas.density
Y_calc[0, :] = gas.Y


# Area slope
if A_in > A_out:
    k = -1
elif A_out > A_in:
    k = 1
else:
    k = 0
dAdx = abs(A_in - A_out) / L
dx = L / n



x_calc = np.linspace(0, L, n+1)


def PFR_ode(x, y, gas, mdot, A_in, dAdx, k):
    rho = y[0]
    T = y[1]
    Y = y[2:]
    gas.TDY = T, rho, Y
    
    #if statement for combustor section
    if x < L/2:
        dAdx = 0
        
    A = A_in + k*x*dAdx
    # velocity from mass flow
    u = mdot / (rho * A)
    # derivative of density and temperature w.r.t x
    # here you would include energy and species balances
    omega = gas.net_production_rates  # kmol/m^3/s
    dYdx = omega * gas.molecular_weights / (rho * u)  # dY/dx

    # --- Energy / temperature evolution ---
    # sum of enthalpy release from reactions
    dTdx = - np.dot(omega, gas.partial_molar_enthalpies) / (rho * u * gas.cp_mass)

    # --- Optional density change from area variation (continuity) ---
    drhodx = - rho / A * dAdx

    
    return np.concatenate(([drhodx, dTdx], dYdx))

y0 = np.concatenate(([rho_calc[0], T_calc[0]], Y_calc[0,:]))
#sol = solve_ivp(PFR_ode, [0, L], y0, args=(gas, mdot, A_in, dAdx, k), t_eval=x_calc)

sol = solve_ivp(PFR_ode, [0, L], y0,
                args=(gas, mdot, A_in, dAdx, k),
                method='BDF',  # stiff solver
                t_eval=x_calc,
                rtol=1e-8,
                atol=1e-12)

rho_calc = sol.y[0,:]
T_calc = sol.y[1,:]
Y_calc = sol.y[2:,:].T

# Compute area, velocity, Mach, pressure
A_calc = A_in + k*x_calc*dAdx
vx_calc = mdot / (rho_calc * A_calc)
R_calc = ct.gas_constant / gas.mean_molecular_weight
P_calc = rho_calc * R_calc * T_calc
M_calc = vx_calc / gas.sound_speed

import matplotlib.pyplot as plt

# -----------------------------
# Plot Temperature along PFR
# -----------------------------
plt.figure()
plt.plot(x_calc, T_calc)
plt.xlabel('x [m]')
plt.ylabel('Temperature [K]')
plt.title('Temperature along PFR')
plt.grid()
plt.show()

# -----------------------------
# Plot Pressure along PFR
# -----------------------------
plt.figure()
plt.plot(x_calc, P_calc/1e5)  # convert to bar
plt.xlabel('x [m]')
plt.ylabel('Pressure [bar]')
plt.title('Pressure along PFR')
plt.grid()
plt.show()

# -----------------------------
# Plot Velocity along PFR
# -----------------------------
plt.figure()
plt.plot(x_calc, vx_calc)
plt.xlabel('x [m]')
plt.ylabel('Velocity [m/s]')
plt.title('Velocity along PFR')
plt.grid()
plt.show()

# -----------------------------
# Plot Mach number along PFR
# -----------------------------
plt.figure()
plt.plot(x_calc, M_calc)
plt.xlabel('x [m]')
plt.ylabel('Mach number')
plt.title('Mach number along PFR')
plt.grid()
plt.show()

# -----------------------------
# Plot selected species
# -----------------------------
species_to_plot = ['CH4','O2','CO2','H2O']
plt.figure()
for s in species_to_plot:
    plt.plot(x_calc, Y_calc[:, gas.species_index(s)], label=s)
plt.xlabel('x [m]')
plt.ylabel('Mass fraction')
plt.title('Species profiles along PFR')
plt.legend()
plt.grid()
plt.show()

