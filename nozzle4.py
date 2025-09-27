# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 19:09:48 2025

@author: trefr
"""
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

T_0 = 1500.0  # inlet temperature [K]
pressure = ct.one_atm  # constant pressure [Pa]
composition_0 = 'H2:2, O2:1, AR:0.1'
length = 1.5e-7  # *approximate* PFR length [m]
u_0 = .006  # inflow velocity [m/s]
area = 1.e-4  # cross-sectional area [m**2]



# input file containing the reaction mechanism
reaction_mechanism = 'h2o2.yaml'




# Resolution: The PFR will be simulated by 'n_steps' time steps or by a chain
# of 'n_steps' stirred reactors.
n_steps = 2000

# import the gas model and set the initial conditions
gas2 = ct.Solution(reaction_mechanism)
gas2.TPX = T_0, pressure*2, composition_0
mass_flow_rate2 = u_0 * gas2.density * area
dz = length / n_steps
r_vol = area * dz


gas3 = ct.Solution(reaction_mechanism)
T_1 = 2000
composition_1 = 'O2:1'
gas3.TPX = T_1, pressure, composition_1
mdot = u_0 * gas3.density * .05e-4
# create a new reactor

r2 = ct.IdealGasReactor(gas2)
r2.volume = r_vol

# create a reservoir to represent the reactor immediately upstream. Note
# that the gas object is set already to the state of the upstream reactor
upstream = ct.Reservoir(gas2, name='upstream')

# create a reservoir for the reactor to exhaust into. The composition of
# this reservoir is irrelevant.
downstream = ct.Reservoir(gas2, name='downstream')

# The mass flow rate into the reactor will be fixed by using a
# MassFlowController object.
m = ct.MassFlowController(upstream, r2, mdot=mass_flow_rate2)

# We need an outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference.
v = ct.PressureController(r2, downstream, primary=m, K=1e-5)



sim2 = ct.ReactorNet([r2])



# define time, space, and other information vectors
z2 = (np.arange(n_steps) + 1) * dz
t_r2 = np.zeros_like(z2)  # residence time in each reactor
u2 = np.zeros_like(z2)
t2 = np.zeros_like(z2)
m2 = np.zeros_like(z2)
states2 = ct.SolutionArray(r2.thermo)
# iterate through the PFR cells
for n in range(n_steps):
    if n < int(n_steps/4):
        # Set the state of the reservoir to match that of the previous reactor
        gas2.TDY = r2.thermo.TDY
        upstream.syncState()
        # integrate the reactor forward in time until steady state is reached
        sim2.reinitialize()
        sim2.advance_to_steady_state()
        # compute velocity and transform into time
        m2[n] =mass_flow_rate2
        u2[n] = mass_flow_rate2 / area / r2.thermo.density
        t_r2[n] = r2.mass / mass_flow_rate2  # residence time in this reactor
        t2[n] = np.sum(t_r2)
        # write output data
        states2.append(r2.thermo.state)
    if n >= int(n_steps/4):
        
        # Set the state of the reservoir to match that of the previous reactor
        
        # integrate the reactor forward in time until steady state is reached
        if n == int(n_steps/4):
            
            # compute velocity and transform into time
            
            #injector section
            #===============================
            upstream2 = ct.Reservoir(gas2)
            mfc = ct.MassFlowController(upstream2, r2, mdot=mdot)
            injector_valve = ct.PressureController(upstream2, r2, primary=mfc, K=1e-5)
            
           
            # write output data
            gas2.TDY = r2.thermo.TDY
            upstream.syncState()
            states2.append(r2.thermo.state)
            sim2.reinitialize()
            sim2.advance_to_steady_state()
            
            m2[n] =mass_flow_rate2 + mdot
            u2[n] = ((mass_flow_rate2 + mdot )/ area) / r2.thermo.density
            t_r2[n] = r2.mass / (mass_flow_rate2 + mdot )  # residence time in this reactor
            t2[n] = np.sum(t_r2)
            
            
        
        if n != int(n_steps/4):
            # compute velocity and transform into time
            
            # write output data
            gas2.TDY = r2.thermo.TDY
            upstream.syncState()
            states2.append(r2.thermo.state)
            sim2.reinitialize()
            sim2.advance_to_steady_state()
            
            m2[n] =mass_flow_rate2 +mdot
            u2[n] = ((mass_flow_rate2 + mdot) / area) / r2.thermo.density
            t_r2[n] = r2.mass / (mass_flow_rate2 + mdot )  # residence time in this reactor
            t2[n] = np.sum(t_r2)
        
        
        
    
    
plt.figure()

plt.plot(z2, states2.T, label='Reactor Chain')
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [K]')
plt.legend(loc=0)
plt.show()
plt.savefig('pfr_T_z.png')

plt.figure()

plt.plot(t2, states2.X[:, gas2.species_index('H2')], label='Reactor Chain')
plt.xlabel('$t$ [s]')
plt.ylabel('$X_{H_2}$ [-]')
plt.legend(loc=0)
plt.show()
plt.savefig('pfr_XH2_t.png')

plt.figure()

plt.plot(z2,m2 , label='Reactor Chain')
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [m dot]')
plt.legend(loc=0)
plt.show()
plt.savefig('pfr_T_z.png')

plt.figure()
plt.plot(z2, states2.X[:, gas2.species_index('O2')], label='O₂ Mole Fraction')
plt.xlabel('$z$ [m]')
plt.ylabel('$X_{O_2}$ [-]')
plt.title('O₂ Profile Along Reactor')
plt.legend()
plt.grid()
plt.show()


