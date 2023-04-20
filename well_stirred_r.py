"""
Modeling of a single well-stirred reactor.
"""

import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import sys
from scipy.integrate import ode
import warnings

warnings.filterwarnings('ignore')

class WellStirredReactor():
    def __init__(self, volume=67.4, res_time=0.03e-3, inlet_temp=298, reactor_temp=2000, heat_loss=0, pressure=ct.one_atm, 
                 inlet_gas='H2:2,O2:1,N2:4', mechanism='h2o2.yaml'):
        """
        Initialize the well-stirred reactor object.
        
        Args:
        - volume (float): Volume of the reactor in m^3.
        - res_time (float): Residence time of the reactor in seconds.
        - inlet_temp (float): Inlet gas temperature in Kelvin.
        - reactor_temp (float): Temperature of the reactor in Kelvin.
        - heat_loss (float): Heat loss in the reactor in Watts.
        - pressure (float): Pressure in the reactor in Pascals.
        - inlet_gas (str): Composition of the gas at the inlet in mole fractions separated by commas.
        - mechanism (str): Mechanism of interest for the gas-phase reactions.
        """
        self.vol = volume   # reactor volume
        self.res_time = res_time # residence time
        self.T_ = inlet_temp     # inlet gas temperature
        self.Q = heat_loss     # Heat loss in reactor
        self.P = pressure  # pressure
        self.T = reactor_temp # Temperature of reactor
        self.X_inlet = inlet_gas # composition of gas at inlet
        self.Mechanism = mechanism # Mechanism of interest
        
        ## Set conditions of gas at inlet
        self.gas_inlet = ct.Solution(self.Mechanism)  # creating gas object using cantera for mechanism
        self.gas_inlet.TPX = self.T_, self.P, self.X_inlet  # setting temperature, pressure and mol fractions of gas into reactor
        self.Y_init = self.gas_inlet.Y # mass fractions of species at gas inlet
        self.h_init = self.gas_inlet.partial_molar_enthalpies/self.gas_inlet.molecular_weights #specific enthalpies of species in gas inlet
        
        self.species_names = self.gas_inlet.species_names # extract names of all species in reactions
        self.species_names.insert(0, 'T') # add Temperature to array of species names
        
        ### set conditions in the reactor, we use equilibrium conditions for reactor
        
        # We create the reactor, and fill it initially with a mixture consisting of the
        # equilibrium products of the inlet mixture. This state corresponds to the state
        # the reactor would reach with infinite residence time, and thus provides a good
        # initial condition from which to reach a steady-state solution on the reacting
        # branch.
        self.gas = ct.Solution('h2o2.yaml')  # creating gas object in reactor
    
        self.X = {key: 1e-16 for key in self.gas.species_names} #set gas in reactor to have no of mol for each species = 1e-10
        # update no of moles for major species products, this approximates final mol fractions of products
        self.X['H2O'] = 1.8  
        self.X['N2'] = 3.76 
        # self.X['CO2'] = 0.97
        
        self.gas.TPX = self.T, self.P ,self.X  # setting temperature, pressure and mol fractions into reactor
        self.Y_init_r = self.gas.Y # storing equilibrium state of reactor
        
    def ReactorOde(self, t, Y):
        
        """
        ODE for the well-stirred reactor model.
        
        Args:
        - t (float): Time in seconds.
        - Y (numpy array): State vector of the reactor including temperature and species mass fractions.
        Returns:
        - dydt (numpy array): Derivatives of the state vector.
        """

        T = Y[0] # Extract temperature from state vector
        Y_species = Y[1:] # Extract species mass fractions from state vector
        
        self.gas.TPY = T, self.P, Y_species # Set the state of the gas at the current time
        
        wdot = self.gas.net_production_rates # Calculate the source term for species mass fractions
        
        # Calculate the derivatives of temperature and species mass fractions according to paper by P. Glarborg
        species_rates = -1/self.res_time * (self.gas.Y - self.Y_init) + wdot * self.gas.molecular_weights/self.gas.density
   
        species_h =  self.gas.partial_molar_enthalpies/self.gas.molecular_weights
        first_term = (1/self.res_time) * np.sum( self.Y_init * (self.h_init - species_h))
        second_term = -np.sum(species_h * wdot*self.gas.molecular_weights/self.gas.density) - (self.Q/self.gas.density*self.vol )
        temp_rate = (first_term + second_term)/self.gas.cp_mass

        return np.concatenate((np.array([temp_rate]),species_rates), axis = 0)
    
    def solve_ode(self,dt = 1e-7):
        
        """
        Integrate the ODE for the well-stirred reactor model up to residence time.
        
        Args:
        - dt: time step for integration
        
        """
        r = ode(self.ReactorOde).set_integrator('vode', method = 'bdf', with_jacobian=True, rtol=1e-15, atol=1e-15) #, method='bdf', with_jacobian=True)
        Y0 = np.concatenate((np.asarray(self.T)[None], self.Y_init_r), axis = 0)
        r.set_initial_value(Y0, 0)
        
        t = []
        ys = []

        while r.t < self.res_time:
            r.integrate(r.t + dt)
            t += [r.t]
            ys += [r.y]

        t_ = np.array(t)
        ys = np.stack(ys)
        
        return (t_,ys)
        
           

    
    
        