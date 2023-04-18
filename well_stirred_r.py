import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import sys
from scipy.integrate import ode
import warnings
warnings.filterwarnings('ignore')

class well_stirred_reactor():
    def __init__(self,volume = 67.4,res_time = 0.03e-3,inlet_temp = 298, reactor_temp = 1000,heat_loss = 0,init_pressure = ct.one_atm,rxn_mol = 'H2:2,O2:1,N2:4'):
        self.vol = volume   # reactor volume
        self.res_time = res_time # residence time
        self.T_ = inlet_temp     # inlet temperature
        self.Q = heat_loss     # Heat loss in reactor
        self.P = init_pressure  # pressure in reactor
        self.X_ = rxn_mol  #mol Fraction
        self.T = reactor_temp
        
        self.gas_inlet = ct.Solution('Gri30.yaml')  # creating gas object using cantera for hydrogen mechanism
        self.gas_inlet.TP = self.T_,self.P  # setting temperature, pressure and mol fractions into reactor
        self.gas_inlet.set_equivalence_ratio(0.5, 'CH4:1.0', 'O2:1.0, N2:3.76')
        self.species_names = self.gas_inlet.species_names
        self.species_names.insert(0,'T')
        self.Y_init = self.gas_inlet.Y # initial mass fraction of species
        self.h_init = self.gas_inlet.partial_molar_enthalpies/self.gas_inlet.molecular_weights
        
        self.gas = ct.Solution('Gri30.yaml')  # creating gas object using cantera for hydrogen mechanism
        self.X = {key: 1e-10 for key in self.gas.species_names}
        self.X['H2O'] = 1.99999
        self.X['N2'] = 3.9999 *2
        self.X['CO2'] = 1
        
        # X = 'H2:2,O2:1,N2:4'
        self.gas.TPX = self.T,self.P ,self.X  # setting temperature, pressure and mol fractions into reactor
        self.Y_init_r = self.gas.Y
        
    def ReactorOde(self, t, Y):
        T = Y[0]
        YY = Y[1:]
        
        self.gas.TPY = T, self.P, YY
        species_rates = -1/self.res_time * (self.gas.Y - self.Y_init) + self.gas.net_production_rates * self.gas.molecular_weights/self.gas.density
        
        species_h =  self.gas.partial_molar_enthalpies/self.gas.molecular_weights
        first_term = (1/self.res_time) * np.sum( self.Y_init * (self.h_init - species_h))
        second_term = -np.sum(species_h * self.gas.net_production_rates*self.gas.molecular_weights/self.gas.density) - (self.Q/self.gas.density*self.vol )
        temp_rate = (first_term + second_term)/self.gas.cp_mass

        return np.concatenate((np.array([temp_rate]),species_rates), axis = 0)
    
    def solve_ode(self,dt = 1e-7,n_steps = 10**20):
        r = ode(self.ReactorOde).set_integrator('vode', method = 'bdf', with_jacobian=True, rtol=1e-15, atol=1e-15) #, method='bdf', with_jacobian=True)
        Y0 = np.concatenate((np.asarray(self.T)[None], self.Y_init_r), axis = 0)
        r.set_initial_value(Y0, 0)
        
        
        t = []
        ys = []
        t += [r.t]
        ys += [r.y]
        for i in range(n_steps):
            r.integrate(r.t + dt)
            t += [r.t]
            ys += [r.y]
            # if r.y[0] > 1600 and misc == False:
            #     igtime = r.t
            #     misc = True
            #   T_old = r.y[0]
            if r.t > self.res_time:
                # print("reaction complete")
                break
            
            # print(i, t[-1])

        t_ = np.array(t)
        
        ys = np.stack(ys)
        return (t_,ys)
        
           

    
    
        