import numpy as np
from elements.simparams import SimulationParameters
from scipy.optimize import fsolve

class Nozzle:
    def __init__(self, entry_diameter: float, throat_diameter: float,
                 convergent_angle: float, super_area_ratio: float):
        self.throat_area = np.pi * (throat_diameter / 2) ** 2
        self.convergent_angle = convergent_angle
        self.entry_area = np.pi * (entry_diameter / 2) ** 2

        self.super_area_ratio = super_area_ratio
        self.exit_area = self.super_area_ratio * self.throat_area


    def initialize_nozzle(self, simulation_parameters: SimulationParameters):
        self.mass_flow_nozzle = 0
        self.exaust_pressure = simulation_parameters.environment.atmospheric_pressure


    def update_exaust(self, gamma: float, chamber_pressure: float, MW_comb_gas: float, 
                      chamber_temperature: float, simulation_parameters: SimulationParameters):
        
        R_specific = simulation_parameters.environment.R / MW_comb_gas

        self.mass_flow_nozzle = self.throat_area * chamber_pressure * \
            np.sqrt(gamma / (R_specific * chamber_temperature)) * \
                (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1))) 

        func = lambda P_exit: abs(self.compute_expansion_ratio(P_exit, chamber_pressure, gamma)) - self.super_area_ratio

        P_exit_initial_guess = chamber_pressure / 2
        self.exaust_pressure = fsolve(func, P_exit_initial_guess, xtol=1e-6, factor=0.1, maxfev=1000)[0]

    def compute_expansion_ratio(self, P_exit, P_chamber, gamma):
        P_exit = np.clip(P_exit, 100, P_chamber * 0.99) # clip to avoid fsolve negative/zero values
        
        return 1 / ((((gamma + 1) / 2) ** (1 / (gamma - 1))) * \
            ((P_exit / P_chamber) ** (1 / gamma)) * \
                np.sqrt(((gamma + 1) / (gamma - 1)) * (1 - (P_exit / P_chamber) ** ((gamma - 1) / gamma))))

