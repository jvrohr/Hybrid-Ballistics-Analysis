from elements.grain import *
from elements.nozzle import *
from utilities.environment import Environment
import numpy as np
from utilities.convert import *
from scipy.interpolate import RegularGridInterpolator

class Chamber:
    def __init__(self, pressure: float, diameter: float, pre_combustor_length: float, 
                 post_combustor_length: float, grain: Grain, nozzle: Nozzle):
        self.transversal_area = np.pi * (diameter ** 2) / 4
        self.post_combustor_length = post_combustor_length
        self.pre_combustor_length = pre_combustor_length
        self.pressure = pressure
        self.grain = grain
        self.nozzle = nozzle


    def initialize_chamber(self, environment: Environment):
        self.gas_mass = 0 # self.get_chamber_internal_volume() * simulation_parameters.environment.GetAirDensity()
        self.gamma = environment.gamma_air
        self.MW_comb_gas = environment.MW_air
        self.combustion_temperature = environment.T

        # Carregar dados
        data = np.load("CEA_results.npz")
        gamma = data["gamma"]
        MW_comb_gas = data["MW_comb_gas"]
        combustion_temperature = data["combustion_temperature"]
        Pchamber_vec = data["Pchamber_vec"]
        OF_vec = data["OF_vec"]

        self.interp_temp = RegularGridInterpolator((Pchamber_vec, OF_vec), combustion_temperature)
        self.interp_gamma = RegularGridInterpolator((Pchamber_vec, OF_vec), gamma)
        self.interp_MW = RegularGridInterpolator((Pchamber_vec, OF_vec), MW_comb_gas)


    def get_chamber_internal_volume(self) -> float:
        return self.transversal_area * (self.post_combustor_length + self.pre_combustor_length) + \
            self.grain.length * self.grain.get_port_trans_area()


    def update_chamber_pressure(self, time_step: float, environment: Environment,
                                oxidizer_name: str, fuel_name: str):
        cea_object = cea.CEA_Obj(oxName=oxidizer_name, fuelName=fuel_name)

        before_pressure = 0
        gas_mass = self.gas_mass
        while(abs(before_pressure - self.pressure) > 0.1):
            gas_mass = self.gas_mass
            before_pressure = self.pressure
            chamber_pressure_psi = convert_pa_2_psia(self.pressure)

            # self.MW_comb_gas, self.gamma = cea_object.get_Chamber_MolWt_gamma(chamber_pressure_psi, 
            #                                                                 self.grain.instant_OF) # [g/mol | -]
            # self.combustion_temperature = cea_object.get_Tcomb(chamber_pressure_psi, self.grain.instant_OF)
            ## converting from american units to SI
            # self.combustion_temperature = rankine_2_kelvin(self.combustion_temperature)

            self.MW_comb_gas = self.interp_MW([chamber_pressure_psi, self.grain.instant_OF])[0]
            self.gamma = self.interp_gamma([chamber_pressure_psi, self.grain.instant_OF])[0]
            self.combustion_temperature = self.interp_temp([chamber_pressure_psi, self.grain.instant_OF])[0]
            
            self.MW_comb_gas = g_mol_2_kg(self.MW_comb_gas) # kg/mol

            gas_mass_variation = (self.grain.instant_mass_generation_rate - self.nozzle.mass_flow_nozzle) * time_step
            gas_mass = gas_mass + gas_mass_variation
            density_comb_gas = gas_mass / self.get_chamber_internal_volume()

            self.pressure = (environment.R / self.MW_comb_gas) * \
                self.combustion_temperature * density_comb_gas

        self.gas_mass = gas_mass
