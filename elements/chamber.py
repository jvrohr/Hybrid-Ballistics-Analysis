from elements.grain import *
from elements.nozzle import *
from elements.simparams import SimulationParameters
import numpy as np
from utilities.convert import *

class Chamber:
    def __init__(self, pressure: float, diameter: float, pre_combustor_length: float, 
                 post_combustor_length: float, grain: Grain, nozzle: Nozzle):
        self.transversal_area = np.pi * (diameter ** 2) / 4
        self.post_combustor_length = post_combustor_length
        self.pre_combustor_length = pre_combustor_length
        self.pressure = pressure
        self.grain = grain
        self.nozzle = nozzle


    def initialize_chamber(self, simulation_parameters: SimulationParameters):
        self.gas_mass = 0 # self.get_chamber_internal_volume() * simulation_parameters.environment.GetAirDensity()
        self.gamma = simulation_parameters.environment.gamma_air
        self.MW_comb_gas = simulation_parameters.environment.MW_air
        self.combustion_temperature = simulation_parameters.environment.T


    def get_chamber_internal_volume(self) -> float:
        return self.transversal_area * (self.post_combustor_length + self.pre_combustor_length) + \
            self.grain.length * self.grain.get_port_trans_area()


    def update_chamber_pressure(self, simulation_parameters: SimulationParameters, 
                                oxidizer_name: str, fuel_name: str):
        deltat = simulation_parameters.time_step
        cea_object = cea.CEA_Obj(oxName=oxidizer_name, fuelName=fuel_name)

        before_pressure = 0
        gas_mass = self.gas_mass
        while(abs(before_pressure - self.pressure) > 0.1):
            gas_mass = self.gas_mass
            before_pressure = self.pressure
            chamber_pressure_psi = convert_pa_2_psia(self.pressure)

            self.MW_comb_gas, self.gamma = cea_object.get_Chamber_MolWt_gamma(chamber_pressure_psi, 
                                                                            self.grain.instant_OF) # [g/mol | -]
            self.combustion_temperature = cea_object.get_Tcomb(chamber_pressure_psi, self.grain.instant_OF)
            
            ## converting from american units to SI
            self.combustion_temperature = rankine_2_kelvin(self.combustion_temperature)
            self.MW_comb_gas = g_mol_2_kg(self.MW_comb_gas) # kg/mol

            gas_mass_variation = (self.grain.instant_mass_generation_rate - self.nozzle.mass_flow_nozzle) * deltat
            gas_mass = gas_mass + gas_mass_variation
            density_comb_gas = gas_mass / self.get_chamber_internal_volume()

            self.pressure = (simulation_parameters.environment.R / self.MW_comb_gas) * \
                self.combustion_temperature * density_comb_gas

        self.gas_mass = gas_mass
