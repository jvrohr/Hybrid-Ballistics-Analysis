import numpy as np
from elements.paraffin import Paraffin

class Grain:
    def __init__(self, material: Paraffin, length: float, internal_diameter: float, external_diameter: float):
        self.material = material                    # Usually a Paraffin object
        self.internal_diameter = internal_diameter    # Grain hole diameter (assumes circular profile) [m]
        self.external_diameter = external_diameter    # Grain external diameter [m]
        self.length = length                        # Grain length [m]

        self.volume_variation = 0
        self.fuel_mass_flow = 0
        self.instant_OF = 0
        self.instant_mass_generation_rate = 0


    def get_port_trans_area(self) -> float:
        return np.pi * (self.internal_diameter / 2) ** 2


    def update_fuel_regression(self, deltat: float, oxidizer_mass_flow: float):
        before_area = self.get_port_trans_area()
        oxidizer_flux = oxidizer_mass_flow / before_area                            # [kg/m^2s]
        
        regretion_rate = 1e-3 * self.material.burn_coefficient * \
            (oxidizer_flux ** self.material.burn_exponent)                          # [m/s]
        
        ## update inner diameter
        self.internal_diameter = self.internal_diameter + deltat * regretion_rate*2 # [m]
        after_area = self.get_port_trans_area()

        self.fuel_mass_flow = np.pi * self.internal_diameter * self.length * \
            regretion_rate * self.material.density                                  # [kg/s]

        self.volume_variation = (before_area - after_area) * self.length

        self.instant_OF = oxidizer_mass_flow / self.fuel_mass_flow # [-]
        self.instant_mass_generation_rate = oxidizer_mass_flow + self.fuel_mass_flow # [kg/s]

