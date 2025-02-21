from elements.nitrousoxide import *

class Aluminum:
    A = [4.8, 0.00322, 155.239]

    def get_CP(self, temperature: float):
        return (self.A[0] + self.A[1]*temperature)*self.A[2]

class Tank:
    def __init__(self, material: Aluminum, fluid: NitrousOxide, environment: Environment, volume: float, tank_mass: float, loaded_fluid_mass: float):
        self.material = material
        self.fluid = fluid
        self.volume = volume                    # [m^3]
        self.tank_mass = tank_mass                # [kg]
        self.loaded_fluid_mass = loaded_fluid_mass  # [kg]

        self.total_N2O = self.loaded_fluid_mass/fluid.molecular_mass # [mol]

        vapor_pressure_gaseous = fluid.get_vapor_pressure()
        molar_volume_liquid = fluid.get_molar_volume_liquid()

        denominator = (environment.R*fluid.temperature - vapor_pressure_gaseous*molar_volume_liquid)

        self.quantity_gaseous = vapor_pressure_gaseous*(self.volume - molar_volume_liquid*self.total_N2O) / denominator # [mol]
        self.quantity_liquid = (self.total_N2O*environment.R*fluid.temperature - vapor_pressure_gaseous*self.volume) / denominator # [mol]

    # Specific heat of aluminum tank casing [J/(kg*K)]
    def get_CP_tank_material(self, temperature: float):
        return self.material.get_CP(self.material, temperature=temperature)
    
    # Set Quantity of gaseous N2O molecules [mol]
    def add_molar_quantity_gaseous_variation(self, quantity_gaseous_variation):
        self.quantity_gaseous = self.quantity_gaseous + quantity_gaseous_variation
        
    # Set Quantity of liquid N2O molecules [mol]
    def add_molar_quantity_liquid_variation(self, quantity_liquid_variation):
        self.quantity_liquid = self.quantity_liquid + quantity_liquid_variation
     