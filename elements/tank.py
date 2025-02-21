from elements.nitrousoxide import *

class Aluminum:
    A = [4.8, 0.00322, 155.239]

    def get_CP(self, temperature: float):
        return (self.A[0] + self.A[1]*temperature)*self.A[2]

class Tank:
    def __init__(self, material: Aluminum, fluid: NitrousOxide, environment: Environment, volume: float, tank_mass: float, loaded_fluid_mass: float):
        self.material = material
        self.fluid = fluid
        self.environment = environment
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

    def calculate_tank_derivatives(self, oxidizer_mass_flow: float) -> None:
        fluid = self.fluid

        if(self.quantity_liquid > 0):
            fluid.phase = Phase.LIQUID
        else:
            self.quantity_liquid = 0
            fluid.Z = fluid.pressure * self.volume / (fluid.temperature * self.environment.R * self.quantity_gaseous)
            fluid.phase = Phase.GAS

        fluid.pressure = fluid.Z * self.quantity_gaseous * self.environment.R * fluid.temperature / \
            (self.volume - self.quantity_liquid * fluid.get_molar_volume_liquid())          # [Pa]
        a = self.tank_mass * self.get_CP_tank_material(fluid.temperature) + \
            self.quantity_gaseous * fluid.get_CP_gaseous() + \
                self.quantity_liquid * fluid.get_CP_liquid()                   # [J/K]
        b = fluid.pressure * fluid.get_molar_volume_liquid()                               # [J/mol]
        e = - fluid.get_vaporization_eat() + self.environment.R * fluid.temperature           # [J/mol]
        f = - oxidizer_mass_flow / fluid.molecular_mass                                                         # [mol/s]
        j = - fluid.get_molar_volume_liquid() * fluid.get_vapor_pressure()                   # [J/mol]
        k = (self.volume - self.quantity_liquid * fluid.get_molar_volume_liquid()) * \
            fluid.get_vapor_pressure_deriv_temp()                                           # [J/K]
        m = self.environment.R * fluid.temperature                                           # [J/mol]
        q = self.environment.R * self.quantity_gaseous                                        # [J/K]

        if(fluid.phase == Phase.LIQUID):
            self.delta_N_gaseous = (-f * (-j * a + (q - k) * b)) / (a * (m + j) + (q - k) * (e - b))      # [mol/s]
            self.delta_N_liquid = (-self.delta_N_gaseous * (m * a + (q - k) * e)) / (-j * a + (q - k) * b)  # [mol/s] 
        else:
            self.delta_N_gaseous = -f         # [mol/s]
            self.delta_N_liquid = 0           # [mol/s]
        
        self.delta_temperature = (b * self.delta_N_liquid + e * self.delta_N_gaseous) / a

    def update_oxidizer_blowdown(self, deltat: float, oxidizer_mass_flow: float) -> None:
        self.calculate_tank_derivatives(oxidizer_mass_flow)

        self.fluid.add_temperature_variation(self.delta_temperature * deltat)
        self.add_molar_quantity_gaseous_variation(self.delta_N_gaseous * deltat)
        self.add_molar_quantity_liquid_variation(self.delta_N_liquid * deltat)
     