
class Material:
    def __init__(self, name: str, density: float):
        self.name = name
        self.density = density

class Oxidizer(Material):
    def __init__(self, name, density, temperature, pressure, phase, molecular_mass, Z):
        super().__init__(name, density)
        self.temperature = temperature
        self.pressure = pressure
        self.phase = phase
        self.molecular_mass = molecular_mass
        self.Z = Z

class Fuel(Material):
    def __init__(self, name, density, burn_coefficient, burn_exponent, cea_input_string):
        super().__init__(name, density)
        self.burn_coefficient = burn_coefficient
        self.burn_exponent = burn_exponent
        self.cea_input_string = cea_input_string

# class Gas(Material):
#     def __init__(self, name, density, molecular_mass, critical_temp, critical_pressure):
#         super().__init__(name, density)
#         self.molecular_mass = molecular