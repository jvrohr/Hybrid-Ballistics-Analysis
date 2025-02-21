class Paraffin:
    density = 800              # Paraffin's density [kg/m^3]
    burn_coefficient = 0.132 # 0.098    # Linear term of the paraffin burn model [mm/s]
    burn_exponent = 0.55 # 0.61        # Exponent of the paraffin burn model [-]

class Grain:
    def __init__(self, material: Paraffin, length: float, internal_diameter: float, external_diameter: float):
        self.material = material                    # Usually a Paraffin object
        self.internal_diameter = internal_diameter    # Grain hole diameter (assumes circular profile) [m]
        self.external_diameter = external_diameter    # Grain external diameter [m]
        self.length = length                        # Grain length [m]

        self.volume_variation = 0

    def add_internal_diameter_variation(self, internal_diameter_variation: float):
        self.internal_diameter = self.internal_diameter + internal_diameter_variation
