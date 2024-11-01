class Paraffin:
    density = 800              # Paraffin's density [kg/m^3]
    burnCoefficient = 0.098    # Linear term of the paraffin burn model [mm/s]
    burnExponent = 0.61        # Exponent of the paraffin burn model [-]

class Grain:
    def __init__(self, material: Paraffin, length: float, internalDiameter: float, externalDiameter: float):
        self.material = material                    # Usually a Paraffin object
        self.internalDiameter = internalDiameter    # Grain hole diameter (assumes circular profile) [m]
        self.externalDiameter = externalDiameter    # Grain external diameter [m]
        self.length = length                        # Grain length [m]

        self.volumeVariation = 0

    def AddInternalDiameterVariation(self, internalDiameterVariation: float):
        self.internalDiameter = self.internalDiameter + internalDiameterVariation
