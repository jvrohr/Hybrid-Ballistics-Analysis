from elements.nitrousoxide import *

class Aluminum:
    A = [4.8, 0.00322, 155.239]

    def GetSpecificHeatCP(self, temperature: float):
        return (self.A[0] + self.A[1]*temperature)*self.A[2]

class Tank:
    def __init__(self, material: Aluminum, fluid: NitrousOxide, environment: Environment, volume: float, tankMass: float, loadedFluidMass: float):
        self.material = material
        self.fluid = fluid
        self.volume = volume                    # [m^3]
        self.tankMass = tankMass                # [kg]
        self.loadedFluidMass = loadedFluidMass  # [kg]

        self.totalN2O = self.loadedFluidMass/fluid.molecularMass # [mol]

        vaporPressureGaseous = fluid.GetVaporPressure()
        molarVolumeLiquid = fluid.GetMolarVolumeLiquid()

        denominator = (environment.R*fluid.temperature - vaporPressureGaseous*molarVolumeLiquid)

        self.quantityGaseous = vaporPressureGaseous*(self.volume - molarVolumeLiquid*self.totalN2O) / denominator # [mol]
        self.quantityLiquid = (self.totalN2O*environment.R*fluid.temperature - vaporPressureGaseous*self.volume) / denominator # [mol]

    # Specific heat of aluminum tank casing [J/(kg*K)]
    def GetSpecificHeatCPofTankMaterial(self, temperature: float):
        return self.material.GetSpecificHeatCP(self.material, temperature=temperature)
    
    # Set Quantity of gaseous N2O molecules [mol]
    def AddMolarQuantityGaseousVariation(self, quantityGaseousVariation):
        self.quantityGaseous = self.quantityGaseous + quantityGaseousVariation
        
    # Set Quantity of liquid N2O molecules [mol]
    def AddMolarQuantityLiquidVariation(self, quantityLiquidVariation):
        self.quantityLiquid = self.quantityLiquid + quantityLiquidVariation
     