from elements.nitrousoxide import *

class Injector:
    def __init__(self, dischargeCoefficient: float, holeArea: float, nbHoles: int):
        self.dischargeCoeffient = dischargeCoefficient
        self.holeArea = holeArea
        self.nbHoles = nbHoles
        self.totalInjectorArea = self.nbHoles*holeArea

    # parameter used in the injection model
    def InjectionKParameter(self, fluid: NitrousOxide, downstreamPressure: float) -> float:
        if fluid.GetVaporPressure() - downstreamPressure <= 0:
            return 0
        else:
            return np.sqrt((fluid.pressure - downstreamPressure)/(fluid.GetVaporPressure() - downstreamPressure))

    # Temperature post injection [K]
    def DownstreamTemperature(self, fluid: NitrousOxide, phase: Phase, downstreamPressure: float) -> float:
        gamma = fluid.GetSpecificHeatsRatio(phase)
        return fluid.temperature/((fluid.pressure/downstreamPressure)**((gamma - 1)/gamma))

    # Mass of N2O liquid or gaseous flowing through the injector [kg/s]
    def GetMassFlow(self, downstreamPressure: float, phase: Phase, fluid: NitrousOxide):
        if downstreamPressure >= fluid.pressure:
            return 0
        else:
            if phase == Phase.GAS:
                return self.dischargeCoeffient*self.totalInjectorArea*np.sqrt(fluid.GetVaporDensity(fluid.temperature)*(fluid.pressure - downstreamPressure))
            elif phase == Phase.LIQUID:
                auxVariable = 1/(1 + self.InjectionKParameter(fluid, downstreamPressure))
                part1 = (1 - auxVariable)*np.sqrt(2*fluid.GetLiquidDensity()*(fluid.pressure - downstreamPressure))
                part2 = auxVariable*fluid.GetVaporDensity(self.DownstreamTemperature(fluid, phase, downstreamPressure))
                part3 = np.sqrt(2*fluid.GetSpecificHeatCP(phase)*(fluid.temperature - self.DownstreamTemperature(fluid, phase, downstreamPressure)))
                return self.dischargeCoeffient*self.totalInjectorArea*(part1 + part2*part3)
