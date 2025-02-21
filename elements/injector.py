from elements.nitrousoxide import *
from CoolProp.CoolProp import PropsSI
from scipy.optimize import minimize_scalar
import numpy as np

class Injector:
    def __init__(self, dischargeCoefficient: float, holeDiameter: float, nbHoles: int):
        self.dischargeCoeffient = dischargeCoefficient
        self.holeArea = np.pi*(holeDiameter/2)**2
        self.nbHoles = nbHoles
        self.totalInjectorArea = self.nbHoles*self.holeArea

        self.oxidizerMassFlow = 0

    # Temperature post injection [K]
    def DownstreamTemperature(self, fluid: NitrousOxide, phase: Phase, downstreamPressure: float) -> float:
        gamma = fluid.GetSpecificHeatsRatio(phase)
        return fluid.temperature/((fluid.pressure/downstreamPressure)**((gamma - 1)/gamma))

    def computeHEM(self, upstreamPressure, downstreamPressure):
        h1 = PropsSI('H', 'P', upstreamPressure, 'Q', 0, "N2O")
        s1 = PropsSI('S', 'P', upstreamPressure, 'Q', 0, "N2O")
        h2 = PropsSI('H', 'P', downstreamPressure, 'S', s1, "N2O")
        dens = PropsSI('D', 'P', downstreamPressure, 'S', s1, "N2O")

        return self.dischargeCoeffient * self.totalInjectorArea * dens * np.sqrt(2 * (h1 - h2))

    def FindHEMMassFlowRate(self, upstreamPressure, downstreamPressure):
        """Finds the mass flow rate and ensures it does not decrease past the critical point"""

        # Find the downstream pressure that maximizes mass flow rate
        result = minimize_scalar(lambda P_d: -self.computeHEM(upstreamPressure, P_d),
                            bounds=(upstreamPressure * 0.1, upstreamPressure),
                            method='bounded')

        
        critical_pressure = result.x
        max_mass_flow = -result.fun  # Since we minimized -mass_flow, invert it

        # If downstream pressure is to the right of the max, return max mass flow
        if downstreamPressure < critical_pressure:
            return max_mass_flow
        else:
            return self.computeHEM(upstreamPressure, downstreamPressure)

    # Mass of N2O liquid or gaseous flowing through the injector [kg/s]
    def GetMassFlow(self, downstreamPressure: float, phase: Phase, fluid: NitrousOxide):
        if downstreamPressure >= fluid.pressure:
            return 0
        else:
            if phase == Phase.GAS:
                return self.dischargeCoeffient*self.totalInjectorArea*np.sqrt(fluid.GetVaporDensity(fluid.temperature)*(fluid.pressure - downstreamPressure))
            elif phase == Phase.LIQUID:
                deltaP = fluid.pressure - downstreamPressure

                dens = PropsSI('D', 'P', downstreamPressure, 'Q', 0, "N2O")
                mdot_spi = self.dischargeCoeffient * self.totalInjectorArea * \
                    np.sqrt(2*dens*deltaP) # fluid.GetLiquidDensity()

                mdot_hem = self.FindHEMMassFlowRate(fluid.pressure, downstreamPressure)

                return (mdot_hem + mdot_spi) / 2
            
