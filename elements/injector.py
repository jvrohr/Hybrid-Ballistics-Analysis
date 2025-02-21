from elements.nitrousoxide import *
from CoolProp.CoolProp import PropsSI
from scipy.optimize import minimize_scalar
import numpy as np

class Injector:
    def __init__(self, discharge_coefficient: float, hole_diameter: float, nb_holes: int):
        self.discharge_coefficient = discharge_coefficient
        self.hole_area = np.pi*(hole_diameter/2)**2
        self.nb_holes = nb_holes
        self.total_injector_area = self.nb_holes*self.hole_area

        self.oxidizer_mass_flow = 0

    def compute_HEM(self, upstream_pressure, downstream_pressure):
        h1 = PropsSI('H', 'P', upstream_pressure, 'Q', 0, "N2O")
        s1 = PropsSI('S', 'P', upstream_pressure, 'Q', 0, "N2O")
        h2 = PropsSI('H', 'P', downstream_pressure, 'S', s1, "N2O")
        dens = PropsSI('D', 'P', downstream_pressure, 'S', s1, "N2O")

        return self.discharge_coefficient * self.total_injector_area * dens * np.sqrt(2 * (h1 - h2))

    def find_HEM_mass_flow_rate(self, upstream_pressure, downstream_pressure):
        """Finds the mass flow rate and ensures it does not decrease past the critical point"""

        # Find the downstream pressure that maximizes mass flow rate
        result = minimize_scalar(lambda P_d: -self.compute_HEM(upstream_pressure, P_d),
                            bounds=(upstream_pressure * 0.1, upstream_pressure),
                            method='bounded')

        
        critical_pressure = result.x
        max_mass_flow = -result.fun  # Since we minimized -mass_flow, invert it

        # If downstream pressure is to the right of the max, return max mass flow
        if downstream_pressure < critical_pressure:
            return max_mass_flow
        else:
            return self.compute_HEM(upstream_pressure, downstream_pressure)

    # Mass of N2O liquid or gaseous flowing through the injector [kg/s]
    def get_mass_flow(self, downstream_pressure: float, phase: Phase, fluid: NitrousOxide):
        if downstream_pressure >= fluid.pressure:
            return 0
        else:
            if phase == Phase.GAS:
                return self.discharge_coefficient*self.total_injector_area*np.sqrt(fluid.get_vapor_density(fluid.temperature)*(fluid.pressure - downstream_pressure))
            elif phase == Phase.LIQUID:
                deltaP = fluid.pressure - downstream_pressure

                dens = PropsSI('D', 'P', downstream_pressure, 'Q', 0, "N2O")
                mdot_spi = self.discharge_coefficient * self.total_injector_area * \
                    np.sqrt(2*dens*deltaP) # fluid.get_liquid_density()

                mdot_hem = self.find_HEM_mass_flow_rate(fluid.pressure, downstream_pressure)

                return (mdot_hem + mdot_spi) / 2
            
