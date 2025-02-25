from materials.oxidizers import *
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
    def update_mass_flow(self, downstream_pressure: float, fluid: Oxidizer, pressure_drop: float):
        ## account for line pressure drop
        upstream_pressure = fluid.pressure - pressure_drop
        if downstream_pressure >= upstream_pressure:
            self.oxidizer_mass_flow = 0
        else:
            if fluid.phase == Phase.GAS:
                gamma = fluid.get_specific_heats_ratio(fluid.phase)
                critical_pressure = downstream_pressure / ((2 / (gamma + 1)) ** (gamma / (gamma - 1)))
                if(upstream_pressure >= critical_pressure):
                    self.oxidizer_mass_flow = self.discharge_coefficient * self.total_injector_area * \
                        np.sqrt(gamma * fluid.density * upstream_pressure * \
                                ((2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1))))
                else:
                    self.oxidizer_mass_flow = self.discharge_coefficient * self.total_injector_area * \
                        fluid.density * np.sqrt(2 * fluid.get_CP(fluid.phase) * fluid.temperature * \
                                                (((downstream_pressure / upstream_pressure) ** (2 / gamma)) - \
                                                 ((downstream_pressure / upstream_pressure) ** ((gamma + 1) / gamma))))
            
            elif fluid.phase == Phase.LIQUID:
                deltaP = upstream_pressure - downstream_pressure

                mdot_spi = self.discharge_coefficient * self.total_injector_area * \
                    np.sqrt(2 * fluid.density * deltaP)

                mdot_hem = self.find_HEM_mass_flow_rate(upstream_pressure, downstream_pressure)

                self.oxidizer_mass_flow = (mdot_hem + mdot_spi) / 2
            
