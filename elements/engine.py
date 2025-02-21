from elements.injector import *
from elements.tank import *
from elements.nozzle import *
from elements.chamber import *

class RocketEngine:
    burn_effectiveness  = 0.85 # Efficiency with which N2O and paraffin mixture will burn

    def __init__(self, injector: Injector, tank: Tank, nozzle: Nozzle, grain: Grain, chamber: Chamber):
        self.injector = injector
        self.tank = tank
        self.nozzle = nozzle
        self.grain = grain
        self.chamber = chamber

        self.instant_m_dot_out = 0
        self.instant_c_star = 0
        self.instant_isp = 0
        self.I = 0
        self.thrust = 0
        self.gas_mass = 0

    # Engine internal volume for burn gases in the combustion chamber [m^3]
    def get_engine_internal_volume(self):
        return self.chamber.get_chamber_internal_volume(self.grain) + self.nozzle.get_convergent_section_internal_volume()