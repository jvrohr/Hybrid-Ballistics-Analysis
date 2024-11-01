from elements.injector import *
from elements.tank import *
from elements.nozzle import *
from elements.chamber import *

class RocketEngine:
    burnEffectiveness  = 0.85 # Efficiency with which N2O and paraffin mixture will burn

    def __init__(self, injector: Injector, tank: Tank, nozzle: Nozzle, grain: Grain, chamber: Chamber):
        self.injector = injector
        self.tank = tank
        self.nozzle = nozzle
        self.grain = grain
        self.chamber = chamber

        self.instantMDotOut = 0
        self.instantCStar = 0
        self.instantIsp = 0
        self.I = 0
        self.thrust = 0
        self.gasMass = 0

    # Engine internal volume for burn gases in the combustion chamber [m^3]
    def GetEngineInternalVolume(self):
        return self.chamber.GetChamberInternalVolume(self.grain) + self.nozzle.GetConvergentSectionInternalVolume()