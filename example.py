from HBA import *
import numpy as np

envi = Environment(293)
prop = NitrousOxide(293)
tank = Tank(prop, envi, 0.00795, 5.8)
injec = Injector(0.66, 0.001, 37)
noz = Nozzle(0.09, 0.027, 45*np.pi/180, 5)
grain = Grain(Paraffin(burn_coefficient = 0.132, burn_exponent = 0.55), 
              0.25, 0.045, 0.086)
chamber = Chamber(envi.atmospheric_pressure, 0.09, 0.075, 0.075, grain, noz)

feed = FeedingSystem(injec, tank, pressure_drop=1e6)

engine = RocketEngine(feed, chamber, 0.9)

simParams = SimulationParameters(envi, 0.002, 12)
sim = SolveSimulation(feed, simParams)

sim.run("blowdown")
sim.plot("blowdown")
