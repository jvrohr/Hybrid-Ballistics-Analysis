from HBA import *
from elements.paraffin import *
import numpy as np

envi = Environment(293)
prop = NitrousOxide(293)
tank = Tank(Aluminum, prop, envi, 0.00795, 2, 5.8)
injec = Injector(0.66, 0.001, 37)
noz = Nozzle(0.09, 0.027, 45*np.pi/180, 5)
grain = Grain(Paraffin(burn_coefficient = 0.132, burn_exponent = 0.55), 
              0.25, 0.045, 0.086)
chamber = Chamber(envi.atmospheric_pressure, 0.09, 0.075, 0.075, grain, noz)

engine = RocketEngine(injec, tank, chamber, 0.9)

simParams = SimulationParameters(envi, 0.002, 12)
sim = SolveSimulation(engine, simParams)

# sim.Run("blowdown")
# sim.plot.plot_results_blowdown()

sim.Run("burn")
sim.plot.plot_results_burn()
