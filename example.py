from HBA import *
import numpy as np

envi = Environment(298)
prop = NitrousOxide(300)
injec = Injector(0.8, np.pi*0.001**2, 20)
tank = Tank(Aluminum, prop, envi, 0.0046, 2, 2.57)
noz = Nozzle(0.1, 0.1, 0.002, 45*np.pi/180, 0.9)
grain = Grain(Paraffin(), 0.5, 0.01, 0.5)
chamber = Chamber(envi.atmosphericPressure, 0.5**2*np.pi, 0.2, 0.2)

engine = RocketEngine(injec, tank, noz, grain, chamber)

simParams = SimulationParameters(envi, 0.01, 5)
sim = SolveSimulation(engine, simParams)

sim.Run("blowdown")
sim.plot.PlotResultsBlowdown()

sim.Run("burn")
sim.plot.PlotResultsBurn()
