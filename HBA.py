from elements.chamber import *
from elements.engine import *
from elements.feedingsystem import *
from elements.grain import *
from elements.injector import *
from elements.nozzle import *
from elements.simobject import *
from elements.simparams import *
from elements.tank import *
from materials.fuels import *
from materials.oxidizers import *
from utilities.convert import *
from utilities.environment import *
from utilities.plot import *
import numpy as np

class SolveSimulation:
    results_dict = {}

    def __init__(self, sim_object: SimulationObject, simulation_parameters: SimulationParameters):
        self.sim_object = sim_object
        self.simulation_parameters = simulation_parameters
        self.ran = "no"

        self.delta_temperature = 0
        self.delta_N_gaseous = 0
        self.delta_N_liquid = 0
        self.delta_pressure = 0

        self.sim_object.initialize(simulation_parameters)


    def run(self, option="burn"):
        self.plot_obj = PlotResults()
        deltat = self.simulation_parameters.time_step
        max_time = self.simulation_parameters.max_time
        self.ran = option

        if self.sim_object.__class__.__name__ == "FeedingSystem":
            iteration_function = self.sim_object.run
        else:
            if option == "burn":
                iteration_function = self.sim_object.run
            elif option == "blowdown":
                iteration_function = self.sim_object.feeding_system.run
            else:
                raise ValueError(option)

        i = 0
        while(self.sim_object.tank.quantity_gaseous > 0 and i * deltat < max_time):
            iteration_function(i * deltat, self.plot_obj.results_dict)
            i = i + 1


    def plot(self, option="burn"):
        if self.sim_object.__class__.__name__ == "FeedingSystem":
            self.plot_obj.plot_results_blowdown()
        else:
            if option == "burn" and self.ran == "burn":
                self.plot_obj.plot_results_burn()
            elif option == "blowdown":
                self.plot_obj.plot_results_blowdown()
            else:
                raise ValueError(option)