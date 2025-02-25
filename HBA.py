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
import os
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

        self.check_CEA_file()


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


    def check_CEA_file(self):
        if self.sim_object.__class__.__name__ == "RocketEngine":
            oxidizer_name = self.sim_object.feeding_system.tank.fluid.name
            fuel_name = self.sim_object.chamber.grain.material.name
            name = f"CEA_results_{oxidizer_name}_{fuel_name}.npz"
            
            materials_dir = os.path.join(os.getcwd(), "materials", "CEA_data")
            file_path = os.path.join(materials_dir, name)

            if os.path.exists(file_path):
                self.sim_object.chamber.initialize_CEA_files(file_path)
            else:
                print(f"[INFO]: Creating CEA file for {oxidizer_name} and {fuel_name} in \t{file_path}")
                os.makedirs(os.path.dirname(file_path))
                self.sim_object.chamber.create_tables_CEA(oxidizer_name, fuel_name, file_path)
                self.sim_object.chamber.initialize_CEA_files(file_path)
                print("[INFO]: CEA file saved.")

