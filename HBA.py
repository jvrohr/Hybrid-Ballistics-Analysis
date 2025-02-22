import numpy as np
from elements.engine import *
from utilities.convert import *
from utilities.plot import *
from elements.simparams import SimulationParameters

paraffin_input_string = """
    fuel paraffin(S)  C 73.0   H 124.0     wt%=100.00
    h,cal=-444694.0724016     t(k)=298.15   rho=1.001766
    """
cea.add_new_fuel('Paraffin', paraffin_input_string)

class SolveSimulation:
    results_dict = {}

    def __init__(self, rocket_engine: RocketEngine, simulation_parameters: SimulationParameters):
        self.rocket_engine = rocket_engine
        self.simulation_parameters = simulation_parameters

        self.delta_temperature = 0
        self.delta_N_gaseous = 0
        self.delta_N_liquid = 0
        self.delta_pressure = 0

        self.rocket_engine.chamber.initialize_chamber(simulation_parameters)
        self.rocket_engine.chamber.nozzle.initialize_nozzle(simulation_parameters)

    def Run(self, option="burn"):
        self.plot = PlotResults()
        deltat = self.simulation_parameters.time_step
        time = self.simulation_parameters.total_time

        if option == "burn":
            iteration_function = self.run_burn_iteration
        elif option == "blowdown":
            iteration_function = self.run_blow_down_iteration
        else:
            raise ValueError(option)

        for i in np.arange(0, time, deltat):
            iteration_function(i)


    def run_blow_down_iteration(self, time: float):
        self.rocket_engine.injector.update_mass_flow(self.rocket_engine.chamber.pressure,
                                                     self.rocket_engine.tank.fluid)
        self.rocket_engine.tank.update_oxidizer_blowdown(self.simulation_parameters.time_step, 
                                                         self.rocket_engine.injector.oxidizer_mass_flow)
        self.save_results_blowdown(time)


    def run_burn_iteration(self, time: float):
        self.rocket_engine.injector.update_mass_flow(self.rocket_engine.chamber.pressure,
                                                     self.rocket_engine.tank.fluid)
        
        self.rocket_engine.tank.update_oxidizer_blowdown(self.simulation_parameters.time_step, 
                                                         self.rocket_engine.injector.oxidizer_mass_flow)

        self.rocket_engine.chamber.grain.update_fuel_regression(self.simulation_parameters.time_step, 
                                                         self.rocket_engine.injector.oxidizer_mass_flow)

        self.rocket_engine.chamber.update_chamber_pressure(self.simulation_parameters)

        self.rocket_engine.chamber.nozzle.update_exaust(self.rocket_engine.chamber.gamma,
                                                self.rocket_engine.chamber.pressure,
                                                self.rocket_engine.chamber.MW_comb_gas,
                                                self.rocket_engine.chamber.combustion_temperature,
                                                self.simulation_parameters)

        self.rocket_engine.update_thrust(self.simulation_parameters)

        self.save_results_burn(time)

    def save_results_blowdown(self, time: float):
        tank = self.rocket_engine.tank
        injector = self.rocket_engine.injector
        
        self.plot.results_dict["Time"].append(time)
        self.plot.results_dict["Temperature"].append(tank.fluid.temperature)
        self.plot.results_dict["Quantity Gas"].append(tank.quantity_gaseous)
        self.plot.results_dict["Quantity Liquid"].append(tank.quantity_liquid)
        self.plot.results_dict["Pressure Tank"].append(tank.fluid.pressure)
        self.plot.results_dict["Oxidizer Mass Flow"].append(injector.oxidizer_mass_flow)
        self.plot.results_dict["Pressure Chamber"].append(self.rocket_engine.chamber.pressure)

    def save_results_burn(self, time):
        self.plot.results_dict["Time"].append(time)
        self.plot.results_dict["Thrust"].append(self.rocket_engine.thrust)
        self.plot.results_dict["Pressure Chamber"].append(self.rocket_engine.chamber.pressure)
        self.plot.results_dict["Pressure Tank"].append(self.rocket_engine.tank.fluid.pressure)
