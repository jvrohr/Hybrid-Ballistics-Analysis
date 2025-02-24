from elements.injector import *
from elements.tank import *
from elements.simparams import SimulationParameters
from elements.simobject import SimulationObject

class FeedingSystem(SimulationObject):
    def __init__(self, injector: Injector, tank: Tank, pressure_drop = 0.0):
        self.injector = injector
        self.tank = tank
        self.pressure_drop = pressure_drop


    def initialize(self, simulation_parameters: SimulationParameters):
        super().initialize(simulation_parameters)
        self.consumed_oxidizer = 0


    def update_feed(self):
        self.consumed_oxidizer = self.consumed_oxidizer + \
            self.injector.oxidizer_mass_flow * self.time_step

    
    def run(self, time: float, results_dict: dict):
        self.injector.update_mass_flow(self.environment.atmospheric_pressure, self.tank.fluid, 
                                       self.pressure_drop)
        self.tank.update_oxidizer_blowdown(self.time_step, self.injector.oxidizer_mass_flow)
        self.update_feed()
        self.save_results(time, results_dict)


    def save_results(self, time: float, results_dict: dict):
        tank = self.tank
        injector = self.injector
        
        results_dict["Time"].append(time)
        results_dict["Temperature"].append(tank.fluid.temperature)
        results_dict["Quantity Gas"].append(tank.quantity_gaseous)
        results_dict["Quantity Liquid"].append(tank.quantity_liquid)
        results_dict["Pressure Tank"].append(tank.fluid.pressure)
        results_dict["Oxidizer Mass Flow"].append(injector.oxidizer_mass_flow)
        results_dict["Oxidizer Mass"].append((tank.quantity_gaseous + tank.quantity_liquid) * \
                                             tank.fluid.molecular_mass)
