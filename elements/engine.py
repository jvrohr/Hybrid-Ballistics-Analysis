from elements.chamber import *
from elements.simparams import SimulationParameters
from elements.feedingsystem import FeedingSystem
from elements.simobject import SimulationObject

class RocketEngine(SimulationObject):
    def __init__(self, feeding_system: FeedingSystem, chamber: Chamber, burn_effectiveness: float):
        self.feeding_system = feeding_system
        self.injector = feeding_system.injector
        self.tank = feeding_system.tank
        self.chamber = chamber
        self.burn_effectiveness = burn_effectiveness

        self.instant_isp = 0
        self.total_isp = 0
        self.impulse = 0
        self.thrust = 0
        self.thrust_coefficient = 0
        self.consumed_oxidizer = 0
        self.consumed_fuel = 0


    def initialize(self, simulation_parameters: SimulationParameters):
        super().initialize(simulation_parameters)
        self.chamber.initialize_chamber(self.environment)
        self.chamber.grain.initialize_grain()
        self.chamber.nozzle.initialize_nozzle(self.environment)
        self.feeding_system.initialize(simulation_parameters)


    def update_thrust(self):
        deltat = self.time_step
        P_ambient = self.environment.atmospheric_pressure
        g0 = self.environment.gravity_acceleration
        gamma = self.chamber.gamma
        P_exit = self.chamber.nozzle.exaust_pressure
        P_chamber = self.chamber.pressure
        epsilon = self.chamber.nozzle.exit_area

        self.thrust_coefficient = np.sqrt(((2 * gamma ** 2) / (gamma - 1)) * \
                                          ((2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1))) * \
                                            (1 - (P_exit / P_chamber) ** ((gamma - 1) / gamma))) + \
                                                (P_exit - P_ambient) * epsilon / P_chamber
        
        self.thrust = self.thrust_coefficient * P_chamber * \
            self.chamber.nozzle.throat_area * self.burn_effectiveness

        self.impulse = self.thrust * deltat
        self.consumed_fuel = self.chamber.grain.fuel_mass_flow * deltat
        self.consumed_oxidizer = self.injector.oxidizer_mass_flow * deltat
        
        self.instant_isp = self.thrust / \
            ((self.consumed_fuel + self.consumed_oxidizer) * g0)
        self.total_isp = self.impulse / (self.consumed_fuel + self.consumed_oxidizer)


    def run(self, time: float, results_dict: dict):
        self.injector.update_mass_flow(self.chamber.pressure, self.tank.fluid, 
                                       self.feeding_system.pressure_drop)
        
        self.tank.update_oxidizer_blowdown(self.time_step, self.injector.oxidizer_mass_flow)

        self.chamber.grain.update_fuel_regression(self.time_step, self.injector.oxidizer_mass_flow)

        self.chamber.update_chamber_pressure(self.time_step, self.environment, 
                                             self.tank.fluid.name, self.chamber.grain.material.name)

        self.chamber.nozzle.update_exaust(self.chamber.gamma, self.chamber.pressure, 
                                          self.chamber.MW_comb_gas, self.chamber.combustion_temperature,
                                          self.environment)

        self.update_thrust()

        self.save_results(time, results_dict)


    def save_results(self, time, results_dict: dict):
        results_dict["Time"].append(time)
        results_dict["Thrust"].append(self.thrust)
        results_dict["Pressure Chamber"].append(self.chamber.pressure)
        results_dict["Pressure Tank"].append(self.tank.fluid.pressure)
        results_dict["O/F"].append(self.chamber.grain.instant_OF)
