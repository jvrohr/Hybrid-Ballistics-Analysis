from elements.injector import *
from elements.tank import *
from elements.chamber import *
from elements.simparams import SimulationParameters

class RocketEngine:
    def __init__(self, injector: Injector, tank: Tank, chamber: Chamber, burn_effectiveness: float):
        self.injector = injector
        self.tank = tank
        self.chamber = chamber
        self.burn_effectiveness = burn_effectiveness

        self.instant_isp = 0
        self.total_isp = 0
        self.impulse = 0
        self.thrust = 0
        self.thrust_coefficient = 0
        self.consumed_oxidizer = 0
        self.consumed_fuel = 0


    def update_thrust(self, sim_parameters: SimulationParameters):
        deltat = sim_parameters.time_step
        P_ambient = sim_parameters.environment.atmospheric_pressure
        g0 = sim_parameters.environment.gravity_acceleration
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

