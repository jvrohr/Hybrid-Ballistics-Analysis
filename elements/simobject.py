from elements.simparams import SimulationParameters

class SimulationObject:
    def __init__(self):
        pass

    def initialize(self, simulation_parameters: SimulationParameters):
        self.time_step = simulation_parameters.time_step
        self.environment = simulation_parameters.environment

    def run(self, time: float, results_dict: dict):
        pass

    def save_results(self, time: float, results_dict: dict):
        pass
