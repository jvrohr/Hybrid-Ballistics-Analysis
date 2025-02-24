from utilities.environment import Environment

class SimulationParameters:
    def __init__(self, environment: Environment, time_step: float, total_time: float):
        self.environment = environment
        self.time_step = time_step
        self.total_time = total_time