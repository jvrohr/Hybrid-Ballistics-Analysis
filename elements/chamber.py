from elements.grain import *
import numpy as np

class Chamber:
    def __init__(self, pressure: float, transversal_area: float, pre_combustor_length: float, post_combustor_length: float):
        self.transversal_area = transversal_area
        self.post_combustor_length = post_combustor_length
        self.pre_combustor_length = pre_combustor_length
        self.pressure = pressure

        self.instant_mass_generation_rate = 0
        self.instant_OF = 0

    def add_chamber_pressure_variation(self, pressure_variation):
        self.pressure = self.pressure + pressure_variation

    def get_chamber_internal_volume(self, grain: Grain) -> float:
        return self.transversal_area * (self.post_combustor_length + self.pre_combustor_length) + grain.length * np.pi * (grain.internal_diameter ** 2) / 4
