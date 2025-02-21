import numpy as np

class Nozzle:
    def __init__(self, entry_area: float, exit_area: float, throat_area: float, convergent_angle: float, discharge_coefficient: float):
        self.exit_area = exit_area
        self.throat_area = throat_area
        self.convergent_angle = convergent_angle
        self.entry_area = entry_area
        self.discharge_coefficient = discharge_coefficient

        self.super_area_ratio = self.exit_area/self.throat_area
        self.mass_flow_nozzle = 0

    def get_convergent_section_internal_volume(self) -> float:
        chamber_radius = np.sqrt(self.entry_area/np.pi)
        throat_radius = np.sqrt(self.throat_area/np.pi)
        return np.pi * (chamber_radius - throat_radius) * (chamber_radius**2 + chamber_radius * throat_radius + throat_radius**2) / (3 * np.tan(self.convergent_angle))
