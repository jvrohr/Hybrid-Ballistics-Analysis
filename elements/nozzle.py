import numpy as np

class Nozzle:
    def __init__(self, entryArea: float, exitArea: float, throatArea: float, convergentAngle: float, dischargeCoefficient: float):
        self.exitArea = exitArea
        self.throatArea = throatArea
        self.convergentAngle = convergentAngle
        self.entryArea = entryArea
        self.dischargeCoefficient = dischargeCoefficient

        self.superAreaRatio = self.exitArea/self.throatArea
        self.massFlowNozzle = 0

    def GetConvergentSectionInternalVolume(self) -> float:
        chamberRadius = np.sqrt(self.entryArea/np.pi)
        throatRadius = np.sqrt(self.throatArea/np.pi)
        return np.pi * (chamberRadius - throatRadius) * (chamberRadius**2 + chamberRadius * throatRadius + throatRadius**2) / (3 * np.tan(self.convergentAngle))
