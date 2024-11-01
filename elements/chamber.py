from elements.grain import *
import numpy as np

class Chamber:
    def __init__(self, pressure: float, transversalArea: float, preCombustorLength: float, postCombustorLength: float):
        self.transversalArea = transversalArea
        self.postCombustorLength = postCombustorLength
        self.preCombustorLength = preCombustorLength
        self.pressure = pressure

        self.instantMassGenerationRate = 0
        self.instantOF = 0

    def AddChamberPressureVariation(self, pressureVariation):
        self.pressure = self.pressure + pressureVariation

    def GetChamberInternalVolume(self, grain: Grain) -> float:
        return self.transversalArea * (self.postCombustorLength + self.preCombustorLength) + grain.length * np.pi * (grain.internalDiameter ** 2) / 4
