import numpy as np
from typing import List
import rocketcea.cea_obj as cea
from elements import *

paraffinInputString = """
    fuel paraffin(S)  C 73.0   H 124.0     wt%=100.00
    h,cal=-444694.0724016     t(k)=298.15   rho=1.001766
    """
cea.add_new_fuel('Paraffin', paraffinInputString)

class SimulationParameters:
    def __init__(self, environment: Environment, timeStep: float, totalTime: float):
        self.environment = environment
        self.timeStep = timeStep
        self.totalTime = totalTime

class SolveSimulation:
    resultsDict = {}

    def __init__(self, rocketEngine: RocketEngine, simulationParameters: SimulationParameters):
        self.rocketEngine = rocketEngine
        self.simulationParameters = simulationParameters

        self.deltaTemperature = 0
        self.deltaNGaseous = 0
        self.deltaNLiquid = 0
        self.deltaPressure = 0

    def SaveResultsBlowdown(self, time: float):
        tank = self.rocketEngine.tank
        
        self.resultsDict["Time"].append(time)
        self.resultsDict["Temperature"].append(tank.fluid.temperature)
        self.resultsDict["Quantity Gas"].append(tank.quantityGaseous)
        self.resultsDict["Quantity Liquid"].append(tank.quantityLiquid)

    def RunBlowDown(self):
        deltat = self.simulationParameters.timeStep
        time = self.simulationParameters.totalTime

        # Initialize the results dict structure
        self.resultsDict["Time"] = []
        self.resultsDict["Temperature"] = []
        self.resultsDict["Quantity Gas"] = []
        self.resultsDict["Quantity Liquid"] = []

        tank = self.rocketEngine.tank

        for i in np.arange(0, time, deltat):
            self.UpdateOxidizerBlowdown()
            self.SaveResultsBlowdown(i)

    def UpdateOxidizerBlowdown(self):
        tank = self.rocketEngine.tank
        deltat = self.simulationParameters.timeStep

        self.CalculateTankDerivatives()

        tank.fluid.AddTemperatureVariation(self.deltaTemperature*deltat)
        tank.AddMolarQuantityGaseousVariation(self.deltaNGaseous*deltat)
        tank.AddMolarQuantityLiquidVariation(self.deltaNLiquid*deltat)

    def RunBurn(self):
        deltat = self.simulationParameters.timeStep
        time = self.simulationParameters.totalTime

        tank = self.rocketEngine.tank

        for i in np.arange(0, time, deltat):
            self.UpdateOxidizerBlowdown()

            self.UpdateFuelRegression()
            self.RunCEA()
            self.UpdateChamberPressure()
            self.SaveResultsBurn(i)

    def UpdateChamberPressure(self):
        self.rocketEngine.chamber.pressure = self.simulationParameters.timeStep*self.CalculateChamberPressureDerivative()

    def CalculateChamberPressureDerivative(self):
        nozzle = self.rocketEngine.nozzle
        chamber = self.rocketEngine.chamber
        grain = self.rocketEngine.grain

        nozzle.massFlowNozzle = chamber.pressure*nozzle.throatArea/self.rocketEngine.instantCStar

        massGain = chamber.instantMassGenerationRate - nozzle.massFlowNozzle

        self.deltaPressure = chamber.pressure*(massGain/self.rocketEngine.gasMass - grain.volumeVariation/self.rocketEngine.GetEngineInternalVolume())

    def ConvertPa2Psia(pressurePa):
        return 0.000145038*pressurePa

    def ConvertFts2Ms(velocity):
        return velocity / 3.2808398950131

    def RunCEA(self) -> float:
        chamber = self.rocketEngine.chamber
        nozzle = self.rocketEngine.nozzle

        ceaObject = cea.CEA_Obj(oxName='N2O', fuelName='Paraffin')

        atmosphericPressurePsi = self.ConvertPa2Psia(Environment.atmosphericPressure)
        chamberPressurePsi = self.ConvertPa2Psia(chamber.pressure)

        self.rocketEngine.instantCStar = self.ConvertFts2Ms(ceaObject.get_Cstar(Pc=chamberPressurePsi, \
                                  MR=chamber.instantOF))
        Cf = ceaObject.get_PambCf(Pamb=atmosphericPressurePsi, Pc=chamberPressurePsi, \
                                  MR=chamber.instantOF, eps=nozzle.superAreaRatio)[1]

        # Pe = m_out * Cstar / At
        # Force = Cf * Pe * At
        # I = Force * tstep + I
        # Cstar = cstar(OF)

    def UpdateFuelRegression(self):
        tank = self.rocketEngine.tank
        fluid = tank.fluid
        injector = self.rocketEngine.injector
        chamber = self.rocketEngine.chamber
        grain = self.rocketEngine.grain
        deltat = self.simulationParameters.timeStep

        if(tank.quantityLiquid >= 0):
            phase = Phase.LIQUID
        else:
            phase = Phase.GAS

        oxidizerMassFlow = injector.GetMassFlow(fluid.pressure, chamber.pressure, phase, fluid)                 # [kg/s]
        oxidizerFlux = oxidizerMassFlow / (np.pi * (grain.internalDiameter / 2) ** 2)                           # [kg/m^2s]
        regretionRate = 1e-3 * grain.material.burnCoefficient * (oxidizerFlux ** grain.material.burnExponent)   # [m/s]
        grain.AddInternalDiameterVariation(deltat * regretionRate*2)                                            # [m]
        fuelMassFlow = np.pi * grain.internalDiameter * grain.length * regretionRate * grain.material.density   # [kg/s]

        grain.volumeVariation = np.pi*((grain.internalDiameter/2 + regretionRate*deltat)**2 - (grain.internalDiameter/2)**2)

        chamber.instantOF = oxidizerMassFlow / fuelMassFlow # [-]
        chamber.instantMassGenerationRate = oxidizerMassFlow + fuelMassFlow # [kg/s]
        
    # Returns
    # 0 -> derivative of nitrous oxide temperature
    # 1 -> derivative of the geseous molar quantity
    # 2 -> derivative of the liquid molar quantity
    def CalculateTankDerivatives(self) -> None:
        tank = self.rocketEngine.tank
        fluid = tank.fluid
        injector = self.rocketEngine.injector
        chamber = self.rocketEngine.chamber

        if(tank.quantityLiquid >= 0):
            phase = Phase.LIQUID
        else:
            phase = Phase.GAS

        P = tank.quantityGaseous * Environment.R * fluid.temperature / \
            (tank.volume - tank.quantityLiquid * fluid.GetMolarVolumeLiquid())          # [J/m^3]
        a = tank.tankMass * tank.GetSpecificHeatCPofTankMaterial(fluid.temperature) + \
            tank.quantityGaseous * fluid.GetSpecificHeatCPGaseous() + \
                tank.quantityLiquid * fluid.GetSpecificHeatCPLiquid()                   # [J/K]
        b = P * fluid.GetMolarVolumeLiquid()                                            # [J/mol]
        e = - fluid.GetVaporizationHeat() + Environment.R * fluid.temperature           # [J/mol]
        f = - injector.GetMassFlow(fluid.pressure, chamber.pressure, phase, fluid) / \
            fluid.molecularMass                                                         # [mol/s]
        j = - fluid.GetMolarVolumeLiquid() * fluid.GetVaporPressure()                   # [J/mol]
        k = (tank.volume - tank.quantityLiquid * fluid.GetMolarVolumeLiquid()) * \
            fluid.GetVaporPressureDerivTemp()                                           # [J/K]
        m = Environment.R * fluid.temperature                                           # [J/mol]
        q = Environment.R * tank.quantityGaseous                                        # [J/K]

        self.deltaNGaseous = (-f * (-j * a + (q - k) * b)) / (a * (m + j) + (q - k) * (e - b))      # [mol/s]
        self.deltaNLiquid = (-self.deltaNGaseous * (m * a + (q - k) * e)) / (-j * a + (q - k) * b)  # [mol/s]
        self.deltaTemperature = (b * self.deltaNLiquid + e * self.deltaNGaseous) / a                # [K/s]

def SaveResultsBurn(self):
    pass
