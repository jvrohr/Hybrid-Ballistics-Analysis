import numpy as np
from typing import List
import rocketcea.cea_obj as cea
from elements import *

class SimulationParameters:
    def __init__(self, environment: Environment, timeStep: float):
        self.environment = environment
        self.timeStep = timeStep

class SolveSimulation:
    resultsDict = {}

    def __init__(self, rocketEngine: RocketEngine, simulationParameters: SimulationParameters):
        self.rocketEngine = rocketEngine
        self.simulationParameters = simulationParameters

        self.deltaTemperature = 0
        self.deltaNGaseous = 0
        self.deltaNLiquid = 0

    def SaveResultsBlowdown(self, time: float):
        tank = self.rocketEngine.tank
        
        self.resultsDict["Time"].append(time)
        self.resultsDict["Temperature"].append(tank.fluid.temperature)
        self.resultsDict["Quantity Gas"].append(tank.quantityGaseous)
        self.resultsDict["Quantity Liquid"].append(tank.quantityLiquid)

    def RunBlowDown(self):
        deltat = 0.01
        time = 5

        # Initialize the results dict structure
        self.resultsDict["Time"] = []
        self.resultsDict["Temperature"] = []
        self.resultsDict["Quantity Gas"] = []
        self.resultsDict["Quantity Liquid"] = []

        tank = self.rocketEngine.tank

        for i in np.arange(0, time, deltat):
            self.UpdateBlowdown()
            self.SaveResultsBlowdown(i)

    def UpdateBlowdown(self):
        tank = self.rocketEngine.tank
        deltat = self.simulationParameters.timeStep

        self.CalculateTankDerivatives()

        tank.fluid.AddTemperatureVariation(self.deltaTemperature*deltat)
        tank.AddMolarQuantityGaseousVariation(self.deltaNGaseous*deltat)
        tank.AddMolarQuantityLiquidVariation(self.deltaNLiquid*deltat)

    def RunBurn(self):
        deltat = 0.01
        time = 5

        tank = self.rocketEngine.tank

        for i in np.arange(0, time, deltat):
            self.UpdateBlowdown()

            self.UpdateBurn()
            self.RunCEA()
            self.SaveResultsBurn(i)

    def AddParaffin():
        paraffinInputString = """
        fuel paraffin(S)  C 73.0   H 124.0     wt%=100.00
        h,cal=-444694.0724016     t(k)=298.15   rho=1.001766
        """
        cea.add_new_fuel('Paraffin', paraffinInputString)

    def ConvertPa2Psia(pressurePa):
        return 0.000145038*pressurePa

    def RunCEA(self) -> float:
        self.AddParaffin()

        ceaObject = cea.CEA_Obj(oxName='N2O', fuelName='Paraffin')

        atmosphericPressurePsi = self.ConvertPa2Psia(Environment.atmosphericPressure)
        chamberPressurePsi = self.ConvertPa2Psia(self.rocketEngine.chamber.pressure)

        Cf = ceaObject.get_PambCf(Pamb=atmosphericPressurePsi, Pc=chamberPressurePsi, MR=1, eps=self.rocketEngine.nozzle.superAreaRatio)[1]

        # Pe = m_out * Cstar / At
        # Force = Cf * Pe * At
        # I = Force * tstep + I
        # Cstar = cstar(OF)

    def UpdateBurn(self) -> List[float]:
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
        grain.AddInternalDiameterVariation(deltat * 2 * regretionRate)                                   # [m]
        fuelMassFlow = np.pi * grain.internalDiameter * grain.length * regretionRate * grain.material.density   # [kg/s]

        return oxidizerMassFlow + fuelMassFlow, oxidizerMassFlow / fuelMassFlow # Total Mass Flow, OF ratio [kg/s, -]
        
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
