import numpy as np
from typing import List
import rocketcea.cea_obj as cea
from elements import *
import matplotlib.pyplot as plt

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

        # Initialize results dict structure
        self.resultsDict["Time"] = []
        self.resultsDict["Thrust"] = []
        self.resultsDict["Pc"] = []
        self.resultsDict["Pt"] = []

        for i in np.arange(0, time, deltat):
            self.UpdateOxidizerBlowdown()

            self.UpdateFuelRegression()
            self.RunCEA()
            self.UpdateChamberPressure()
            self.SaveResultsBurn(i)

    def UpdateChamberPressure(self):
        self.CalculateChamberPressureDerivative()

        deltat = self.simulationParameters.timeStep
        self.rocketEngine.chamber.AddChamberPressureVariation(deltat*self.deltaPressure)

    def CalculateChamberPressureDerivative(self):
        nozzle = self.rocketEngine.nozzle
        chamber = self.rocketEngine.chamber
        grain = self.rocketEngine.grain
        environment = self.simulationParameters.environment

        if(chamber.pressure - environment.atmosphericPressure > 0):
            nozzle.massFlowNozzle = nozzle.dischargeCoefficient*chamber.pressure*nozzle.throatArea/self.rocketEngine.instantCStar
        else:
            nozzle.massFlowNozzle = 0

        massGain = chamber.instantMassGenerationRate - nozzle.massFlowNozzle

        self.deltaPressure = chamber.pressure*(massGain/self.rocketEngine.gasMass - grain.volumeVariation/self.rocketEngine.GetEngineInternalVolume())

    def ConvertPa2Psia(self, pressurePa: float):
        return 0.000145038*pressurePa

    def ConvertFts2Ms(self, velocity: float):
        return velocity / 3.2808398950131
    
    def ConvertLbmFt32Kgm3(self, density: float):
        return density*16.01846

    def RunCEA(self) -> float:
        chamber = self.rocketEngine.chamber
        nozzle = self.rocketEngine.nozzle

        ceaObject = cea.CEA_Obj(oxName='N2O', fuelName='Paraffin')

        atmosphericPressurePsi = self.ConvertPa2Psia(Environment.atmosphericPressure)
        chamberPressurePsi = self.ConvertPa2Psia(chamber.pressure)

        instantCStar =ceaObject.get_Cstar(Pc=chamberPressurePsi, \
                                  MR=chamber.instantOF)
        self.rocketEngine.instantCStar = self.ConvertFts2Ms(instantCStar)

        Cf = ceaObject.get_PambCf(Pamb=atmosphericPressurePsi, Pc=chamberPressurePsi, \
                                  MR=chamber.instantOF, eps=nozzle.superAreaRatio)[1]

        combustionGasDensity = ceaObject.get_Densities(\
            Pc=chamberPressurePsi, MR=chamber.instantOF, eps=nozzle.superAreaRatio)[0]
        combustionGasDensity = self.ConvertLbmFt32Kgm3(combustionGasDensity)

        self.rocketEngine.gasMass = self.rocketEngine.GetEngineInternalVolume()*combustionGasDensity
        self.rocketEngine.thrust = Cf*self.rocketEngine.chamber.pressure*self.rocketEngine.nozzle.throatArea

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

        oxidizerMassFlow = injector.GetMassFlow(chamber.pressure, phase, fluid)                                 # [kg/s]
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
        environment = self.simulationParameters.environment

        if(tank.quantityLiquid > 0):
            phase = Phase.LIQUID
        else:
            tank.quantityLiquid = 0
            fluid.Z = fluid.pressure * tank.volume / (fluid.temperature * environment.R * tank.quantityGaseous)
            phase = Phase.GAS

        fluid.pressure = fluid.Z * tank.quantityGaseous * environment.R * fluid.temperature / \
            (tank.volume - tank.quantityLiquid * fluid.GetMolarVolumeLiquid())          # [Pa]
        a = tank.tankMass * tank.GetSpecificHeatCPofTankMaterial(fluid.temperature) + \
            tank.quantityGaseous * fluid.GetSpecificHeatCPGaseous() + \
                tank.quantityLiquid * fluid.GetSpecificHeatCPLiquid()                   # [J/K]
        b = fluid.pressure * fluid.GetMolarVolumeLiquid()                               # [J/mol]
        e = - fluid.GetVaporizationHeat() + environment.R * fluid.temperature           # [J/mol]
        f = - injector.GetMassFlow(chamber.pressure, phase, fluid) / \
            fluid.molecularMass                                                         # [mol/s]
        j = - fluid.GetMolarVolumeLiquid() * fluid.GetVaporPressure()                   # [J/mol]
        k = (tank.volume - tank.quantityLiquid * fluid.GetMolarVolumeLiquid()) * \
            fluid.GetVaporPressureDerivTemp()                                           # [J/K]
        m = environment.R * fluid.temperature                                           # [J/mol]
        q = environment.R * tank.quantityGaseous                                        # [J/K]

        if(phase == Phase.LIQUID):
            self.deltaNGaseous = (-f * (-j * a + (q - k) * b)) / (a * (m + j) + (q - k) * (e - b))      # [mol/s]
            self.deltaNLiquid = (-self.deltaNGaseous * (m * a + (q - k) * e)) / (-j * a + (q - k) * b)  # [mol/s] 
        else:
            self.deltaNGaseous = -f         # [mol/s]
            self.deltaNLiquid = 0           # [mol/s]
        
        self.deltaTemperature = (b * self.deltaNLiquid + e * self.deltaNGaseous) / a

    def SaveResultsBurn(self, time):
        self.resultsDict["Time"].append(time)
        self.resultsDict["Thrust"].append(self.rocketEngine.thrust)
        self.resultsDict["Pc"].append(self.rocketEngine.chamber.pressure)
        self.resultsDict["Pt"].append(self.rocketEngine.tank.fluid.pressure)

    def PlotResultsBurn(self):
        time = self.resultsDict["Time"]
        self.PlotData(time, self.resultsDict["Thrust"], "Thrust [N]", "")
        self.PlotData(time, self.resultsDict["Pc"], "Chamber Pressure [Pa]", "")
        self.PlotData(time, self.resultsDict["Pt"], "Tank Pressure [Pa]", "")

    def PlotData(self, X, Y, ylabel, title, xlabel="Time [s]", show=True):
        plt.plot(X, Y)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid()
        plt.title(title)
        if show:
            plt.show()