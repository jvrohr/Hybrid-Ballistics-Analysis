import numpy as np
import rocketcea.cea_obj as cea
from elements.engine import *
from utilities.convert import *
from utilities.plot import *

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

    def Run(self, option="burn"):
        self.plot = PlotResults()
        deltat = self.simulationParameters.timeStep
        time = self.simulationParameters.totalTime

        if option == "burn":
            iterationFunction = self.RunBurnIteration
        elif option == "blowdown":
            iterationFunction = self.RunBlowDownIteration
        else:
            raise ValueError(option)

        for i in np.arange(0, time, deltat):
            iterationFunction(i)

    def RunBlowDownIteration(self, time: float):
        self.UpdateOxidizerBlowdown()
        self.SaveResultsBlowdown(time)

    def SaveResultsBlowdown(self, time: float):
        tank = self.rocketEngine.tank
        
        self.plot.resultsDict["Time"].append(time)
        self.plot.resultsDict["Temperature"].append(tank.fluid.temperature)
        self.plot.resultsDict["Quantity Gas"].append(tank.quantityGaseous)
        self.plot.resultsDict["Quantity Liquid"].append(tank.quantityLiquid)
        self.plot.resultsDict["Pressure Tank"].append(tank.fluid.pressure)

    def RunBurnIteration(self, time: float):
        self.UpdateOxidizerBlowdown()
        self.UpdateFuelRegression()
        self.RunCEA()
        self.UpdateChamberPressure()
        self.SaveResultsBurn(time)

    def UpdateOxidizerBlowdown(self):
        tank = self.rocketEngine.tank
        deltat = self.simulationParameters.timeStep

        self.CalculateTankDerivatives()

        tank.fluid.AddTemperatureVariation(self.deltaTemperature*deltat)
        tank.AddMolarQuantityGaseousVariation(self.deltaNGaseous*deltat)
        tank.AddMolarQuantityLiquidVariation(self.deltaNLiquid*deltat)

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

    def RunCEA(self) -> float:
        chamber = self.rocketEngine.chamber
        nozzle = self.rocketEngine.nozzle

        ceaObject = cea.CEA_Obj(oxName='N2O', fuelName='Paraffin')

        atmosphericPressurePsi = ConvertPa2Psia(Environment.atmosphericPressure)
        chamberPressurePsi = ConvertPa2Psia(chamber.pressure)

        instantCStar = ceaObject.get_Cstar(Pc=chamberPressurePsi, \
                                  MR=chamber.instantOF)
        self.rocketEngine.instantCStar = ConvertFts2Ms(instantCStar)

        Cf = ceaObject.get_PambCf(Pamb=atmosphericPressurePsi, Pc=chamberPressurePsi, \
                                  MR=chamber.instantOF, eps=nozzle.superAreaRatio)[1]

        combustionGasDensity = ceaObject.get_Densities(\
            Pc=chamberPressurePsi, MR=chamber.instantOF, eps=nozzle.superAreaRatio)[0]
        combustionGasDensity = ConvertLbmFt32Kgm3(combustionGasDensity)

        self.rocketEngine.gasMass = self.rocketEngine.GetEngineInternalVolume()*combustionGasDensity
        self.rocketEngine.thrust = Cf*self.rocketEngine.chamber.pressure*self.rocketEngine.nozzle.throatArea        

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

    def SaveResultsBurn(self, time):
        self.plot.resultsDict["Time"].append(time)
        self.plot.resultsDict["Thrust"].append(self.rocketEngine.thrust)
        self.plot.resultsDict["Pressure Chamber"].append(self.rocketEngine.chamber.pressure)
        self.plot.resultsDict["Pressure Tank"].append(self.rocketEngine.tank.fluid.pressure)
