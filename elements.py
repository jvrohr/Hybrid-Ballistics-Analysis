import numpy as np
import scipy as scp
from enum import Enum
from typing import List

class Phase(Enum):
    LIQUID = "liquid"
    GAS = "gas"

class NitrousOxide:
    # N2O Properties' coefficients
    G = [96.512, -4045, -12.277, 2.886e-5, 2]               # Vapour pressure coefficients [Pa] (valid for 182.3 K - 309.57 K)
    J = [2.3215e7, 0.384, 0, 0];                            # Vaporization heat coefficients [J/kmol]
    D = [132.632, 0.052187, -0.364923, -1.20233, 0.536141]  # Specific heat at constant pressure for gaseous N2O [J/(kmol*K)]
    E = [2.49973, 0.023454, -3.80136, 13.0945, -14.5180]    # Specific heat at constant pressure for liquid N2O [J/(kmol*K)]
    Q = [2.781, 0.27244, 309.57, 0.2882]                    # Specific molar volume coefficients for liquid N2O [m^3/kmol]
    F = [-1.009, -6.28792, 7.50332, -7.90463, -0.629427]    # Vapor Density of N2O [kg/m^3]
    H = [1.72328, -0.83950, 0.51060, -0.10412]              # Liquid N2O Density [kg/m^3]

    molecularMass = 44.013          # N2O molecular mass [kg/kmol]
    criticalPressure = 7251000      # Critical Pressure of N2O [Pa]
    criticalTemperature = 309.57    # Critical Temperature of N2O [K]
    criticalDensity = 452           # Critical Density of N2O [kg/m^3]

    def __init__(self, temperature: float):
        self.temperature = temperature  # Current approximately uniform N2O temperature [K]

        self.pressure = self.GetVaporPressure()
        self.temperatureRatio = self.temperature/self.criticalTemperature

    def SetTemperature(self, temperature: float) -> None:
        self.temperature = temperature
        self.pressure = self.GetVaporPressure()
        self.temperatureRatio = self.temperature/self.criticalTemperature

    def SetPressure(self, pressure: float) -> None:
        self.pressure = pressure

    # Molar volume of liquid N2O [m^3/mol]
    def GetMolarVolumeLiquid(self) -> float:
        return 1e-3*self.Q[1]**(1+(1-self.temperature/self.Q[2])**self.Q[3])/self.Q[0]

    # Vapor Pressure of Saturated N2O [Pa]
    def GetVaporPressure(self) -> float:
        return np.exp(self.G[0] + self.G[1]/self.temperature + self.G[2]*np.log(self.temperature) + self.G[3]*(self.temperature**self.G[4]))

    # Vapor Pressure Temperature Derivative of N2O [Pa/T]
    def GetVaporPressureDerivTemp(self) -> float:
        return (-self.G[1]/(self.temperature**2) + self.G[2]/self.temperature + self.G[3]*self.G[4]*self.temperature**(self.G[4]-1)) * np.exp(self.G[0] + self.G[1]/self.temperature + self.G[2]*np.log(self.temperature) + self.G[3]*self.temperature**self.G[4])

    # Gaseous N2O Specific Heat for constant pressure [J/(mol*K)]
    def GetSpecificHeatCPGaseous(self) -> float:
        auxVariable = 1 - self.temperatureRatio
        return self.molecularMass*(self.D[0]*(1 + self.D[1]*(auxVariable**(-2/3)) + self.D[2]*(auxVariable**(-1/3)) + self.D[3]*(auxVariable**(1/3)) + self.D[4]*(auxVariable**(2/3))))

    # Liquid N2O Specific Heat for constant pressure [J/(mol*K)]
    def GetSpecificHeatCPLiquid(self) -> float:
        auxVariable = 1 - self.temperatureRatio
        return self.molecularMass*(self.E[0]*(1 + self.E[1]/auxVariable + self.E[2]*auxVariable + self.E[3]*(auxVariable**2) + self.E[4]*(auxVariable**3)))

    # Specific Heat for constant pressure [J/(mol*K)]
    def GetSpecificHeatCP(self, phase: Phase) -> float:
        if phase == Phase.LIQUID:
            return self.GetSpecificHeatCPLiquid()
        elif phase == Phase.GAS:
            return self.GetSpecificHeatCPGaseous()
        else:
            raise ValueError(phase)

    # Vaporization Heat of N2O [J/kmol]
    def GetVaporizationHeat(self) -> float:
        return self.J[0]*(1 - self.temperatureRatio)**(self.J[1] + self.J[2]*self.temperatureRatio + self.J[3]*self.temperatureRatio**2)

    # Density of gaseous/vapor N2O [kg/m^3]
    def GetVaporDensity(self, temperature: float) -> float:
        auxVariable = 1/(temperature/self.criticalTemperature) - 1
        return self.criticalDensity*(np.exp(self.F[0]*(auxVariable**(1/3)) + self.F[1]*(auxVariable**(2/3)) + self.F[2]*auxVariable + self.F[3]*(auxVariable**(4/3)) + self.F[4]*(auxVariable**(5/3))))

    # Density of liquid N2O [kg/m^3]
    def GetLiquidDensity(self) -> float:
        auxVariable = 1 - self.temperatureRatio
        return self.criticalDensity*(np.exp(self.H[0]*(auxVariable**(1/3)) + self.H[1]*(auxVariable**(2/3)) + self.H[2]*auxVariable + self.H[3]*(auxVariable**(4/3))))

    # Adiabatic expansion coefficient, specific heats ratio Cv/Cp [-]
    def GetSpecificHeatsRatio(self, phase: Phase):
        return 1/(1 - Environment.R/self.GetSpecificHeatCP(phase))

class Injector:
    def __init__(self, dischargeCoefficient: float, holeArea: float, nbHoles: int):
        self.dischargeCoeffient = dischargeCoefficient
        self.holeArea = holeArea
        self.nbHoles = nbHoles
        self.totalInjectorArea = self.nbHoles*holeArea

    # parameter used in the injection model
    def InjectionKParameter(self, fluid: NitrousOxide, downstreamPressure: float, upstreamPressure: float) -> float:
        if fluid.GetVaporPressure() - upstreamPressure <= 0:
            return 0
        else:
            return np.sqrt((upstreamPressure - downstreamPressure)/(fluid.GetVaporPressure() - downstreamPressure))

    # Temperature post injection [K]
    def DownstremTemperature(self, fluid: NitrousOxide, phase: Phase, upstreamPressure: float, downstreamPressure: float) -> float:
        gamma = fluid.GetSpecificHeatsRatio(phase)
        return fluid.temperature/((upstreamPressure/downstreamPressure)**((gamma - 1)/gamma))

    def GetFluidMassFlow(self, upstreamPressure: float, downstreamPressure: float, phase: Phase, fluid: NitrousOxide):
        if downstreamPressure >= upstreamPressure:
            return 0
        else:
            if phase == Phase.LIQUID:
                return self.dischargeCoeffient*self.totalInjectorArea*np.sqrt(fluid.GetVaporDensity(fluid.temperature)*(upstreamPressure - downstreamPressure))
            elif phase == Phase.GAS:
                auxVariable = 1/(1 + self.InjectionKParameter(fluid, downstreamPressure, upstreamPressure))
                part1 = (1 - auxVariable)*np.sqrt(2*fluid.GetLiquidDensity()*(upstreamPressure - downstreamPressure))
                part2 = auxVariable*fluid.GetVaporDensity(self.DownstremTemperature(fluid, phase, upstreamPressure, downstreamPressure))
                part3 = np.sqrt(2*fluid.GetSpecificHeatCP(phase)*(fluid.temperature - self.DownstremTemperature(fluid, phase, upstreamPressure, downstreamPressure)))
                return self.dischargeCoeffient*self.totalInjectorArea*(part1 + part2*part3)

class Environment:
    R = 8.3143                      # Universal Gas Constant [J/(mol*K)]
    atmosphericPressure = 101325    # Atmospheric Pressure [Pa]
    def __init__(self, temperature: float):
        self.T = temperature # Ambient Temperature [K]

class Aluminum:
    A = [4.8, 0.00322, 155.239]

    def GetSpecificHeatCP(self, temperature: float):
        return (self.A[0] + self.A[1]*temperature)*self.A[2]

class Paraffin:
    density = 800              # Paraffin's density [kg/m^3]
    burnCoefficient = 0.098    # Linear term of the paraffin burn model [mm/s]
    burnExponent = 0.61        # Exponent of the paraffin burn model [-]

class Grain:
    def __init__(self, material, length: float, internalDiameter: float, externalDiameter: float):
        self.material = material                    # Usually a Paraffin object
        self.internalDiameter = internalDiameter    # Grain hole diameter (assumes circular profile) [m]
        self.externalDiameter = externalDiameter    # Grain external diameter [m]
        self.length = length                        # Grain length [m]

class Tank:
    def __init__(self, material: Aluminum, fluid: NitrousOxide, environment: Environment, volume: float, tankMass: float, loadedFluidMass: float):
        self.material = material
        self.fluid = fluid
        self.volume = volume
        self.tankMass = tankMass
        self.loadedFluidMass = loadedFluidMass  

        self.totalN2O = 1e3*self.loadedFluidMass/NitrousOxide.molecularMass # [mol]

        vaporPressureGaseous = fluid.GetVaporPressure()
        molarVolumeLiquid = fluid.GetMolarVolumeLiquid()
        self.quantityGaseous = vaporPressureGaseous*(self.volume - molarVolumeLiquid*self.totalN2O) / (environment.R*fluid.temperature - vaporPressureGaseous*molarVolumeLiquid)
        self.quantityLiquid = (self.totalN2O*environment.R*fluid.temperature - vaporPressureGaseous*self.volume) / (environment.R*fluid.temperature - vaporPressureGaseous*molarVolumeLiquid)

    # Specific heat of aluminum tank casing [J/(kg*K)]
    def GetSpecificHeatCPofTankMaterial(self, temperature: float):
        return self.material.GetSpecificHeatCP(self.material, temperature=temperature)
    
    # Set Quantity of gaseous N2O molecules [mol]
    def SetMolarQuantityGaseous(self, quantityGaseous):
        self.quantityGaseous = quantityGaseous
        
    # Set Quantity of liquid N2O molecules [mol]
    def SetMolarQuantityLiquid(self, quantityLiquid):
        self.quantityLiquid = quantityLiquid
        
class Nozzle:
    def __init__(self, entryArea: float, exitArea: float, throatArea: float, convergentAngle: float):
        self.exitArea = exitArea
        self.throatArea = throatArea
        self.convergentAngle = convergentAngle
        self.entryArea = entryArea

        self.superAreaRatio = self.exitArea/self.throatArea

    def GetConvergentSectionInternalVolume(self) -> float:
        chamberRadius = np.sqrt(self.entryArea/np.pi)
        throatRadius = np.sqrt(self.throatArea/np.pi)
        return np.pi * (chamberRadius - throatRadius) * (chamberRadius**2 + chamberRadius * throatRadius + throatRadius**2) / (3 * np.tan(self.convergentAngle))

class Chamber:
    def __init__(self, pressure: float, transversalArea: float, preCombustorLength: float, postCombustorLength: float):
        self.transversalArea = transversalArea
        self.postCombustorLength = postCombustorLength
        self.preCombustorLength = preCombustorLength
        self.pressure = pressure

    def GetChamberInternalVolume(self, grain: Grain) -> float:
        return self.transversalArea * (self.postCombustorLength + self.preCombustorLength) + grain.length * np.pi * (grain.internalDiameter ** 2) / 4

class RocketEngine:
    burnEffectiveness  = 0.85 # Efficiency with which N2O and paraffin will mixture burn

    def __init__(self, injector: Injector, tank: Tank, nozzle: Nozzle, grain: Grain, chamber: Chamber):
        self.injector = injector
        self.tank = tank
        self.nozzle = nozzle
        self.grain = grain
        self.chamber = chamber

    # Engine internal volume for burn gases in the combustion chamber [m^3]
    def GetEngineInternalVolume(self):
        return self.chamber.GetChamberInternalVolume() + self.nozzle.GetConvergentSectionInternalVolume()
        
class SimulationParameters:
    def __init__(self, environment: Environment):
        self.environment = environment

class SolveSimulation:
    def __init__(self, rocketEngine: RocketEngine, simulationParameters: SimulationParameters):
        self.rocketEngine = rocketEngine
        self.simulationParameters = simulationParameters

    def RunBlowDown(self):
        deltat = 0.01
        time = 5

        for i in np.arange(0, time, deltat):
            [dT, dng, dnl] = self.CalculateTankDerivatives()

            self.rocketEngine.tank.fluid.SetTemperature(self.rocketEngine.tank.fluid.temperature + dT*deltat)
            self.rocketEngine.tank.SetMolarQuantityGaseous(self.rocketEngine.tank.quantityGaseous + dng*deltat)
            self.rocketEngine.tank.SetMolarQuantityLiquid(self.rocketEngine.tank.quantityLiquid + dnl*deltat)

        print(self.rocketEngine.tank.fluid.temperature)
        print(self.rocketEngine.tank.quantityLiquid)
        print(self.rocketEngine.tank.quantityGaseous)


    def RunSimulation(self):
        self.CalculateTankDerivatives()

        self.RunCEA()
        self.UpdateBurn()
        self.SaveResults()

        # run tank sim
        # run cea maybe
        # update burn model / grain size
        blabla = 1

    def RunCEA(self) -> float:
        

    # def UpdateBurn(self) -> float:
        

    # def SaveResults(self) -> float:
        
        
    # Returns 
    # 0 -> derivative of nitrous oxide temperature
    # 1 -> derivative of the geseous molar quantity
    # 2 -> derivative of the liquid molar quantity
    def CalculateTankDerivatives(self) -> List[float]:
        tank = self.rocketEngine.tank
        fluid = tank.fluid
        injector = self.rocketEngine.injector
        chamber = self.rocketEngine.chamber

        if(tank.quantityLiquid >= 0):
            phase = Phase.LIQUID
        else:
            phase = Phase.GAS

        P = tank.quantityGaseous * Environment.R * fluid.temperature / (tank.volume - tank.quantityLiquid * fluid.GetMolarVolumeLiquid())
        a = tank.tankMass * tank.GetSpecificHeatCPofTankMaterial(fluid.temperature) + tank.quantityGaseous * fluid.GetSpecificHeatCPGaseous() + tank.quantityLiquid * fluid.GetSpecificHeatCPLiquid()
        b = P * fluid.GetMolarVolumeLiquid()
        e = - fluid.GetVaporizationHeat() + Environment.R * fluid.temperature
        f = - injector.GetFluidMassFlow(fluid.pressure, chamber.pressure, phase, fluid) / fluid.molecularMass
        j = - fluid.GetMolarVolumeLiquid() * fluid.GetVaporPressure()
        k = (tank.volume - tank.quantityLiquid * fluid.GetMolarVolumeLiquid()) * fluid.GetVaporPressureDerivTemp()
        m = Environment.R * fluid.temperature
        q = Environment.R * tank.quantityGaseous
        Z = (-f * (-j * a + (q - k) * b)) / (a * (m + j) + (q - k) * (e - b))
        W = (-Z * (m * a + (q - k) * e)) / (-j * a + (q - k) * b)

        return [(b * W + e * Z) / a, Z, W]
