import numpy as np
from enum import Enum
from elements.environment import *

class Phase(Enum):
    LIQUID = "liquid"
    GAS = "gas"

class NitrousOxide:
    # N2O Properties' coefficients, ESDU 91022
    G = [96.512, -4045, -12.277, 2.886e-5, 2]                   # Vapour pressure coefficients [Pa] (valid for 182.3 K - 309.57 K)
    J = [-200, 116.043, -917.225, 794.779, -589.587]            # Liquid N2O Specific enthalpy coefficients [kJ/kg]
    L = [-200, 440.055, -459.701, 434.081, -485.338]            # Gaseous N2O Specific enthalpy coefficients [kJ/kg]
    D = [132.632, 0.052187, -0.364923, -1.20233, 0.536141]      # Specific heat at constant pressure for gaseous N2O [kJ/(kg*K)]
    E = [2.49973, 0.023454, -3.80136, 13.0945, -14.5180]        # Specific heat at constant pressure for liquid N2O [kJ/(kg*K)]
    Q = [2.781, 0.27244, 309.57, 0.2882]                        # Specific molar volume coefficients for liquid N2O [m^3/kmol]
    F = [-1.009, -6.28792, 7.50332, -7.90463, -0.629427]        # Vapor Density of N2O [kg/m^3]
    H = [1.72328, -0.83950, 0.51060, -0.10412]                  # Liquid N2O Density [kg/m^3]

    molecularMass = 0.044013        # N2O molecular mass [kg/mol]
    criticalPressure = 7251000      # Critical Pressure of N2O [Pa]
    criticalTemperature = 309.57    # Critical Temperature of N2O [K]
    criticalDensity = 452           # Critical Density of N2O [kg/m^3]

    def __init__(self, temperature: float):
        self.temperature = temperature  # Current approximately uniform N2O temperature [K]

        self.pressure = self.GetVaporPressure()
        self.temperatureRatio = self.temperature/self.criticalTemperature

        self.Z = 1

    def AddTemperatureVariation(self, temperatureVariation: float) -> None:
        self.temperature = self.temperature + temperatureVariation
        self.pressure = self.GetVaporPressure()
        self.temperatureRatio = self.temperature/self.criticalTemperature

    def SetPressure(self, pressure: float) -> None:
        self.pressure = pressure

    # Molar volume of liquid N2O [m^3/mol]
    def GetMolarVolumeLiquid(self) -> float:
        return self.molecularMass/self.GetLiquidDensity()

    # Vapor Pressure of Saturated N2O [Pa]
    def GetVaporPressure(self) -> float:
        return np.exp(self.G[0] + self.G[1]/self.temperature + self.G[2]*np.log(self.temperature) + self.G[3]*(self.temperature**self.G[4]))

    # Vapor Pressure Temperature Derivative of N2O [Pa/Ks]
    def GetVaporPressureDerivTemp(self) -> float:
        return (-self.G[1]/(self.temperature**2) + self.G[2]/self.temperature + self.G[3]*self.G[4]*self.temperature**(self.G[4]-1)) * np.exp(self.G[0] + self.G[1]/self.temperature + self.G[2]*np.log(self.temperature) + self.G[3]*self.temperature**self.G[4])

    # Gaseous N2O Specific Heat for constant pressure [J/(mol*K)]
    def GetSpecificHeatCPGaseous(self) -> float:
        auxVariable = 1 - self.temperatureRatio
        return self.molecularMass*1e3*(self.D[0]*(1 + self.D[1]*(auxVariable**(-2/3)) + self.D[2]*(auxVariable**(-1/3)) + self.D[3]*(auxVariable**(1/3)) + self.D[4]*(auxVariable**(2/3))))

    # Liquid N2O Specific Heat for constant pressure [J/(mol*K)]
    def GetSpecificHeatCPLiquid(self) -> float:
        auxVariable = 1 - self.temperatureRatio
        return self.molecularMass*1e3*(self.E[0]*(1 + self.E[1]/auxVariable + self.E[2]*auxVariable + self.E[3]*(auxVariable**2) + self.E[4]*(auxVariable**3)))

    # Specific Heat for constant pressure [J/(mol*K)]
    def GetSpecificHeatCP(self, phase: Phase) -> float:
        if phase == Phase.LIQUID:
            return self.GetSpecificHeatCPLiquid()
        elif phase == Phase.GAS:
            return self.GetSpecificHeatCPGaseous()
        else:
            raise ValueError(phase)

    # Adiabatic expansion coefficient, specific heats ratio Cp/Cv [-]
    def GetSpecificHeatsRatio(self, phase: Phase):
        return 1/(1 - Environment.R/self.GetSpecificHeatCP(phase))

    # Vaporization Heat of N2O [J/mol]
    def GetVaporizationHeat(self) -> float:
        return self.GetSpecificEnthalpyGaseous() - self.GetSpecificEnthalpyLiquid()
    
    # Specific Enthalpy Liquid [J/mol]
    def GetSpecificEnthalpyLiquid(self):
        auxVariable = 1 - self.temperatureRatio
        return self.molecularMass*1e3*(self.J[0] + self.J[1]*auxVariable**(1/3) + self.J[2]*auxVariable**(2/3) + self.J[3]*auxVariable + self.J[4]*auxVariable**(4/3))

    # Specific Enthalpy Gaseous [J/mol]
    def GetSpecificEnthalpyGaseous(self):
        auxVariable = 1 - self.temperatureRatio
        return self.molecularMass*1e3*(self.L[0] + self.L[1]*auxVariable**(1/3) + self.L[2]*auxVariable**(2/3) + self.L[3]*auxVariable + self.L[4]*auxVariable**(4/3))

    # Density of gaseous/vapor N2O [kg/m^3]
    def GetVaporDensity(self, temperature: float) -> float:
        auxVariable = 1/(temperature/self.criticalTemperature) - 1
        return self.criticalDensity*(np.exp(self.F[0]*(auxVariable**(1/3)) + self.F[1]*(auxVariable**(2/3)) + self.F[2]*auxVariable + self.F[3]*(auxVariable**(4/3)) + self.F[4]*(auxVariable**(5/3))))

    # Density of liquid N2O [kg/m^3]
    def GetLiquidDensity(self) -> float:
        auxVariable = 1 - self.temperatureRatio
        return self.criticalDensity*(np.exp(self.H[0]*(auxVariable**(1/3)) + self.H[1]*(auxVariable**(2/3)) + self.H[2]*auxVariable + self.H[3]*(auxVariable**(4/3))))
