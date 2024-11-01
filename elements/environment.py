class Environment:
    R = 8.3143                      # Universal Gas Constant [J/(mol*K)]
    atmosphericPressure = 101325    # Atmospheric Pressure [Pa]
    MWair = 0.0289652               # Molecular weight air [kg/mol]

    def __init__(self, temperature: float, pressure=atmosphericPressure):
        self.T = temperature # Ambient Temperature [K]
        self.pressure = pressure
    
    def GetAirDensity(self):
        return self.pressure*self.MWair/(self.R*self.T)
