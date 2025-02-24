class Environment:
    R = 8.3143                       # Universal Gas Constant [J/(mol*K)]
    atmospheric_pressure = 101325    # Atmospheric Pressure [Pa]
    MW_air = 0.0289652               # Molecular weight air [kg/mol]
    gamma_air = 1.4                  # Specific Heat Ratio Air [-]
    gravity_acceleration = 9.81      # Gravity Acceleration [m/s^2]

    def __init__(self, temperature: float, pressure=atmospheric_pressure):
        self.T = temperature # Ambient Temperature [K]
        self.pressure = pressure
    
    def GetAirDensity(self):
        return self.pressure * self.MW_air / (self.R * self.T)
