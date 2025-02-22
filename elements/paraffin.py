class Paraffin:
    density = 800              # Paraffin's density [kg/m^3]
    
    def __init__(self, burn_coefficient = 0.098, burn_exponent = 0.61):
        self.burn_coefficient = burn_coefficient
        self.burn_exponent = burn_exponent