from materials.material import Fuel

class Paraffin(Fuel):
    paraffin_input_string = """
    fuel paraffin(S)  C 73.0   H 124.0     wt%=100.00
    h,cal=-444694.0724016     t(k)=298.15   rho=1.001766
    """
    def __init__(self, density = 800, burn_coefficient = 0.098, burn_exponent = 0.61, 
                 cea_input_string = paraffin_input_string):
        super().__init__("Paraffin", density, burn_coefficient, burn_exponent, cea_input_string)

