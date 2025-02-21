def convert_pa_to_psia(pressure_pa: float) -> float:
    """
    Convert pressure from Pascal to PSI.
    
    :param pressure_pa: Pressure in Pascal [Pa]
    :return: Pressure in PSI [psi]
    """
    return 0.000145038 * pressure_pa

def convert_fts_to_ms(velocity: float) -> float:
    """
    Convert velocity from feet per second to meters per second.
    
    :param velocity: Velocity in feet per second [ft/s]
    :return: Velocity in meters per second [m/s]
    """
    return velocity / 3.2808398950131

def convert_lbmft3_to_kgm3(density: float) -> float:
    """
    Convert density from pounds per cubic foot to kilograms per cubic meter.
    
    :param density: Density in pounds per cubic foot [lbm/ft^3]
    :return: Density in kilograms per cubic meter [kg/m^3]
    """
    return density * 16.01846