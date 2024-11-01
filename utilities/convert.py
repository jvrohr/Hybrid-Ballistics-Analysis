def ConvertPa2Psia(pressurePa: float):
        return 0.000145038*pressurePa

def ConvertFts2Ms(velocity: float):
    return velocity / 3.2808398950131

def ConvertLbmFt32Kgm3(density: float):
    return density*16.01846