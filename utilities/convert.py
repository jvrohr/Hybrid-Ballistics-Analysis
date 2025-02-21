def convert_pa_2_psia(pressurePa: float):
        return 0.000145038*pressurePa

def convert_fts_2_ms(velocity: float):
    return velocity / 3.2808398950131

def convert_lbmft3_2_kgm3(density: float):
    return density*16.01846