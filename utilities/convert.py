def convert_pa_2_psia(pressurePa: float):
        return 0.000145038*pressurePa

def convert_fts_2_ms(velocity: float):
    return velocity / 3.2808398950131

def convert_lbmft3_2_kgm3(density: float):
    return density*16.01846

def rankine_2_kelvin(rankine):
    return rankine / 1.8

def g_mol_2_kg(g_mol):
    return g_mol / 1000