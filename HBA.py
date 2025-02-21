import numpy as np
import rocketcea.cea_obj as cea
from elements.engine import *
from utilities.convert import *
from utilities.plot import *

paraffin_input_string = """
    fuel paraffin(S)  C 73.0   H 124.0     wt%=100.00
    h,cal=-444694.0724016     t(k)=298.15   rho=1.001766
    """
cea.add_new_fuel('Paraffin', paraffin_input_string)


class SimulationParameters:
    def __init__(self, environment: Environment, time_step: float, total_time: float):
        self.environment = environment
        self.time_step = time_step
        self.total_time = total_time


class SolveSimulation:
    results_dict = {}

    def __init__(self, rocket_engine: RocketEngine, simulation_parameters: SimulationParameters):
        self.rocket_engine = rocket_engine
        self.simulation_parameters = simulation_parameters

        self.delta_temperature = 0
        self.delta_N_gaseous = 0
        self.delta_N_liquid = 0
        self.delta_pressure = 0

    def Run(self, option="burn"):
        self.plot = PlotResults()
        deltat = self.simulation_parameters.time_step
        time = self.simulation_parameters.total_time

        if option == "burn":
            iteration_function = self.run_burn_iteration
        elif option == "blowdown":
            iteration_function = self.run_blow_down_iteration
        else:
            raise ValueError(option)

        for i in np.arange(0, time, deltat):
            iteration_function(i)

    def run_blow_down_iteration(self, time: float):
        self.update_oxidizer_blowdown()
        self.save_results_blowdown(time)

    def run_burn_iteration(self, time: float):
        self.update_oxidizer_blowdown()
        self.update_fuel_regression()
        self.run_CEA()
        self.update_chamber_pressure()
        self.save_results_burn(time)

    def save_results_blowdown(self, time: float):
        tank = self.rocket_engine.tank
        injector = self.rocket_engine.injector
        
        self.plot.results_dict["Time"].append(time)
        self.plot.results_dict["Temperature"].append(tank.fluid.temperature)
        self.plot.results_dict["Quantity Gas"].append(tank.quantity_gaseous)
        self.plot.results_dict["Quantity Liquid"].append(tank.quantity_liquid)
        self.plot.results_dict["Pressure Tank"].append(tank.fluid.pressure)
        self.plot.results_dict["Oxidizer Mass Flow"].append(injector.oxidizer_mass_flow)
        self.plot.results_dict["Pressure Chamber"].append(self.rocket_engine.chamber.pressure)

    def update_oxidizer_blowdown(self):
        tank = self.rocket_engine.tank
        deltat = self.simulation_parameters.time_step

        self.calculate_tank_derivatives()

        tank.fluid.add_temperature_variation(self.delta_temperature*deltat)
        tank.add_molar_quantity_gaseous_variation(self.delta_N_gaseous*deltat)
        tank.add_molar_quantity_liquid_variation(self.delta_N_liquid*deltat)

    def calculate_tank_derivatives(self) -> None:
        tank = self.rocket_engine.tank
        fluid = tank.fluid
        injector = self.rocket_engine.injector
        chamber = self.rocket_engine.chamber
        environment = self.simulation_parameters.environment

        if(tank.quantity_liquid > 0):
            phase = Phase.LIQUID
        else:
            tank.quantity_liquid = 0
            fluid.Z = fluid.pressure * tank.volume / (fluid.temperature * environment.R * tank.quantity_gaseous)
            phase = Phase.GAS

        injector.oxidizer_mass_flow = injector.get_mass_flow(chamber.pressure, phase, fluid)
        fluid.pressure = fluid.Z * tank.quantity_gaseous * environment.R * fluid.temperature / \
            (tank.volume - tank.quantity_liquid * fluid.get_molar_volume_liquid())          # [Pa]
        a = tank.tank_mass * tank.get_CP_tank_material(fluid.temperature) + \
            tank.quantity_gaseous * fluid.get_CP_gaseous() + \
                tank.quantity_liquid * fluid.get_CP_liquid()                   # [J/K]
        b = fluid.pressure * fluid.get_molar_volume_liquid()                               # [J/mol]
        e = - fluid.get_vaporization_eat() + environment.R * fluid.temperature           # [J/mol]
        f = - injector.oxidizer_mass_flow / fluid.molecular_mass                                                         # [mol/s]
        j = - fluid.get_molar_volume_liquid() * fluid.get_vapor_pressure()                   # [J/mol]
        k = (tank.volume - tank.quantity_liquid * fluid.get_molar_volume_liquid()) * \
            fluid.get_vapor_pressure_deriv_temp()                                           # [J/K]
        m = environment.R * fluid.temperature                                           # [J/mol]
        q = environment.R * tank.quantity_gaseous                                        # [J/K]

        if(phase == Phase.LIQUID):
            self.delta_N_gaseous = (-f * (-j * a + (q - k) * b)) / (a * (m + j) + (q - k) * (e - b))      # [mol/s]
            self.delta_N_liquid = (-self.delta_N_gaseous * (m * a + (q - k) * e)) / (-j * a + (q - k) * b)  # [mol/s] 
        else:
            self.delta_N_gaseous = -f         # [mol/s]
            self.delta_N_liquid = 0           # [mol/s]
        
        self.delta_temperature = (b * self.delta_N_liquid + e * self.delta_N_gaseous) / a

    def update_fuel_regression(self):
        tank = self.rocket_engine.tank
        injector = self.rocket_engine.injector
        chamber = self.rocket_engine.chamber
        grain = self.rocket_engine.grain
        deltat = self.simulation_parameters.time_step

        if(tank.quantity_liquid >= 0):
            phase = Phase.LIQUID
        else:
            phase = Phase.GAS

        oxidizerFlux = injector.oxidizer_mass_flow / (np.pi * (grain.internal_diameter / 2) ** 2)                           # [kg/m^2s]
        regretion_rate = 1e-3 * grain.material.burn_coefficient * (oxidizerFlux ** grain.material.burn_exponent)   # [m/s]
        grain.add_internal_diameter_variation(deltat * regretion_rate*2)                                            # [m]
        fuel_mass_flow = np.pi * grain.internal_diameter * grain.length * regretion_rate * grain.material.density   # [kg/s]

        grain.volume_variation = np.pi*((grain.internal_diameter/2 + regretion_rate*deltat)**2 - (grain.internal_diameter/2)**2)

        chamber.instant_OF = injector.oxidizer_mass_flow / fuel_mass_flow # [-]
        chamber.instant_mass_generation_rate = injector.oxidizer_mass_flow + fuel_mass_flow # [kg/s]

    def run_CEA(self) -> float:
        chamber = self.rocket_engine.chamber
        nozzle = self.rocket_engine.nozzle

        cea_object = cea.CEA_Obj(oxName='N2O', fuelName='Paraffin')

        atmospheric_pressure_psi = convert_pa_2_psia(Environment.atmospheric_pressure)
        chamber_pressure_psi = convert_pa_2_psia(chamber.pressure)

        instant_c_star = cea_object.get_Cstar(Pc=chamber_pressure_psi, \
                                  MR=chamber.instant_OF)
        self.rocket_engine.instant_c_star = convert_fts_2_ms(instant_c_star)

        Cf = cea_object.get_PambCf(Pamb=atmospheric_pressure_psi, Pc=chamber_pressure_psi, \
                                  MR=chamber.instant_OF, eps=nozzle.super_area_ratio)[1]

        combustion_gas_density = cea_object.get_Densities(\
            Pc=chamber_pressure_psi, MR=chamber.instant_OF, eps=nozzle.super_area_ratio)[0]
        combustion_gas_density = convert_lbmft3_2_kgm3(combustion_gas_density)

        self.rocket_engine.gas_mass = self.rocket_engine.get_engine_internal_volume()*combustion_gas_density
        self.rocket_engine.thrust = Cf*self.rocket_engine.chamber.pressure*self.rocket_engine.nozzle.throat_area        

    def update_chamber_pressure(self):
        self.calculate_chamber_pressure_derivative()

        deltat = self.simulation_parameters.time_step
        self.rocket_engine.chamber.add_chamber_pressure_variation(deltat*self.delta_pressure)

    def calculate_chamber_pressure_derivative(self):
        nozzle = self.rocket_engine.nozzle
        chamber = self.rocket_engine.chamber
        grain = self.rocket_engine.grain
        environment = self.simulation_parameters.environment

        if(chamber.pressure - environment.atmospheric_pressure > 0):
            nozzle.mass_flow_nozzle = nozzle.discharge_coefficient*chamber.pressure*nozzle.throat_area/self.rocket_engine.instant_c_star
        else:
            nozzle.mass_flow_nozzle = 0

        mass_gain = chamber.instant_mass_generation_rate - nozzle.mass_flow_nozzle

        self.delta_pressure = chamber.pressure*(mass_gain/self.rocket_engine.gas_mass - grain.volume_variation/self.rocket_engine.get_engine_internal_volume())

    def save_results_burn(self, time):
        self.plot.results_dict["Time"].append(time)
        self.plot.results_dict["Thrust"].append(self.rocket_engine.thrust)
        self.plot.results_dict["Pressure Chamber"].append(self.rocket_engine.chamber.pressure)
        self.plot.results_dict["Pressure Tank"].append(self.rocket_engine.tank.fluid.pressure)
