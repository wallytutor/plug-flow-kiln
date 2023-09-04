# -*- coding: utf-8 -*-
from dataclasses import dataclass
from pathlib import Path
from subprocess import run
from IPython import embed
from majordome.unit_conversion import FlowUnits as fu
from scipy.integrate import cumtrapz
from scipy.optimize import root
import codecs
import json
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Add mechanism to path.
ct.add_directory("../../src/octave/+kiln/data")

CASE_REFERENCE = "cases/reference.json"
# CASE_REFERENCE = "cases/reference-coating.json"
""" Case containing reference parameters for simulation. """

MECHANISM = "1S_CH4_MP1.yaml"
""" Mechanism for which computations are carried out. """

SPECIES = ["CH4", "O2", "CO2", "H2O", "AR", "N2"]
""" Order of species reporting for model compatibility. """

OCTAVE = "octave-cli-7.3.0.exe"
# OCTAVE = "octave-cli-8.1.0.exe"
""" Octave executable used for running case. """

# SIMULATION_KIND = "EBU"
SIMULATION_KIND = "precombusted"
""" Simulation type (for combustion). """

AIR_GAS_RATIO = 1.0
""" Reference equivalence ratio [kg/kg]. """

FUEL = "CH4: 0.95, N2: 0.05"
""" Assumed fuel composition [mole fractions]. """

OXID = "N2: 0.78, O2: 0.21, AR: 0.01, H2O: 0.006"
""" Assumed oxidizer composition [mole fractions]. """

TEMP_INJECTION = 300.0
""" Injection temperature for premix. [K]. """


@dataclass
class Results:
    """ Structure for returning preprocessed data. """
    equivalent_temperature: float = None
    mixture_flow_rate: float = None
    secondary_flow_rate: float = None
    equivalent_composition: list = None


def get_gas(comp=None):
    """ Get gas at standard state conditions. """
    gas = ct.Solution(MECHANISM)

    match comp:
        case None:
            gas.TP = 273.15, ct.one_atm
        case _:
            gas.TPX = 273.15, ct.one_atm, comp

    return gas


def get_mixture_mass_flow(qdot_nat, phi, fuel, oxid):
    """ Compute oxidizer flow for given parameters.
    
    This function follows the logic:
    - Create a "pure" fuel and air solutions.
    - Use fuel to convert normal to mass flow rate.
    - Compute premixed gas under given equivalence ratio.
    - Find a target composition and methane content.
    - Create "quantities" of fuel an air for optimization.
    - Change air mass until air-gas mixture matches target methane.
    - Verify mass averaged mixture matches equivalence ratio.
    - Return mixture with computed "primary air mass".
    """
    # Create solutions with given compositions.
    mix_fuel = get_gas(fuel)
    mix_oxid = get_gas(oxid)

    # Get natural gas mass flow rate [kg/s].
    mw = mix_fuel.mean_molecular_weight
    mdot_nat = fu.normal_flow_to_mass_flow(qdot_nat, mw)

    # Compute composition of premixed air-gas.
    mix_phi = get_gas()
    mix_phi.set_equivalence_ratio(phi, fuel, oxid)
    
    # Get premixed dictionary and target CH4 content.
    mix = mix_phi.mole_fraction_dict()
    x_ch4 = mix["CH4"]

    # Create quantities for mixture optimization.
    q_fuel = ct.Quantity(mix_fuel, mass=mdot_nat)
    q_oxid = ct.Quantity(mix_oxid, mass=mdot_nat)

    def mix_gases(mdot_air):
        """ Helper for objective and validation. """
        q_oxid.mass = mdot_air
        return q_fuel + q_oxid

    def objective(mdot_air):
        """ Residual based on methane fraction. """
        xm = mix_gases(mdot_air).mole_fraction_dict()
        return (x_ch4 - xm["CH4"])**2

    # Check if optimization reportedly converged.
    if not (sol := root(objective, [mdot_nat])).success:
        raise ValueError("Could not find mixture flow!")

    # Retrieve air flow rate and mixture composition.
    mdot_air = sol.x[0]
    xm = mix_gases(mdot_air).mole_fraction_dict()

    # Validate final mixture matches equivalence ratio.
    for k, v in mix.items():
        if not np.allclose(v, xm[k]):
            raise RuntimeError(f"Optimization failed for {k}")

    return mix_gases(mdot_air)


def get_secondary_air_flow(mdot_air_gas, air, mix, x_o2):
    """ Get secondary air flow to match measured oxygen.
    
    This function follows the logic:
    - Create "quantities" of air and premixed gas for optimization.
    - Equilibrate gas to find adiabatic flame temperature.
    - Change air mass until air-premix mixture matches target oxygen.
    - Return mixture with computed "secondary air mass".
    """
    # Create quantities for mixture optimization.
    q_air = ct.Quantity(air, mass=mdot_air_gas)
    q_mix = ct.Quantity(mix, mass=mdot_air_gas)

    # Find adiabatic flame conditions.
    q_mix.equilibrate("HP")

    def mix_gases(mdot_leak):
        """ Helper for objective and validation. """
        q_air.mass = mdot_leak
        return q_mix + q_air
    
    def objective(mdot_leak):
        """ Residual based on methane fraction. """
        xm = mix_gases(mdot_leak).mole_fraction_dict()
        return (x_o2 - xm["O2"])**2

    # Check if optimization reportedly converged.
    if not (sol := root(objective, [0.0])).success:
        raise ValueError("Could not find secondary flow!")
    
    mdot_leak = sol.x[0]
    inlet = mix_gases(mdot_leak)

    return inlet


def report_preprocess(mdot_premix, mdot_leak, inlet):
    """ Report preprocessing results and return composition. """
    qdot = inlet.mass * inlet.cp_mass * inlet.T / 1_000_000

    print(
        f"\nEquivalent temperature .... {inlet.T:.3f} K"
        f"\nMixture flow rate ......... {mdot_premix:.6f} kg/s"
        f"\nAir leak rate ............. {mdot_leak:.6f} kg/s"
        f"\nEnergy input .............. {qdot:.2f} MW"
        f"\n\nEquivalent composition:"
    )

    arr = []
    for species in SPECIES:
        fraction = inlet.X[inlet.species_index(species)]
        print(f"\t{species:10s} ... {fraction:.12e}")
        arr.append(fraction)

    return arr


def get_command(octave_persist):
    """ Get command for running simulation. """
    cmd = [OCTAVE, "simulate.m"]

    if octave_persist:
        cmd.insert(1, "-q")
        cmd.insert(2, "--persist")

    return cmd


def load_results(results):
    """ Load simulation results from JSON file. """
    with codecs.open(results, encoding="unicode_escape") as fp:
        text = fp.read()
        data = json.loads(text)
        df = pd.DataFrame(data)
    return df


def plot_fluxes(z, qdot, q):
    """ Plot heat losses through kiln shell. """
    plt.close("all")
    plt.style.use("seaborn-white")

    ax1 = plt.subplot(111)
    ax2 = ax1.twinx()

    ax1.plot(z, qdot, "r-")
    ax1.grid(linestyle=":")
    ax1.set_xlabel("Coordinate [m]")
    ax1.set_ylabel("Heat flux [kW/m]", fontdict={"color": "r"})

    ax2.plot(z, q, "b-")
    ax2.set_xlabel("Coordinate [m]")
    ax2.set_ylabel("Total loss [MW]", fontdict={"color": "b"})

    plt.tight_layout()
    plt.show()


def plot_shell_temperature(sim, shell):
    """ Plot shell temperature for validation. """
    plt.close("all")
    plt.style.use("seaborn-white")

    if shell is not None:
        z1 = shell["z"].to_numpy()
        z1 = max(z1) - z1
        T_min = shell["T_min"].to_numpy()
        T_max = shell["T_max"].to_numpy()

        plt.scatter(z1, T_max, c="r")
        plt.scatter(z1, T_min, c="b")

    z2 = sim["z"].to_numpy()
    z2 = max(z2) - z2
    T_sim = sim["T_sh"].to_numpy() - 273.15

    plt.plot(z2, T_sim)
    plt.grid(linestyle=":")
    plt.xlabel("Coordinate [m]")
    plt.ylabel("Temperature [°C]")
    plt.tight_layout()
    plt.show()


def inner_preprocess_data(
        qdot_nat,
        phi,
        fuel,
        oxid,
        x_o2,
        temp_air_gas,
        temp_air_leak,
    ):
    """ Find flow rates for setting up simulation.
    
    The only known quantities for energy input in the kiln are the
    gas flow rate and equivalence ratio. From these values we need
    to find out the required primary air flow. Because there is a
    secondary air leak at known temperature and oxygen level at the
    tail of kiln is known, we need to compute what additional air
    flow is compatible with measurements. 
    """
    ###############################################################
    # COMPUTE PRIMARY AIR MIXTURE
    ###############################################################

    # Get appropriated air-gas mixture at `phi`.
    premix = get_mixture_mass_flow(qdot_nat, phi, fuel, oxid)

    # Premixed equivalent injection flow rate [kg/s].
    mdot_premix = premix.mass

    ###############################################################
    # COMPUTE SECONDARY AIR LEAK AND EQUIVALENT INLET
    ###############################################################

    # Create burner premixed gas and air objects.
    mix = get_gas(premix.X)
    air = get_gas(oxid)

    # Set known temperatures to solutions.
    mix.TP = temp_air_gas, ct.one_atm
    air.TP = temp_air_leak, ct.one_atm

    # Compute equivalent inlet gas to match O2 measurement.
    inlet = get_secondary_air_flow(mdot_premix, air, mix, x_o2)

    # Secondary air flow rate [kg/s].
    mdot_leak = inlet.mass - mdot_premix

    ###############################################################
    # COMPUTE EQUIVALENT "COLD" INLET FOR EBU/MAK SETUP
    ###############################################################

    # Reset state of premixed air-gas.
    mix.TPX = temp_air_gas, ct.one_atm, premix.X

    # Create quantities for mixture optimization.
    q_air = ct.Quantity(air, mass=mdot_leak)
    q_mix = ct.Quantity(mix, mass=mdot_premix)

    # Mix quantities without equilibrating premixed air-gas.
    inlet_cold = q_mix + q_air

    ###############################################################
    # REPORT RESULTS
    ###############################################################

    print("\n *** EQUILIBRATED GAS")
    arr_hot = report_preprocess(mdot_premix, mdot_leak, inlet)

    print("\n *** COLD GAS")
    arr_ebu = report_preprocess(mdot_premix, mdot_leak, inlet_cold)

    ###############################################################
    # RETURN RESULTS
    ###############################################################

    # Select which preprocessed data to use.
    match SIMULATION_KIND:
        case "precombusted":
            results = Results(
                equivalent_temperature=inlet.T,
                mixture_flow_rate=mdot_premix,
                secondary_flow_rate=mdot_leak,
                equivalent_composition=arr_hot
            )
        case _:
            results = Results(
                equivalent_temperature=inlet_cold.T,
                mixture_flow_rate=mdot_premix,
                secondary_flow_rate=mdot_leak,
                equivalent_composition=arr_ebu
            )

    return results


def preprocess_data(product):
    """ Manage data preprocessing and case generation. """
    # Load preset for preprocessing.
    with open(f"products/{product}.json",) as fp:
        preset = json.load(fp)

    # Reference natural gas flow rate [Nm³/h].
    qdot_nat = preset["qdot_nat"]

    # Residual oxigen at fumes outlet [mole fraction]
    x_o2 = preset["x_o2"]

    # Secondary air temperature [K].
    temp_air_leak = preset["temp_air_leak"] + 273.15

    # Get structure with preprocessed data.
    results = inner_preprocess_data(
        qdot_nat, AIR_GAS_RATIO, FUEL, OXID, x_o2,
        TEMP_INJECTION, temp_air_leak)
    
    # Load kiln reference conditions/geometry.
    with open(CASE_REFERENCE) as fp:
        reference = json.load(fp)

    # Load selected bed material data file.
    with open(preset["material_file"]) as fp:
        bed_data = json.load(fp)

    # To provide external ProjectData.m file.
    reference["project_data"] = preset["project_data"]

    # NOTE: this will modify directly the `reference` data without
    # need to make `reference["operation"] = operation` later!
    operation = reference["operation"]
    operation["rotation_rate"] = preset["rotation_rate"]
    operation["bed_feed_rate"] = preset["bed_feed_rate"]
    operation["gas_flow_rate"] = results.mixture_flow_rate
    operation["gas_leak_rate"] = results.secondary_flow_rate
    operation["gas_feed_temperature"] = results.equivalent_temperature
    operation["gas_feed_composition"] = results.equivalent_composition

    hypothesis = reference["hypothesis"]

    match SIMULATION_KIND:
        case "precombusted":
            hypothesis["gas_equilibrated"] = True
        case "EBU":
            hypothesis["gas_equilibrated"] = False
            hypothesis["gas_kinetics"] = "EBU"
        case "MAK":
            hypothesis["gas_equilibrated"] = False
            hypothesis["gas_kinetics"] = "MAK"
        case _:
            raise ValueError(f"Uknown kind {SIMULATION_KIND}")
    
    reference["materials"]["bed_data"] = bed_data

    with open("cases/simulate.json", "w") as fp:
        json.dump(reference, fp, indent=4)


def postprocess_data(product):
    """ Manage standard data post-processing. """
    # Load data for postprocessing.
    sim = load_results("results/simulate.json")

    # Handle shell data availability.
    shell = None
    shell_file = Path(f"products/{product}_shell.csv")
    if shell_file.exists():
        shell = pd.read_csv(shell_file)

    # Nodes coordinates [m].
    z = sim["z"].to_numpy()

    # Environment convective HTC [W/(m.K)].
    h_env = sim["h_env"].to_numpy()

    # Steel shell emissivity [-].
    e_env = 0.79

    # Environment temperature [K].
    T_env = sim["T_env"].to_numpy()

    # Shell temperature [K].
    T_sh = sim["T_sh"].to_numpy()

    # Shell perimeter [m].
    P_env = 2 * np.pi * sim["R_sh"].to_numpy()

    # Compute contributions from different terms.
    con = h_env * (T_sh - T_env)
    rad = e_env * ct.stefan_boltzmann * (T_sh**4 - T_env**4)

    # Compute heat flux density [kW/m].
    qdot = P_env * (con + rad) / 1000

    # Integrate data.
    q = cumtrapz(qdot, z, initial=0.0) / 1000.0

    plot_fluxes(z, qdot, q)
    plot_shell_temperature(sim, shell)

    # Estimation of energy consumption by bed ~1.1 MW.
    # dh = 191.1 / 0.1019613
    # Tf = sim.T_b[0]
    # bed_heat = sim.mdot_b * dh / 1000000

    bal_b = sim["bal_b"].sum()
    bal_g = sim["bal_g"].sum()
    bal_c = sim["bal_c"].sum()
    check = bal_g - bal_c - bal_b

    if not np.allclose(check, 0.0):
        print("Simulation balance fails!")

    # Retrieve matrix of mole fractions.
    X = sim[[f"X_{c}" for c in SPECIES]].to_numpy()

    # Create solution array with computed states.
    gas = ct.Solution(MECHANISM)
    sol = ct.SolutionArray(gas, shape=(sim.shape[0],))
    sol.TPX = sim.T_g, ct.one_atm, X

    # Compute gas enthalpy change.
    h = sol.enthalpy_mass
    u = sol.int_energy_mass
    dh = h[-1] - h[0]
    du = u[-1] - u[0]

    embed(colors="Linux")


def simulate(product, octave_persist):
    """ Manage process simulation and embed processing. """
    preprocess_data(product)
    run(get_command(octave_persist))
    postprocess_data(product)


if __name__ == "__main__":
    # Select product.
    # product = "c2_6"
    product = "h3_0"
    # product = "c2_6_b"

    simulate(product, octave_persist=True)
