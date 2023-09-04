# -*- coding: utf-8 -*-
from pathlib import Path
from IPython import embed
from scipy.optimize import root
import cantera as ct


def find_air_leak(
        mechanism: str | Path,
        mdot_air_gas: float,
        temp_air_gas: float,
        temp_air_leak: float,
        air_gas_equiv: float,
        fuel: str | dict[str, float],
        oxid: str | dict[str, float],
        air: str | dict[str, float],
        x_o2: float
    ) -> ct.Quantity:
    """ Find air leak flow explaining measured oxygen content. """
    gas1 = ct.Solution(mechanism)
    gas1.TP = temp_air_gas, ct.one_atm
    gas1.set_equivalence_ratio(air_gas_equiv, fuel, oxid)

    gas2 = ct.Solution(mechanism)
    gas2.TPX = temp_air_leak, ct.one_atm, air

    # Original mixture for source gas.
    q0 = ct.Quantity(gas1, mass=mdot_air_gas)

    # Original mixture for combustion fumes gas.
    q1 = ct.Quantity(gas1, mass=mdot_air_gas)
    q1.equilibrate("HP")

    idx_o2 = q1.species_index("O2")

    def objective(f2):
        """ Minimize target of equilibrium composition. """
        qm = q1 + ct.Quantity(gas2, mass=f2)
        return x_o2 - qm.X[idx_o2]

    if not (sol := root(objective, [mdot_air_gas])).success:
        raise ValueError("Could not find required flow!")

    q2 = ct.Quantity(gas2, mass=sol.x[0])

    # This is the equivalent pre-equilibrated gas.
    qm = q1 + q2

    # This is the equivalent source gas.
    qn = q0 + q2

    return qm, qn


def pretty_report(qm, qn, mdot_air_gas):
    """ Formatted report of results. """
    print(
        f"\n *** EQUILIBRATED GAS"
        f"\nEquivalent temperature .... {qm.T:.2f} K"
        f"\nAir leak rate ............. {qm.mass-mdot_air_gas:.4f} kg/s"
        f"\nEquivalent composition .... {qm.mole_fraction_dict()}"
    )

    print(
        f"\n *** SOURCE GAS"
        f"\nEquivalent temperature .... {qn.T:.2f} K"
        f"\nAir leak rate ............. {qn.mass-mdot_air_gas:.4f} kg/s"
        f"\nEquivalent composition .... {qn.mole_fraction_dict()}"
    )

    print("\n *** END ***\n")


if __name__ == "__main__":
    here = Path(__file__).resolve().parents[1]
    mechanism = here / "octave/+kiln/data/1S_CH4_MP1.yaml"
    # mechanism = "gri30.yaml"

    x_o2 = 0.096
    mdot_air_gas = 0.115
    temp_air_gas = 300.0
    air_gas_equiv = 1.0
    fuel = "CH4: 0.90 N2: 0.10"
    oxid = "N2: 0.78 O2: 0.21 AR: 0.01 H2O: 0.006"

    temp_air_leak = temp_air_gas
    air = oxid

    qm, qn = find_air_leak(
        mechanism,
        mdot_air_gas,
        temp_air_gas,
        temp_air_leak,
        air_gas_equiv,
        fuel,
        oxid,
        air,
        x_o2
    )

    pretty_report(qm, qn, mdot_air_gas)

    embed(colors="Linux")
