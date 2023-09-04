# -*- coding: utf-8 -*-
import cantera as ct
import pandas as pd

ct.add_directory("../../src/+kiln/data")

fuel = "CH4: 1.00"
oxid = "N2: 0.78 O2: 0.21 AR: 0.01 H2O: 0.006"


def density_mass_ref(X):
    """ Gas density at reference state [kg/m³]. """
    gas = ct.Solution("1S_CH4_MP1.yaml")
    gas.TPX = 298.15, ct.one_atm, X
    return gas.density_mass


def get_combusted_gas(mdot_fuel, mdot_oxid):
    """ Compute equivalent combustion gas. """
    gas_fuel = ct.Solution("1S_CH4_MP1.yaml")
    gas_oxid = ct.Solution("1S_CH4_MP1.yaml")

    gas_fuel.TPX = 298.15, ct.one_atm, fuel
    gas_oxid.TPX = 298.15, ct.one_atm, oxid

    q1 = ct.Quantity(gas_fuel, mass=mdot_fuel)
    q2 = ct.Quantity(gas_oxid, mass=mdot_oxid)

    qm = q1 + q2
    qm.equilibrate("HP")

    return qm


if __name__ == "__main__":
    df = pd.read_csv("data/data-barr.csv")

    df["mdot_fuel"] = df["Fuel flow rate (L/s)"]
    df["mdot_fuel"] *= density_mass_ref(fuel) * (1/1000)

    df["mdot_oxid"] = df["Primary air flow rate (L/s)"]
    df["mdot_oxid"] += df["Secondary air flow rate (L/s)"]
    df["mdot_oxid"] *= density_mass_ref(oxid) * (1/1000)

    for k, row in df.iterrows():
        mdot_fuel = row["mdot_fuel"]
        mdot_oxid = row["mdot_oxid"]
        qm = get_combusted_gas(mdot_fuel, mdot_oxid)
        X = repr(qm.X)[7:-2].replace("\n", "").replace(" ", "")
        df.loc[k, ("Tf", "Xf", "mdot")] = qm.T, X, qm.mass

    # TODO automate generation of input files.
    print(df)
