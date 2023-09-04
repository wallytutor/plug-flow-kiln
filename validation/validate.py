# -*- coding: utf-8 -*-
from pathlib import Path
from scipy.io import savemat
import cantera as ct
import numpy as np


def main():
    """ Generate validation data for package. """
    root = Path(__file__).resolve().parents[1]
    mech = root / "src/+kiln/data/1S_CH4_MP1.yaml"

    n_points = 100

    fuel = "CH4:0.90, N2:0.10"
    oxid = "N2:0.78, O2:0.21, AR:0.01"
    T = np.linspace(300, 1500, n_points)
    phi = np.linspace(0.0, 2.0, n_points)

    gas = ct.Solution(mech)
    sol = ct.SolutionArray(gas, shape=T.shape)

    sol.TP = T, ct.one_atm
    sol.set_equivalence_ratio(phi, fuel, oxid)
    sol.equilibrate("HP")

    qdot = -1 * sol.heat_release_rate
    h = sol.partial_molar_enthalpies / sol.molecular_weights
    wdotk = sol.net_production_rates * sol.molecular_weights

    data = np.vstack((
        phi,                       # 1
        sol.T,                     # 2
        sol.cp_mass,               # 3
        sol.enthalpy_mass,         # 4
        qdot,                      # 5
        sol.mean_molecular_weight, # 6
        sol.density_mass,          # 7
        sol.viscosity,             # 8
        sol.thermal_conductivity,  # 9
        sol.X.T,                   # 10-15
        sol.Y.T,                   # 16-21
        h.T,                       # 22-27
        wdotk.T                    # 28-33
    )).T

    savemat("validate.mat", {"data": data})

if __name__ == "__main__":
    main()
