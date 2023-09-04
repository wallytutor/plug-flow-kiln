# -*- coding: utf-8 -*-
from scipy.integrate import cumtrapz
import codecs
import json
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_results():
    """ Load simulation results from JSON file. """
    with codecs.open("results.json", encoding="unicode_escape") as fp:
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
    ax2.set_ylabel("Total loss [kW]", fontdict={"color": "b"})

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Load data for postprocessing.
    df = load_results()

    # Nodes coordinates [m].
    z = df["z"].to_numpy()

    # Environment convective HTC [W/(m.K)].
    h_env = df["h_env"].to_numpy()

    # Steel shell emissivity [-].
    e_env = 0.79

    # Environment temperature [K].
    T_env = df["T_env"].to_numpy()

    # Shell temperature [K].
    T_sh = df["T_sh"].to_numpy()

    # Shell perimeter [m].
    P_env = 2 * np.pi * df["R_sh"].to_numpy()

    # Compute contributions from different terms.
    con = h_env * (T_sh - T_env)
    rad = e_env * ct.stefan_boltzmann * (T_sh**4 - T_env**4)

    # Compute heat flux density [kW/m].
    qdot = P_env * (con + rad) / 1000

    # Integrate data.
    q = cumtrapz(qdot, z, initial=0.0)

    plot_fluxes(z, qdot, q)
