# -*- coding: utf-8 -*-
import json
import matplotlib.pyplot as plt
import pandas as pd


with open("results/results-barr.json") as fp:
    df = pd.DataFrame(json.load(fp))

ref = pd.read_csv("data/data-barr.csv")
ref = ref.loc[ref["Experiment ID"] == "T4"]

Tg_off_wall = ["Tg_off_wall 1", "Tg_off_wall 2", "Tg_off_wall 3",
               "Tg_off_wall 4", "Tg_off_wall 5", "Tg_off_wall 6",
               "Tg_off_wall 7", "Tg_off_wall 8", "Tg_off_wall 9"]

Tg_off_bed = ["Tg_off_bed 1", "Tg_off_bed 2", "Tg_off_bed 3",
              "Tg_off_bed 4", "Tg_off_bed 5", "Tg_off_bed 6",
              "Tg_off_bed 7", "Tg_off_bed 8", "Tg_off_bed 9"]

Ts = ["Ts1", "Ts2", "Ts3", "Ts4", "Ts5", "Ts6",
      "Ts7", "Ts8", "Ts9", "Ts10", "Ts11", "Ts12"]

Tw = ["Tw1", "Tw2", "Tw3", "Tw4", "Tw5", "Tw6", "Tw7", "Tw8"]

zg = [0.11, 0.89, 2.15, 2.51, 2.85, 3.20, 3.91, 4.44, 4.95]

zs = [0.11, 0.87, 1.44, 2.14, 2.50, 2.85,
      3.20, 3.91, 4.48, 4.95, 5.25, 5.50]

zw = [1.33, 2.32, 2.67, 3.03, 3.38, 3.83, 4.40, 4.78]

Tg1 = ref[Tg_off_wall].T.to_numpy().ravel() - 273.15
Tg2 = ref[Tg_off_bed].T.to_numpy().ravel() - 273.15
Ts1 = ref[Ts].T.to_numpy().ravel() - 273.15
Tw1 = ref[Tw].T.to_numpy().ravel() - 273.15

z = df["z"].to_numpy()
dz = z[-1] / z.shape[0]

q_cgw = df["q_cgw"][::-1] / (1000.0 * dz)
q_cgb = df["q_cgb"][::-1] / (1000.0 * dz)
q_cwb = df["q_cwb"][::-1] / (1000.0 * dz)
q_rgw = df["q_rgw"][::-1] / (1000.0 * dz)
q_rgb = df["q_rgb"][::-1] / (1000.0 * dz)
q_rwb = df["q_rwb"][::-1] / (1000.0 * dz)

T_g = df["T_g"][::-1] - 273.15
T_b = df["T_b"][::-1] - 273.15
T_w = df["T_wg"][::-1] - 273.15

plt.close("all")
plt.style.use("classic")
plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.plot(z, T_g, "k", label=r"Freeboard")
plt.plot(z, T_b, "r", label=r"Bed")
plt.plot(z, T_w, "b", label=r"Wall")
plt.scatter(zg, Tg1, color="k", label="_none_")
plt.scatter(zg, Tg2, color="k", label="_none_")
plt.scatter(zs, Ts1, color="r", label="_none_")
plt.scatter(zw, Tw1, color="b", label="_none_")
plt.grid(linestyle=":")
plt.legend(loc=4, frameon=True, framealpha=1.0)
plt.xlabel("Position [m]")
plt.ylabel("Temperature [°C]")
plt.xlim(0, 5.5)
plt.ylim(0, 1000)

plt.subplot(1, 2, 2)
plt.plot(z, q_cgw, "g", label=r"$\dot{Q}_{CGW}$")
plt.plot(z, q_cgb, "b", label=r"$\dot{Q}_{CGB}$")
plt.plot(z, q_cwb, "r", label=r"$\dot{Q}_{CWB}$")
plt.plot(z, q_rgw, "m", label=r"$\dot{Q}_{RGW}$")
plt.plot(z, q_rgb, "k", label=r"$\dot{Q}_{RGB}$")
plt.plot(z, q_rwb, "c", label=r"$\dot{Q}_{RWB}$")
plt.grid(linestyle=":")
plt.legend(loc=1, frameon=True, framealpha=1.0)
plt.xlabel("Position [m]")
plt.ylabel("Heat flux [kW]")
plt.xlim(0, 5.5)
plt.ylim(-1, 6)

plt.tight_layout()
plt.savefig("../../docs/media/validation-barr", dpi=300)
