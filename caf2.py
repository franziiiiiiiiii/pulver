import os
from os import mkdir
from pathlib import Path

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat

def gerade_fit(x, a, b):
    return a * x + b

# PATHING
name = "CaF2"
data_path = Path("./data")
fig_path = Path("./figures")

caf2_fig_path = fig_path / name
caf2_path = data_path / name

if not os.path.exists(data_path):
    mkdir(data_path)
if not os.path.exists(fig_path):
    mkdir(fig_path)
if not os.path.exists(caf2_fig_path):
    mkdir(caf2_fig_path)

# DATA INIT CaF2
data_a = pd.read_csv(caf2_path / "CaF2_a.asc", skiprows=1, names = ["T", "a", "delta-a"], sep=r"\s+")
data_bckg1 = pd.read_csv(caf2_path / "CaF2_bckg1.asc", skiprows=1, names=["T", "bckg1", "delta-bckg1"], sep=r"\s+")
data_bckg2 = pd.read_csv(caf2_path / "CaF2_bckg2.asc", skiprows=1, names=["T", "bckg2", "delta-bckg2"], sep=r"\s+")
data_gw = pd.read_csv(caf2_path / "CaF2_GW.asc", skiprows=1, names=["T", "gw", "delta-gw"], sep=r"\s+")
data_ly = pd.read_csv(caf2_path / "CaF2_LY.asc", skiprows=1, names=["T", "ly", "delta-ly"], sep=r"\s+")
data_shift = pd.read_csv(caf2_path / "CaF2_shift.asc", skiprows=1, names=["T", "shift", "delta-shift"], sep=r"\s+")

# PLOTTING CaF2

# a
param_opt, param_cov = curve_fit(gerade_fit, xdata=data_a["T"] + 273.15, ydata=data_a["a"], sigma=data_a["delta-a"], absolute_sigma=True)
param_err = np.sqrt(np.diag(param_cov))

alpha_therm = ufloat(param_opt[0], param_err[0]) / ufloat(param_opt[1], param_err[1])

x = np.linspace(start=20 + 273.15, stop=750 + 273.15, num=10)

plt.figure(figsize=[20, 4])
plt.xlabel("T in K")
plt.ylabel("a in \u00C5")
plt.errorbar(x=data_a["T"] + 273.15, y=data_a["a"], yerr=data_a["delta-a"], fmt="o", capsize=5, label="a")
plt.plot(x, gerade_fit(x, param_opt[0], param_opt[1]), label="fit")
plt.legend(loc="upper left")

# Place text at T=100 and a=5.48 (adjust based on your actual data range)
plt.text(600, 5.46, f"$\\alpha_{{therm}}$ = {alpha_therm:.3e} / K", fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

plt.savefig(caf2_fig_path / "CaF2_a.png", dpi=300, bbox_inches='tight')
plt.close()

# a take 6 values
take = 7
param_opt, param_cov = curve_fit(gerade_fit, xdata=data_a["T"][0:take] + 273.15, ydata=data_a["a"][0:take], sigma=data_a["delta-a"][0:take], absolute_sigma=True)
param_err = np.sqrt(np.diag(param_cov))

alpha_therm = ufloat(param_opt[0], param_err[0]) / ufloat(param_opt[1], param_err[1])

stop = 20 + (take * 10) + 273.15
x = np.linspace(start=30 + 273.15, stop=stop, num=10)
plt.figure(figsize=[20, 4])
plt.xlabel("T in K")
plt.ylabel("a in \u00C5")
plt.errorbar(x=data_a["T"][0:take] + 273.15, y=data_a["a"][0:take], yerr=data_a["delta-a"][0:take], fmt="o", capsize=5, label="a")
plt.plot(x, gerade_fit(x, param_opt[0], param_opt[1]), label="fit")
plt.legend(loc="upper left")

# Place text at T=100 and a=5.48 (adjust based on your actual data range)
plt.text(stop - 10, 5.462, f"$\\alpha_{{therm}}$ = {alpha_therm:.3e} / K", fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

plt.savefig(caf2_fig_path / "CaF2_a_first5.png", dpi=300, bbox_inches='tight')
plt.close()

# shift
plt.figure(figsize=[20, 4])
plt.xlabel("T in K")
plt.ylabel("$\u03B8_{0}$")
plt.errorbar(x=data_shift["T"] + 273.15, y=data_shift["shift"], yerr=data_shift["delta-shift"], fmt="o", capsize=5, label="a")

median_val = np.median(data_shift["shift"])
plt.axhline(median_val, color="k", linestyle="--")

plt.savefig(caf2_fig_path / "CaF2_shift.png", dpi=300, bbox_inches='tight')
plt.close()

# gw ly
plt.figure(figsize=[20, 4])
plt.xlabel("T in K")
plt.ylabel("GW- und LY-Anteile")
plt.errorbar(x=data_gw["T"] + 273.15, y=data_gw["gw"], yerr=data_gw["delta-gw"], fmt="o", capsize=5, label="GW")
plt.errorbar(x=data_ly["T"] + 273.15, y=data_ly["ly"], yerr=data_ly["delta-ly"], fmt="o", capsize=5, label="LY")
plt.legend(loc="center left")


plt.savefig(caf2_fig_path / "CaF2_gw_ly.png", dpi=300, bbox_inches='tight')
plt.close()

# bckg
plt.figure(figsize=[20, 4])
plt.xlabel("T in K")
plt.ylabel("Background")
plt.errorbar(x=data_bckg1["T"] + 273.15, y=data_bckg1["bckg1"], yerr=data_bckg1["delta-bckg1"], fmt="o", capsize=5, label="Bckg1")
plt.errorbar(x=data_bckg2["T"] + 273.15, y=data_bckg2["bckg2"], yerr=data_bckg2["delta-bckg2"], fmt="o", capsize=5, label="Bckg2")
plt.legend(loc="center left")


plt.savefig(caf2_fig_path / "CaF2_bckg.png", dpi=300, bbox_inches='tight')
plt.close()
