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
name = "KNO3"
data_path = Path("./data")
fig_path = Path("./figures")

kno3_fig_path = fig_path / name
kno3_alpha_path = Path("./data/KNO3-alpha")
kno3_beta_path = Path("./data/KNO3-beta")
kno3_gamma_path = Path("./data/KNO3-gamma")

if not os.path.exists(data_path):
    mkdir(data_path)
if not os.path.exists(fig_path):
    mkdir(fig_path)
if not os.path.exists(kno3_fig_path):
    mkdir(kno3_fig_path)

data_alpha_a = pd.read_csv(kno3_alpha_path / "KNO3_Cyclic_alpha_a.asc", skiprows=1, names = ["T", "a", "delta-a"], sep=r"\s+")
data_alpha_c = pd.read_csv(kno3_alpha_path / "KNO3_Cyclic_alpha_c.asc", skiprows=1, names = ["T", "c", "delta-c"], sep=r"\s+")
data_beta_a = pd.read_csv(kno3_beta_path / "KNO3_cyclic-beta-a.asc", skiprows=1, names = ["T", "a", "delta-a"], sep=r"\s+")
data_beta_c = pd.read_csv(kno3_beta_path / "KNO3_cyclic-beta_c.asc", skiprows=1, names = ["T", "c", "delta-c"], sep=r"\s+")
data_gamma_a = pd.read_csv(kno3_gamma_path / "KNO3_cyclic-gamma-a.asc", skiprows=1, names = ["T", "a", "delta-a"], sep=r"\s+")
data_gamma_c = pd.read_csv(kno3_gamma_path / "KNO3_cyclic-gamma-c.asc", skiprows=1, names = ["T", "c", "delta-c"], sep=r"\s+")

# alpha
## a
a_param_opt, a_param_cov = curve_fit(gerade_fit, xdata=data_alpha_a["T"] + 273.15, ydata=data_alpha_a["a"], sigma=data_alpha_a["delta-a"], absolute_sigma=True)
a_param_err = np.sqrt(np.diag(a_param_cov))

a_alpha_therm = ufloat(a_param_opt[0], a_param_err[0]) / ufloat(a_param_opt[1], a_param_err[1])

## c
c_param_opt, c_param_cov = curve_fit(gerade_fit, xdata=data_alpha_c["T"] + 273.15, ydata=data_alpha_c["c"], sigma=data_alpha_c["delta-c"], absolute_sigma=True)
c_param_err = np.sqrt(np.diag(c_param_cov))

c_alpha_term = ufloat(c_param_opt[0], c_param_err[0]) / ufloat(c_param_opt[1], c_param_err[1])

## plot a
x = np.linspace(start=17 + 273.15, stop=135 + 273.15, num=10)

plt.figure(figsize=[10, 4])
plt.xlabel("T in K")
plt.ylabel("a in \u00C5")

plt.errorbar(x=data_alpha_a["T"] + 273.15, y=data_alpha_a["a"], yerr=data_alpha_a["delta-a"], fmt="o", capsize=5, label="a")
plt.plot(x, gerade_fit(x, a_param_opt[0], a_param_opt[1]), label="a fit")
# Place text at T=100 and a=5.48 (adjust based on your actual data range)
plt.text(350, 5.414, f"$\\alpha_{{therm\ of\ a}}$ = {a_alpha_therm:.3e} / K", fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

plt.legend(loc="upper left")
plt.savefig(kno3_fig_path / "KNO3_alpha_a.png", dpi=300, bbox_inches='tight')
plt.close()

## plot c
plt.figure(figsize=[10, 4])
plt.xlabel("T in K")
plt.ylabel("c in \u00C5")

plt.errorbar(x=data_alpha_c["T"] + 273.15, y=data_alpha_c["c"], yerr=data_alpha_c["delta-c"], fmt="o", capsize=5, label="c")
plt.plot(x, gerade_fit(x, c_param_opt[0], c_param_opt[1]), label="c fit")
# Place text at T=100 and a=5.48 (adjust based on your actual data range)
plt.text(350, 6.44, f"$\\alpha_{{therm\ of\ c}}$ = {a_alpha_therm:.3e} / K", fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

plt.legend(loc="upper left")
plt.savefig(kno3_fig_path / "KNO3_alpha_c.png", dpi=300, bbox_inches='tight')
plt.close()

# beta
## a
a_param_opt, a_param_cov = curve_fit(gerade_fit, xdata=data_beta_a["T"] + 273.15, ydata=data_beta_a["a"], sigma=data_beta_a["delta-a"], absolute_sigma=True)
a_param_err = np.sqrt(np.diag(a_param_cov))

a_alpha_therm = ufloat(a_param_opt[0], a_param_err[0]) / ufloat(a_param_opt[1], a_param_err[1])

## c
c_param_opt, c_param_cov = curve_fit(gerade_fit, xdata=data_beta_c["T"] + 273.15, ydata=data_beta_c["c"], sigma=data_beta_c["delta-c"], absolute_sigma=True)
c_param_err = np.sqrt(np.diag(c_param_cov))

c_alpha_term = ufloat(c_param_opt[0], c_param_err[0]) / ufloat(c_param_opt[1], c_param_err[1])

## plot a
x = np.linspace(start=125 + 273.15, stop=147 + 273.15, num=10)

plt.figure(figsize=[10, 4])
plt.xlabel("T in K")
plt.ylabel("a in \u00C5")

plt.errorbar(x=data_beta_a["T"] + 273.15, y=data_beta_a["a"], yerr=data_beta_a["delta-a"], fmt="o", capsize=5, label="a")
plt.plot(x, gerade_fit(x, a_param_opt[0], a_param_opt[1]), label="a fit")
# Place text at T=100 and a=5.48 (adjust based on your actual data range)
plt.text(410, 5.42, f"$\\alpha_{{therm\ of\ a}}$ = {a_alpha_therm:.3e} / K", fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

plt.legend(loc="upper left")
plt.savefig(kno3_fig_path / "KNO3_beta_a.png", dpi=300, bbox_inches='tight')
plt.close()

## plot c
plt.figure(figsize=[10, 4])
plt.xlabel("T in K")
plt.ylabel("c in \u00C5")

plt.errorbar(x=data_beta_c["T"] + 273.15, y=data_beta_c["c"], yerr=data_beta_c["delta-c"], fmt="o", capsize=5, label="c")
plt.plot(x, gerade_fit(x, c_param_opt[0], c_param_opt[1]), label="c fit")
# Place text at T=100 and a=5.48 (adjust based on your actual data range)
plt.text(410, 9.71, f"$\\alpha_{{therm\ of\ c}}$ = {a_alpha_therm:.3e} / K", fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

plt.legend(loc="upper left")
plt.savefig(kno3_fig_path / "KNO3_beta_c.png", dpi=300, bbox_inches='tight')
plt.close()

# gamma
## a
a_param_opt, a_param_cov = curve_fit(gerade_fit, xdata=data_gamma_a["T"] + 273.15, ydata=data_gamma_a["a"], sigma=data_gamma_a["delta-a"], absolute_sigma=True)
a_param_err = np.sqrt(np.diag(a_param_cov))

a_alpha_therm = ufloat(a_param_opt[0], a_param_err[0]) / ufloat(a_param_opt[1], a_param_err[1])

## c
c_param_opt, c_param_cov = curve_fit(gerade_fit, xdata=data_gamma_c["T"] + 273.15, ydata=data_gamma_c["c"], sigma=data_gamma_c["delta-c"], absolute_sigma=True)
c_param_err = np.sqrt(np.diag(c_param_cov))

c_alpha_term = ufloat(c_param_opt[0], c_param_err[0]) / ufloat(c_param_opt[1], c_param_err[1])

## plot a
x = np.linspace(start=70 + 273.15, stop=123 + 273.15, num=10)

plt.figure(figsize=[10, 4])
plt.xlabel("T in K")
plt.ylabel("a in \u00C5")

plt.errorbar(x=data_gamma_a["T"] + 273.15, y=data_gamma_a["a"], yerr=data_gamma_a["delta-a"], fmt="o", capsize=5, label="a")
plt.plot(x, gerade_fit(x, a_param_opt[0], a_param_opt[1]), label="a fit")
# Place text at T=100 and a=5.48 (adjust based on your actual data range)
plt.text(370, 5.46, f"$\\alpha_{{therm\ of\ a}}$ = {a_alpha_therm:.3e} / K", fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

plt.legend(loc="upper left")
plt.savefig(kno3_fig_path / "KNO3_gamma_a.png", dpi=300, bbox_inches='tight')
plt.close()

## plot c
plt.figure(figsize=[10, 4])
plt.xlabel("T in K")
plt.ylabel("c in \u00C5")

plt.errorbar(x=data_gamma_c["T"] + 273.15, y=data_gamma_c["c"], yerr=data_gamma_c["delta-c"], fmt="o", capsize=5, label="c")
plt.plot(x, gerade_fit(x, c_param_opt[0], c_param_opt[1]), label="c fit")
# Place text at T=100 and a=5.48 (adjust based on your actual data range)
plt.text(370, 9.05, f"$\\alpha_{{therm\ of\ c}}$ = {a_alpha_therm:.3e} / K", fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

plt.legend(loc="upper left")
plt.savefig(kno3_fig_path / "KNO3_gamma_c.png", dpi=300, bbox_inches='tight')
plt.close()
