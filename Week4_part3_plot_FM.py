# PRA2003 - Week 4 Deliverable (Part 3, figure): Figure 1 of the README
# Draws the asymmetry A = WT - mutant (bacteria per event) of the 6
# wild-type/mutant pairs, with 1-sigma error bars. The numbers are not
# calculated here: they are read from three_sigma_results_FM.csv, made by
# Week4_part3_asymmetry_FM.r (the uncertainty is the spread of the 10
# sub-samples, see that script).
#
# HOW TO RUN:
#   python Week4_part3_plot_FM.py      (after Week4_part3_asymmetry_FM.r)
#
# NEEDS:    Python 3.9 or newer with pandas, numpy and matplotlib
# PRODUCES: asymmetry_plot_FM.png
# RUN TIME: a few seconds.
#
# The differences range from about 0.03 (E. coli) to 0.00004 (Salmonella)
# bacteria per event. On one shared axis the rare strains would sit on top of
# 0, so every strain gets its own panel with its own horizontal scale.

import os
import sys

import matplotlib
matplotlib.use("Agg")   #Draw to a file only, so no screen is needed
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

# ---- Settings ----
input_csv = "three_sigma_results_FM.csv"
plot_file = "asymmetry_plot_FM.png"
threshold = 3

pair_codes = [211, 321, 2212, 3122, 3312, 3334]
pair_names = ["E. coli", "B. subtilis", "P. aeruginosa",
              "S. pneumoniae", "M. tuberculosis", "Salmonella"]

colour_significant = "#2a78d6"
colour_not_significant = "#8b8a85"
colour_axis = "#52514e"

os.chdir(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isfile(input_csv):
    sys.exit(f"Can't find '{input_csv}' - run Week4_part3_asymmetry_FM.r first "
             f"(or keep the included CSV file in this folder).")

results = pd.read_csv(input_csv).set_index("code")
missing = [c for c in pair_codes if c not in results.index]
if missing:
    sys.exit(f"{input_csv} has no row for code(s) {missing}.")

fig, axes = plt.subplots(len(pair_codes), 1, figsize=(8, 7.5))
for ax, code, name in zip(axes, pair_codes, pair_names):
    row = results.loc[code]
    n_sigma = row["n_sigma"]
    #Every panel in units of 10^power bacteria per event, so its tick labels are short
    power = int(np.floor(np.log10(abs(row["difference"]) + row["uncertainty"])))
    a = row["difference"] / 10.0 ** power
    sigma = row["uncertainty"] / 10.0 ** power
    colour = colour_significant if n_sigma > threshold else colour_not_significant

    ax.axvline(0, color=colour_axis, linestyle="--", linewidth=1, zorder=1)
    ax.errorbar(a, 0, xerr=sigma, fmt="o", color=colour, ecolor=colour,
                elinewidth=2, capsize=4, markersize=8, zorder=3)
    ax.annotate(f"{n_sigma:.2f} sigma", (a + sigma, 0), xytext=(6, 0),
                textcoords="offset points", va="center", fontsize=9)

    #Range: the error bar and 0, with room on the right for the n sigma label
    low = min(0, a - sigma)
    high = max(0, a + sigma)
    span = high - low
    ax.set_xlim(low - 0.1 * span, high + 0.45 * span)
    ax.set_ylim(-1, 1)
    ax.set_yticks([])
    ax.set_ylabel(f"{name}\n(\u00d7 $10^{{{power}}}$)", rotation=0, ha="right", va="center", fontsize=11)
    ax.grid(True, axis="x", color="#e8e7e2", zorder=0)
    ax.set_axisbelow(True)
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)

axes[-1].set_xlabel("Mean asymmetry A = WT - mutant\n[average bacteria per event, times the factor under each name]")
fig.suptitle("Wild type minus mutant, 5M events (error bars: 1 sigma from the spread of 10 sub-samples)",
             fontsize=12)
key = [Line2D([], [], marker="o", linestyle="none", color=colour_significant,
              label=f"> {threshold} sigma: significant"),
       Line2D([], [], marker="o", linestyle="none", color=colour_not_significant,
              label=f"<= {threshold} sigma: not significant"),
       Line2D([], [], color=colour_axis, linestyle="--", label="A = 0")]
fig.legend(handles=key, loc="lower center", ncol=3, frameon=False, fontsize=9)
fig.tight_layout(rect=(0, 0.04, 1, 0.97))
fig.savefig(plot_file, dpi=200)
plt.close(fig)
print(f"Saved figure to {plot_file}")
