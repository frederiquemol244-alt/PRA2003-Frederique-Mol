# PRA2003 - Week 5 Deliverable (Part 3 of 3): the group's momentum figures (line plot and fit)
# Makes the two figures of the group for question 3, with the group's method.
# It is the Python version of the group's two R plotting scripts
#   "PlotAsymmetrypermomentum_Function momentum.R"  -> line plot
#   "UncertaintyMomentum.R"                         -> straight-line fit
# with the same 5 momentum bins, the same fit and the same layout (one panel
# per strain, 6 panels per figure). One thing is different: the y-axis is the
# AVERAGE difference per event, not the total number N_WT - N_mutant. The counting that the group does in
# "Asymmetry_Funtionofmomentum.R" is done here by Week5_momentum_counting_FM.py.
#
# HOW TO RUN:
#   python Week5_group_figures_FM.py
#   momentum_sub_sample_results_FM.csv must be in the same folder as this
#   script (it is included, or made by Week5_momentum_counting_FM.py).
#
# NEEDS:    Python 3.9 or newer with pandas, numpy and matplotlib
#           (pip install pandas numpy matplotlib)
# PRODUCES: momentum_lineplot_FM.png          (line plot, 6 panels)
#           momentum_fit_FM.png               (straight-line fit, 6 panels)
#           momentum_difference_results_FM.csv (30 rows: 6 pairs x 5 bins)
#           momentum_fit_results_FM.csv        (6 rows: one fit per pair)
# RUN TIME: a few seconds.
#
# DEFINITIONS (per strain and momentum bin):
#   per set     D_i   = (N_WT - N_mutant) / non-empty events of set i: the
#                       average number of WT per event minus the average
#                       number of mutants per event in that set. (The group's
#                       R scripts use the total number N_WT - N_mutant; the
#                       y-axis has to be an average, so it is divided by the
#                       number of events here.)
#   average     D     = mean of the 10 values D_i, in bacteria per event
#   uncertainty sigma = standard deviation of the 10 values D_i (ddof = 1,
#                       like sd() in R), no division by sqrt(10)
#   significance      = D / sigma; |D / sigma| >= 3 means the point is
#                       significant (filled and labelled in the line plot)
#
# STRAIGHT-LINE FIT (per strain, through the 5 bins):
#   D = intercept + slope * bin      (bin = 1, 2, 3, 4, 5)
#   weighted least squares with weight w = 1 / sigma^2, so bins with a small
#   uncertainty count more:
#     S   = sum(w),   Sx  = sum(w x),   Sxx = sum(w x^2)
#     Sy  = sum(w y), Sxy = sum(w x y), Delta = S*Sxx - Sx^2
#     slope     = (S*Sxy - Sx*Sy) / Delta     +- sqrt(S   / Delta)
#     intercept = (Sxx*Sy - Sx*Sxy) / Delta   +- sqrt(Sxx / Delta)
#   slope / uncertainty = significance of the slope. If |slope| >= 3 sigma,
#   D changes with momentum. chi2 / ndf tells whether a straight line
#   describes the points well (around 1 = good fit).
#
# CONSTANCY TEST (from "Asymmetry_Funtionofmomentum.R"): chi-square of the 5
# bins against one constant D (the weighted average of the bins), 4 degrees
# of freedom. p < 0.003 (3 sigma) means D is not the same in every bin.
#
# DIFFERENCE WITH THE R SCRIPTS: output-Set10.txt stores (pT, phi, eta)
# instead of (px, py, pz). The counting script uses |p| = pT * cosh(eta) for
# that file; the R script uses sqrt(px^2 + py^2 + pz^2) for every file.

import os
import sys

import matplotlib
matplotlib.use("Agg")   #Draw to files only, so no screen is needed
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

# ---- Settings ----
input_csv = "momentum_sub_sample_results_FM.csv"
results_csv = "momentum_difference_results_FM.csv"
fit_csv = "momentum_fit_results_FM.csv"
lineplot_png = "momentum_lineplot_FM.png"
fit_png = "momentum_fit_FM.png"
threshold = 3            #Significant if at least this many sigma from zero
p_value_limit = 0.003    #Constancy test: 3 sigma

pair_codes = [211, 321, 2212, 3122, 3312, 3334]
pair_names = ["E. coli", "B. subtilis", "P. aeruginosa",
              "S. pneumoniae", "M. tuberculosis", "Salmonella"]
momentum_unit = r"$10^{-20}$ kg m/s"

colour_line = "#2F5E9E"        #Points and connecting line (the group's blue)
colour_fit = "#C0392B"         #Fitted line (the group's red)
colour_error_bar = "#4D4D4D"
colour_zero = "#7F7F7F"
figure_size = (12.9, 7.9)      #The group's 1800 x 1100 pixels at 140 dpi

#Protection: settings must line up 1-to-1
if len(pair_codes) != len(pair_names):
    sys.exit("pair_codes and pair_names must be the same length.")


def chi2_p_value(chi2, ndf):
    """P(chi2 >= observed) for an even number of degrees of freedom (no scipy needed)."""
    if ndf <= 0 or ndf % 2 != 0:
        return np.nan
    half = chi2 / 2
    return float(np.exp(-half) * sum(half ** i / np.prod(range(1, i + 1)) for i in range(ndf // 2)))


# ---- Read the per-file counts ----
os.chdir(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isfile(input_csv):
    sys.exit(f"Can't find '{input_csv}' - run Week5_momentum_counting_FM.py first "
             f"(or keep the included CSV file in this folder).")

counts = pd.read_csv(input_csv)
file_names = list(dict.fromkeys(counts["sub_sample"]))      #In the order of the CSV
bin_labels = list(dict.fromkeys(counts["momentum_bin"]))    #In the order of the CSV
n_sets = len(file_names)
bins = np.arange(1, len(bin_labels) + 1)                    #Bin number 1 ... 5, the x of the fit
print(f"Read {len(counts)} rows from {input_csv} ({n_sets} sets, {len(bin_labels)} momentum bins).\n")

#Average difference WT - mutant per event, per set, pair and bin
counts["D"] = (counts["count_X"] - counts["count_minusX"]) / counts["n_events"]

# ---- Average and uncertainty per pair and bin, then the fit per pair ----
rows = []
fit_rows = []
for code, name in zip(pair_codes, pair_names):
    y = np.zeros(len(bin_labels))
    sigma = np.zeros(len(bin_labels))
    for k, label in enumerate(bin_labels):
        cell = counts[(counts["pair_code"] == code) & (counts["momentum_bin"] == label)]
        if len(cell) != n_sets:
            sys.exit(f"{name} {label}: expected {n_sets} sets, found {len(cell)}.")
        d_per_set = cell["D"].to_numpy(dtype=float)
        y[k] = d_per_set.mean()
        sigma[k] = d_per_set.std(ddof=1)
        rows.append({"pair_code": code, "pair_name": name, "bin": k + 1, "momentum_bin": label,
                     "N_WT_total": int(cell["count_X"].sum()),
                     "N_mutant_total": int(cell["count_minusX"].sum()),
                     "D": y[k], "uncertainty": sigma[k]})

    #The fit and the significance divide by sigma
    if np.any(sigma <= 0):
        sys.exit(f"{name}: a bin has the same D in every set (sigma = 0), so it cannot be weighted.")

    #Weighted straight-line fit: D = intercept + slope * bin
    w = 1 / sigma ** 2
    S = w.sum()
    Sx = (w * bins).sum()
    Sxx = (w * bins ** 2).sum()
    Sy = (w * y).sum()
    Sxy = (w * bins * y).sum()
    Delta = S * Sxx - Sx ** 2

    slope = (S * Sxy - Sx * Sy) / Delta
    slope_err = np.sqrt(S / Delta)
    intercept = (Sxx * Sy - Sx * Sxy) / Delta
    intercept_err = np.sqrt(Sxx / Delta)
    chi2 = (((y - (intercept + slope * bins)) / sigma) ** 2).sum()
    ndf = len(bins) - 2

    #Constancy test: chi-square against the weighted average D of the bins
    d_constant = Sy / S
    chi2_constant = (((y - d_constant) / sigma) ** 2).sum()
    ndf_constant = len(bins) - 1
    p_value = chi2_p_value(chi2_constant, ndf_constant)

    fit_rows.append({"pair_code": code, "pair_name": name,
                     "slope_per_bin": slope, "slope_uncertainty": slope_err,
                     "slope_significance": slope / slope_err,
                     "slope_significant": bool(abs(slope / slope_err) >= threshold),
                     "intercept": intercept, "intercept_uncertainty": intercept_err,
                     "chi2": chi2, "ndf": ndf,
                     "constant_D": d_constant, "chi2_constant": chi2_constant,
                     "ndf_constant": ndf_constant, "p_value_constant": p_value,
                     "depends_on_momentum": bool(p_value < p_value_limit)})

results = pd.DataFrame(rows)
results["significance"] = results["D"] / results["uncertainty"]
results["significant"] = results["significance"].abs() >= threshold
fits = pd.DataFrame(fit_rows)

# ---- Save (rounded only here, all calculations use the exact values) ----
saved = results.copy()
saved["D"] = [float(f"{v:.4g}") for v in saved["D"]]                     #4 significant digits
saved["uncertainty"] = [float(f"{v:.3g}") for v in saved["uncertainty"]]
saved["significance"] = saved["significance"].round(2)
saved.to_csv(results_csv, index=False)

saved_fits = fits.copy()
for column in ["slope_per_bin", "slope_uncertainty", "intercept", "intercept_uncertainty", "constant_D"]:
    saved_fits[column] = [float(f"{v:.4g}") for v in saved_fits[column]]   #4 significant digits
for column in ["slope_significance", "chi2", "chi2_constant"]:
    saved_fits[column] = saved_fits[column].round(2)
saved_fits["p_value_constant"] = [float(f"{p:.3g}") for p in saved_fits["p_value_constant"]]
saved_fits.to_csv(fit_csv, index=False)

# ---- Print readable tables ----
with pd.option_context("display.width", 200, "display.max_rows", 200):
    print(f"D = average WT per event - average mutant per event, per momentum bin: mean of the {n_sets} sets,")
    print(f"uncertainty = standard deviation of the {n_sets} sets (1 sigma):\n")
    print(saved[["pair_name", "momentum_bin", "D", "uncertainty", "significance", "significant"]]
          .to_string(index=False))
    print("\nStraight-line fit  D = intercept + slope * bin  (bin = 1 ... 5):\n")
    print(saved_fits[["pair_name", "slope_per_bin", "slope_uncertainty", "slope_significance",
                      "slope_significant", "intercept", "intercept_uncertainty", "chi2", "ndf"]]
          .to_string(index=False))
    print(f"\nIs D the same in every bin? (p < {p_value_limit} = it is not, at 3 sigma):\n")
    print(saved_fits[["pair_name", "constant_D", "chi2_constant", "ndf_constant",
                      "p_value_constant", "depends_on_momentum"]]
          .to_string(index=False, formatters={"p_value_constant": "{:.3g}".format}))


# ---- The parts both figures share: axes, zero line and error bars ----
def unit_of(values):
    """Power of ten to show a strain in, so the axis has numbers like 2.5 and not 0.0025."""
    return 10.0 ** np.floor(np.log10(np.max(np.abs(values))))


def unit_text(unit):
    return rf"10$^{{{int(round(np.log10(unit)))}}}$"


def draw_panel(ax, y, sigma, y_low, y_high, unit):
    ax.axhline(0, color=colour_zero, linestyle="--", linewidth=1, zorder=1)
    ax.errorbar(bins, y, yerr=sigma, fmt="none", ecolor=colour_error_bar,
                elinewidth=1.2, capsize=4, zorder=2)
    ax.set_xlim(0.7, len(bins) + 0.3)
    ax.set_ylim(y_low, y_high)
    ax.set_xticks(bins)
    ax.set_xticklabels(bin_labels, fontsize=8)
    ax.set_xlabel(f"|p| bin ({momentum_unit})")
    ax.set_ylabel("Mean difference D = WT - mutant\n[" + unit_text(unit) + " bacteria per event]")


def add_legend(fig, handles, labels):
    """One legend for all 6 panels, in the bottom margin."""
    fig.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False, fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 1), h_pad=2.0, w_pad=2.0)


# ---- Figure 1: line plot of D against the momentum bin ----
fig, axes = plt.subplots(2, 3, figsize=figure_size)
for ax, code, name in zip(axes.flat, pair_codes, pair_names):
    pair = results[results["pair_code"] == code]
    unit = unit_of(pair["D"].to_numpy() + np.sign(pair["D"].to_numpy()) * pair["uncertainty"].to_numpy())
    y = pair["D"].to_numpy() / unit
    sigma = pair["uncertainty"].to_numpy() / unit
    z = pair["significance"].to_numpy()
    significant = pair["significant"].to_numpy()

    #y range: all error bars and the zero line, with room above for the sigma labels
    low = min(0, np.min(y - sigma))
    high = max(0, np.max(y + sigma))
    span = high - low
    draw_panel(ax, y, sigma, low - 0.05 * span, high + 0.15 * span, unit)

    ax.plot(bins, y, color=colour_line, linewidth=2, zorder=3)
    #Open point = not significant, filled point = at least 3 sigma from zero
    ax.plot(bins[~significant], y[~significant], "o", markerfacecolor="white",
            markeredgecolor=colour_line, markeredgewidth=2, markersize=8, zorder=4)
    ax.plot(bins[significant], y[significant], "o", color=colour_line,
            markeredgewidth=2, markersize=8, zorder=4)
    for xi, top, zi in zip(bins[significant], (y + sigma)[significant], z[significant]):
        ax.annotate(f"{zi:.1f} sigma", (xi, top), xytext=(0, 4), textcoords="offset points",
                    ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_title(name, fontstyle="italic", fontsize=13)

add_legend(fig,
           [Line2D([], [], marker="o", linestyle="none", markerfacecolor="white",
                   markeredgecolor=colour_line, markeredgewidth=2, markersize=8),
            Line2D([], [], marker="o", linestyle="none", color=colour_line, markersize=8),
            Line2D([], [], color=colour_error_bar, linewidth=1.2),
            Line2D([], [], color=colour_zero, linestyle="--", linewidth=1)],
           ["Not significant", f"Significant (>= {threshold} sigma from zero)",
            f"+/- 1 SD over the {n_sets} sets", "D = 0 (no difference)"])
fig.savefig(lineplot_png, dpi=200)
plt.close(fig)
print(f"\nSaved figure to {lineplot_png}")

# ---- Figure 2: the same points with the straight-line fit ----
x_line = np.array([0.7, len(bins) + 0.3])   #The line is drawn over the whole panel
fig, axes = plt.subplots(2, 3, figsize=figure_size)
for ax, code, name in zip(axes.flat, pair_codes, pair_names):
    pair = results[results["pair_code"] == code]
    fit = fits[fits["pair_code"] == code].iloc[0]
    unit = unit_of(pair["D"].to_numpy() + np.sign(pair["D"].to_numpy()) * pair["uncertainty"].to_numpy())
    y = pair["D"].to_numpy() / unit
    sigma = pair["uncertainty"].to_numpy() / unit
    y_line = (fit["intercept"] + fit["slope_per_bin"] * x_line) / unit

    #y range: all error bars, the zero line and the fitted line
    low = min(0, np.min(y - sigma), y_line.min())
    high = max(0, np.max(y + sigma), y_line.max())
    span = high - low
    draw_panel(ax, y, sigma, low - 0.04 * span, high + 0.04 * span, unit)

    ax.plot(x_line, y_line, color=colour_fit, linewidth=2, zorder=3)
    ax.plot(bins, y, "o", markerfacecolor="white", markeredgecolor=colour_line,
            markeredgewidth=2, markersize=8, zorder=4)
    #Two title lines: the strain, and the slope with its significance
    ax.set_title(name, fontstyle="italic", fontsize=13, pad=22)
    ax.text(0.5, 1.03, f"slope = {fit['slope_per_bin'] / unit:.2f} +/- {fit['slope_uncertainty'] / unit:.2f} "
                       f"per bin ({fit['slope_significance']:.1f} sigma)",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=10)

add_legend(fig,
           [Line2D([], [], marker="o", linestyle="none", markerfacecolor="white",
                   markeredgecolor=colour_line, markeredgewidth=2, markersize=8),
            Line2D([], [], color=colour_error_bar, linewidth=1.2),
            Line2D([], [], color=colour_fit, linewidth=2),
            Line2D([], [], color=colour_zero, linestyle="--", linewidth=1)],
           [f"Average of the {n_sets} sets", "+/- 1 SD", "Straight-line fit", "D = 0 (no difference)"])
fig.savefig(fit_png, dpi=200)
plt.close(fig)
print(f"Saved figure to {fit_png}")
print(f"Saved {results_csv} and {fit_csv}")
