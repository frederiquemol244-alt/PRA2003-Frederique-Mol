# PRA2003 - Week 5 Deliverable (Part 2 of 3): asymmetry as a function of momentum
# Reads the per-file counts per momentum bin (from Part 1) and, for each of
# the 6 wild-type/mutant pairs and each momentum bin, calculates the
# asymmetry with its statistical uncertainty from the sub-sampling method.
# Then it draws one figure with a scatter plot per pair (asymmetry against
# momentum, with error bars), all 6 together.
#
# HOW TO RUN:
#   python Week5_asymmetry_vs_momentum_FM.py
#   momentum_sub_sample_results_FM.csv must be in the same folder as this
#   script (it is included, or made by Week5_momentum_counting_FM.py).
#
# NEEDS:    Python 3.9 or newer with pandas, numpy and matplotlib
#           (pip install pandas numpy matplotlib)
# PRODUCES: asymmetry_vs_momentum_FM.csv (30 rows: 6 pairs x 5 bins),
#           asymmetry_momentum_all_FM.png (the 6 scatter plots in one figure)
# RUN TIME: a few seconds.
#
# DEFINITIONS (per momentum bin):
#   per file    d_i   = (count_X - count_minusX) / n_events of file i: the
#                       average number of X per event minus the average
#                       number of -X per event in that file (a plain
#                       difference, no division by the sum)
#   asymmetry   A     = the AVERAGE of the 10 values d_i, in bacteria per
#                       event. This is what is on the y-axis: an average per
#                       event, not a total number of bacteria.
#   uncertainty sigma = standard deviation of the 10 values d_i (ddof = 1,
#                       like sd() in R), no division by sqrt(10), as in the
#                       lecture and in week 4
#   n_sigma           = |A| / sigma
#   X and -X are counted in the same events, so taking the difference inside
#   each file also takes their correlation into account.
#
# MOMENTUM BINS: the 5 bins used by the group (unit 10^-20 kg m/s):
#   0-0.654, 0.654-1.32, 1.32-2.68, 2.68-5.88 and > 5.88. Each bin holds
#   roughly the same number of bacteria, and every bacterium is in one bin.
#   The bins have very different widths and the last one has no upper edge,
#   so in the plots they are drawn next to each other at equal distance, with
#   the momentum range of each bin written on the axis.
#
# Decision rule (fixed before looking at the numbers, same as week 4):
#   n_sigma > 3  -> significant asymmetry in that bin
#   n_sigma <= 3 -> not significant (consistent with zero)

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
output_csv = "asymmetry_vs_momentum_FM.csv"
week4_pairs_csv = "three_sigma_results_FM.csv"   #Optional: used to check the totals
threshold = 3

pair_codes = [211, 321, 2212, 3122, 3312, 3334]
pair_names = ["E. coli", "B. subtilis", "P. aeruginosa",
              "S. pneumoniae", "M. tuberculosis", "Salmonella"]
file_parts = ["ecoli", "bsubtilis", "paeruginosa",
              "spneumoniae", "mtuberculosis", "salmonella"]

#Group bins [lower, upper), unit 10^-20 kg m/s; the last bin is open ended
bin_labels = ["0-0.654", "0.654-1.32", "1.32-2.68", "2.68-5.88", ">5.88"]
bin_low = [0, 0.654, 1.32, 2.68, 5.88]
bin_high = [0.654, 1.32, 2.68, 5.88, np.inf]
momentum_unit = r"$10^{-20}$ kg m/s"

colour_significant = "#2a78d6"
colour_not_significant = "#8b8a85"
colour_error_bar = "#52514e"

#Protection: settings must line up 1-to-1
if not (len(pair_codes) == len(pair_names) == len(file_parts)):
    sys.exit("pair_codes, pair_names and file_parts must be the same length.")
if not (len(bin_labels) == len(bin_low) == len(bin_high)):
    sys.exit("The bin settings must all have the same length.")


def signif(x, digits):
    """Round to a number of significant digits (like signif() in R)."""
    x = np.asarray(x, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        magnitude = np.where(x == 0, 0, np.floor(np.log10(np.abs(x))))
    factor = 10.0 ** (digits - 1 - magnitude)
    return np.round(x * factor) / factor


# ---- Read the per-file counts ----
os.chdir(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isfile(input_csv):
    sys.exit(f"Can't find '{input_csv}' - run Week5_momentum_counting_FM.py first "
             f"(or keep the included CSV file in this folder).")

counts = pd.read_csv(input_csv)
file_names = list(dict.fromkeys(counts["sub_sample"]))   #In the order of the CSV
print(f"Read {len(counts)} rows from {input_csv} ({len(file_names)} sub-samples).\n")

# ---- Asymmetry, uncertainty and significance per pair and bin ----
rows = []
for code, name in zip(pair_codes, pair_names):
    for k, label in enumerate(bin_labels):
        cell = counts[(counts["pair_code"] == code) & (counts["momentum_bin"] == label)]
        if len(cell) != len(file_names):
            sys.exit(f"{name} {label}: expected {len(file_names)} files, found {len(cell)}.")

        count_x = cell["count_X"].to_numpy(dtype=np.int64)
        count_minus_x = cell["count_minusX"].to_numpy(dtype=np.int64)
        n_events = cell["n_events"].to_numpy(dtype=np.int64)

        #One difference per file; a bin with 0 of both gives d_i = 0 and is kept
        d_per_file = (count_x - count_minus_x) / n_events
        #Central value: the average of the 10 sub-samples
        difference = d_per_file.mean()
        sigma = np.std(d_per_file, ddof=1)
        #If all 10 files give the same d_i, sigma = 0 and n_sigma is undefined
        n_sigma = abs(difference) / sigma if sigma > 0 else np.nan

        rows.append({"pair_code": code, "pair_name": name, "momentum_bin": label,
                     "bin_low": bin_low[k], "bin_high": bin_high[k],
                     "difference": difference, "uncertainty": sigma, "n_sigma": n_sigma,
                     "significant": bool(n_sigma > threshold)})

results = pd.DataFrame(rows)

# ---- Save (rounded only here, all calculations use the exact values) ----
saved = results.copy()
saved["difference"] = signif(saved["difference"], 4)   #4 significant digits, so rare strains do not round to 0
saved["uncertainty"] = signif(saved["uncertainty"], 3)
saved["n_sigma"] = saved["n_sigma"].round(2)
saved.to_csv(output_csv, index=False)

# ---- Print a readable table ----
table = saved[["pair_name", "momentum_bin", "difference", "uncertainty", "n_sigma", "significant"]]
with pd.option_context("display.width", 120, "display.max_rows", 200):
    print("Mean asymmetry A = average of the 10 sub-samples of avg(X) - avg(-X), per momentum bin [bacteria per event],")
    print(f"uncertainty = spread of the {len(file_names)} sub-samples (1 sigma):\n")
    print(table.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

summary = [f"{name} {int(results.loc[results['pair_name'] == name, 'significant'].sum())} of {len(bin_labels)}"
           for name in pair_names]
print(f"\nBins above {threshold} sigma: " + ", ".join(summary))

# ---- Check: the 5 bins add up to the total difference of week 4 ----
if os.path.isfile(week4_pairs_csv):
    week4 = pd.read_csv(week4_pairs_csv).set_index("code")
    print(f"\nCheck against {week4_pairs_csv} (sum of the 5 bins vs total difference):")
    for code, name in zip(pair_codes, pair_names):
        bin_sum = results.loc[results["pair_code"] == code, "difference"].sum()
        reference = week4.loc[code, "difference"]
        #The week 4 value is rounded to 4 significant digits
        agrees = abs(bin_sum - reference) <= 1e-3 * abs(reference)
        print(f"  {name:16s} {bin_sum:.6g} vs {reference:.6g}  {'OK' if agrees else 'DIFFERENT'}")
        if not agrees:
            sys.exit(f"The momentum bins of {name} do not add up to the week 4 total.")

# ---- Does the asymmetry depend on momentum? ----
#The bins contain different numbers of bacteria, so A per bin follows the bin
#population even if the mutation had the same effect at every momentum.
#Null hypothesis: the imbalance is the same for every bacterium, whatever its
#momentum. Then A in bin k is the total A times that bin's share of the
#bacteria: expected_k = A_total * N_k / N_total, with N_k = count_X + count_minusX.
#chi2 = sum over the 5 bins of ((A_k - expected_k) / sigma_k)^2,
#with 4 degrees of freedom (5 bins, A_total taken from the data).


def chi2_p_value(chi2, ndf):
    """P(chi2 >= observed) for an even number of degrees of freedom (no scipy needed)."""
    if ndf % 2 != 0:
        return np.nan
    half = chi2 / 2
    return float(np.exp(-half) * sum(half ** i / np.prod(range(1, i + 1)) for i in range(ndf // 2)))


print("\nMomentum dependence: A per bin compared with 'same imbalance at every momentum'")
for code, name in zip(pair_codes, pair_names):
    pair = results[results["pair_code"] == code]
    pair_counts = counts[counts["pair_code"] == code]
    n_per_bin = np.array([(pair_counts.loc[pair_counts["momentum_bin"] == b, "count_X"] +
                           pair_counts.loc[pair_counts["momentum_bin"] == b, "count_minusX"]).sum()
                          for b in pair["momentum_bin"]])
    expected = pair["difference"].sum() * n_per_bin / n_per_bin.sum()
    chi2 = float((((pair["difference"] - expected) / pair["uncertainty"]) ** 2).sum())
    ndf = len(pair) - 1
    print(f"  {name:16s} chi2 / ndf = {chi2:.1f} / {ndf}, p = {chi2_p_value(chi2, ndf):.3g}")


# ---- Scatter plot of the asymmetry against momentum for one pair ----
#The bins are drawn at equal distance (position 1 to 5); the tick labels give
#the momentum range of each bin
x_centre = np.arange(1, len(bin_labels) + 1)
x_half_width = 0.4   #Horizontal bar: marks the bin (its range is on the axis)
tick_labels = [f"{lo:g} to {hi:g}" if np.isfinite(hi) else f"> {lo:g}" for lo, hi in zip(bin_low, bin_high)]


def draw_pair(ax, pair_rows, title):
    y = pair_rows["difference"].to_numpy()
    sigma = pair_rows["uncertainty"].to_numpy()
    n_sigma = pair_rows["n_sigma"].to_numpy()
    colours = [colour_significant if s else colour_not_significant for s in pair_rows["significant"]]

    ax.axhline(0, color=colour_error_bar, linestyle="--", linewidth=1, zorder=1)
    #Vertical bars: 1 sigma; horizontal bars: mark the momentum bin
    ax.errorbar(x_centre, y, yerr=sigma, xerr=x_half_width, fmt="none",
                ecolor=colour_error_bar, elinewidth=1.2, capsize=3, zorder=2)
    ax.scatter(x_centre, y, s=60, c=colours, edgecolors="black", linewidths=0.8, zorder=3)

    #n_sigma next to every point: above the error bar if A >= 0, below it if A < 0
    for xi, yi, si, ni in zip(x_centre, y, sigma, n_sigma):
        #Short label ("2.3σ"), so the labels of neighbouring bins do not overlap
        #Near the threshold 1 decimal is not enough (2.996 would read "3.0" while grey)
        digits = 2 if round(ni, 1) == threshold else 1
        label = "n/a" if np.isnan(ni) else f"{ni:.{digits}f}\u03c3"
        above = yi >= 0
        ax.annotate(label, (xi, yi + si if above else yi - si),
                    xytext=(0, 4 if above else -4), textcoords="offset points",
                    ha="center", va="bottom" if above else "top", fontsize=8)

    #y range: all error bars and the zero line, plus padding, plus room for the labels
    low = min(0, np.min(y - sigma))
    high = max(0, np.max(y + sigma))
    span = high - low
    has_label_below = np.any(y < 0)
    ax.set_ylim(low - span * (0.12 + (0.12 if has_label_below else 0)),
                high + span * (0.12 + (0.12 if np.any(y >= 0) else 0)))

    ax.set_xlim(0.5, len(bin_labels) + 0.5)
    ax.set_xticks(x_centre)
    ax.set_xticklabels(tick_labels, fontsize=8)
    ax.set_title(f"{title} (WT vs mutant)")
    ax.set_xlabel(f"Momentum |p| bin [{momentum_unit}]")
    ax.set_ylabel("Mean asymmetry A = WT - mutant\n[average bacteria per event]")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)   #Factor shown as x10^-2, not 1e-2
    ax.grid(True, color="#e8e7e2", zorder=0)
    ax.set_axisbelow(True)


# ---- All 6 pairs in one 2 x 3 figure ----
fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
for ax, code, name in zip(axes.flat, pair_codes, pair_names):
    draw_pair(ax, results[results["pair_code"] == code], name)
fig.suptitle("Asymmetry (WT minus mutant) as a function of momentum, "
             f"{len(file_names)} sub-samples (error bars: 1 sigma)", fontsize=13)
key = [Line2D([], [], marker="o", linestyle="none", markersize=8, markerfacecolor=colour_significant,
              markeredgecolor="black", label=f"> {threshold} sigma: significant"),
       Line2D([], [], marker="o", linestyle="none", markersize=8, markerfacecolor=colour_not_significant,
              markeredgecolor="black", label=f"<= {threshold} sigma: not significant"),
       Line2D([], [], color=colour_error_bar, marker="|", markersize=10,
              label=f"error bars: 1 sigma from the spread of the {len(file_names)} sub-samples; "
                    "horizontal: the bin")]
fig.legend(handles=key, loc="lower center", ncol=3, frameon=False, fontsize=10)
fig.tight_layout(rect=(0, 0.05, 1, 0.96))
fig.savefig("asymmetry_momentum_all_FM.png", dpi=200)
plt.close(fig)
print("Saved figure to asymmetry_momentum_all_FM.png")
