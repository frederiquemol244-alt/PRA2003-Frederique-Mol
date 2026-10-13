# PRA2003 - Week 7: Salmonella asymmetry vs momentum in equal-count bins
# The 5 group bins hold equal numbers of ALL bacteria together, which are
# mostly E. coli. Salmonella has higher momenta, so the highest group bin
# contains 14 times more Salmonella than the lowest, and the error bars in
# Figure 2 grow from bin to bin for that reason alone.
#
# This script makes the Salmonella plot with 10 bins that each hold the same
# number of Salmonella (the deciles of the Salmonella momenta). The
# statistical uncertainty is then expected to be the same in every bin, so a
# real momentum dependence would stand out against a flat band.
#
# HOW TO RUN:
#   python Week7_salmonella_momentum_FM.py
#
# NEEDS:    Python 3.9 or newer with pandas, numpy and matplotlib
#           sub_sample_results_FM.csv (week 4: events and counts per file)
#           salmonella_momenta_FM.csv, or (to remake it) the raw data files
#           output-Set1.txt ... output-Set10.txt
# PRODUCES: salmonella_momenta_FM.csv (momentum of every Salmonella, made
#           once from the raw data, about 1.5 minutes)
#           salmonella_momentum_FM.csv (asymmetry and uncertainty per bin)
#           salmonella_momentum_FM.png (Figure 5)
#
# METHOD (the same sub-sampling method as in weeks 4 and 5):
#   per bin and file:  A_i = (N_WT - N_mutant) / non-empty events of the file
#   central value:     A = the average of the 10 values A_i (an average per
#                      event, not a total number of bacteria)
#   uncertainty:       standard deviation of the 10 values A_i (no division
#                      by sqrt(10))
#   For comparison the plot also shows the uncertainty expected from
#   counting alone: sqrt(N_WT + N_mutant) per file / events. Salmonella is
#   never found twice in one event, so the two should agree (README
#   section 6).
#
# Momentum: |p| = sqrt(px^2 + py^2 + pz^2). output-Set10.txt stores
# (pT, phi, eta) instead, so there |p| = pT * cosh(eta), as in week 5.

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Settings ----
file_names = [f"output-Set{i}.txt" for i in range(1, 11)]
wt_id, mutant_id = 3334, -3334
n_bins = 10
chunk_size = 2_000_000
pt_phi_eta_files = {"output-Set10.txt"}
week4_csv = "sub_sample_results_FM.csv"
momenta_csv = "salmonella_momenta_FM.csv"
output_csv = "salmonella_momentum_FM.csv"
output_png = "salmonella_momentum_FM.png"

script_folder = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_folder)
data_folders = []
folder = script_folder
for _ in range(4):
    data_folders += [folder, os.path.join(folder, "data")]
    folder = os.path.dirname(folder)


def find_data_file(filename):
    """Return the full path of a data file, or None if it is in no data folder."""
    for folder in data_folders:
        path = os.path.join(folder, filename)
        if os.path.isfile(path):
            return path
    return None


def read_salmonella(path, filename):
    """Return a table (sub_sample, id, p) with every Salmonella in one file."""
    parts = []
    reader = pd.read_csv(path, sep=r"\s+", header=None, names=["a", "b", "c", "d"],
                         dtype=np.float64, chunksize=chunk_size, engine="c")
    for chunk in reader:
        rows = chunk[chunk["d"].isin([wt_id, mutant_id])]
        if len(rows) == 0:
            continue
        if filename in pt_phi_eta_files:
            p = rows["a"] * np.cosh(rows["c"])                       #a = pT, c = eta
        else:
            p = np.sqrt(rows["a"] ** 2 + rows["b"] ** 2 + rows["c"] ** 2)
        parts.append(pd.DataFrame({"sub_sample": filename,
                                   "id": rows["d"].astype(int).to_numpy(),
                                   "p": p.to_numpy()}))
    return pd.concat(parts, ignore_index=True)


# ---- Step 1: the momentum of every Salmonella (from the CSV if it exists) ----
if not os.path.isfile(week4_csv):
    sys.exit(f"Can't find '{week4_csv}' - keep it in the same folder as this script.")
week4 = pd.read_csv(week4_csv)

if os.path.isfile(momenta_csv):
    momenta = pd.read_csv(momenta_csv)
    print(f"Read {len(momenta)} Salmonella from {momenta_csv}")
else:
    missing = [f for f in file_names if find_data_file(f) is None]
    if missing:
        sys.exit(f"Can't find '{momenta_csv}' or these data files: {', '.join(missing)}")
    parts = []
    for filename in file_names:
        print(f"Reading {filename} ...")
        parts.append(read_salmonella(find_data_file(filename), filename))
    momenta = pd.concat(parts, ignore_index=True)
    momenta.to_csv(momenta_csv, index=False, float_format="%.6g")
    print(f"Saved {len(momenta)} rows to {momenta_csv}")

# ---- Check: the number of Salmonella per file equals week 4 ----
events = {}
for filename in file_names:
    in_file = week4[week4["sub_sample"] == filename]
    events[filename] = int(in_file["n_events"].iloc[0])
    for code in (wt_id, mutant_id):
        expected = int(in_file.loc[in_file["code"] == code, "total_count"].iloc[0])
        found = int(((momenta["sub_sample"] == filename) & (momenta["id"] == code)).sum())
        if found != expected:
            sys.exit(f"{filename}, ID {code}: found {found} bacteria, week 4 has {expected}.")
print("Check passed: the counts of all 10 files equal week 4.")
total_events = sum(events.values())

# ---- Step 2: 10 bins with the same number of Salmonella in each ----
#Deciles of the momentum of all Salmonella (WT and mutant together, all 10
#files), rounded to 3 significant digits. WT and mutant are used together, so
#the choice of edges cannot create a difference between them.
quantiles = np.quantile(momenta["p"], np.arange(1, n_bins) / n_bins)
inner_edges = [float(f"{q:.3g}") for q in quantiles]
edges = [0.0] + inner_edges + [np.inf]
labels = [f"{edges[i]:g} to {edges[i + 1]:g}" for i in range(n_bins - 1)] + [f"> {edges[-2]:g}"]
momenta["bin"] = np.searchsorted(edges, momenta["p"], side="right") - 1

# ---- Step 3: asymmetry and uncertainty per bin ----
rows = []
for b in range(n_bins):
    in_bin = momenta[momenta["bin"] == b]
    a_per_file = []
    n_per_file = []
    for filename in file_names:
        f = in_bin[in_bin["sub_sample"] == filename]
        n_wt = int((f["id"] == wt_id).sum())
        n_mutant = int((f["id"] == mutant_id).sum())
        a_per_file.append((n_wt - n_mutant) / events[filename])
        n_per_file.append(n_wt + n_mutant)
    n_wt_total = int((in_bin["id"] == wt_id).sum())
    n_mutant_total = int((in_bin["id"] == mutant_id).sum())
    asymmetry = float(np.mean(a_per_file))                          #Average of the 10 sub-samples
    sigma = float(np.std(a_per_file, ddof=1))                       #Sub-sampling: spread of the 10 files
    sigma_counting = np.sqrt(np.mean(n_per_file)) / (total_events / len(file_names))
    rows.append({"bin": b + 1, "momentum_bin": labels[b],
                 "N_WT": n_wt_total, "N_mutant": n_mutant_total,
                 "asymmetry": asymmetry, "uncertainty": sigma,
                 "n_sigma": abs(asymmetry) / sigma,
                 "counting_uncertainty": sigma_counting})
result = pd.DataFrame(rows)

# ---- Step 4: is the asymmetry the same in every bin? ----
weights = 1 / result["uncertainty"] ** 2
constant = float((weights * result["asymmetry"]).sum() / weights.sum())
chi2 = float((weights * (result["asymmetry"] - constant) ** 2).sum())
ndf = n_bins - 1

result.to_csv(output_csv, index=False, float_format="%.6g")
scale = 1e5
print("\nSalmonella, asymmetry A = WT - mutant per momentum bin [1e-5 bacteria per event]")
for _, r in result.iterrows():
    print(f"  {r['momentum_bin']:>14s}  N = {r['N_WT'] + r['N_mutant']:5d}  "
          f"A = {r['asymmetry'] * scale:+5.2f} +- {r['uncertainty'] * scale:4.2f} "
          f"({r['n_sigma']:.1f} sigma)   counting alone: {r['counting_uncertainty'] * scale:4.2f}")
print(f"Constant A: chi2 / ndf = {chi2:.1f} / {ndf}")
print(f"Saved {output_csv}")

# ---- Step 5: the plot ----
x = np.arange(n_bins)
fig, ax = plt.subplots(figsize=(11, 5.6))
expected = float(result["counting_uncertainty"].mean()) * scale
ax.axhspan(-expected, expected, color="#D9E6F5", zorder=0,
           label="expected from counting alone: 0 ± 1 sigma")
ax.axhline(0, color="#444444", linestyle="--", linewidth=1, label="A = 0 (no asymmetry)")
significant = result["n_sigma"] > 3
for mask, color, label in [(~significant, "#888780", "<= 3 sigma: not significant"),
                           (significant, "#2B74D8", "> 3 sigma: significant")]:
    if mask.any():
        ax.errorbar(x[mask], result.loc[mask, "asymmetry"] * scale,
                    yerr=result.loc[mask, "uncertainty"] * scale,
                    fmt="o", color=color, markeredgecolor="black", markersize=8,
                    ecolor="#444444", elinewidth=1.3, capsize=4, zorder=3,
                    label=f"{label}; error bars: 1 sigma from the spread of the 10 sub-samples")
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=9)
ax.set_xlabel("Momentum |p| bin [10$^{-20}$ kg m/s], each bin holds 1/10 of the Salmonella")
ax.set_ylabel("Mean asymmetry A = WT - mutant\n[10$^{-5}$ average bacteria per event]")
ax.set_title("Salmonella (WT vs mutant): asymmetry as a function of momentum, 10 sub-samples")
ax.grid(True, color="#E6E4DD", linewidth=0.8)
ax.set_axisbelow(True)
limit = 1.25 * float((result["asymmetry"].abs() + result["uncertainty"]).max()) * scale
ax.set_ylim(-limit, limit)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=1, frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(output_png, dpi=200)
print(f"Saved {output_png}")
