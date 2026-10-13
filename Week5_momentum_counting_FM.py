# PRA2003 - Week 5 Deliverable (Part 1 of 3): counting per momentum bin
# Reads the 10 data files output-Set1.txt through output-Set10.txt
# (output-Set0.txt is a test file with 1 event and is excluded) and counts,
# for every file, every wild-type/mutant pair and every momentum bin, how
# many bacteria with ID X and with ID -X there are.
#
# MOMENTUM BINS: the 5 bins used by the group (unit 10^-20 kg m/s):
#   0-0.654, 0.654-1.32, 1.32-2.68, 2.68-5.88 and > 5.88
# (lower edge included, upper edge excluded; the last bin has no upper edge).
# The edges are the quintiles (1/5, 2/5, 3/5, 4/5) of the momentum of the
# bacteria of the 6 pairs in the first 1,000,000 lines of output-Set1.txt,
# rounded to 3 significant digits, so each bin holds roughly the same number
# of bacteria. They are fixed here, so every file uses the same bins and the
# result does not depend on which file is read first. Every bacterium falls
# in exactly one bin. The result is saved in
# one small CSV, so the plots (Part 2) can be remade in seconds without the
# raw data (~800 MB per file).
#
# HOW TO RUN:
#   python Week5_momentum_counting_FM.py
#   The data files are looked for in the folders listed in 'data_folders'
#   below (first the folder of this script). sub_sample_results_FM.csv and
#   asymmetry_vs_momentum_Set1_FM.csv (from weeks 3 and 4) must be in the
#   same folder as this script: they are used for the built-in checks.
#
# NEEDS:    Python 3.9 or newer with pandas and numpy
#           (pip install pandas numpy matplotlib)
# PRODUCES: momentum_sub_sample_results_FM.csv (300 rows: 10 files x 6 pairs
#           x 5 momentum bins)
# RUN TIME: about 2 minutes on a laptop (10 to 15 seconds per file).
#
# WHAT IT DOES:
#   Each event in a file looks like this:
#     - a header line: event_id  number_of_particles     (2 fields)
#     - followed by that many rows of: px py pz ID      (4 fields)
#
#   Every line can be classified on its own (2 fields = header, 4 fields =
#   bacterium), so no information has to be carried from one line to the
#   next and the file can safely be read in chunks. Only running totals
#   (12 IDs x 9 momentum cells) are kept, never one entry per bacterium.
#
#   Momentum: |p| = sqrt(px^2 + py^2 + pz^2).
#   Exception: output-Set10.txt stores (pT, phi, eta) (transverse momentum,
#   angle, pseudorapidity) in the first three columns instead of
#   (px, py, pz). For that file |p| = pT * cosh(eta). The script checks the
#   column ranges of every file and stops if a file does not have the format
#   listed in 'pt_phi_eta_files'.
#
#   Only events with at least 1 bacterium are counted as events (an empty
#   event contains no measurement), the same rule as in week 4.
#
# BUILT-IN CHECKS (the script stops with an error if one fails):
#   1. n_events of every file equals n_events in sub_sample_results_FM.csv
#   2. for every file and ID, the sum over the 5 momentum bins equals
#      total_count in sub_sample_results_FM.csv
#   3. the counts of output-Set1.txt reproduce the week 3 result
#      asymmetry_vs_momentum_Set1_FM.csv (columns N_normal and N_mutant).
#      Week 3 used other bins (edges 0, 1, 2, 4, 8). To make this check
#      possible the bacteria are counted on a finer grid that contains the
#      edges of both sets of bins; the cells are then added up both ways.
#
# If the script is stopped halfway, the files already done are kept in
# momentum_sub_sample_results_FM_partial.csv, and the next run skips them.

import os
import sys
import time

import numpy as np
import pandas as pd

# ---- Settings ----
file_names = [f"output-Set{i}.txt" for i in range(1, 11)]
chunk_size = 2_000_000                          #Number of lines read at a time
output_csv = "momentum_sub_sample_results_FM.csv"
partial_csv = "momentum_sub_sample_results_FM_partial.csv"   #Progress, so a crash does not lose everything
#Group bins [lower, upper), unit 10^-20 kg m/s; the last bin is open ended
edges = [0, 0.654, 1.32, 2.68, 5.88, np.inf]
bin_labels = ["0-0.654", "0.654-1.32", "1.32-2.68", "2.68-5.88", ">5.88"]

#Week 3 bins (asymmetry_vs_momentum_Set1_FM.csv), only used for check 3
week3_edges = [0, 1, 2, 4, 8, np.inf]
week3_labels = ["[0,1)", "[1,2)", "[2,4)", "[4,8)", "[8,Inf)"]

#Files whose first three columns are (pT, phi, eta) instead of (px, py, pz)
pt_phi_eta_files = {"output-Set10.txt"}

#Reference results from weeks 3 and 4, used for the checks
week4_csv = "sub_sample_results_FM.csv"
week3_set1_csv = "asymmetry_vs_momentum_Set1_FM.csv"

#Folders searched for the data files, in this order: this script's folder,
#then up to 3 folders above it (each also with a data/ subfolder), so the
#~8 GB of raw data does not have to be copied next to the scripts
script_folder = os.path.dirname(os.path.abspath(__file__))
data_folders = []
folder = script_folder
for _ in range(4):
    data_folders += [folder, os.path.join(folder, "data")]
    folder = os.path.dirname(folder)

# ---- The 6 wild-type/mutant pairs: ID X (wild type) and -X (mutant) ----
pair_codes = [211, 321, 2212, 3122, 3312, 3334]
pair_names = ["E. coli", "B. subtilis", "P. aeruginosa",
              "S. pneumoniae", "M. tuberculosis", "Salmonella"]

#Protection: codes, names and bins must line up 1-to-1
if len(pair_codes) != len(pair_names):
    sys.exit("pair_codes and pair_names must be the same length.")
if len(edges) - 1 != len(bin_labels) or any(np.diff(edges) <= 0):
    sys.exit("There must be exactly one bin label per momentum bin, and the edges must increase.")

#All 12 IDs, sorted, so np.searchsorted can turn an ID into an index 0..11
all_ids = np.array(sorted(pair_codes + [-c for c in pair_codes]), dtype=np.int64)
n_bins = len(edges) - 1

#Fine grid with the edges of both sets of bins. Every fine cell lies inside
#exactly one group bin and exactly one week 3 bin, so counts on the fine grid
#can be added up to either set.
fine_edges = sorted(set(edges) | set(week3_edges))
n_cells = len(fine_edges) - 1
cell_to_bin = np.searchsorted(edges, fine_edges[:-1], side="right") - 1
cell_to_week3 = np.searchsorted(week3_edges, fine_edges[:-1], side="right") - 1


def add_up(cell_counts, cell_to_target, n_target):
    """Add the counts of the fine cells (last axis) up to wider bins."""
    out = np.zeros(cell_counts.shape[:-1] + (n_target,), dtype=np.int64)
    for cell, target in enumerate(cell_to_target):
        out[..., target] += cell_counts[..., cell]
    return out


def find_data_file(filename):
    """Return the full path of a data file, or None if it is in no data folder."""
    for folder in data_folders:
        path = os.path.join(folder, filename)
        if os.path.isfile(path):
            return os.path.normpath(path)
    return None


def count_file(path, filename):
    """Count the 12 IDs per fine momentum cell in one file. Returns (counts, n_events)."""
    counts = np.zeros(len(all_ids) * n_cells, dtype=np.int64)   #Running totals, index = id_index * n_cells + bin
    n_events = 0
    malformed_lines = 0
    #Column ranges, used to confirm which momentum format the file has
    min_first_column = np.inf
    max_abs_second_column = 0.0
    uses_pt_phi_eta = filename in pt_phi_eta_files

    reader = pd.read_csv(path, sep=r"\s+", header=None, names=["a", "b", "c", "d"],
                         dtype=np.float64, chunksize=chunk_size, engine="c")
    for chunk in reader:
        #A header line has 2 fields, so c and d are empty (NaN)
        is_header = chunk["a"].notna() & chunk["b"].notna() & chunk["c"].isna() & chunk["d"].isna()
        #A bacterium line has 4 numbers
        is_particle = chunk[["a", "b", "c", "d"]].notna().all(axis=1)
        #Anything else (3 fields, a value that is not a number, ...) is malformed: skip it, but count it
        malformed_lines += int((~is_header & ~is_particle).sum())

        #Only events with at least 1 bacterium count (column b is nParticles)
        n_events += int((chunk.loc[is_header, "b"] > 0).sum())

        particles = chunk.loc[is_particle]
        if len(particles) == 0:
            continue
        min_first_column = min(min_first_column, particles["a"].min())
        max_abs_second_column = max(max_abs_second_column, particles["b"].abs().max())

        #Keep only the 12 IDs of the 6 pairs
        ids = particles["d"].to_numpy().astype(np.int64)
        id_index = np.searchsorted(all_ids, ids)
        id_index[id_index == len(all_ids)] = 0          #Out of range: will fail the equality test below
        keep = all_ids[id_index] == ids
        a = particles["a"].to_numpy()[keep]
        b = particles["b"].to_numpy()[keep]
        c = particles["c"].to_numpy()[keep]

        if uses_pt_phi_eta:
            p = a * np.cosh(c)                           #a = pT, c = eta
        else:
            p = np.sqrt(a ** 2 + b ** 2 + c ** 2)        #a, b, c = px, py, pz

        #side="right" puts p == lower edge in that cell, so every cell is [lower, upper)
        bin_index = np.searchsorted(fine_edges, p, side="right") - 1
        if np.any(bin_index < 0) or np.any(bin_index >= n_cells):
            sys.exit(f"{filename}: a momentum fell outside all bins - check the bin edges.")

        counts += np.bincount(id_index[keep] * n_cells + bin_index, minlength=len(counts))

    if malformed_lines > 0:
        print(f"  Note: {malformed_lines} line(s) were skipped (neither a valid header nor a valid particle row).")
    else:
        print("  Malformed lines skipped: 0")

    #Protection against dividing by zero later if the file had no valid events
    if n_events <= 0:
        sys.exit(f"{filename}: no events found - check the file path/format.")

    #Check that the file really has the momentum format the settings assume.
    #(pT, phi, eta): pT is never negative and phi is an angle in [-pi, pi].
    looks_like_pt_phi_eta = min_first_column >= 0 and max_abs_second_column <= np.pi + 1e-6
    if looks_like_pt_phi_eta != uses_pt_phi_eta:
        found = "(pT, phi, eta)" if looks_like_pt_phi_eta else "(px, py, pz)"
        sys.exit(f"{filename}: the columns look like {found}, which does not match "
                 f"'pt_phi_eta_files' in the settings. Update the settings.")

    return counts.reshape(len(all_ids), n_cells), n_events


def to_rows(filename, counts, n_events):
    """Turn the counts of one file into 30 rows (6 pairs x 5 group bins)."""
    counts = add_up(counts, cell_to_bin, n_bins)
    rows = []
    for code, name in zip(pair_codes, pair_names):
        x = counts[np.searchsorted(all_ids, code)]
        minus_x = counts[np.searchsorted(all_ids, -code)]
        for k in range(n_bins):
            rows.append({"sub_sample": filename, "pair_code": code, "pair_name": name,
                         "momentum_bin": bin_labels[k], "n_events": n_events,
                         "count_X": int(x[k]), "count_minusX": int(minus_x[k])})
    return rows


def check_file(file_rows, filename, week4):
    """Built-in checks 1 and 2 for one file. Stops the script if one fails."""
    # Check 1: same number of non-empty events as in week 4
    reference_events = week4.loc[week4["sub_sample"] == filename, "n_events"].unique()
    if len(reference_events) != 1 or file_rows["n_events"].iloc[0] != reference_events[0]:
        sys.exit(f"CHECK 1 FAILED for {filename}: n_events = {file_rows['n_events'].iloc[0]}, "
                 f"{week4_csv} has {reference_events}.")

    # Check 2: the 5 momentum bins add up to the week 4 total of every ID
    for code in pair_codes:
        rows = file_rows[file_rows["pair_code"] == code]
        for signed_code, column in [(code, "count_X"), (-code, "count_minusX")]:
            reference = week4.loc[(week4["sub_sample"] == filename) & (week4["code"] == signed_code),
                                  "total_count"]
            expected = int(reference.iloc[0]) if len(reference) else 0   #No row = never seen in that file
            found = int(rows[column].sum())
            if found != expected:
                sys.exit(f"CHECK 2 FAILED for {filename}, ID {signed_code}: momentum bins add up to "
                         f"{found}, {week4_csv} has {expected} (difference {found - expected}).")


def check_week3(counts, week3_set1):
    """Built-in check 3: the counts of output-Set1.txt, added up to the week 3
    bins, must reproduce the week 3 result. Stops the script if it fails."""
    counts = add_up(counts, cell_to_week3, len(week3_labels))
    for _, ref in week3_set1.iterrows():
        code = pair_codes[pair_names.index(ref["Pair"])]
        k = week3_labels.index(ref["Momentum"])
        found = (int(counts[np.searchsorted(all_ids, code), k]),
                 int(counts[np.searchsorted(all_ids, -code), k]))
        expected = (int(ref["N_normal"]), int(ref["N_mutant"]))
        if found != expected:
            sys.exit(f"CHECK 3 FAILED for {ref['Pair']} {ref['Momentum']}: found {found}, "
                     f"{week3_set1_csv} has {expected}.")


# ---- Fail early: every data file and reference CSV must be present ----
os.chdir(script_folder)   #The CSV files are read and written next to this script

missing_csvs = [f for f in [week4_csv, week3_set1_csv] if not os.path.isfile(f)]
if missing_csvs:
    sys.exit("Can't find " + ", ".join(missing_csvs) + ". These week 3/4 results are needed for "
             "the built-in checks; keep them in the same folder as this script.")

data_paths = {f: find_data_file(f) for f in file_names}
missing_files = [f for f, path in data_paths.items() if path is None]
if missing_files:
    sys.exit("Can't find these data files:\n  " + "\n  ".join(missing_files) +
             "\nDownload them from the course SURFdrive and put them in one of these folders:\n  " +
             "\n  ".join(os.path.normpath(d) for d in data_folders))

if not isinstance(chunk_size, int) or chunk_size <= 0:
    sys.exit("chunk_size must be a positive whole number.")

week4 = pd.read_csv(week4_csv)
week3_set1 = pd.read_csv(week3_set1_csv)
#In the week 3 file Salmonella was called by its species name
week3_set1["Pair"] = week3_set1["Pair"].replace({"S. enterica": "Salmonella"})

# ---- Resume: files already in the partial CSV are not counted again ----
all_rows = []
done_files = set()
if os.path.isfile(partial_csv):
    previous = pd.read_csv(partial_csv)
    done_files = set(previous["sub_sample"])
    all_rows = previous.to_dict("records")
    print(f"Resuming: {len(done_files)} file(s) already counted in {partial_csv}.")

# ---- Count every file ----
start_time = time.time()
for filename in file_names:
    if filename in done_files:
        print(f"=== {filename}: already done, skipped ===")
        continue
    print(f"=== Counting {filename} ({data_paths[filename]}) ===")
    file_start = time.time()
    counts, n_events = count_file(data_paths[filename], filename)
    file_rows = to_rows(filename, counts, n_events)
    check_file(pd.DataFrame(file_rows), filename, week4)
    if filename == "output-Set1.txt":
        check_week3(counts, week3_set1)
    all_rows.extend(file_rows)
    pd.DataFrame(all_rows).to_csv(partial_csv, index=False)   #Save progress after every file
    print(f"  {n_events} events, checks passed, {time.time() - file_start:.0f} s "
          f"(total {time.time() - start_time:.0f} s)")

# ---- Final checks on the complete table, then save ----
results = pd.DataFrame(all_rows)
results["file_order"] = results["sub_sample"].map({f: i for i, f in enumerate(file_names)})
results["pair_order"] = results["pair_code"].map({c: i for i, c in enumerate(pair_codes)})
results["bin_order"] = results["momentum_bin"].map({b: i for i, b in enumerate(bin_labels)})
results = results.sort_values(["file_order", "pair_order", "bin_order"]).drop(
    columns=["file_order", "pair_order", "bin_order"])

#Also re-checks files taken over from a previous (partial) run. Check 3 needs
#the fine cells, so it is only done when output-Set1.txt is counted (above).
for filename in file_names:
    check_file(results[results["sub_sample"] == filename], filename, week4)

expected_rows = len(file_names) * len(pair_codes) * n_bins
if len(results) != expected_rows:
    sys.exit(f"Expected {expected_rows} rows, found {len(results)}.")

results.to_csv(output_csv, index=False)
if os.path.isfile(partial_csv):
    os.remove(partial_csv)   #Complete result saved, progress file no longer needed

n_events_total = results.drop_duplicates("sub_sample")["n_events"].sum()
print(f"\nAll 3 checks passed for all {len(file_names)} files ({n_events_total} non-empty events).")
print(f"Saved {len(results)} rows to {output_csv}")
print(f"Total time elapsed: {(time.time() - start_time) / 60:.1f} min")
