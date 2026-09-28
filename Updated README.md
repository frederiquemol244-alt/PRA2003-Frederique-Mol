# PRA2003 – Week 4: Wild-type vs mutant bacteria in 5 million events


## Contents

- [1. Context and scientific question](#1-context-and-scientific-question)
- [2. Files in this submission](#2-files-in-this-submission)
- [3. Installation](#3-installation)
- [4. Usage: how to run the code](#4-usage-how-to-run-the-code)
- [5. The data](#5-the-data)
- [6. Method](#6-method)
- [7. Results](#7-results)
- [8. Conclusion](#8-conclusion)
- [9. Troubleshooting](#9-troubleshooting)

---

## 1. Context and scientific question

Each bacterial strain in this data comes as a normal (wild-type) form with a
positive ID and a mutant form with the same ID and a minus sign, for example
an antibiotic-resistant strain. A mutation can change how common a strain
is, so a difference in counts between the wild type and its mutant would
show that the mutation has an effect.

### Scientific questions

These are the goals set at the start of the project, in the week 1 README
(`README_Frederique.md`):

> Goal to answer the following questions:
> - What are the average counts of each bacterial stain and
> their statistical uncertainties?
> - Is there any asymmetry between the normal and the
> mutant strain?
>  Quantify!!!
> - Is there any asymmetry as a function of their momentum?
>  Quantify!!!

| Week 1 question | Where it is answered |
|---|---|
| 1. Average counts and their statistical uncertainties | **This week**, using the full 5M-event sample: [section 7.1](#71-average-count-per-event-5m-events-10-sub-samples) |
| 2. Asymmetry between the normal and the mutant strain | **This week**, quantified for every pair: sections [7.2](#72-asymmetry-wild-type-vs-mutant) and [7.3](#73-asymmetry-the-other-12-particleantiparticle-pairs) |
| 3. Asymmetry as a function of momentum | The week 3 deliverable, using one file (`output-Set1.txt`) |

For this week the questions are made precise as follows:

1. **What is the average number of each bacterial strain per event, and
   what is its statistical uncertainty?**
2. **Is there an asymmetry between each wild type and its mutant?** More
   generally, is there an asymmetry for every pair of IDs X and -X? How
   large is it, and is it statistically significant?

**Hypothesis:** if a mutation has no effect on how common a strain is, then
the wild type and the mutant have the same average per event (asymmetry
A = 0). An asymmetry is only claimed when the data disagree with A = 0 by
**more than 3σ**.

---

## 2. Files in this submission

| File | Type | Purpose |
|---|---|---|
| `README.md` | documentation | This file |
| `README_Frederique.md` | documentation | The original week 1 README with the project goals |
| `Week 4.r` | script, part 1 | Counts every ID in each of the 10 data files |
| `week 4.2.r` | script, part 2 | Combines the 10 files with the sub-sampling method |
| `Week 4 part 3.r` | script, part 3 | Asymmetry test for every pair X vs -X, and the figure |
| `Week 4 limitations.r` | script, optional | Evidence for the limitations (paired test, standard error of the mean) |
| `sub_sample_results.csv` | data, output of part 1 | Count and average per ID, for each file |
| `final_results.csv` | data, output of part 2 | Final average ± uncertainty per ID |
| `three_sigma_results.csv` | data, output of part 3 | Difference, n σ and asymmetry for all 18 pairs |
| `asymmetry_plot.png` | figure, output of part 3 | Asymmetry of the 6 wild-type/mutant pairs |

The raw data files (about 8 GB) are **not** included because they are too
large. See [The data](#5-the-data) for how to get them. Parts 2 and 3 can
be run **without** the raw data, because they only need the included CSV
files.

---

## 3. Installation

### Requirements

| What | Version / size | Notes |
|---|---|---|
| **R** | 4.0 or newer (tested with **R 4.6.1**) | Free download from [cran.r-project.org](https://cran.r-project.org) |
| **R packages** | none | Only base R is used, so nothing needs to be installed with `install.packages()` |
| **Disk space** | about 8 GB | Only needed for the raw data files for part 1 |
| **Memory (RAM)** | no special requirement | Part 1 reads the files in chunks, so it never loads a full file into memory |
| *Optional:* an editor | RStudio, or VS Code with the R extension | You can also run everything from a terminal |

### Step by step

1. **Install R** from [cran.r-project.org](https://cran.r-project.org): choose
   your operating system and follow the installer.
2. **Check that R works.** Open a terminal (Mac: *Terminal*, Windows:
   *Command Prompt*) and type:
   ```
   Rscript --version
   ```
   It should print something like `Rscript (R) version 4.6.1`. On Windows,
   if the command is not found, use the full path, for example
   `"C:\Program Files\R\R-4.6.1\bin\Rscript.exe"`, or run the scripts from
   RStudio instead.
3. **Put all submission files in one folder**, for example `PRA2003/`.
4. **Only if you want to rerun part 1:** download the data files (see
   [The data](#5-the-data)) into the **same folder**.

Your folder should then look like this:

```
PRA2003/
├── README.md
├── README_Frederique.md      ← week 1 README (project goals)
├── Week 4.r                  ← part 1
├── week 4.2.r                ← part 2
├── Week 4 part 3.r           ← part 3
├── Week 4 limitations.r      ← optional: evidence for the limitations
├── sub_sample_results.csv
├── final_results.csv
├── three_sigma_results.csv
├── asymmetry_plot.png
└── output-Set1.txt … output-Set10.txt   ← raw data, only needed for part 1
```

---

## 4. Usage: how to run the code

The scripts must be run **in order**, because each part uses the output of
the one before:

```
 output-Set1..10.txt ──► Week 4.r ──► sub_sample_results.csv
                                          │
                                          ▼
                        week 4.2.r ──► final_results.csv
                                          │
                                          ▼
                   Week 4 part 3.r ──► three_sigma_results.csv + asymmetry_plot.png
```

### Option A: terminal (recommended)

Go to the folder and run the three scripts:

```bash
cd path/to/PRA2003
Rscript "Week 4.r"          # part 1: about 21 min (about 2 min per file). Needs the raw data.
Rscript "week 4.2.r"        # part 2: a few seconds
Rscript "Week 4 part 3.r"   # part 3: a few seconds
```

The quotes are needed because the file names contain spaces.

### Option B: RStudio

1. Open the script in RStudio.
2. Choose **Session → Set Working Directory → To Source File Location**.
3. Click **Source**, or press `Ctrl + Shift + S` (on a Mac: `Cmd + Shift + S`).
4. Repeat for the next script.

### Option C: VS Code

1. Open the `PRA2003` folder in VS Code, with the R extension installed.
2. Open a script and press **Run Source** (the ▶ button).

### Quick check without the raw data

Because the output of part 1 is included, you can reproduce **all results
in a few seconds** without downloading 8 GB:

```bash
Rscript "week 4.2.r"
Rscript "Week 4 part 3.r"
```

### Expected output

The last script, part 3, ends with:

```
Biology strains: 2 of 6 pairs show an asymmetry above 3 sigma: Pseudomonas aeruginosa WT, Streptococcus pneumoniae
All pairs: 2 of 18 pairs show an asymmetry above 3 sigma.
Saved results to three_sigma_results.csv
Saved figure to asymmetry_plot.png
```

**Tested:** the full pipeline (parts 1, 2 and 3) was run from start to
finish on macOS (Apple M2, 16 GB RAM) with R 4.6.1 and reproduces every
number in this README.

### Optional: check the limitations

After part 3, this script prints the evidence tables from
[Limitations](#limitations). It does not write any files.

```bash
Rscript "Week 4 limitations.r"
```

It ends with:

```
Paired test: 3 of 6 pairs above 3 sigma: E. coli WT, Pseudomonas aeruginosa WT, Streptococcus pneumoniae
Standard error of the mean: 2 of 6 pairs above 3 sigma: Pseudomonas aeruginosa WT, Streptococcus pneumoniae
```

---

## 5. The data

### Input: the raw data files

- **Source:** the 10 files `output-Set1.txt` … `output-Set10.txt`, downloaded
  from the course's SURFdrive link (see the course page).
- **Size:** 500,000 events per file, 5,000,000 in total, about 800 MB per file.
- **Not used:** `output-Set0.txt` is a small test file.

Each event starts with a **header line** (`eventNumber nParticles`), followed
by **one line per bacterium** (`px py pz ID`, where px, py and pz are the
momentum components). For example, from `output-Set1.txt`:

```
1 33                                    ← header: event 1 contains 33 bacteria
0.729758 0.610603 -24.4335 -321         ← one bacterium: momentum px py pz, ID = -321
-0.0978169 -0.336023 -4.641 3222
0.539257 -0.0158658 -1.39 -3122
...
19 0                                    ← an empty event, with no bacteria
```

**Empty events** (header `N 0`) contain no measurement, so they are **not
counted**. That leaves **4,617,993 events** (461,650 – 462,025 per file).
Lines that are neither a header nor a bacterium are skipped and reported.

### Output: the included CSV files

**`sub_sample_results.csv`**: one row per ID per file (output of part 1)

| Column | Meaning |
|---|---|
| `code` | ID |
| `name` | Strain name ("Unknown" for IDs that are not one of the 12 strains) |
| `n_events` | Number of non-empty events in this file |
| `total_count` | How often this ID occurs in this file |
| `average_per_event` | `total_count / n_events` |
| `uncertainty` | Poisson uncertainty for this file: `sqrt(total_count) / n_events` |
| `sub_sample` | Which data file the row comes from |

**`final_results.csv`**: one row per ID (output of part 2)

| Column | Meaning |
|---|---|
| `code`, `name` | ID and strain name |
| `total_count` | Total count in all 10 files |
| `average_per_event` | Average per event over the full sample |
| `uncertainty` | Statistical uncertainty: the spread between the 10 sub-samples |

**`three_sigma_results.csv`**: one row per pair X vs -X (output of part 3)

| Column | Meaning |
|---|---|
| `code`, `name`, `partner_name` | The pair |
| `biology_strain` | `TRUE` for the 6 wild-type/mutant pairs |
| `avg_code`, `avg_minus_code` | Average per event of X and -X |
| `difference`, `uncertainty` | X - (-X), and its uncertainty |
| `n_sigma` | Size of the difference in units of its uncertainty |
| `asymmetry_percent`, `asymmetry_uncertainty_percent` | A = (X - (-X)) / (X + (-X)), in % |
| `significant` | `TRUE` if n σ > 3 |

---

## 6. Method

### Part 1: counting (`Week 4.r`)

- **Reading:** each file is read in **chunks of 1,000,000 lines**, so memory
  use stays low.
- **Counting:** each line is classified as a header (2 fields) or a bacterium
  (4 fields). The IDs in each chunk are tallied with `table()` and added to a
  running total. Only the counts are kept, not every individual bacterium.
- **Output per ID and file:** the script records the number of events, the
  total count, the average per event `N / n_events`, and the Poisson
  uncertainty `sqrt(N) / n_events`.
- **Defensive checks:**
  - It stops with a clear message if a file is missing or contains no
    events.
  - It counts and reports malformed lines.
  - It makes sure every ID gets the right strain name.

### Part 2: sub-sampling method (`week 4.2.r`)

Each of the 10 files is treated as **one independent sub-sample**:

- **Average per event** = total count in all 10 files ÷ total number of
  events in all 10 files. This is the exact average over the full sample.
- **Statistical uncertainty** = **standard deviation of the 10 per-file
  averages**. This measures directly how much the result varies from one
  sub-sample to another.

All values are calculated from the exact counts, not from rounded numbers.
If an ID does not appear in a file, that file counts as 0 for that ID.

### Part 3: asymmetry test (`Week 4 part 3.r`)

For every ID X for which -X also occurs, the script calculates the
following. There are 18 such pairs: the 6 biology strains and 12 others.

```
difference   = avg_X - avg_-X
sigma        = sqrt(sigma_X² + sigma_-X²)
n_sigma      = |difference| / sigma

asymmetry A  = (avg_X - avg_-X) / (avg_X + avg_-X)
sigma_A      = 2 / (a + b)² · sqrt(b²·sigma_a² + a²·sigma_b²)     (error propagation)
```

**What is compared:** for each pair, the average count of the wild type (X)
is compared with the average count of its mutant (-X). The asymmetry A
gives the size of that difference as a percentage. `n_sigma` is the
**Z-score** of the comparison. It says how many standard deviations the
wild type and the mutant are apart, and so whether the asymmetry is larger
than random fluctuations can explain.

**Decision rule, set before looking at the results:**

- **n σ > 3: significant asymmetry.** A pure statistical fluctuation this
  large happens by chance only **0.27%** of the time.
- **n σ ≤ 3: no significant asymmetry.** The data are consistent with
  A = 0. In that case, A + 3·σ_A gives an **upper limit** on how large a real
  asymmetry could be.

---

## 7. Results

### 7.1 Average count per event (5M events, 10 sub-samples)

| ID | Strain | Total count | Average / event | Uncertainty |
|---:|---|---:|---:|---:|
| 211 | E. coli WT | 92,126,688 | 19.9495 | ± 0.0327 |
| -211 | E. coli mutant | 91,977,542 | 19.9172 | ± 0.0319 |
| 321 | Bacillus subtilis WT | 11,587,227 | 2.50915 | ± 0.00477 |
| -321 | Bacillus subtilis mutant | 11,560,946 | 2.50346 | ± 0.00550 |
| 2212 | Pseudomonas aeruginosa WT | 5,578,693 | 1.20803 | ± 0.00190 |
| -2212 | P. aeruginosa antibiotic-resistant | 5,468,447 | 1.18416 | ± 0.00241 |
| 3122 | Streptococcus pneumoniae | 1,277,330 | 0.276599 | ± 0.00107 |
| -3122 | Capsule-deficient S. pneumoniae | 1,254,690 | 0.271696 | ± 0.000985 |
| 3312 | Mycobacterium tuberculosis | 182,139 | 0.0394412 | ± 0.000284 |
| -3312 | Drug-resistant M. tuberculosis | 180,104 | 0.0390005 | ± 0.000402 |
| 3334 | Salmonella enterica | 5,482 | 0.00118710 | ± 0.0000417 |
| -3334 | Salmonella mutant | 5,318 | 0.00115158 | ± 0.0000508 |

The data also contain 26 other IDs that are not one of the 12 strains. They
are listed as "Unknown" in `final_results.csv`.

### 7.2 Asymmetry: wild type vs mutant

![Asymmetry of the 6 wild-type/mutant pairs](asymmetry_plot.png)

*Figure 1: Asymmetry A = (WT - mutant)/(WT + mutant) for each strain, with
1σ error bars. The dashed line at A = 0 means no asymmetry. Blue points
(P. aeruginosa and S. pneumoniae) are more than 3σ away from 0. The gray
points are consistent with 0.*

| Pair | X - (-X) | n σ | Asymmetry A | Asymmetric? |
|---|---:|---:|---:|:---:|
| 211 vs -211 (E. coli) | 0.0323 ± 0.0457 | 0.71 | (0.08 ± 0.12)% | no |
| 321 vs -321 (B. subtilis) | 0.00569 ± 0.00728 | 0.78 | (0.11 ± 0.15)% | no |
| 2212 vs -2212 (P. aeruginosa) | 0.02387 ± 0.00307 | **7.78** | **(1.00 ± 0.13)%** | **yes** |
| 3122 vs -3122 (S. pneumoniae) | 0.00490 ± 0.00145 | **3.37** | **(0.89 ± 0.27)%** | **yes** |
| 3312 vs -3312 (M. tuberculosis) | 0.00044 ± 0.00049 | 0.90 | (0.56 ± 0.63)% | no |
| 3334 vs -3334 (Salmonella) | 0.000036 ± 0.000066 | 0.54 | (1.5 ± 2.8)% | no |

The averages of X and -X are in table 7.1.

#### The answer for each pair, and why

- **211 vs -211, E. coli: no asymmetry.**
  - The WT is 0.0323 per event higher, but the uncertainty on that
    difference is 0.0457, so the difference is only **0.71σ**.
  - A difference this size is expected from statistical fluctuations alone.
  - A = (0.08 ± 0.12)% is consistent with zero, so any real asymmetry is
    smaller than about **0.4%**.
- **321 vs -321, Bacillus subtilis: no asymmetry.**
  - The difference is **0.78σ**, and A = (0.11 ± 0.15)% is consistent with
    zero.
  - Any real asymmetry is smaller than about **0.5%**.
- **2212 vs -2212, Pseudomonas aeruginosa: yes, a clear asymmetry.**
  - The wild type is 0.0239 per event more common than the
    antibiotic-resistant strain.
  - That is **7.8σ**, far above the 3σ threshold, so it cannot be a
    statistical fluctuation.
  - A = **(1.00 ± 0.13)%**.
  - The wild type is also higher in **all 10** sub-samples.
- **3122 vs -3122, Streptococcus pneumoniae: yes, an asymmetry.**
  - The wild type is more common than the capsule-deficient strain by
    **3.37σ**, which is above the 3σ threshold.
  - A = **(0.89 ± 0.27)%**.
  - The wild type is higher in **all 10** sub-samples.
  - Because 3.37σ is only just above 3σ, this result is weaker than the one
    for Pseudomonas.
- **3312 vs -3312, Mycobacterium tuberculosis: no asymmetry.**
  - The difference is **0.90σ**, and A = (0.56 ± 0.63)% is consistent with
    zero.
  - With about 180,000 counts per strain, the data can only rule out
    asymmetries larger than about **2.5%**.
- **3334 vs -3334, Salmonella: no asymmetry.**
  - The difference is **0.54σ**.
  - With only about 5,400 counts per strain, A = (1.5 ± 2.8)% is very
    imprecise, so an asymmetry smaller than about **10%** cannot be ruled
    out.

### 7.3 Asymmetry: the other 12 particle/antiparticle pairs

The same test applied to the 12 other pairs of IDs in the data:

| Pair | X - (-X) | n σ | Asymmetry A | Asymmetric? |
|---|---:|---:|---:|:---:|
| 3212 vs -3212 | 0.00173 ± 0.00097 | 1.78 | (0.57 ± 0.32)% | no |
| 3222 vs -3222 | 0.00173 ± 0.00083 | 2.07 | (0.58 ± 0.28)% | no |
| 3112 vs -3112 | 0.00202 ± 0.00080 | 2.52 | (0.68 ± 0.27)% | no (hint) |
| 3322 vs -3322 | 0.00026 ± 0.00035 | 0.76 | (0.34 ± 0.45)% | no |
| 431 vs -431 | -0.00005 ± 0.00021 | 0.26 | (-0.23 ± 0.90)% | no |
| 531 vs -531 | 0.00002 ± 0.00007 | 0.32 | (1.3 ± 3.9)% | no |
| 4232 vs -4232 | 0.000003 ± 0.000042 | 0.07 | (0.2 ± 3.6)% | no |
| 4132 vs -4132 | 0.000008 ± 0.000045 | 0.17 | (0.7 ± 3.9)% | no |
| 5132 vs -5132 | 0.000014 ± 0.000020 | 0.70 | (14 ± 19)% | no |
| 5232 vs -5232 | -0.000002 ± 0.000017 | 0.13 | (-2 ± 16)% | no |
| 4332 vs -4332 | 0.000004 ± 0.000010 | 0.46 | (13 ± 31)% | no |
| 5332 vs -5332 | 0.0000006 ± 0.0000032 | 0.20 | (18 ± 86)% | no |

- **None of these 12 pairs reaches 3σ,** so none shows a significant
  asymmetry.
- **3112 vs -3112** is the closest, at 2.5σ. That is a hint, but not enough
  to claim an asymmetry.
- **The rarest pairs** (5132, 5232, 4332 and 5332, with 7–260 counts each)
  have such large uncertainties that they cannot show anything either way.

---

## 8. Conclusion

| Question | Answer |
|---|---|
| **1. Average count per event and its uncertainty** | Measured for all 12 strains with the sub-sampling method (table 7.1). The averages range from 19.95 ± 0.03 per event (E. coli WT) to 0.00115 ± 0.00005 per event (Salmonella mutant). |
| **2. Asymmetry between wild type and mutant** | **Yes for 2 of the 6 strains.** P. aeruginosa has A = (1.00 ± 0.13)%, 7.8σ, and S. pneumoniae has A = (0.89 ± 0.27)%, 3.4σ. In both, the **wild type is about 1% more common** than the mutant. The other 4 strains, and all 12 other ID pairs, show no significant asymmetry. |

### Limitations

The results above use the sub-sampling method from the lecture. Two of its
assumptions were checked with `Week 4 limitations.r`, which reproduces
every number in the tables below.

#### Limitation 1: the wild type and mutant are not independent

The test in part 3 adds the uncertainties of X and -X as if the two were
independent. In reality, their per-file averages **rise and fall together**:
a file with more E. coli WT also has more E. coli mutants. Because of this,
the combined uncertainty is overestimated.

The better test is a **paired comparison**: calculate the difference
X - (-X) inside each file, then take the spread of those 10 differences.

| Strain | WT higher in | Correlation X vs -X | σ, independent (part 3) | n σ | σ, paired | n σ, paired |
|---|:---:|:---:|---:|---:|---:|---:|
| E. coli | 10 / 10 files | **0.99** | 0.0457 | 0.71 | 0.00452 | **7.15** |
| B. subtilis | 9 / 10 files | 0.80 | 0.00728 | 0.78 | 0.00329 | 1.73 |
| P. aeruginosa | 10 / 10 files | 0.41 | 0.00307 | **7.78** | 0.00237 | **10.06** |
| S. pneumoniae | 10 / 10 files | 0.84 | 0.00145 | **3.37** | 0.000583 | **8.41** |
| M. tuberculosis | 9 / 10 files | 0.02 | 0.000492 | 0.90 | 0.000488 | 0.90 |
| Salmonella | 6 / 10 files | -0.09 | 0.0000657 | 0.54 | 0.0000685 | 0.52 |

**What this shows:**

- **E. coli changes result.** X and -X are almost perfectly correlated
  (0.99), so part 3 overestimates the uncertainty by about 10×. The paired
  test gives **7.15σ**, and the WT is higher in all 10 files.
- **B. subtilis stays below 3σ** (1.73σ), even with the paired test.
- **P. aeruginosa and S. pneumoniae become stronger** (10.06σ and 8.41σ).
- **M. tuberculosis and Salmonella do not change,** because X and -X are
  uncorrelated for these strains (correlation of about 0).

#### Limitation 2: the spread of one file vs the uncertainty of the mean

Part 2 uses the **standard deviation of the 10 per-file averages** as the
uncertainty. That is the uncertainty of **one** file. The final average
uses all 10 files, so its uncertainty (the standard error of the mean) is
smaller by a factor √10 ≈ 3.16.

| Strain | σ, sd of 10 files (part 3) | n σ | σ, sd / √10 | n σ, sd / √10 |
|---|---:|---:|---:|---:|
| E. coli | 0.0457 | 0.71 | 0.0145 | 2.24 |
| B. subtilis | 0.00728 | 0.78 | 0.00230 | 2.47 |
| P. aeruginosa | 0.00307 | **7.78** | 0.000971 | **24.59** |
| S. pneumoniae | 0.00145 | **3.37** | 0.000459 | **10.69** |
| M. tuberculosis | 0.000492 | 0.90 | 0.000156 | 2.83 |
| Salmonella | 0.0000657 | 0.54 | 0.0000208 | 1.71 |

**What this shows:** all n σ values grow by 3.16×, but the conclusion stays
the same. P. aeruginosa and S. pneumoniae are the only pairs above 3σ.
M. tuberculosis (2.83σ) gets close.

#### Summary of the limitations

| Method | Pairs above 3σ |
|---|---|
| Sub-sampling, independent (**used in this README**) | P. aeruginosa, S. pneumoniae |
| Paired comparison per file | E. coli, P. aeruginosa, S. pneumoniae |
| Standard error of the mean (sd / √10) | P. aeruginosa, S. pneumoniae |

- **The main result is robust:** P. aeruginosa and S. pneumoniae are above
  3σ with every method.
- **The E. coli result depends on the method.** Its "no asymmetry" should be
  read as "not shown with the conservative method", not as "no effect".
- The sub-sampling method was kept for the main results because it is the
  method from the lecture and gives the most conservative uncertainty.

---

## 9. Troubleshooting

| Problem | Solution |
|---|---|
| `Can't find the file at 'output-Set1.txt'` | The data files are not in the working directory. Put them in the same folder as the scripts, or in RStudio use *Session → Set Working Directory → To Source File Location*. |
| `Can't find 'sub_sample_results.csv' - run Part 1 first.` | Run the scripts in order, or make sure the included CSV file is in the same folder. |
| `Rscript: command not found` | R is not on your PATH. Use the full path to `Rscript`, or run the scripts from RStudio. |
| `cannot open file 'Week': No such file or directory` | You forgot the quotes. Use `Rscript "Week 4.r"`. |
| Part 1 seems to hang | It is still working, since each file takes about 2 minutes. You can uncomment the `cat("Processing chunk"...)` line in `Week 4.r` to see its progress. |
