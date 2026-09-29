# PRA2003 — Wild-type vs Mutant Bacteria in 5 Million Events (Biology Track)

## Contents

1. [Overview](#1-overview)
2. [Data](#2-data)
3. [Method](#3-method)
4. [Main result](#4-main-result)
5. [Is there an asymmetry between wild type and mutant?](#5-is-there-an-asymmetry-between-wild-type-and-mutant)
6. [Consistency check: Poisson vs sub-sampling](#6-consistency-check-poisson-vs-sub-sampling)
7. [Asymmetry as a function of momentum](#7-asymmetry-as-a-function-of-momentum)
8. [Discussion and limitations](#8-discussion-and-limitations)
9. [How to install and reproduce](#9-how-to-install-and-reproduce)
10. [Files in this submission](#10-files-in-this-submission)

---

## 1. Overview

Each bacterial strain in the data has a **wild-type** form (positive ID) and
a **mutant** form (the same ID with a minus sign), for example an
antibiotic-resistant strain. If a mutation changes how common a strain is,
the wild type and its mutant will occur with different frequencies.

**Why this matters:** mutations such as antibiotic resistance often carry a
*fitness cost*, which makes the mutant slightly less successful. Measuring
whether mutants are less common than their wild type is a direct way to
look for such a cost. Because the expected effect is small, a large sample
and a careful uncertainty estimate are needed.

**Questions** (from the week 1 README, `README_week1_FM.md`):

| Question | Where it is answered |
|---|---|
| 1. What is the **average number of each strain per event**, and its statistical uncertainty? | [Section 4](#4-main-result) |
| 2. Is there an **asymmetry between each wild type and its mutant**? | [Section 5](#5-is-there-an-asymmetry-between-wild-type-and-mutant) |
| 3. Does the asymmetry depend on **momentum**? | [Section 7](#7-asymmetry-as-a-function-of-momentum) (week 3) |

**Hypothesis:** if a mutation has no effect, the wild type and mutant have
the same average (asymmetry A = 0). An asymmetry is only claimed above
**3σ**.

---

## 2. Data

| Property | Value |
|---|---|
| Files used | `output-Set1.txt` … `output-Set10.txt`, from the course SURFdrive (not included: ~800 MB each) |
| Not used | `output-Set0.txt`, a test file with 1 event |
| Events | 500,000 per file, **5,000,000** in total |
| Non-empty events (used) | **4,617,993** in total, 461,650 – 462,025 per file |
| Different IDs found | 38: the 12 strains below + 26 other IDs |

Each file is plain text, made of repeating blocks:

```
<eventNumber> <nParticles>          ← header line
<px> <py> <pz> <ID>                  ← one line per bacterium (momentum + ID)
...
```

Events with 0 bacteria (e.g. a header `19 0`) contain no measurement, so
they are **not counted** in the number of events.

### Bacterial ID → strain

| ID | Strain | ID | Strain |
|---:|---|---:|---|
| 211 | *E. coli* WT | -211 | *E. coli* mutant |
| 321 | *Bacillus subtilis* WT | -321 | *B. subtilis* mutant |
| 2212 | *Pseudomonas aeruginosa* WT | -2212 | *P. aeruginosa* antibiotic-resistant |
| 3122 | *Streptococcus pneumoniae* | -3122 | Capsule-deficient *S. pneumoniae* |
| 3312 | *Mycobacterium tuberculosis* | -3312 | Drug-resistant *M. tuberculosis* |
| 3334 | *Salmonella enterica* | -3334 | *Salmonella* mutant |

---

## 3. Method

**Step 1: counting per file** (`Week4_part1_counting_FM.r`). Each of the 10 files is read
in chunks of 1,000,000 lines, so a file is never loaded into memory at
once. The IDs in each chunk are tallied with `table()`. For each ID and
file the script saves the count, the number of non-empty events, the
average per event and the Poisson uncertainty `sqrt(count) / events`.

**Step 2: combining the 10 sub-samples** (`Week4_part2_subsampling_FM.r`), with the
**sub-sampling method** from the lecture:

| Quantity | How it is calculated |
|---|---|
| **Central value** | total count ÷ total events over all 10 files (the per-file averages, weighted by their number of events) |
| **Statistical uncertainty** | standard deviation of the 10 per-file averages, used directly (no division by √10), as in the lecture |

**Step 3: asymmetry test** (`Week4_part3_asymmetry_FM.r`) for every pair X vs -X. As
on the lecture slides, the difference is calculated **in each sub-sample
separately**, and the spread of those 10 results is its uncertainty:

| Quantity | Central value (all 10 files) | Uncertainty (sub-sampling) |
|---|---|---|
| Difference | avg_X − avg_-X | standard deviation of the 10 per-file differences (no division by √10) |
| Asymmetry A | (avg_X − avg_-X) / (avg_X + avg_-X) | standard deviation of the 10 per-file values of A |
| n σ | \|difference\| / uncertainty | |

The central difference over all 10 files equals the mean of the 10
per-file differences, weighted by their number of events. (The unweighted
mean gives the same values to the precision shown.)

X and -X are counted in the **same events**, so their numbers rise and fall
together from file to file. Taking the difference inside each file accounts
for this correlation.

All values are calculated from the exact counts, not from rounded per-file
averages, because rounding to a few decimals noticeably changes the result
for the rarer strains.

**Decision rule, fixed in advance:** n σ > 3 means a significant asymmetry
(lecture slide 11: a 1-in-740 chance of a statistical fluctuation).
Otherwise, A + 3σ_A is an upper limit on any real asymmetry.

---

## 4. Main result

The number of bacteria of each strain per event, from the entire 5M-event
sample (**4,617,993 non-empty events**), with the statistical uncertainty
from the sub-sampling method. Unit: bacteria per event.

| ID | Strain | Total count | Average / event | Statistical uncertainty |
|---:|---|---:|---:|---:|
| 211 | *E. coli* WT | 92,126,688 | 19.9495 | ± 0.0327 |
| -211 | *E. coli* mutant | 91,977,542 | 19.9172 | ± 0.0319 |
| 321 | *Bacillus subtilis* WT | 11,587,227 | 2.50915 | ± 0.00477 |
| -321 | *B. subtilis* mutant | 11,560,946 | 2.50346 | ± 0.00550 |
| 2212 | *Pseudomonas aeruginosa* WT | 5,578,693 | 1.20803 | ± 0.00190 |
| -2212 | *P. aeruginosa* antibiotic-resistant | 5,468,447 | 1.18416 | ± 0.00241 |
| 3122 | *Streptococcus pneumoniae* | 1,277,330 | 0.276599 | ± 0.00107 |
| -3122 | Capsule-deficient *S. pneumoniae* | 1,254,690 | 0.271696 | ± 0.000985 |
| 3312 | *Mycobacterium tuberculosis* | 182,139 | 0.0394412 | ± 0.000284 |
| -3312 | Drug-resistant *M. tuberculosis* | 180,104 | 0.0390005 | ± 0.000402 |
| 3334 | *Salmonella enterica* | 5,482 | 0.00118710 | ± 0.0000417 |
| -3334 | *Salmonella* mutant | 5,318 | 0.00115158 | ± 0.0000508 |

The 26 other IDs are listed as "Unknown" in `final_results_FM.csv`.

**Systematic uncertainty:** not evaluated. The data contain no information
about the measuring device (such as detection efficiency), so only the
statistical uncertainty is reported.

---

## 5. Is there an asymmetry between wild type and mutant?

![Asymmetry of the 6 wild-type/mutant pairs](asymmetry_plot_FM.png)

*Figure 1: Asymmetry A = (WT − mutant) / (WT + mutant) for each strain, with
1σ error bars. Blue points are more than 3σ from A = 0 (dashed line).*

| Pair | X − (−X) | n σ | Asymmetry A | Asymmetric? | WT higher in | Upper limit on A |
|---|---:|---:|---:|:---:|:---:|---:|
| 211 vs -211 (E. coli) | 0.0323 ± 0.0045 | **7.15** | **(0.08 ± 0.01)%** | **yes** | 10 / 10 files | – |
| 321 vs -321 (B. subtilis) | 0.00569 ± 0.00329 | 1.73 | (0.11 ± 0.07)% | no | 9 / 10 files | 0.3% |
| 2212 vs -2212 (P. aeruginosa) | 0.02387 ± 0.00237 | **10.06** | **(1.00 ± 0.10)%** | **yes** | 10 / 10 files | – |
| 3122 vs -3122 (S. pneumoniae) | 0.00490 ± 0.00058 | **8.41** | **(0.89 ± 0.11)%** | **yes** | 10 / 10 files | – |
| 3312 vs -3312 (M. tuberculosis) | 0.00044 ± 0.00049 | 0.90 | (0.56 ± 0.62)% | no | 9 / 10 files | 2.4% |
| 3334 vs -3334 (Salmonella) | 0.000036 ± 0.000069 | 0.52 | (1.5 ± 2.9)% | no | 6 / 10 files | 10% |

**Conclusion: 3 of the 6 strains show a significant asymmetry**, and in all
three the wild type is the more common form:

- ***P. aeruginosa*** (10.1σ) and ***S. pneumoniae*** (8.4σ): the wild type
  is about **1%** more common than the resistant / capsule-deficient form.
- ***E. coli*** (7.2σ): a very small asymmetry of 0.08%, significant only
  because *E. coli* is so common that it is measured very precisely.
- ***B. subtilis*, *M. tuberculosis* and *Salmonella***: consistent with
  A = 0. The upper limit shows how large a real asymmetry could still be.

### The other 12 pairs of IDs

The same test for the 12 other pairs X vs -X in the data:

| Pair | X − (−X) | n σ | Asymmetry A | Asymmetric? |
|---|---:|---:|---:|:---:|
| 3212 vs -3212 | 0.00173 ± 0.00076 | 2.26 | (0.57 ± 0.25)% | no |
| 3222 vs -3222 | 0.00173 ± 0.00059 | 2.93 | (0.58 ± 0.20)% | no |
| 3112 vs -3112 | 0.00202 ± 0.00046 | **4.43** | **(0.68 ± 0.15)%** | **yes** |
| 3322 vs -3322 | 0.00026 ± 0.00030 | 0.88 | (0.34 ± 0.39)% | no |
| 431 vs -431 | -0.00005 ± 0.00017 | 0.32 | (-0.23 ± 0.74)% | no |
| 531 vs -531 | 0.00002 ± 0.00008 | 0.32 | (1.3 ± 3.9)% | no |
| 4232 vs -4232 | 0.000003 ± 0.000033 | 0.09 | (0.2 ± 2.7)% | no |
| 4132 vs -4132 | 0.000008 ± 0.000041 | 0.19 | (0.7 ± 3.5)% | no |
| 5132 vs -5132 | 0.000014 ± 0.000020 | 0.70 | (14 ± 19)% | no |
| 5232 vs -5232 | -0.000002 ± 0.000016 | 0.14 | (-2 ± 15)% | no |
| 4332 vs -4332 | 0.000004 ± 0.000012 | 0.37 | (13 ± 36)% | no |
| 5332 vs -5332 | 0.0000006 ± 0.0000035 | 0.18 | (18 ± 86)% | no |

- **1 of the 12 pairs is above 3σ:** 3112, with A = (0.68 ± 0.15)%. 3222
  (2.93σ) is just below the threshold.
- **The rarest pairs** (5132, 5232, 4332, 5332) have too few counts to show
  anything either way.

---

## 6. Consistency check: Poisson vs sub-sampling

Step 1 also gives each file's **Poisson uncertainty**, `sqrt(count) /
events`. If the counts were pure Poisson noise, the 10 per-file averages
would scatter by about that much. The table compares that prediction with
the real spread between the files (the uncertainty used in section 4).

| Strain | Real spread (10 files) | Poisson, one file | Ratio |
|---|---:|---:|---:|
| *E. coli* WT | 0.0327 | 0.00657 | **5.0×** |
| *E. coli* mutant | 0.0319 | 0.00657 | **4.9×** |
| *B. subtilis* WT | 0.00477 | 0.00233 | **2.0×** |
| *B. subtilis* mutant | 0.00550 | 0.00233 | **2.4×** |
| *P. aeruginosa* WT | 0.00190 | 0.00162 | 1.2× |
| *P. aeruginosa* resistant | 0.00241 | 0.00160 | 1.5× |
| *S. pneumoniae* | 0.00107 | 0.000774 | 1.4× |
| Capsule-deficient *S. pneumoniae* | 0.000985 | 0.000767 | 1.3× |
| *M. tuberculosis* | 0.000284 | 0.000292 | 1.0× |
| Drug-resistant *M. tuberculosis* | 0.000402 | 0.000291 | 1.4× |
| *Salmonella enterica* | 0.0000417 | 0.0000507 | 0.8× |
| *Salmonella* mutant | 0.0000508 | 0.0000499 | 1.0× |

**Conclusion:** for rare strains the ratio is about 1, so Poisson works.
For common strains the real spread is up to **5× larger** than Poisson
predicts. The number of bacteria per event varies a lot, so counts within
one event are not independent. This is why the sub-sampling method is used:
it measures the spread directly instead of assuming Poisson statistics. The
same effect makes X and -X move together, which is why the asymmetry test
takes the difference inside each file.

---

## 7. Asymmetry as a function of momentum

This was the week 3 deliverable, using one file (`output-Set1.txt`) and
Poisson uncertainties. Bacteria were grouped by momentum
|p| = √(px² + py² + pz²). The results are in
`asymmetry_vs_momentum_Set1_FM.csv`. For the two strains with the clearest
asymmetry:

| Momentum bin | *P. aeruginosa*: A | n σ | *S. pneumoniae*: A | n σ |
|---|---:|---:|---:|---:|
| 0 – 1 | 0.66% | 2.7 | 0.63% | 1.1 |
| 1 – 2 | 1.35% | 6.3 | 0.81% | 1.7 |
| 2 – 4 | 0.89% | 4.4 | 0.32% | 0.8 |
| 4 – 8 | 0.99% | 4.6 | 0.78% | 1.8 |
| ≥ 8 | 1.47% | 7.6 | 1.43% | 3.7 |

The wild type is more common at **every momentum**, with the largest
asymmetry at high momentum. This uses only 1 of the 10 files and Poisson
uncertainties, which section 6 shows are too small for common strains, so
it is an indication rather than a final result.

---

## 8. Discussion and limitations

### What the asymmetry means

| Observation | Value |
|---|---|
| Significant asymmetries (strains) | 1.00%, 0.89% and 0.08% |
| Strains where the wild type is more common | **6 of 6** |
| All pairs where X is more common than −X | **16 of 18** |
| Chance of 16 of 18 if there were no effect (sign test) | p = 0.0013 |

- **The wild type is the more common form** in every strain. A possible
  explanation is the fitness cost of the mutation (antibiotic resistance,
  loss of the capsule), which makes the mutant slightly less common.
- **The same direction appears almost everywhere.** With no real effect,
  each sign would be a coin flip, so 16 of 18 in the same direction is very
  unlikely by chance. The non-significant pairs may still have a small
  asymmetry in the same direction, too small to measure one by one.
- **The effect is small**: at most about 1%, which is why 5 million events
  are needed to see it.

### Limitations

The asymmetry test was compared with the more common error-propagation
method using `Week4_limitations_FM.r`:

- **Propagation** (σ = √(σ_X² + σ_-X²)) treats X and -X as independent.
- **But X and -X are correlated:** they are counted in the same events, so
  they rise and fall together from file to file. For the common strains the
  correlation is strong, so propagation **overestimates** σ.

| Strain | Correlation X vs -X | n σ (used: per-file differences) | n σ (propagation) |
|---|:---:|---:|---:|
| E. coli | 0.99 | **7.15** | 0.71 |
| B. subtilis | 0.80 | 1.73 | 0.78 |
| P. aeruginosa | 0.41 | **10.06** | **7.78** |
| S. pneumoniae | 0.84 | **8.41** | **3.37** |
| M. tuberculosis | 0.02 | 0.90 | 0.90 |
| Salmonella | -0.09 | 0.52 | 0.54 |

- ***P. aeruginosa* and *S. pneumoniae* are robust:** above 3σ with both
  methods.
- ***E. coli* depends on the correlation:** it is only significant when the
  correlation between X and -X is taken into account, which is why the
  per-file differences are used.
- **Other limitations:** no systematic uncertainty is included
  (section 4), and the momentum result uses only one file (section 7).

---

## 9. How to install and reproduce

### Requirements

| Requirement | Details |
|---|---|
| **R** | 4.0 or newer (tested with **R 4.6.1**), free from [cran.r-project.org](https://cran.r-project.org) |
| **R packages** | None: only base R is used |
| **Disk space** | About 8 GB, only to rerun step 1 on the raw data |
| **Memory** | No special requirement: step 1 reads the files in chunks |

### Installation

1. **Install R** from [cran.r-project.org](https://cran.r-project.org) and
   follow the installer for your operating system.
2. **Check that it works:** open a terminal and type `Rscript --version`.
   It should print `Rscript (R) version 4.x`. On Windows, if the command is
   not found, use the full path (e.g.
   `"C:\Program Files\R\R-4.6.1\bin\Rscript.exe"`) or use RStudio.
3. **Put all submission files in one folder** (see
   [section 10](#10-files-in-this-submission)).
4. **Only to rerun step 1:** download `output-Set1.txt` … `output-Set10.txt`
   from the course SURFdrive into the same folder.

### Running the analysis

Run the scripts **in order** from the project folder. Each step uses the
output of the one before.

| Step | Script | Needs | Produces | Time |
|:---:|---|---|---|---|
| 1 | `Week4_part1_counting_FM.r` | `output-Set1.txt` … `output-Set10.txt` | `sub_sample_results_FM.csv` | ~21 min |
| 2 | `Week4_part2_subsampling_FM.r` | `sub_sample_results_FM.csv` | `final_results_FM.csv` | seconds |
| 3 | `Week4_part3_asymmetry_FM.r` | `final_results_FM.csv`, `sub_sample_results_FM.csv` | `three_sigma_results_FM.csv`, `asymmetry_plot_FM.png` | seconds |
| opt. | `Week4_limitations_FM.r` | the three CSV files | printed table only | seconds |

```bash
cd path/to/PRA2003
Rscript Week4_part1_counting_FM.r      # step 1 (skip this if you do not have the raw data)
Rscript Week4_part2_subsampling_FM.r   # step 2
Rscript Week4_part3_asymmetry_FM.r     # step 3
Rscript Week4_limitations_FM.r         # optional
```

- **No raw data?** Skip step 1: the included `sub_sample_results_FM.csv`
  reproduces all results in a few seconds.
- **RStudio:** open a script, choose *Session → Set Working Directory →
  To Source File Location*, then click **Source**. Repeat for each step.

**Expected output** at the end of step 3:

```
Biology strains: 3 of 6 pairs show an asymmetry above 3 sigma: E. coli WT, Pseudomonas aeruginosa WT, Streptococcus pneumoniae
All pairs: 4 of 18 pairs show an asymmetry above 3 sigma.
```

**Tested:** the full pipeline (steps 1–3) was run from start to finish on
macOS with R 4.6.1 and reproduces every number in this README.

### Troubleshooting

| Problem | Solution |
|---|---|
| `Can't find the file at 'output-Set1.txt'` | Put the data files in the same folder as the scripts, and set that folder as the working directory. |
| `Can't find 'sub_sample_results_FM.csv'` | Run the steps in order, or keep the included CSV files in the folder. |
| `Rscript: command not found` | Use the full path to `Rscript`, or run the scripts from RStudio. |
| `cannot open file 'Week4_...'` | Check that your terminal is in the submission folder (`cd path/to/PRA2003`). |
| Step 1 seems to hang | Each file takes about 2 minutes. Uncomment the `cat("Processing chunk"...)` line in `Week4_part1_counting_FM.r` to see progress. |

---

## 10. Files in this submission

| File | Purpose |
|---|---|
| `Week4_part1_counting_FM.r` | Step 1: counts every ID in each of the 10 data files |
| `Week4_part2_subsampling_FM.r` | Step 2: combines the 10 sub-samples into the final averages and uncertainties |
| `Week4_part3_asymmetry_FM.r` | Step 3: asymmetry test for every pair X vs -X, and Figure 1 |
| `Week4_limitations_FM.r` | Optional: compares the asymmetry test with error propagation |
| `sub_sample_results_FM.csv` | Per-file results: count, events, average and Poisson uncertainty per ID (371 rows) |
| `final_results_FM.csv` | Final average ± uncertainty per ID (38 IDs) |
| `three_sigma_results_FM.csv` | Difference, n σ and asymmetry for all 18 pairs |
| `asymmetry_plot_FM.png` | Figure 1 |
| `asymmetry_vs_momentum_Set1_FM.csv` | Week 3 result used in section 7 |
| `README_week1_FM.md` | The week 1 README with the project goals |
| `README_FM.md` | This file |
