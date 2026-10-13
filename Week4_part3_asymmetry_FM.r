# PRA2003 - Week 4 Deliverable (Part 3): asymmetry test for every pair of IDs
# Reads the final results over all 10 files (from Part 2) and the per-file
# results (from Part 1). For every pair of IDs X and -X that both occur in
# the data (normal strain vs mutant, or particle vs antiparticle), it tests
# whether their average count per event is significantly different.
#
# The asymmetry is the plain difference, in bacteria per event (no division):
#   A_i = avg_X - avg_minusX in file i
#
# Sub-sampling method: the difference is calculated in each of the 10 files
# separately. The mean asymmetry is the average of those 10 results, and
# their spread (standard deviation) is its uncertainty. X and -X come from
# the same events, so this also takes their correlation into account.
#   A          = mean(difference in file 1 ... file 10)
#   sigma      = sd(difference in file 1 ... file 10)
#   n_sigma    = |A| / sigma
#
# Decision rule (fixed before looking at the numbers):
#   n_sigma > 3  -> significant asymmetry
#   n_sigma <= 3 -> no significant asymmetry (consistent with zero)

input_csv <- "final_results_FM.csv"
sub_sample_csv <- "sub_sample_results_FM.csv"
output_csv <- "three_sigma_results_FM.csv"
threshold <- 3

if (!file.exists(input_csv)) {
  stop(paste0("Can't find '", input_csv, "' - run Part 2 first."))
}
if (!file.exists(sub_sample_csv)) {
  stop(paste0("Can't find '", sub_sample_csv, "' - run Part 1 first."))
}

final_results <- read.csv(input_csv)
sub_samples <- read.csv(sub_sample_csv)

# ---- Find every pair X / -X present in the data ----
#The 6 biology-track strains come first, then the other pairs by abundance.
strain_codes <- c(211, 321, 2212, 3122, 3312, 3334)
positive_codes <- final_results$code[final_results$code > 0]
paired_codes <- positive_codes[-positive_codes %in% final_results$code]
other_codes <- setdiff(paired_codes, strain_codes)
pair_codes <- c(strain_codes, other_codes)

#Protection: every biology-track strain needs its mutant partner in the results
if (!all(strain_codes %in% paired_codes)) {
  stop("A wild-type or mutant strain code is missing from final_results_FM.csv.")
}

x <- final_results[match(pair_codes, final_results$code), ]
anti <- final_results[match(-pair_codes, final_results$code), ]

a <- x$average_per_event
b <- anti$average_per_event


# ---- Uncertainties from the spread of the 10 sub-samples ----
events_per_file <- tapply(sub_samples$n_events, sub_samples$sub_sample, function(v) v[1])
file_names <- names(events_per_file)

#Exact per-file averages; a code missing from a file counts as 0 there
per_file_average <- function(code) {
  counts <- setNames(rep(0, length(file_names)), file_names)
  rows <- sub_samples[sub_samples$code == code, ]
  counts[rows$sub_sample] <- rows$total_count
  counts / events_per_file
}

difference <- numeric(length(pair_codes))   #The mean asymmetry A, in bacteria per event
sigma <- numeric(length(pair_codes))
for (i in seq_along(pair_codes)) {
  file_a <- per_file_average(pair_codes[i])
  file_b <- per_file_average(-pair_codes[i])
  difference[i] <- mean(file_a - file_b)   #Average of the 10 sub-samples
  sigma[i] <- sd(file_a - file_b)          #Spread of the 10 sub-samples
}
n_sigma <- abs(difference) / sigma

results <- data.frame(
  code = pair_codes,
  name = x$name,
  partner_name = anti$name,
  biology_strain = pair_codes %in% strain_codes,
  avg_code = a,
  avg_minus_code = b,
  difference = signif(difference, 4),
  uncertainty = signif(sigma, 3),
  n_sigma = round(n_sigma, 2),
  significant = n_sigma > threshold
)

options(scipen = 999)
print(results[, c("code", "name", "difference", "uncertainty", "n_sigma", "significant")],
      row.names = FALSE)

strains <- results[results$biology_strain, ]
cat("\nBiology strains:", sum(strains$significant), "of", nrow(strains),
    "pairs show an asymmetry above", threshold, "sigma:",
    paste(strains$name[strains$significant], collapse = ", "), "\n")
cat("All pairs:", sum(results$significant), "of", nrow(results),
    "pairs show an asymmetry above", threshold, "sigma.\n")

write.csv(results, output_csv, row.names = FALSE)
cat("\nSaved results to", output_csv, "\n")

# ---- Figure ----
#Figure 1 (asymmetry_plot_FM.png) is drawn in Python from the CSV above,
#because the course asks for all plots in Python: run Week4_part3_plot_FM.py.
