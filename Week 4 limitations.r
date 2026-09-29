# PRA2003 - Week 4 Deliverable (optional): check of the asymmetry method
# Reads the results of Parts 1, 2 and 3 and, for the 6 wild-type/mutant
# pairs, compares the method used in Part 3 (spread of the per-file
# differences) with error propagation: combine the uncertainties of X and
# -X from Part 2 as if they were independent,
#   sigma = sqrt(sigma_X^2 + sigma_minusX^2).
# This ignores that X and -X rise and fall together from file to file,
# which is measured here with the correlation of their per-file averages.
#
# Nothing is written to disk: the script only prints the table.

files_needed <- c("sub_sample_results.csv", "final_results.csv",
                  "three_sigma_results.csv")
for (f in files_needed) {
  if (!file.exists(f)) {
    stop(paste0("Can't find '", f, "' - run Parts 1, 2 and 3 first."))
  }
}

sub_samples <- read.csv("sub_sample_results.csv")
final_results <- read.csv("final_results.csv")
three_sigma <- read.csv("three_sigma_results.csv")
strains <- three_sigma[three_sigma$biology_strain, ]
threshold <- 3

rows <- list()
for (i in seq_len(nrow(strains))) {
  code <- strains$code[i]
  wt <- sub_samples[sub_samples$code == code, ]
  mutant <- sub_samples[sub_samples$code == -code, ]
  both <- merge(wt, mutant, by = "sub_sample", suffixes = c("_x", "_minus"))

  #Exact per-file averages, not the rounded ones in the CSV
  avg_x <- both$total_count_x / both$n_events_x
  avg_minus <- both$total_count_minus / both$n_events_minus

  sigma_x <- final_results$uncertainty[final_results$code == code]
  sigma_minus <- final_results$uncertainty[final_results$code == -code]
  sigma_propagated <- sqrt(sigma_x^2 + sigma_minus^2)

  rows[[i]] <- data.frame(
    code = code,
    name = strains$name[i],
    correlation = round(cor(avg_x, avg_minus), 2),
    n_sigma_used = strains$n_sigma[i],
    n_sigma_propagated = round(abs(strains$difference[i]) / sigma_propagated, 2)
  )
}
checks <- do.call(rbind, rows)

options(scipen = 999)
print(checks, row.names = FALSE)

for (method in c("n_sigma_used", "n_sigma_propagated")) {
  above <- checks$name[checks[[method]] > threshold]
  cat(method, ":", length(above), "of", nrow(checks), "pairs above", threshold,
      "sigma:", paste(above, collapse = ", "), "\n")
}
