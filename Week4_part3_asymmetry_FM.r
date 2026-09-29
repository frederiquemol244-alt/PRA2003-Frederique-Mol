# PRA2003 - Week 4 Deliverable (Part 3): asymmetry test for every pair of IDs
# Reads the final results over all 10 files (from Part 2) and the per-file
# results (from Part 1). For every pair of IDs X and -X that both occur in
# the data (normal strain vs mutant, or particle vs antiparticle), it tests
# whether their average count per event is significantly different.
#
# Central values come from the full sample (Part 2):
#   difference = avg_X - avg_minusX
#   asymmetry  A = (avg_X - avg_minusX) / (avg_X + avg_minusX)
#
# Uncertainties use the sub-sampling method: the difference and the
# asymmetry are calculated in each of the 10 files separately, and the
# spread (standard deviation) of those 10 results is the uncertainty.
# X and -X come from the same events, so this also takes their correlation
# into account.
#   sigma      = sd(difference in file 1 ... file 10)
#   sigma_A    = sd(A in file 1 ... file 10)
#   n_sigma    = |difference| / sigma
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

difference <- a - b
asymmetry <- (a - b) / (a + b)

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

sigma <- numeric(length(pair_codes))
sigma_asymmetry <- numeric(length(pair_codes))
for (i in seq_along(pair_codes)) {
  file_a <- per_file_average(pair_codes[i])
  file_b <- per_file_average(-pair_codes[i])
  sigma[i] <- sd(file_a - file_b)
  #A is undefined in a file where neither code occurs; such files are skipped
  sigma_asymmetry[i] <- sd((file_a - file_b) / (file_a + file_b), na.rm = TRUE)
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
  asymmetry_percent = round(100 * asymmetry, 3),
  asymmetry_uncertainty_percent = round(100 * sigma_asymmetry, 3),
  significant = n_sigma > threshold
)

options(scipen = 999)
print(results[, c("code", "name", "difference", "uncertainty", "n_sigma",
                  "asymmetry_percent", "asymmetry_uncertainty_percent", "significant")],
      row.names = FALSE)

strains <- results[results$biology_strain, ]
cat("\nBiology strains:", sum(strains$significant), "of", nrow(strains),
    "pairs show an asymmetry above", threshold, "sigma:",
    paste(strains$name[strains$significant], collapse = ", "), "\n")
cat("All pairs:", sum(results$significant), "of", nrow(results),
    "pairs show an asymmetry above", threshold, "sigma.\n")

write.csv(results, output_csv, row.names = FALSE)
cat("\nSaved results to", output_csv, "\n")

# ---- Figure: asymmetry of the 6 biology pairs, with 1-sigma error bars ----
plot_file <- "asymmetry_plot_FM.png"
short_names <- c("E. coli", "B. subtilis", "P. aeruginosa",
                 "S. pneumoniae", "M. tuberculosis", "Salmonella")
y <- rev(seq_len(nrow(strains)))   #First strain at the top
A <- strains$asymmetry_percent
err <- strains$asymmetry_uncertainty_percent
point_colour <- ifelse(strains$significant, "#2a78d6", "#8b8a85")
x_range <- c(-2, 6)   #Asymmetry range shown, in %

png(plot_file, width = 1600, height = 900, res = 200)
par(mar = c(4.5, 9, 3, 1), family = "sans", col.axis = "#52514e", fg = "#52514e")
plot(NA, xlim = x_range, ylim = c(0.5, nrow(strains) + 0.5),
     yaxt = "n", xlab = "Asymmetry A = (WT - mutant) / (WT + mutant)  [%]", ylab = "",
     main = "Wild type vs mutant asymmetry (5M events, error bars = 1 sigma)",
     cex.main = 0.95, font.main = 1, bty = "n")
abline(v = seq(x_range[1], x_range[2], 1), col = "#e8e7e2")
abline(v = 0, col = "#52514e", lty = 2)
axis(2, at = y, labels = short_names, las = 1, tick = FALSE)
#Error bars are clipped to the plot range
arrows(pmax(A - err, x_range[1]), y, pmin(A + err, x_range[2]), y, angle = 90, code = 3,
       length = 0.04, col = point_colour, lwd = 2)
points(A, y, pch = 19, cex = 1.3, col = point_colour)
text(pmin(A + err, x_range[2]), y, paste0(format(strains$n_sigma, nsmall = 2), " sigma"),
     pos = 4, cex = 0.75, col = "#0b0b0b", xpd = TRUE)
legend("topright", legend = c("> 3 sigma: significant", "<= 3 sigma: not significant"),
       col = c("#2a78d6", "#8b8a85"), pch = 19, bty = "n", cex = 0.75)
invisible(dev.off())
cat("Saved figure to", plot_file, "\n")
