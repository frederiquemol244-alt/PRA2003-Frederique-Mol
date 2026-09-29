# PRA2003 - Week 4 Deliverable (Part 2 of 3)
# Reads the combined per-sub-sample results (from Part 1) and applies the
# sub-sampling method. For each particle code:
#   - the central value is the average over the full sample:
#       (total count in all 10 files) / (total events in all 10 files)
#     This is exact, unbiased, and needs no weights.
#   - the statistical uncertainty is the standard deviation of the 10
#     sub-sample averages (the "spread" method described in the slides)
#
# All averages are recalculated from the exact total_count and n_events
# columns, not from the rounded average_per_event column, so rare codes
# keep their precision.
#
# A code that never appears in a file does not get a row for that file
# in Part 1. Such a file is counted here as a sub-sample with 0 counts,
# so that every code is averaged over all 10 sub-samples.

input_csv <- "sub_sample_results_FM.csv"
output_csv <- "final_results_FM.csv"

if (!file.exists(input_csv)) {
  stop(paste0("Can't find '", input_csv, "' - run Part 1 first."))
}

combined <- read.csv(input_csv)

if (!("n_events" %in% names(combined))) {
  stop("sub_sample_results_FM.csv has no n_events column - rerun Part 1 (Week4_part1_counting_FM.r).")
}

# ---- Number of events in each sub-sample (the same for every code in that file) ----
events_per_file <- tapply(combined$n_events, combined$sub_sample, function(x) x[1])
file_names <- names(events_per_file)

final_result <- do.call(rbind, lapply(split(combined, combined$code), function(group) {

  #Counts per file, with 0 for files where this code never appeared
  counts <- setNames(rep(0, length(file_names)), file_names)
  counts[group$sub_sample] <- group$total_count

  sub_sample_averages <- counts / events_per_file

  data.frame(
    code = group$code[1],
    name = group$name[1],
    total_count = sum(counts),
    average_per_event = signif(sum(counts) / sum(events_per_file), 6),
    uncertainty = signif(sd(sub_sample_averages), 3)   #Spread across the sub-samples
  )
}))

final_result <- final_result[order(-final_result$average_per_event), ]

options(scipen = 999)   #Disables scientific notation in the printout
cat("Final result: full-sample average, uncertainty from spread across",
    length(file_names), "sub-samples (", sum(events_per_file), "events ):\n\n")
print(final_result, row.names = FALSE)

write.csv(final_result, output_csv, row.names = FALSE)
cat("\nSaved final results to", output_csv, "\n")
