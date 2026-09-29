# PRA2003 - Week 4 Deliverable (Part 1 of 3)
# Runs the Week 3 analysis on data files output-Set1.txt through
# output-Set10.txt (output-Set0.txt was a test file and is excluded),
# and saves each file's per-code average-per-event into one combined CSV.
# This CSV is then used by the second script to calculate the overall
# result using the sub-sampling method: average the 10 sub-sample
# means, and use their standard deviation as the uncertainty.
#
# HOW TO RUN:
#   Make sure output-Set1.txt through output-Set10.txt are sitting
#   in the same folder as this script.
#
#   If you're in VS Code: Just press Run Code (Ctrl + Alt + N)
#
# WHAT IT DOES:
#   Each event in the file looks like this:
#     - a header line: event_id  number_of_particles
#     - followed by that many rows of: px py pz code
#
#   The file is read in chunks rather than all at once. Instead of
#   storing every particle code in memory (25 million numbers, growing
#   the vector every chunk), each chunk is tallied on its own with
#   table(), and only the running COUNTS per code are kept - a handful
#   of numbers, not one per particle. This is both faster and far
#   lighter on memory.
#
#   For each code we calculate:
#     - the average number of that particle per event
#     - the statistical uncertainty, using Poisson statistics:
#       for a total count N, sigma = sqrt(N), so the uncertainty
#       on the average (N / n_events) is sqrt(N) / n_events
#
#   Known biology-track particle codes are given readable names;
#   any other code encountered is labeled "Unknown".

# ---- Settings ----
chunk_size <- 1000000                     #Number of lines read at a time
output_csv <- "sub_sample_results_FM.csv"    #Where the combined per-file results are saved

# ---- Known strain codes and their names ----
known_codes <- c(211, -211, 321, -321, 2212, -2212,
                  3122, -3122, 3312, -3312, 3334, -3334)
known_names <- c("E. coli WT", "E. coli mutant",
                  "Bacillus subtilis WT", "Bacillus subtilis mutant",
                  "Pseudomonas aeruginosa WT", "Pseudomonas aeruginosa antibiotic-resistant",
                  "Streptococcus pneumoniae", "Capsule-deficient S. pneumoniae",
                  "Mycobacterium tuberculosis", "Drug-resistant M. tuberculosis",
                  "Salmonella enterica", "Salmonella mutant")

#Protection: known_codes and known_names must line up 1-to-1, or the
#name lookup later would silently attach the wrong name to a code.
if (length(known_codes) != length(known_names)) {
  stop("known_codes and known_names must be the same length.")
}

analyse_all_particles <- function(filename, chunk_size) {
  
  # ---- Initialize everything used below, before the loop starts ----
  start_time <- Sys.time()
  n_events <- 0
  running_counts <- setNames(numeric(0), character(0))  #code (as text) -> running total count
  malformed_lines <- 0
  chunk_number <- 0
  
  cat("Started at:", format(start_time), "\n\n")
  
  #Fail early with a clear message
  if (!file.exists(filename)) {
    stop(paste0("Can't find the file at '", filename, "'.\n",
                "Double check the data file is there, or update 'filename'."))
  }
  if (!is.numeric(chunk_size) || chunk_size <= 0) {
    stop("chunk_size must be a positive number.")
  }
  
  con <- file(filename, "r")
  on.exit(close(con))   #Closes the connection even if something errors out mid-loop
  
  repeat {
    chunk <- readLines(con, n = chunk_size)
    if (length(chunk) == 0) {
      break
    }
    chunk_number <- chunk_number + 1
    
    #The next line can be used for debugging, to see if the system is running correctly, or if it is stuck somewhere
    #cat("Processing chunk", chunk_number, "-", format(Sys.time()), "\n")
    
    split_lines <- strsplit(trimws(chunk), "\\s+")
    line_lengths <- lengths(split_lines)
    
    #Anything that isn't a 2-field header or a 4-field particle row is
    #malformed (blank line, corrupted row, etc.) - skip it, but count it.
    malformed_lines <- malformed_lines + sum(!(line_lengths %in% c(2, 4)))
    
    #Only events with at least 1 particle are counted; empty events
    #(header "N 0") contain no measurement. To count them too, use instead:
    #n_events <- n_events + sum(line_lengths == 2)
    header_lines <- split_lines[line_lengths == 2]
    n_events <- n_events + sum(sapply(header_lines, function(h) as.numeric(h[2]) > 0), na.rm = TRUE)
    
    particle_lines <- split_lines[line_lengths == 4]
    if (length(particle_lines) > 0) {
      particle_matrix <- matrix(as.numeric(unlist(particle_lines)), ncol = 4, byrow = TRUE)
      chunk_codes <- particle_matrix[, 4]
      chunk_codes <- chunk_codes[!is.na(chunk_codes)]  #Drop any code that failed to parse as a number
      
      #Tally just this chunk, then fold it into the running total - this
      #keeps memory use to "one entry per unique code", not per particle.
      chunk_table <- table(chunk_codes)
      new_codes <- setdiff(names(chunk_table), names(running_counts))
      if (length(new_codes) > 0) {
        running_counts[new_codes] <- 0
      }
      running_counts[names(chunk_table)] <- running_counts[names(chunk_table)] + as.numeric(chunk_table)
    }
  }
  
  if (malformed_lines > 0) {
    cat("\nNote:", malformed_lines, "line(s) were skipped (neither a valid header nor a valid particle row).\n")
  }
  
  #Protection against dividing by zero if the file had no valid events.
  if (n_events == 0) {
    stop("No events found in file - check the file path/format.")
  }
  
  average_per_event <- running_counts / n_events
  uncertainty <- sqrt(running_counts) / n_events
  
  codes_numeric <- as.numeric(names(running_counts))
  matched_index <- match(codes_numeric, known_codes)
  particle_name <- ifelse(is.na(matched_index), "Unknown", known_names[matched_index])
  
  results <- data.frame(
    code = codes_numeric,
    name = particle_name,
    n_events = n_events,   #Saved so Part 2 can recompute exact averages without rounding
    total_count = as.numeric(running_counts),
    #5 significant digits, so rare codes do not round to 0
    average_per_event = signif(as.numeric(average_per_event), 5),
    uncertainty = signif(as.numeric(uncertainty), 5)
  )
  
  results <- results[order(-results$total_count), ]   #Sort by total count, descending (highest first)
  
  options(scipen = 999)   #Disables scientific notation for this session
  print(results)
  
  end_time <- Sys.time()
  cat("Done! Processed", n_events, "events across", chunk_number, "chunk(s).\n")
  cat("Total time elapsed:", format(end_time - start_time), "\n\n")

  return(results)
}

# ---- Run the function on all 10 files (output-Set1.txt through output-Set10.txt) ----
file_names <- paste0("output-Set", 1:10, ".txt")

all_sub_sample_results <- list()

for (filename in file_names) {
  cat("=== Analysing", filename, "===\n")
  sub_sample_result <- analyse_all_particles(filename, chunk_size)
  sub_sample_result$sub_sample <- filename   #Tag which file this result came from
  all_sub_sample_results[[filename]] <- sub_sample_result
}

# ---- Combine all 10 results into one table and save to CSV ----
combined_results <- do.call(rbind, all_sub_sample_results)
write.csv(combined_results, output_csv, row.names = FALSE)
cat("\nSaved combined per-sub-sample results to", output_csv, "\n")