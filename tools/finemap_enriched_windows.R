# Load tidyverse for data manipulation and I/O helpers (readr, dplyr, tibble, etc.). 
library(tidyverse)
set.seed(5) # Global seed for reproducibility of any stochastic steps in this script. 

# Parse command-line arguments: 
# args[1] = path to counts TSV with fixed column order, args[2] = output directory, args[3] = output stem. 
args = commandArgs(trailingOnly=TRUE)

# Read sliding-window count data with explicit column names and types to avoid guessing issues. 
count_data = read_tsv(args[1], col_names = c("chr","start","end","name","score","strand","window_n","input","clip"), col_types = c("ciiciciii"))
output_directory = args[2]
output_stem = args[3]
window_size = 75 # Smoothing and peak selection window size (nt). 

# Ensure the output directory exists. 
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

# If there is no data, emit an empty file and stop early. 
if(nrow(count_data) == 0) {
        file.create(paste0(output_directory, "/", output_stem, ".finemapped_windows.bed.gz"))
        quit()
}

# Note: some R versions have memory issues on large joins; computations below avoid unnecessary joins where possible. 

# Aggregate <window_size> nt sliding windows to smooth coverage by convolving counts and averaging positions. 
# The convolution with a length-window_size boxcar computes running sums, then clip/input are rounded to integers. 
smoothed_data = count_data %>% group_by(window_n,strand,chr) %>% 
        do(tibble(
                pos = as.integer(round(convolve(.$end, rep(1,window_size), type = "filter") / window_size)), 
                input_sum = convolve(.$input, rep(1,window_size), type = "filter"), 
                clip_sum = convolve(.$clip, rep(1,window_size), type = "filter"))
        ) %>% 
        mutate(clip_sum = round(clip_sum), input_sum = round(input_sum))

# Helper to find indices of local maxima in a vector using sign changes in the first derivative. 
# Returns positions of increasing-to-decreasing transitions, with edge handling for ties at the start. 
# Source reference: https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima. 
localMaxima = function(x) {
  if(length(x) == 1) {
    return(1)
  }
  # Use -Inf alternative when x is numeric; here an integer sentinel is used for speed. 
  y = diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y = cumsum(rle(y)$lengths)
  y = y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y = y[-1]
  }
  y
}

# Select finemapped candidate windows within each (chr, strand, window_n) track. 
# 1) Compute a heuristic enrichment score per smoothed position. 
# 2) Identify the top position, break ties randomly, and mark a ±window_size claim region. 
# 3) Keep local maxima above the median enrichment to seed iterative selection. 
set.seed(0) # Seed for tie-breaking when multiple positions share the maximum. 
candidate_sites = smoothed_data %>% 
        mutate(window_n = window_n, pos = pos,enrichment_heuristic =log2((clip_sum+20) / (input_sum+20))) %>% 
        group_by(window_n,strand,chr) %>% 
        mutate(enrichment_median = median(enrichment_heuristic)) %>%
        mutate(top_pos = enrichment_heuristic == max(enrichment_heuristic)) %>% 
        mutate(pos_selected = ifelse(sum(top_pos) > 1, sample(pos[top_pos], 1), pos[top_pos])) %>% 
        mutate(finemapped = pos == pos_selected, rank = ifelse(pos == pos_selected, 1, NA)) %>%
        # mutate(principal_window = (pos >= pos_principal - 12) & (pos <= pos_principal + 12), 
        mutate(claimed = ((pos >= pos_selected - window_size) & (pos <= pos_selected + window_size))) %>% 
        mutate(local_maximum = row_number() %in% localMaxima(enrichment_heuristic) | finemapped) %>%
        filter(enrichment_heuristic > enrichment_median, local_maximum) 

# Iteratively add additional non-overlapping local maxima, ranking by descending enrichment and claiming a ±window_size region each round. 
round = 1
finemapped_data = candidate_sites %>% filter(finemapped)
while(nrow(filter(candidate_sites, !claimed)) > 0) {
        round = round + 1 
        candidate_sites = candidate_sites %>% filter(!claimed) %>%
                mutate(top_pos = enrichment_heuristic == max(enrichment_heuristic)) %>%
                mutate(pos_selected = ifelse(sum(top_pos) == 1, pos[top_pos], sample(pos[top_pos],1))) %>%
                mutate(finemapped = pos == pos_selected, rank = replace(rank, pos == pos_selected, round)) %>%
                mutate(claimed = (pos >= pos_selected - window_size) & (pos <= pos_selected + window_size)) 

        finemapped_data = bind_rows(finemapped_data, filter(candidate_sites, finemapped))
}

# Convert ranked finemapped positions into fixed-width BED windows, attach counts and a simple CDF of ranks, and write sorted output. 
finemapped_data %>% mutate(cdf = 1 - (rank - 1) / max(rank)) %>% ungroup %>% filter(finemapped) %>%
        transmute(
                chr, 
                start = pos - (window_size + 1) / 2, 
                end = pos + (window_size -1) / 2, 
                name = paste0("MW:", window_n, ":", rank), 
                score = 0, 
                strand, 
                input = replace(input_sum, input_sum == 0, 0), 
                clip = replace(clip_sum, clip_sum == 0, 0),
                cdf
        ) %>% ungroup %>% arrange(chr, start) %>%
        write_tsv(paste0(output_directory, "/", output_stem, ".finemapped_windows.bed.gz"),col_names = FALSE)
