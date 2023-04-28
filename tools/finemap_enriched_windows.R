library(tidyverse)
set.seed(5)

args = commandArgs(trailingOnly=TRUE)

count_data = read_tsv(args[1], col_names = c("chr","start","end","name","score","strand","window_n","input","clip"), col_types = c("ciiciciii"))
output_directory = args[2]
output_stem = args[3]
window_size = 75

dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

if(nrow(count_data) == 0) {
        file.create(paste0(output_directory, "/", output_stem, ".finemapped_windows.bed.gz"))
        quit()
}

# memory issues on joining in certain R versions

# aggregate <window_size> nt sliding windows to smooth out coverage
smoothed_data = count_data %>% group_by(window_n,strand,chr) %>% 
        do(tibble(
                pos = as.integer(round(convolve(.$end, rep(1,window_size), type = "filter") / window_size)), 
                input_sum = convolve(.$input, rep(1,window_size), type = "filter"), 
                clip_sum = convolve(.$clip, rep(1,window_size), type = "filter"))
        ) %>% 
        mutate(clip_sum = round(clip_sum), input_sum = round(input_sum))

# function to find local maxima
# https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
localMaxima = function(x) {
  if(length(x) == 1) {
    return(1)
  }
  # Use -Inf instead if x is numeric (non-integer)
  y = diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y = cumsum(rle(y)$lengths)
  y = y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y = y[-1]
  }
  y
}

# pick the most enriched <window_size nt> window per enriched feature and retain other local maxima above the median enrichment
# mark overlapping windows as claimed
set.seed(0)
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

# while nonoverlapping local maxima remain, pick additional <window_size> nt windows in order of enrichment level
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
