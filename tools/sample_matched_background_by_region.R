library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

enriched_windows = read_tsv(args[1])
all_windows = read_tsv(args[2])
fixed_window_size = as.integer(args[3])
output_directory = args[4]
output_stem = args[5]

dir.create(paste0(output_directory, "/variable/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_directory, "/fixed/"), showWarnings = FALSE, recursive = TRUE)

if(nrow(enriched_windows) == 0) {
	file.create(paste0(output_directory, "/variable/", output_stem, ".sampled_variable_windows.bed.gz"))
	file.create(paste0(output_directory, "/fixed/", output_stem, ".sampled_fixed_windows.bed.gz"))
	quit()	
}
enriched_windows$sampling_group = "Other"
enriched_windows[grepl("SS",enriched_windows$feature_types),"sampling_group"] = "SS"
enriched_windows[grepl("EXON_LNCRNA|EXON_SMALL",enriched_windows$feature_types),"sampling_group"] = "EXON_NC"
enriched_windows[grepl("EXON_MRNA",enriched_windows$feature_types),"sampling_group"] = "EXON_MRNA"

all_windows$sampling_group = "Other"
all_windows[grepl("SS",all_windows$feature_types),"sampling_group"] = "SS"
all_windows[grepl("EXON_LNCRNA|EXON_SMALL",all_windows$feature_types),"sampling_group"] = "EXON_NC"
all_windows[grepl("EXON_MRNA",all_windows$feature_types),"sampling_group"] = "EXON_MRNA"

# Aim to select at least 20000 windows for matched background
sampling_factor = 1 + as.integer(20000 / nrow(enriched_windows))

# Randomly sample region matched bins
sampled_data = enriched_windows %>% 
	group_by(sampling_group) %>% count(name="n_enriched") %>% 
	inner_join(all_windows %>% select(name, sampling_group)) %>% group_by(sampling_group) %>% 
	summarize(tibble(name = sample(name[!name %in% enriched_windows$name],pmin(sampling_factor*n_enriched,sum(!name %in% enriched_windows$name)) ))) %>%
	ungroup %>% select(name) %>% inner_join(all_windows,.)

write_tsv(sampled_data[,1:7], paste0(output_directory, "/variable/", output_stem, ".sampled_variable_windows.bed.gz"), col_names = FALSE)

sampled_data[,1:6] %>% rowwise %>% 
	mutate(start = ifelse(end - start == 1, start, sample(seq(start, end + 1), 1)), end = start + fixed_window_size) %>% ungroup %>%
write_tsv(paste0(output_directory, "/fixed/", output_stem, ".sampled_fixed_windows.bed.gz"), col_names = FALSE)
