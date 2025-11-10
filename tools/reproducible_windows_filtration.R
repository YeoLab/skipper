# Load necessary packages. 
library(glue)
library(tidyverse)
library(DescTools)

# Create output directories for figures and data.
dir.create("output/figures/reproducible_enriched_windows/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/reproducible_enriched_windows/", showWarnings = FALSE, recursive = TRUE)

# Define the path and condition names for the data. 
args = commandArgs(trailingOnly=TRUE)
input_window_path = args[1]
input_nt_path = args[2]
prefix = args[3]

# Load blacklist if provided, otherwise create empty tibble.
if(length(args) > 3) {
	blacklist = read_tsv(args[4], col_names = c("chr","start","end","name","score","strand"), col_types = "cddcdc")
} else {
	blacklist = tibble(chr=character(),start=numeric(),end=numeric(),name=character(),score=numeric(),strand=character())
}

# Load in the window data, nucleotide counts data, and the blacklist. 
window_data = read_tsv(input_window_path)
nt_counts = read_tsv(input_nt_path, col_names = c("chr","start","end","name","score","strand","window_n","input","clip"), col_types = c("ciiciciii"))

# Check to make sure it is not an instance of 0 reproducible windows. 
if (nrow(window_data) == 0){
	# Output placeholder "No data" plots
	pdf(paste0('output/figures/reproducible_enriched_windows/', prefix, '.reproducible_enriched_window_counts.linear.pdf'), height = 1, width = 2)
		print(ggplot() + annotate("text", x = 1, y = 1, label = "No data") + theme_void())
	dev.off()
	pdf(paste0('output/figures/reproducible_enriched_windows/', prefix, '.reproducible_enriched_window_counts.log10.pdf'), height = 1, width = 2)
		print(ggplot() + annotate("text", x = 1, y = 1, label = "No data") + theme_void())
	dev.off()

    # Construct an empty dataframe with expected columns
	columns= c("chr","start","end","name","score","strand","gc",
	"gc_bin","chrom","feature_id","feature_bin","feature_type_top","feature_types",
	"gene_name","gene_id","transcript_ids","gene_type_top","transcript_type_top",
	"gene_types","transcript_types", "input_sum","clip_sum","enrichment_n",
	"enrichment_l2or_min","enrichment_l2or_mean","enrichment_l2or_max","p_max","p_min",
	"q_max","q_min") 
	reproducible_enriched_window_data = data.frame(matrix(nrow = 0, ncol = length(columns))) 
	colnames(reproducible_enriched_window_data) = columns

	# Save empty table and exit
	write_tsv(reproducible_enriched_window_data, paste0("output/reproducible_enriched_windows/", prefix, ".reproducible_enriched_windows.tsv.gz"))
	quit()
}	

# Filter out the blacklist
window_data = window_data %>% anti_join(blacklist %>% select(-name)) 

# Calculate the gini-coefficient of each window from the nt counts. 
window_metrics <- nt_counts %>% 
    group_by(window_n) %>% 
    summarize(
    gini     = DescTools::Gini(clip),
    .groups  = "drop"
    )
    
# add the ratios to nt_counts. 
nt_counts = left_join(nt_counts, window_metrics, by = "window_n")

# Create ID columns for both nt_counts and window_data for merging. 
window_data$ID = paste(window_data$chr, window_data$start, sep = "_")
nt_counts$ID = paste(nt_counts$chr, nt_counts$start, sep = "_")

# Merge the ratio information into the window data
window_data = left_join(window_data, dplyr::select(nt_counts, c("ID", "window_n", "gini")), by = "ID")

# Use the ratios to filter the data. 
window_data_filtered <- window_data %>%
  filter(gini < 0.95)

# Save the resulting filtered data. 
write_tsv(window_data_filtered, paste0('output/reproducible_enriched_windows/', prefix, '.reproducible_enriched_windows.tsv.gz'))

# Use the ratios to filter the data. 
discarded_windows <- window_data %>%
  filter(gini >= 0.95)

# Save the resulting filtered data. 
write_tsv(discarded_windows, paste0('output/filtered_out_windows/', prefix, '.filtered_out_windows.tsv.gz'))

# Plot reproducible enriched window counts (linear scale).
pdf(paste0('output/figures/reproducible_enriched_windows/', prefix, '.reproducible_enriched_window_counts.linear.pdf'), height = 1.8, width = 2.2)
window_data_filtered %>% 
	mutate(feature_group = sub("_.*","", feature_type_top)) %>% 
ggplot(aes(feature_type_top, fill = feature_group, group = feature_type_top)) + theme_bw(base_size = 7) + 
	geom_bar() + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	xlab("Type of feature") + ylab("# enriched")
dev.off()

# Plot reproducible enriched window counts (log scale).
pdf(paste0('output/figures/reproducible_enriched_windows/', prefix, '.reproducible_enriched_window_counts.log10.pdf'), height = 1.8, width = 2.2)
window_data_filtered %>% 
	mutate(feature_group = sub("_.*","", feature_type_top)) %>% 
	group_by(feature_group, feature_type_top) %>% count %>%
ggplot(aes(feature_type_top, n, color = feature_group, group = feature_type_top)) + theme_bw(base_size = 7) + 
	geom_point(stroke = 0) + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	xlab("Type of feature") + ylab("# enriched") + scale_y_log10()
dev.off()
