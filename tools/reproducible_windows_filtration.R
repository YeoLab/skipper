# Load necessary packages. 
library(glue)
library(tidyverse)
library(DescTools)
library(data.table)

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
windows = read_tsv(input_window_path)
nucs = read_tsv(input_nt_path, col_names = c("chr","start","end","name","score","strand","window_n","input","clip"), col_types = c("ciiciciii"))

# Check to make sure it is not an instance of 0 reproducible windows. 
if (nrow(windows) == 0){
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
	reproducible_enriched_windows = data.frame(matrix(nrow = 0, ncol = length(columns))) 
	colnames(reproducible_enriched_windows) = columns

	# Save empty table and exit
	write_tsv(reproducible_enriched_windows, paste0("output/reproducible_enriched_windows/", prefix, ".reproducible_enriched_windows.tsv.gz"))
	quit()
}	

# Filter out the blacklist
windows = windows %>% anti_join(blacklist %>% select(-name)) 

# Convert to data.table
dt_windows <- as.data.table(windows)
dt_nucs    <- as.data.table(nucs)

# Set keys (required for foverlaps)
setkey(dt_windows, chr, start, end)
setkey(dt_nucs,    chr, start, end)

# Overlap join: which window contains each nucleotide?
res <- foverlaps(
  dt_nucs,
  dt_windows,
  by.x = c("chr", "start", "end"),
  type = "within",  # nucleotide interval must lie within window
  mult = "first"    # windows don't overlap, so there is at most 1 anyway
)

# Clean up: keep original nuc pos + window info
nucs_with_window <- res[
  , .(
    chr  = chr,
    start = i.start,
    end = i.end,
    name,
    input,
    clip,
    window_start = start,
    window_end   = end
  )
]

# Remove all nucleotides that appear in no windows. 
nucs_with_window_filt <- na.omit(nucs_with_window)

# Calculate the gini-coefficient of each window from the nt counts. 
window_metrics <- nucs_with_window_filt %>% 
    group_by(name) %>% 
    summarize(
    gini     = DescTools::Gini(clip),
    .groups  = "drop"
    )

# Merge the ratio information into the window data
window_data = left_join(windows, window_metrics, by = "name")

# Use the ratios to filter the data. 
windows_filtered <- window_data %>%
  filter(gini < 0.9)

# Save the resulting filtered data. 
write_tsv(windows_filtered, paste0('output/reproducible_enriched_windows/', prefix, '.reproducible_enriched_windows.tsv.gz'))

# Use the ratios to filter the data. 
discarded_windows <- window_data %>%
  filter(gini >= 0.9)

# Save the resulting filtered data. 
write_tsv(discarded_windows, paste0('output/filtered_out_windows/', prefix, '.filtered_out_windows.tsv.gz'))

# Plot reproducible enriched window counts (linear scale).
pdf(paste0('output/figures/reproducible_enriched_windows/', prefix, '.reproducible_enriched_window_counts.linear.pdf'), height = 1.8, width = 2.2)
windows_filtered %>% 
	mutate(feature_group = sub("_.*","", feature_type_top)) %>% 
ggplot(aes(feature_type_top, fill = feature_group, group = feature_type_top)) + theme_bw(base_size = 7) + 
	geom_bar() + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	xlab("Type of feature") + ylab("# enriched")
dev.off()

# Plot reproducible enriched window counts (log scale).
pdf(paste0('output/figures/reproducible_enriched_windows/', prefix, '.reproducible_enriched_window_counts.log10.pdf'), height = 1.8, width = 2.2)
windows_filtered %>% 
	mutate(feature_group = sub("_.*","", feature_type_top)) %>% 
	group_by(feature_group, feature_type_top) %>% count %>%
ggplot(aes(feature_type_top, n, color = feature_group, group = feature_type_top)) + theme_bw(base_size = 7) + 
	geom_point(stroke = 0) + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	xlab("Type of feature") + ylab("# enriched") + scale_y_log10()
dev.off()
