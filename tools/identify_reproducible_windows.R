library(tidyverse)

# Create output directories for figures and data.
dir.create("output/figures/unfiltered_reproducible_enriched_windows/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/unfiltered_reproducible_enriched_windows/", showWarnings = FALSE, recursive = TRUE)

# Command-line arguments:
args = commandArgs(trailingOnly=TRUE)
data_directory = args[1]
prefix = args[2]

# Collect all enriched window files for the given experiment prefix.
enriched_window_files = list.files(path = data_directory, pattern = paste0("^", prefix, "\\..*enriched_windows.tsv.gz"), full.names = TRUE)

# Read in enriched window data:
enriched_window_data = enriched_window_files %>%
    setNames(sub("\\.enriched_windows\\.tsv.gz", "", basename(.))) %>% 
    map(function(x) {
        read_tsv(x, col_types="cddcdcdddcddddcddccccccccc") %>%
            mutate(name = as.character(name))
    }) %>% 
    Filter(function(x) nrow(x) > 0, .) %>%
    bind_rows(.id = "clip_replicate_label")

# Handle case: no enriched windows across all replicates.
if (nrow(enriched_window_data) == 0){
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
	write_tsv(reproducible_enriched_window_data, paste0("output/unfiltered_reproducible_enriched_windows/", prefix, ".unfiltered_reproducible_enriched_windows.tsv.gz"))
	quit()
}

# Handle case: only single-replicate enrichment (no overlap).
if (nrow(enriched_window_data %>% group_by(name) %>% filter(n() > 1)) == 0) {
	# Save structure-matching empty tibble
	reproducible_enriched_window_data = enriched_window_data %>% group_by(across(-c(clip_replicate_label,baseline_l2or,input,clip,enrichment_l2or,pvalue,qvalue))) %>% summarize %>% head(0)
	write_tsv(reproducible_enriched_window_data, paste0("output/unfiltered_reproducible_enriched_windows/", prefix, ".unfiltered_reproducible_enriched_windows.tsv.gz"))
	quit()
}	

# Aggregate reproducible enriched windows across replicates.
reproducible_enriched_window_data = enriched_window_data %>%
	group_by(across(-c(clip_replicate_label,baseline_l2or,input,clip,enrichment_l2or,pvalue,qvalue))) %>%
	summarize(
		input_sum = sum(input),                               # total input reads
		clip_sum = sum(clip),                                 # total CLIP reads
		enrichment_n = sum(qvalue < 0.2),                     # number of replicates enriched (q<0.2)
		enrichment_l2or_min = min(enrichment_l2or),           # min log2 enrichment
		enrichment_l2or_mean = mean(enrichment_l2or),         # mean log2 enrichment
		enrichment_l2or_max = max(enrichment_l2or),           # max log2 enrichment
		p_max = max(pvalue),                                  # worst p-value
		p_min = min(pvalue),                                  # best p-value
		q_max = max(qvalue),                                  # worst q-value
		q_min = min(qvalue)                                   # best q-value
	) %>%
	filter(enrichment_n > 1) %>%                             # require reproducibility across replicates
	arrange(desc(enrichment_l2or_mean))                      # rank by strongest enrichment

# Save reproducible enriched window data.
write_tsv(reproducible_enriched_window_data, paste0("output/unfiltered_reproducible_enriched_windows/", prefix, ".unfiltered_reproducible_enriched_windows.tsv.gz"))
